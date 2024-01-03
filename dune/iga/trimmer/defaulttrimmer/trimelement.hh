// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>

#include <dune/iga/geometrykernel/findintersection.hh>

#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>

namespace Dune::IGANEW::DefaultTrim {

  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(
      const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
      const PatchTrimData& patchTrimData) {
    std::cout << "START " << std::endl;
    using namespace Clipper2Lib;
    auto geo = element.geometry();
    std::array<FieldVector<ctype, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    PathD elementPath;
    for (const auto i : std::views::iota(0, 4))
      elementPath.push_back({corners[i][0], corners[i][1], i});

    PathsD clipPaths;
    PathD tempPath;
    constexpr int N = 10;
    // @todo store param value of sampled points on trimmingCurve
    for (auto loop : patchTrimData.loops()) {
      tempPath.clear();
      for (int i = 100; auto& curve : loop.curves()) {
        for (auto v : Utilities::linspace(curve.domain()[0], N)) {
          auto fV = curve.global({v});
          tempPath.emplace_back(fV[0], fV[1], i++);
        }
        // Add an additional point just outside of the element
        auto localLastPoint = curve.domain()[0][1];
        auto lastPoint      = curve.global(localLastPoint);

        auto scale = element.geometry().volume() / 10;

        auto dx = curve.jacobian({localLastPoint}) * scale;
        tempPath.emplace_back(lastPoint[0] + dx[0], lastPoint[1] + dx[1], i++);

        i += 99 - N;
      }
      clipPaths.push_back(tempPath);
    }

    auto [flag, result] = Impl::clipElementRectangle(elementPath, clipPaths);

    // Create ElementTrimData with exact intersection Points and a geometry representation of the element edges
    ElementTrimData elementTrimData(flag);

    if (flag != ElementTrimFlag::trimmed) return elementTrimData;

    auto nextEntity  = [&](const int i) { return (i + 1) % result.vertices_.size(); };
    auto isNewVertex = [](const auto& vV) { return std::holds_alternative<Impl::ClippingResult::NewVertex>(vV); };
    auto getPt       = [](auto&& vV) { return vV.pt; };
    auto getHostIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::HostVertex>(vV).hostIdx; };
    auto getEdgeIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::NewVertex>(vV).onEdgeIdx; };
    auto getTrimmingCurveZ = [](const auto& vV) {return std::get<Impl::ClippingResult::NewVertex>(vV).trimmingCurveZ; };

    // @todo could be done via Jacobian of the edges
    std::array<FieldVector<ScalarType, dim>, 4> dirs {FieldVector<ScalarType, dim>{1.0, 0.0}, FieldVector<ScalarType, dim>{0, 1},
      FieldVector<ScalarType, dim>{-1, 0}, FieldVector<ScalarType, dim>{0, -1}};

    auto approxSamePoint = [](const Clipper2Lib::PointD& pt1, const Dune::FieldVector<ScalarType, dim>& pt2, const double prec) -> bool {
      return FloatCmp::eq(pt1.x, pt2[0], prec) and FloatCmp::eq(pt1.y, pt2[1], prec);
    };
    auto getTrimmingCurveIdx = [&](auto& vV) -> size_t {
      return static_cast<size_t>(std::floor((getTrimmingCurveZ(vV)/100)-1));
    };

    auto callFindIntersection = [&](const size_t tcIdx, const int edgeIdx, auto ip) -> std::pair<double, FieldVector<ScalarType, dim>> {
      // @todo gerneralize for more than one loop
      auto curvePatchGeo = patchTrimData.loops().front().curves()[tcIdx];
      auto pos = corners[edgeIdx];
      auto dir = dirs[edgeIdx];
      // @todo make better guess by maybe using z Val with lookup-table
      auto guessTParam = FieldVector<ScalarType, dim>{curvePatchGeo.domainMidPoint()[0], 0.5};

      auto [success, tParam, curvePoint] = findIntersectionCurveAndLine(curvePatchGeo, pos, dir, guessTParam);
      if (success == FindIntersectionCurveAndLineResult::sucess)
        return std::make_pair(tParam[0], curvePoint);
      if (success == FindIntersectionCurveAndLineResult::linesParallel) {
        // Check ip against first point of curve
        if (approxSamePoint(ip, curvePatchGeo.corner(0), 1e-4))
          return std::make_pair(curvePatchGeo.domain()[0].front(), curvePatchGeo.corner(0));
        DUNE_THROW(Dune::GridError, "Edge is parallel to curve");

      }

      // Check domain[0] and domain[1]
      auto domain = curvePatchGeo.domain();

      if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].front()})))
        return std::make_pair(domain[0][0], curvePoint);
      if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].back()})))
        return std::make_pair(domain[1][0], curvePoint);

      DUNE_THROW(Dune::GridError, "Couldn't find intersection Point");

    };

    for (const auto i : std::views::iota(0u, result.vertices_.size())) {
      auto vV1 = result.vertices_[i];
      auto vV2 = result.vertices_[nextEntity(i)];

      // First case edge is completly untrimmed
      if (!isNewVertex(vV1) and !isNewVertex(vV2)) {
        elementTrimData.addEdge(Impl::giveEdgeIdx(getHostIdx(vV1), getHostIdx(vV2)));
        continue;
      }

      auto pt1 = std::visit(getPt, vV1);
      auto pt2 = std::visit(getPt, vV2);

      // Second case edge begins on a hostVertes and ends on a newVertex
      if (!isNewVertex(vV1) and isNewVertex(vV2)) {
        // Find intersection point on edge(vV2)
        auto [tParam, curvePoint] = callFindIntersection( getTrimmingCurveIdx(vV2), getEdgeIdx(vV2), pt2);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam << std::endl;
        continue;
      }
      // Third case newVertex - newVertex
      if (isNewVertex(vV1), isNewVertex(vV2)) {
        // Find intersection point on edge(vV2)
        // First intersection point should be already available
        auto [tParam, curvePoint] = callFindIntersection( getTrimmingCurveIdx(vV2), getEdgeIdx(vV2), pt2);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam << std::endl;
        continue;
      }
      // Fourth case edge begins on a newVertex and ends in a HostVertex
      if (isNewVertex(vV1), !isNewVertex(vV2)) {
        // First intersection point should be already available
        continue;
      }
    }

    return elementTrimData;
  }
}  // namespace Dune::IGANEW::DefaultTrim
