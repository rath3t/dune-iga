// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>

#include <dune/iga/geometrykernel/findintersection.hh>
#include <dune/iga/geometrykernel/slicecurve.hh>
#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>

namespace Dune::IGANEW::DefaultTrim {

  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(
      const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
      const PatchTrimData& patchTrimData) {
    std::cout << "START " << std::endl;
    auto geo = element.geometry();
    std::array<FieldVector<double, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    auto [flag, result] = Impl::clipElementRectangle(element, patchTrimData);

    // Create ElementTrimData with exact intersection Points and a geometry representation of the element edges
    ElementTrimData elementTrimData(flag);

    if (flag != ElementTrimFlag::trimmed) return elementTrimData;

    auto nextEntity  = [&](const int i) { return (i + 1) % result.vertices_.size(); };
    auto isNewVertex = [](const auto& vV) { return std::holds_alternative<Impl::ClippingResult::NewVertex>(vV); };
    auto getPt       = [](auto&& vV) { return vV.pt; };
    auto getHostIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::HostVertex>(vV).hostIdx; };
    auto getEdgeIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::NewVertex>(vV).onEdgeIdx; };
    auto getTrimmingCurveZ
        = [](const auto& vV) { return std::get<Impl::ClippingResult::NewVertex>(vV).trimmingCurveZ; };

    // @todo could be done via Jacobian of the edges
    std::array<FieldVector<ScalarType, dim>, 4> dirs{
        FieldVector<ScalarType, dim>{1.0, 0.0}, FieldVector<ScalarType, dim>{0, 1}, FieldVector<ScalarType, dim>{-1, 0},
        FieldVector<ScalarType, dim>{0, -1}};

    auto approxSamePoint
        = [](const Clipper2Lib::PointD& pt1, const Dune::FieldVector<ScalarType, dim>& pt2, const double prec) -> bool {
      return FloatCmp::eq(pt1.x, pt2[0], prec) and FloatCmp::eq(pt1.y, pt2[1], prec);
    };
    auto getTrimmingCurveIdx
        = [&](auto& vV) -> size_t { return static_cast<size_t>(std::floor((getTrimmingCurveZ(vV) / 100) - 1)); };

    auto distance = [](const Clipper2Lib::PointD& p1, const auto& p2) -> double {
      return std::hypot(p1.x - p2[0], p1.y - p2[1]);
    };

    auto findGoodStartingPoint = [&](const auto& curve, const Clipper2Lib::PointD& pt, int N = 10) -> double {
      auto linSpace = Utilities::linspace(curve.domain().front(), N);
      std::vector<double> distances;
      std::ranges::transform(linSpace, std::back_inserter(distances), [&](const auto u) {
        return distance(pt, curve.global(u));
      });
      auto min_idx = std::ranges::distance(distances.begin(), std::ranges::min_element(distances));
      return linSpace[min_idx];
    };

    auto callFindIntersection
        = [&](const size_t tcIdx, const int edgeIdx, auto ip) -> std::pair<double, FieldVector<ScalarType, dim>> {
      // @todo gerneralize for more than one loop
      auto curvePatchGeo = patchTrimData.loops().front().curves()[tcIdx];
      auto pos           = corners[edgeIdx];
      auto dir           = dirs[edgeIdx];
      // @todo make better guess by maybe using z Val with lookup-table
      auto guessTParam = FieldVector<ScalarType, 2>{findGoodStartingPoint(curvePatchGeo, ip), 0.5};

      auto [success, tParam, curvePoint] = findIntersectionCurveAndLine(curvePatchGeo, pos, dir, guessTParam);
      if (success == FindIntersectionCurveAndLineResult::sucess) return std::make_pair(tParam[0], curvePoint);
      if (success == FindIntersectionCurveAndLineResult::linesParallel) {
        DUNE_THROW(Dune::GridError, "Couldn't find intersection Point, lines are parallel");
      }

      assert(approxSamePoint(ip, curvePoint, 1e-4) && "Intersection Point and Curve Point are not the same");

      auto domain = curvePatchGeo.domain();
      if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].front()})))
        return std::make_pair(domain[0][0], curvePoint);
      if (FloatCmp::eq(curvePoint, curvePatchGeo.global({domain[0].back()})))
        return std::make_pair(domain[0].back(), curvePoint);

      DUNE_THROW(Dune::GridError, "Couldn't find intersection Point");
    };

    auto createHostGeometry = [&](auto& vertex1, auto& vertex2) -> TrimmingCurve {
      const std::array<std::vector<double>, 1> knotSpans = {{{0, 0, 1, 1}}};

      const std::vector<typename TrimmingCurve::ControlPointType> controlPoints
          = {{{.p = vertex1, .w = 1}, {.p = vertex2, .w = 1}}};
      auto controlNet = NURBSPatchData<1, 2>::ControlPointNetType(controlPoints);
      typename TrimmingCurve::PatchData patchData(knotSpans, controlNet, {1});

      return GeometryKernel::NURBSPatch(patchData);
    };

    auto createTrimmingCurveSlice = [&](auto& curveIdx, double t1, double t2) -> TrimmingCurve {
      return sliceCurve(patchTrimData.loops().front().curves()[curveIdx], {t1, t2});
    };

    ///
    // Here the actual code is starting
    ///

    // State
    std::vector<FieldVector<ScalarType, dim>> foundVertices;
    bool isOnNewEdge       = false;
    size_t currentCurveIdx = std::numeric_limits<size_t>::infinity();
    double currentT        = std::numeric_limits<double>::infinity();

    for (const auto i : std::views::iota(0u, result.vertices_.size())) {
      auto vV1 = result.vertices_[i];
      auto vV2 = result.vertices_[nextEntity(i)];

      auto pt1 = std::visit(getPt, vV1);
      auto pt2 = std::visit(getPt, vV2);

      // First case edge is completly untrimmed
      if (!isNewVertex(vV1) and !isNewVertex(vV2)) {
        elementTrimData.addEdge(Impl::giveEdgeIdx(getHostIdx(vV1), getHostIdx(vV2)));
        foundVertices.push_back({pt2.x, pt2.y});
        continue;
      }

      // Second case edge begins on a hostVertes and ends on a newVertex
      if (!isNewVertex(vV1) and isNewVertex(vV2)) {
        currentCurveIdx           = getTrimmingCurveIdx(vV2);
        auto [tParam, curvePoint] = callFindIntersection(currentCurveIdx, getEdgeIdx(vV2), pt2);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam
                  << std::endl;

        FieldVector<ScalarType, dim> p
            = foundVertices.empty() ? FieldVector<ScalarType, dim>{pt1.x, pt1.y} : foundVertices.back();
        auto trimmedEdge = createHostGeometry(p, curvePoint);
        elementTrimData.addEdgeHostNew(getEdgeIdx(vV2), trimmedEdge, curvePoint);
        foundVertices.push_back(curvePoint);

        isOnNewEdge = true;
        currentT    = tParam;
        continue;
      }
      // Third case newVertex - newVertex
      if (isNewVertex(vV1), isNewVertex(vV2)) {
        assert(getTrimmingCurveIdx(vV2) == currentCurveIdx);

        auto [tParam, curvePoint] = callFindIntersection(currentCurveIdx, getEdgeIdx(vV2), pt2);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam
                  << std::endl;

        auto elementTrimmingCurve = createTrimmingCurveSlice(currentCurveIdx, currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        foundVertices.push_back(curvePoint);

        isOnNewEdge = false;
        currentT    = std::numeric_limits<double>::infinity();

        continue;
      }
      // Fourth case edge begins on a newVertex and ends in a HostVertex
      if (isNewVertex(vV1), !isNewVertex(vV2)) {
        auto v2 = Dune::FieldVector<ScalarType, dim>{pt2.x, pt2.y};

        FieldVector<ScalarType, dim> p;
        if (foundVertices.empty())
          std::tie(std::ignore, p) = callFindIntersection(getTrimmingCurveIdx(vV1), getEdgeIdx(vV1), pt1);
        else
          p = foundVertices.back();

        auto trimmedEdge = createHostGeometry(p, v2);
        elementTrimData.addEdgeNewHost(getEdgeIdx(vV1), trimmedEdge, getHostIdx(vV2));

        foundVertices.push_back({pt2.x, pt2.y});
        continue;
      }
    }

    return elementTrimData;
  }
}  // namespace Dune::IGANEW::DefaultTrim
