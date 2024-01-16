// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>

#include <dune/iga/geometrykernel/findintersection.hh>
#include <dune/iga/geometrykernel/slicecurve.hh>
#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/trimutils.hh>

namespace Dune::IGANEW::DefaultTrim {

  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(const auto& element, const PatchTrimData& patchTrimData) {
    // std::cout << "START " << std::endl;
    auto geo = element.geometry();
    std::array<FieldVector<double, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    auto [flag, result] = Util::clipElementRectangle(element, patchTrimData);

    // Create ElementTrimData with exact intersection Points and a geometry representation of the element edges
    ElementTrimData elementTrimData(flag, element);

    if (flag != ElementTrimFlag::trimmed) return elementTrimData;

    /*
     * Visitors to get vertex properties
     */

    auto nextEntity  = [&](const int i) { return (i + 1) % result.vertices_.size(); };
    auto isNewVertex = [](const auto& vV) { return std::holds_alternative<Util::ClippingResult::NewVertex>(vV); };
    auto isHostVertex = [](const auto& vV) { return std::holds_alternative<Util::ClippingResult::HostVertex>(vV); };
    auto isInsideVertex = [](const auto& vV) { return std::holds_alternative<Util::ClippingResult::InsideVertex>(vV); };
    auto getPt       = [](auto&& vV) { return vV.pt; };
    auto getHostIdx  = [](const auto& vV) { return std::get<Util::ClippingResult::HostVertex>(vV).hostIdx; };
    auto getEdgeIdx  = [](const auto& vV) { return std::get<Util::ClippingResult::NewVertex>(vV).onEdgeIdx; };
    auto getCurveI = [](const auto& vV) {return std::get<Util::ClippingResult::InsideVertex>(vV).curveIdxI; };
    auto getCurveJ = [](const auto& vV) {return std::get<Util::ClippingResult::InsideVertex>(vV).curveIdxJ; };
    auto getLoopIdx = [](const auto& vV) {return std::get<Util::ClippingResult::InsideVertex>(vV).loopIdx; };
    auto getTrimmingCurveIdx = [&](auto& vV) -> std::pair<size_t, size_t> {
      return patchTrimData.getIndices(std::get<Util::ClippingResult::NewVertex>(vV).trimmingCurveZ);
    };

    ///
    // Here the actual code is starting
    ///

    // Major todo: Create fallback to straight line if alog is not able to find correct Trimming CurveIdx
    // Also todo: Make a second algo that just connects the points as lines (as template?)

    // State
    std::vector<FieldVector<ScalarType, dim>> foundVertices;
    bool isOnNewEdge          = false;
    std::pair currentCurveIdx = {std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity()};
    double currentT           = std::numeric_limits<double>::infinity();

    for (const auto i : std::views::iota(0u, result.vertices_.size())) {
      auto vV1 = result.vertices_[i];
      auto vV2 = result.vertices_[nextEntity(i)];

      auto pt1 = std::visit(getPt, vV1);
      auto pt2 = std::visit(getPt, vV2);

      // First case edge is completly untrimmed
      if (isHostVertex(vV1) and isHostVertex(vV2)) {
        elementTrimData.addEdge(Util::giveEdgeIdx(getHostIdx(vV1), getHostIdx(vV2)));
        foundVertices.push_back({pt2.x, pt2.y});
        continue;
      }

      // Second case edge begins on a hostVertes and ends on a newVertex
      if (isHostVertex(vV1) and isNewVertex(vV2)) {
        currentCurveIdx = getTrimmingCurveIdx(vV2);
        auto [tParam, curvePoint]
                = Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), getEdgeIdx(vV2), pt2, corners);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam
                  << std::endl;

        FieldVector<ScalarType, dim> p
            = foundVertices.empty() ? FieldVector<ScalarType, dim>{pt1.x, pt1.y} : foundVertices.back();
        auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, curvePoint);
        elementTrimData.addEdgeHostNew(getEdgeIdx(vV2), trimmedEdge, curvePoint);
        foundVertices.push_back(curvePoint);

        isOnNewEdge = true;
        currentT    = tParam;
        continue;
      }
      // Third case newVertex - newVertex
      if (isNewVertex(vV1) and isNewVertex(vV2)) {
        if (foundVertices.empty()) {
          currentCurveIdx = getTrimmingCurveIdx(vV1);
          FieldVector<ScalarType, dim> curvePoint;
          std::tie(currentT, curvePoint)
              = Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), getEdgeIdx(vV1), pt1, corners);
          std::cout << "Found: " << curvePoint << " From Clipping: " << pt1.x << " " << pt1.y << " t: " << currentT
                    << std::endl;
          foundVertices.push_back(curvePoint);
        }
        assert(getTrimmingCurveIdx(vV2) == currentCurveIdx);

        auto [tParam, curvePoint]
        = Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), getEdgeIdx(vV2), pt2, corners);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam
                  << std::endl;

        if (currentT > tParam) {
          assert(getEdgeIdx(vV1) == getEdgeIdx(vV2));
          auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(foundVertices.back(), curvePoint);
          elementTrimData.addEdgeNewNewOnHost(getEdgeIdx(vV2), trimmedEdge, curvePoint);
          isOnNewEdge = true;
          currentT    = tParam;
        } else {
          auto elementTrimmingCurve
              = Util::createTrimmingCurveSlice(patchTrimData.getCurve(currentCurveIdx), currentT, tParam);
          elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
          isOnNewEdge = false;
          currentT    = std::numeric_limits<double>::infinity();
        }
        foundVertices.push_back(curvePoint);

        continue;
      }
      // Fourth case edge begins on a newVertex and ends in a HostVertex
      if (isNewVertex(vV1), isHostVertex(vV2)) {
        auto v2 = Dune::FieldVector<ScalarType, dim>{pt2.x, pt2.y};

        FieldVector<ScalarType, dim> p;
        if (foundVertices.empty())
          std::tie(std::ignore, p) = Util::callFindIntersection(patchTrimData.getCurve(getTrimmingCurveIdx(vV1)),
                                                                getEdgeIdx(vV1), pt1, corners);
        else
          p = foundVertices.back();

        auto trimmedEdge = Util::createHostGeometry<TrimmingCurve>(p, v2);
        elementTrimData.addEdgeNewHost(getEdgeIdx(vV1), trimmedEdge, getHostIdx(vV2));

        foundVertices.push_back({pt2.x, pt2.y});
        continue;
      }
      // Additional cases to cover inside vertices
      if (isNewVertex(vV1), isInsideVertex(vV2)) {
        if (foundVertices.empty()) {
          currentCurveIdx = getTrimmingCurveIdx(vV1);
          FieldVector<ScalarType, dim> curvePoint;
          std::tie(currentT, curvePoint)
              = Util::callFindIntersection(patchTrimData.getCurve(currentCurveIdx), getEdgeIdx(vV1), pt1, corners);
          std::cout << "Found: " << curvePoint << " From Clipping: " << pt1.x << " " << pt1.y << " t: " << currentT
                    << std::endl;
          foundVertices.push_back(curvePoint);
        }
        const auto& curve = patchTrimData.loops()[getLoopIdx(vV2)].curves()[getCurveI(vV2)];
        double tParam = curve.domain().front().back();
        auto elementTrimmingCurve = Util::createTrimmingCurveSlice(curve, currentT, tParam);
        auto curvePoint = curve.corner(1);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        foundVertices.push_back(curvePoint);
        continue;
      }
      if (isInsideVertex(vV1), isNewVertex(vV2)) {
        const auto& curve = patchTrimData.loops()[getLoopIdx(vV1)].curves()[getCurveJ(vV1)];
        currentT = curve.domain().front().front();
        auto [tParam, curvePoint]
            = Util::callFindIntersection(curve, getEdgeIdx(vV2), pt2, corners);
        std::cout << "Found: " << curvePoint << " From Clipping: " << pt2.x << " " << pt2.y << " t: " << tParam
                  << std::endl;
        auto elementTrimmingCurve
              = Util::createTrimmingCurveSlice(curve, currentT, tParam);
        elementTrimData.addEdgeNewNew(elementTrimmingCurve, curvePoint);
        foundVertices.push_back(curvePoint);
        continue;
      }
    }
    elementTrimData.finalize();
    return elementTrimData;
  }
}  // namespace Dune::IGANEW::DefaultTrim
