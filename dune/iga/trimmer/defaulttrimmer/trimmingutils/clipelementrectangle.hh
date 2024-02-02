// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>

#include "cliputils.hh"

#include <dune/common/float_cmp.hh>

#include <dune/iga/geometrykernel/findintersection.hh>

namespace Dune::IGANEW::DefaultTrim::Util {

  auto clipElementRectangle(const auto& element, const auto& patchTrimData)
      -> std::tuple<ElementTrimFlag, ClippingResult> {
    using namespace Clipper2Lib;

    auto eleGeo = element.geometry();
    static constexpr int numberOfCorners = 4;

    PathD eleRect;
    for (const auto vertexIndex : Dune::range(numberOfCorners)) {
      auto corner = eleGeo.corner(vertexIndexMapping[vertexIndex]);
      eleRect.emplace_back(corner[0], corner[1], vertexIndex);
    }

    const auto& trimmingCurves = patchTrimData.clipperLoops();

    // @todo check z value if same as eleRect -> full
    auto isFullElement
        = [&](const auto& clippedEdges) { return FloatCmp::eq(Area(clippedEdges), Area(eleRect), 1e-8); };

    // First determine if element is trimmed
    const auto intersectResult = Intersect({eleRect}, trimmingCurves, FillRule::NonZero, 8);

    if (intersectResult.empty()) return std::make_tuple(ElementTrimFlag::empty, ClippingResult{eleRect});
    if (isFullElement(intersectResult.front())) {
      assert(intersectResult.size()==1 && "If the element is full, there should be only one path");
      return std::make_tuple(ElementTrimFlag::full, ClippingResult{eleRect});
    }

    // Element is trimmed
    ClipperD clipper(8);
    ClippingResult result(eleRect);

    auto checkParallel = [&](const auto& curve, const int edgeIndex) -> bool {
      return curve.affine()
             and std::get<0>(findIntersectionLinearCurveAndLine(curve, eleGeo.corner(vertexIndexMapping[edgeIndex]),
                                                                edgeDirections[edgeIndex],
                                                                {curve.domainMidPoint()[0], 0.5}))
                     == IntersectionCurveAndLine::parallel;
    };

    clipper.SetZCallback(
        [&](const PointD& rectangleFirstPoint, const PointD& rectangleSecondPoint, const PointD& firstTrimmingCurvePoint, const PointD& secondTrimmingCurvePoint,
          const PointD& intersectionPoint) {

          // We are only interested in intersections with the edges
          // if (rectangleFirstPoint.z > 3) return;
          assert(rectangleFirstPoint.z < 4 && "Intersection should only occur with the rectangle edges");

          auto [isHost, idx] = isCornerVertex(intersectionPoint, eleRect);
          if (isHost) {
            result.addOriginalVertex(idx);
            return;
          }

          const auto edgeIdx = giveEdgeIdx(rectangleFirstPoint.z, rectangleSecondPoint.z);
          assert(firstTrimmingCurvePoint.z==secondTrimmingCurvePoint.z &&
            "The z values of the trimming curve points should have a difference of 1, since they should belong to the same curve");
          const auto curveZ  = secondTrimmingCurvePoint.z;

          // Now check that we don't have a parallel intersection
          if (checkParallel(patchTrimData.getCurve(curveZ), edgeIdx)) return;

          result.addNewVertex(edgeIdx, intersectionPoint, curveZ);
        });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    clipper.AddClip(trimmingCurves);
    clipper.AddSubject({eleRect});

    clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);

    for (const auto& pt : resultClosedPaths.front()) {
      if (auto [isHost, idx] = isCornerVertex(pt, eleRect); isHost) result.addOriginalVertex(idx);
    }

    // Now we have to check for points inside the elements
    for (auto i : Dune::range(patchTrimData.loops().size())) {
      for (const auto& pt : patchTrimData.getPointsInPatch(i)) {
        if (const PointD cp{pt.pt[0], pt.pt[1]}; PointInPolygon(cp, eleRect) == PointInPolygonResult::IsInside) {
          result.addInsideVertex(cp, pt.curveIdxI, pt.curveIdxJ, i);
        }
      }
    }

    // Add startpoints of trimming curves to the vertices if they are on one of the element edges
    for (auto cI = 0; cI < patchTrimData.loops().front().curves().size(); ++cI) {
      const auto& curve = patchTrimData.loops().front().curves()[cI];

      for (int i = 0; i < 2; ++i) {
        auto pt = curve.corner(i);
        PointD ptClipper{pt[0], pt[1]};

        auto [isHost, idx] = isCornerVertex(ptClipper, eleRect);
        if (isHost) {
          result.addOriginalVertex(idx);
          continue;
        }

        for (auto e = 0; e < edgeLookUp.size(); ++e) {
          if (const auto& edgeIdx = edgeLookUp[e];
              isPointOnLine(pt, eleGeo.corner(edgeIdx.front()), eleGeo.corner(edgeIdx.back()))
              && !checkParallel(curve, e)) {
            result.addNewVertex(e, ptClipper, patchTrimData.getZValue(cI,0));
            break;
          }
        }
      }
    }

    result.finish(patchTrimData);

    // While developing
    result.report();

    return std::make_tuple(ElementTrimFlag::trimmed, result);
  }
}  // namespace Dune::IGANEW::DefaultTrim::Util
