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

    auto eleGeo                          = element.geometry();
    static constexpr int numberOfCorners = 4;

    PathD eleRect;
    for (const auto vertexIndex : Dune::range(numberOfCorners)) {
      auto corner = eleGeo.corner(vertexIndexMapping[vertexIndex]);
      eleRect.emplace_back(corner[0], corner[1], vertexIndex);
    }

    const auto& trimmingCurves = patchTrimData.clipperLoops();

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

    clipper.SetZCallback([&](const PointD& rectangleFirstPoint, const PointD& rectangleSecondPoint,
                             const PointD& firstTrimmingCurvePoint, const PointD& secondTrimmingCurvePoint,
                             const PointD& intersectionPoint) {
      assert(rectangleFirstPoint.z < 4 && "Intersection should only occur with the rectangle edges");
      assert(rectangleSecondPoint.z < 4 && "Intersection should only occur with the rectangle edges");

      const auto [isHost, idx] = isCornerVertex(intersectionPoint, eleRect);
      if (isHost) {
        result.addOriginalVertex(idx);
        return;
      }

      const auto edgeIdx = giveEdgeIdx(rectangleFirstPoint.z, rectangleSecondPoint.z);
      assert(patchTrimData.getIndices(firstTrimmingCurvePoint.z).first
                 == patchTrimData.getIndices(secondTrimmingCurvePoint.z).first
             && "The points of the trimming curves should be on the same loop");
      assert(firstTrimmingCurvePoint.z==secondTrimmingCurvePoint.z or secondTrimmingCurvePoint.z== intersectionPoint.z or firstTrimmingCurvePoint.z== intersectionPoint.z &&
  "The indices of the trimming curves should be the same or the intersection point should be on one of the two second curves");
      const auto vertexZValue = secondTrimmingCurvePoint.z;

      // Now check that we don't have a parallel intersection
      if (checkParallel(patchTrimData.getCurve(vertexZValue), edgeIdx)) {
        return;
      }

      result.addNewVertex(edgeIdx, intersectionPoint, vertexZValue);
    });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    clipper.AddClip(trimmingCurves);
    clipper.AddSubject({eleRect});

    clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);
    assert(resultOpenPaths.empty() && "There should be no open paths in the clipping result");

    auto isFullElement
        = [&](const auto& clippedEdges) {
          return clippedEdges.size()==4 and (clippedEdges[0].z < 4 and clippedEdges[1].z <4 and clippedEdges[2].z <4 and clippedEdges[3].z <4);
        };

    if (resultClosedPaths.empty()) return std::make_tuple(ElementTrimFlag::empty, ClippingResult{eleRect});
    if (isFullElement(resultClosedPaths.front())) {
      assert(resultClosedPaths.size() == 1 && "If the element is full, there should be only one path");
      return std::make_tuple(ElementTrimFlag::full, ClippingResult{eleRect});
    }

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
