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

    auto eleGeo                      = element.geometry();

    PathD eleRect;
    for (const auto i : std::views::iota(0u, 4u)) {
      auto corner = eleGeo.corner(vIdxMapping[i]);
      eleRect.emplace_back(corner[0], corner[1], i);
    }

    const auto& trimmingCurves = patchTrimData.clipperLoops();

    // @todo check z value if same as eleRect -> full
    auto isFullElement
        = [&](const auto& clippedEdges) { return FloatCmp::eq(Area(clippedEdges), Area(eleRect), 1e-8); };

    // First determine if element is trimmed
    const auto intersectResult = Intersect({eleRect}, trimmingCurves, FillRule::NonZero, 8);

    if (intersectResult.empty()) return std::make_tuple(ElementTrimFlag::empty, ClippingResult{eleRect});
    if (isFullElement(intersectResult.front())) {
      return std::make_tuple(ElementTrimFlag::full, ClippingResult{eleRect});
    }

    // Element is trimmed
    ClipperD clipper(8);
    ClippingResult result(eleRect);

    auto checkParallel = [&](const auto& curve, const int edgeIndex) -> bool {
      return curve.degree().front() == 1
             and std::get<0>(findIntersectionLinearCurveAndLine(curve, eleGeo.corner(vIdxMapping[edgeIndex]),
                                                                edgeDirections[edgeIndex],
                                                                {curve.domainMidPoint()[0], 0.5}))
                     == FindIntersectionCurveAndLineResult::linesParallel;
    };

    clipper.SetZCallback(
        [&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top, const PointD& pt) {
          std::cout << "Z-Callback\n";
          std::cout << e1bot << " " << e1top << std::endl;
          std::cout << e2bot << " " << e2top << std::endl;
          std::cout << pt << std::endl;

          // We are only interested in intersections with the edges
          if (e1bot.z > 3) return;

          if (auto [isHost, idx] = isHostVertex(pt, eleRect); isHost) {
            result.addOriginalVertex(idx);
            return;
          }

          const auto edgeIdx  = giveEdgeIdx(e1bot.z, e1top.z);
          const auto curveZ   = e2bot.z;

          // Now check that we don't have a parallel intersection
          if (checkParallel(patchTrimData.getCurve(curveZ), edgeIdx))
            return;

          result.addNewVertex(edgeIdx, pt, curveZ);
        });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    clipper.AddClip(trimmingCurves);
    clipper.AddSubject({eleRect});

    clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);

    for (const auto& pt : resultClosedPaths.front()) {
      if (auto [isHost, idx] = isHostVertex(pt, eleRect); isHost) result.addOriginalVertex(idx);
    }

    // Add startpoints of trimming curves to the vertices if they are on one of the element edges
    for (auto cI = 0; cI < patchTrimData.loops().front().curves().size(); ++cI) {
      const auto& curve = patchTrimData.loops().front().curves()[cI];

      for (int i = 0; i < 2; ++i) {
        auto pt = curve.corner(i);
        PointD ptClipper{pt[0], pt[1]};

        if (auto [isHost, idx] = isHostVertex(ptClipper, eleRect); isHost)
          result.addOriginalVertex(idx);

        for (auto e = 0; e < edgeLookUp.size(); ++e) {
          if (const auto& edgeIdx = edgeLookUp[e];
            isPointOnLine(pt, eleGeo.corner(edgeIdx.front()), eleGeo.corner(edgeIdx.back()))
              && !checkParallel(curve, e)) {
            // \todo this needs to be done more generically
            const size_t newIdx = i == 0 ? (cI + 1) * patchTrimData.getSplitter() : (cI + 1) * patchTrimData.getSplitter() + patchTrimData.getSplitter() - 1;
            result.addNewVertex(e, ptClipper, newIdx);
            break;
          }
        }
      }
    }

    result.finish();

    // While developing
    result.report();

    return std::make_tuple(ElementTrimFlag::trimmed, result);
  }
}  // namespace Dune::IGANEW::DefaultTrim::Util
