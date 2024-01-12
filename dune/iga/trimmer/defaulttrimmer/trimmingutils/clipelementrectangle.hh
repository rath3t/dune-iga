// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>

#include "cliputils.hh"

#include <dune/common/float_cmp.hh>

#include <dune/iga/geometrykernel/findintersection.hh>

namespace Dune::IGANEW::DefaultTrim::Util {

  // static int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }
  constexpr std::array<std::array<int, 2>, 4> edgeLookUp{std::array{0, 1}, {1, 3}, {3, 2}, {2, 0}};

  inline auto toCurveIdx = [](const size_t z) { return static_cast<size_t>(std::floor((z / 100) - 1)); };

  auto clipElementRectangle(const auto& element, const auto& patchTrimData)
      -> std::tuple<ElementTrimFlag, ClippingResult> {
    using namespace Clipper2Lib;

    auto eleGeo                      = element.geometry();
    constexpr std::array vIdxMapping = {0u, 1u, 3u, 2u};

    PathD eleRect;
    for (const auto i : std::views::iota(0u, 4u)) {
      auto corner = eleGeo.corner(vIdxMapping[i]);
      eleRect.emplace_back(corner[0], corner[1], i);
    }

    PathsD trimmingCurves;
    PathD tempPath;
    constexpr int N = 80;
    // @todo store param value of sampled points on trimmingCurve
    for (auto loop : patchTrimData.loops()) {
      tempPath.clear();
      for (int i = 100; auto& curve : loop.curves()) {
        for (auto v : Utilities::linspace(curve.domain()[0], N)) {
          auto fV = curve.global({v});
          tempPath.emplace_back(fV[0], fV[1], i++);
        }
        i += 100 - N;
      }
      trimmingCurves.push_back(tempPath);
    }

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
          const auto curveIdx = toCurveIdx(curveZ);

          // Now check that we don't have a parallel intersection
          if (checkParallel(patchTrimData.loops().front().curves()[curveIdx], edgeIdx)) return;

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

        if (auto [isHost, idx] = isHostVertex(ptClipper, eleRect); isHost) result.addOriginalVertex(idx);

        for (auto e = 0; e < edgeLookUp.size(); ++e) {
          const auto& edgeIdx = edgeLookUp[e];
          if (isPointOnLine(pt, eleGeo.corner(edgeIdx.front()), eleGeo.corner(edgeIdx.back()))
              && !checkParallel(curve, e)) {
            result.addNewVertex(e, ptClipper, (cI + 1) * 100);
            break;
          }
        }
      }
    }

    result.finish();

    std::cout << "Vertices found\n";
    struct Visitor {
      void operator()(const ClippingResult::HostVertex& v) const {
        std::cout << "Pt: " << v.pt << " Host Idx: " << v.hostIdx << std::endl;
      }
      void operator()(const ClippingResult::NewVertex& v) const {
        std::cout << "Edge: " << v.onEdgeIdx << " Pt: " << v.pt << " On TC: " << v.trimmingCurveZ << std::endl;
      }
    };
    for (auto& vV : result.vertices_)
      std::visit(Visitor{}, vV);

    return std::make_tuple(ElementTrimFlag::trimmed, result);
  }
}  // namespace Dune::IGANEW::DefaultTrim::Util
