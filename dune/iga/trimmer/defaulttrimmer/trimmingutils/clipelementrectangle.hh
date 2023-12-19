// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>
#include <dune/common/float_cmp.hh>

namespace Dune::IGANEW::DefaultTrim {

  // enum class error_type { clipperNotSucessfull, malformedCurve };
  enum class TrimmingFlag {
    full, empty, trimmed
  };
  static int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }

  struct ClippingResult {
    std::array<bool, 4> edgesVisited{false, false, false, false};

    struct HostVertex {
      std::size_t originalIdx{};
      Clipper2Lib::PointD pt{};
    };

    struct NewVertex {
      int onEdgeIdx{};
      bool goesIn{};
      std::size_t boundaryCurveIdx{};
      Clipper2Lib::PointD pt{};
    };

    using VertexVariant = std::variant<HostVertex, NewVertex>;
    std::vector<VertexVariant> vertices{};

    void addOriginalVertex(const size_t originalIdx, const Clipper2Lib::PointD& pt) {
      vertices.emplace_back(HostVertex(originalIdx, pt));
    }

    void addNewVertex(const int edgeIdx, const bool goesIn, const std::size_t boundaryCurveIdx,
                      const Clipper2Lib::PointD& pt) {
      edgesVisited[edgeIdx] = true;
      vertices.emplace_back(NewVertex(edgeIdx, goesIn, boundaryCurveIdx, pt));
    }

    void finish(const Clipper2Lib::PathD& eleRect) {
      // For now
      if (vertices.empty()) {
        for (const auto i : std::views::iota(0u, 4u))
          addOriginalVertex(i, eleRect[i]);
        return;
      }

      // For now only 2 Intersections, on different edges
      assert(std::get<1>(vertices[0]).goesIn and vertices.size() == 2 and not std::get<1>(vertices[1]).goesIn);

      const auto idx1 = std::get<1>(vertices[0]).onEdgeIdx;
      auto idx2       = std::get<1>(vertices[1]).onEdgeIdx;

      while (idx2 != idx1) {
        idx2 = nextEntityIdx(idx2, 1);
        addOriginalVertex(idx2, eleRect[idx2]);
      }

      // Now sort counter-clockwise
    }
  };

  inline int giveEdgeIdx(const std::size_t e1, const std::size_t e2) {
    assert(e1 < 4 and e2 < 4);
    if ((e1 == 0 and e2 == 1) or (e1 == 1 and e2 == 0)) return 0;
    if ((e1 == 1 and e2 == 2) or (e1 == 2 and e2 == 1)) return 1;
    if ((e1 == 2 and e2 == 3) or (e1 == 3 and e2 == 2)) return 2;
    if ((e1 == 3 and e2 == 0) or (e1 == 0 and e2 == 3)) return 3;
    assert(false);
  }

  inline bool goesIn(const std::size_t e1, const std::size_t e2) {
    if (giveEdgeIdx(e1, e2) < 3) return e1 > e2;
    return e1 < e2;
  }

  inline auto clipElementRectangle(Clipper2Lib::PathD& eleRect, Clipper2Lib::PathsD& trimmingCurves) -> std::tuple<TrimmingFlag, ClippingResult> {
    using namespace Clipper2Lib;

    auto isFullElement = [&](const auto& clippedEdges) {
      return FloatCmp::eq(Area(clippedEdges), Area(eleRect));
    };


    // First determine if element is trimmed
    const auto intersectResult = Intersect({eleRect}, trimmingCurves, FillRule::NonZero, 8);
    if (intersectResult.empty())
      return std::make_tuple(TrimmingFlag::empty, ClippingResult{});
    if (isFullElement(intersectResult.front()))
      return std::make_tuple(TrimmingFlag::full, ClippingResult{});


    ClipperD clipper(8);

    auto isOnHostVertex = [&](const auto& point)-> bool {
      auto it =  std::ranges::find_if(eleRect, [&](const auto& elementPoint) {
        return FloatCmp::eq(elementPoint.x, point.x) and FloatCmp::eq(elementPoint.y, point.y);
      });
      return it != eleRect.end();
    };

    ClippingResult result{};
    clipper.SetZCallback(
        [&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top, const PointD& pt) {
          // if (e1bot.z > 3 or e1top.z > 3)
          //   return;
          std::cout << "New Intersection x: " << pt.x << " y: " << pt.y << std::endl;
          std::cout << e1bot.z << " " << e1top.z << std::endl;
          std::cout << e2bot.z << " " << e2top.z << std::endl;

          result.addNewVertex(giveEdgeIdx(e1bot.z, e1top.z), goesIn(e1bot.z, e1top.z), -1, pt);
        });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    for (const auto i : std::views::iota(0, 4)) {
      std::cout << "Edge " << i << std::endl;

      clipper.Clear();
      clipper.AddClip(trimmingCurves);

      clipper.AddOpenSubject(PathsD{{eleRect[i], eleRect[nextEntityIdx(i, 1)]}});
      clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);
      std::cout << resultClosedPaths;
    }
    // result.finish(eleRect);

    std::cout << resultClosedPaths;

    std::cout << "Vertices found\n";

    struct Visitor {
      void operator()(const ClippingResult::HostVertex& v) const {
        std::cout << "Original Idx: " << v.originalIdx << std::endl;
      }
      void operator()(const ClippingResult::NewVertex& v) const {
        std::cout << "Edge: " << v.onEdgeIdx << ", goes in " << v.goesIn << " Pt: " << v.pt << std::endl;
      }
    };

    for (auto& vV : result.vertices)
      std::visit(Visitor{}, vV);

    return std::make_tuple(TrimmingFlag::trimmed, result);
  }
}  // namespace Dune::IGANEW::DefaultTrim
