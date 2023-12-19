// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>

namespace Dune::IGANEW::DefaultTrim {

  enum class error_type { clipperNotSucessfull, malformedCurve };

  struct ClippingResult {
    std::size_t nVertices{};
    std::size_t newVertices = 4;
    std::array<bool, 4> edgesVisited{false, false, false, false};

    struct HostVertex {
      std::size_t idx{};
      std::size_t originalIdx{};
      Clipper2Lib::PointD pt{};
    };

    struct NewVertex {
      std::size_t idx{};
      int onEdgeIdx{};
      bool goesIn{};
      std::size_t boundaryCurveIdx{};
      Clipper2Lib::PointD pt{};
    };

    static int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }

    using VertexVariant = std::variant<HostVertex, NewVertex>;
    std::vector<VertexVariant> vertices{};

    void addOriginalVertex(const size_t originalIdx, const Clipper2Lib::PointD& pt) {
      vertices.emplace_back(HostVertex(nVertices++, originalIdx, pt));
    }

    void addNewVertex(const int edgeIdx, const bool goesIn, const std::size_t boundaryCurveIdx,
                      const Clipper2Lib::PointD& pt) {
      edgesVisited[edgeIdx] = true;
      vertices.emplace_back(NewVertex(nVertices++, edgeIdx, goesIn, boundaryCurveIdx, pt));
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
    __builtin_unreachable();
  }

  inline bool goesIn(const std::size_t e1, const std::size_t e2) {
    if (giveEdgeIdx(e1, e2) < 3) return e1 > e2;
    return e1 < e2;
  }

  inline auto clipElementRectangle(Clipper2Lib::PathD& eleRect, Clipper2Lib::PathsD& trimmingCurves) {
    using namespace Clipper2Lib;

    ClipperD clipper(8);
    // Counters
    std::size_t counter     = 0;
    std::size_t loopCounter = 0;
    std::vector<size_t> counterBreaks{};

    // Give points in eleRect the z-Val from 0 to 4 counter-clockwise
    std::ranges::for_each(eleRect, [&](auto& point) { point.z = counter++; });
    counterBreaks.push_back(counter);

    for (auto& trimmingCurve : trimmingCurves) {
      std::ranges::for_each(trimmingCurve, [&](auto& point) { point.z = counter++; });
      counterBreaks.push_back(counter);
      loopCounter++;
    }

    ClippingResult result{};

    clipper.SetZCallback(
        [&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top, const PointD& pt) {
          std::cout << "New Intersection x: " << pt.x << " y: " << pt.y << std::endl;
          std::cout << e1bot.z << " " << e1top.z << std::endl;
          std::cout << e2bot.z << " " << e2top.z << std::endl;

          // result.addNewVertex(giveEdgeIdx(e1bot.z, e1top.z), goesIn(e1bot.z, e1top.z), e2bot.z, pt);
        });

    clipper.AddClip(trimmingCurves);
    clipper.AddSubject({eleRect});
    //clipper.AddOpenSubject(PathsD{{eleRect[3], eleRect[0]}});

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    clipper.Execute(ClipType::Intersection, FillRule::NonZero, trimmingCurves);
    // result.finish(eleRect);

    std::cout << "Vertices found\n";

    struct Visitor {
      void operator()(const ClippingResult::HostVertex& v) const {
        std::cout << v.idx << " Original Idx: " << v.originalIdx << std::endl;
      }

      void operator()(const ClippingResult::NewVertex& v) const {
        std::cout << v.idx << " Edge: " << v.onEdgeIdx << ", goes in " << v.goesIn << std::endl;
      }
    };

    for (auto& vV : result.vertices)
      std::visit(Visitor{}, vV);

    return result;
  }
}  // namespace Dune::IGANEW::DefaultTrim
