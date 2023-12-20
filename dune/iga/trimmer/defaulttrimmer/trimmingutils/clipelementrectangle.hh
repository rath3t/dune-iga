// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <matplot/matplot.h>
#include <variant>

#include <dune/common/float_cmp.hh>

namespace Dune::IGANEW::DefaultTrim::Impl {

  // enum class error_type { clipperNotSucessfull, malformedCurve };
  enum class TrimmingFlag { full, empty, trimmed };
  static int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }


  struct ClippingResult {
    struct HostVertex {
      Clipper2Lib::PointD pt{};
    };

    struct NewVertex {
      int onEdgeIdx{};
      std::size_t boundaryCurveIdx{};
      Clipper2Lib::PointD pt{};
    };
    using VertexVariant = std::variant<HostVertex, NewVertex>;

    auto isAlreadyThere(const auto& pt) {
      const auto it = std::ranges::find_if(vertices_, [&](const VertexVariant& vertexVariant) {
        auto isSame = false;
        std::visit([&](auto&& vertex) { isSame = FloatCmp::eq(vertex.pt.x, pt.x) and FloatCmp::eq(vertex.pt.y, pt.y); },
                   vertexVariant);
        return isSame;
      });
      return it != vertices_.end();
    }

    void addOriginalVertex(const Clipper2Lib::PointD& pt) {
      if (not isAlreadyThere(pt)) vertices_.emplace_back(HostVertex(pt));
    }

    void addNewVertex(const int edgeIdx, const std::size_t boundaryCurveIdx,
                      const Clipper2Lib::PointD& pt) {
      if (not isAlreadyThere(pt)) vertices_.emplace_back(NewVertex(edgeIdx, boundaryCurveIdx, pt));
    }

    void finish() {
      // Sort the points in counter clockwise manner such that the first point is in the lower left (courtesy to ChatGPT)
      auto getPoint = [](auto&& vertexVariant) {
        return vertexVariant.pt;
      };

      const auto minVertex = *std::ranges::min_element(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
        const auto a = std::visit(getPoint, aV);
        const auto b = std::visit(getPoint, bV);
        return (a.y < b.y) || (a.y == b.y && a.x < b.x);
      });

      auto polarAngle = [](const auto& p, const auto& reference) -> double {
        return atan2(p.y - reference.y, p.x - reference.x);
      };

      // Sort the points based on polar angle with respect to the reference point
      auto comparePoints  = [&](const VertexVariant& aV, const VertexVariant& bV, const VertexVariant& referenceV) -> bool {
        const auto a = std::visit(getPoint, aV);
        const auto b = std::visit(getPoint, bV);
        const auto reference = std::visit(getPoint, referenceV);
        const double angleA = polarAngle(a, reference);
        const double angleB = polarAngle(b, reference);

        if (angleA != angleB)
          return angleA < angleB;

        //If two points have the same angle, sort based on distance to the reference point
        return (a.x - reference.x) * (a.x - reference.x) + (a.y - reference.y) * (a.y - reference.y) <
               (b.x - reference.x) * (b.x - reference.x) + (b.y - reference.y) * (b.y - reference.y);

      };
      std::ranges::sort(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
        return comparePoints(aV, bV, minVertex);
      });

    }
    std::vector<VertexVariant> vertices_{};
  };

  inline int giveEdgeIdx(const std::size_t e1, const std::size_t e2) {
    assert(e1 < 4 and e2 < 4);
    if ((e1 == 0 and e2 == 1) or (e1 == 1 and e2 == 0)) return 0;
    if ((e1 == 1 and e2 == 2) or (e1 == 2 and e2 == 1)) return 1;
    if ((e1 == 2 and e2 == 3) or (e1 == 3 and e2 == 2)) return 2;
    if ((e1 == 3 and e2 == 0) or (e1 == 0 and e2 == 3)) return 3;
    assert(false);
  }

  inline bool goesIn(const std::size_t e1, const std::size_t e2) { return e1 > e2; }

  inline auto clipElementRectangle(Clipper2Lib::PathD& eleRect, const Clipper2Lib::PathsD& trimmingCurves)
      -> std::tuple<TrimmingFlag, ClippingResult> {
    using namespace Clipper2Lib;

    auto isFullElement = [&](const auto& clippedEdges) { return FloatCmp::eq(Area(clippedEdges), Area(eleRect)); };

    // First determine if element is trimmed
    const auto intersectResult = Intersect({eleRect}, trimmingCurves, FillRule::NonZero, 2);
    if (intersectResult.empty()) return std::make_tuple(TrimmingFlag::empty, ClippingResult{});
    if (isFullElement(intersectResult.front())) {
      matplot::rectangle(eleRect[0].x, eleRect[1].y, eleRect[1].x - eleRect[0].x, eleRect[3].y - eleRect[0].y);
      return std::make_tuple(TrimmingFlag::full, ClippingResult{});
    }

    // Element is trimmed
    ClipperD clipper(8);
    ClippingResult result{};

    clipper.SetZCallback(
        [&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top, const PointD& pt) {
          result.addNewVertex(giveEdgeIdx(e1bot.z, e1top.z), pt.z, pt);
        });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    for (const auto i : std::views::iota(0, 4)) {
      // std::cout << "Edge " << i << std::endl;

      clipper.Clear();
      clipper.AddClip(trimmingCurves);

      clipper.AddOpenSubject(PathsD{{eleRect[i], eleRect[nextEntityIdx(i, 1)]}});
      clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);

      if (not resultOpenPaths.empty())
        for (const auto& p : resultOpenPaths.front())
          result.addOriginalVertex(p);
    }

    auto isPointInTrimmingCurves = [&](const auto& pt1) -> bool {
      return std::ranges::any_of(trimmingCurves, [&](const auto& curve) {
        auto it =  std::ranges::find_if(
                 curve,
                 [&](const auto& pt2) { return FloatCmp::eq(pt1.x, pt2.x) and FloatCmp::eq(pt1.y, pt2.y); });
        return it != curve.end();
      });
    };

    // check if any element vertex is a starting point of one of the trimming curves as the algorihtm cannot find them
    for (const auto& cP : eleRect)
      if (isPointInTrimmingCurves(cP))
        result.addOriginalVertex(cP);

    result.finish();

    std::cout << "Vertices found\n";

    struct Visitor {
      static void plotEllipse(const PointD& pt) {
        constexpr auto w = 0.025;
        const auto c = matplot::ellipse(pt.x - (w / 2), pt.y - (w / 2), w, w);
        c->color("blue");
      }

      void operator()(const ClippingResult::HostVertex& v) const {
        plotEllipse(v.pt);
        std::cout << "Pt: " << v.pt << std::endl;
      }
      void operator()(const ClippingResult::NewVertex& v) const {
        plotEllipse(v.pt);
        std::cout << "Edge: " << v.onEdgeIdx << " Pt: " << v.pt
                  << " On BC: " << v.boundaryCurveIdx << std::endl;
      }
    };

    for (auto& vV : result.vertices_)
      std::visit(Visitor{}, vV);

    matplot::hold("on");
    matplot::axis(matplot::square);
    matplot::rectangle(eleRect[0].x, eleRect[1].y, eleRect[1].x - eleRect[0].x, eleRect[3].y - eleRect[0].y);

    return std::make_tuple(TrimmingFlag::trimmed, result);
  }
}  // namespace Dune::IGANEW::DefaultTrim
