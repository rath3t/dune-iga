// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>

#include <dune/common/float_cmp.hh>
#include <dune/iga/geometrykernel/findintersection.hh>

namespace Dune::IGANEW::DefaultTrim::Impl {
  // enum class error_type { clipperNotSucessfull, malformedCurve };
  static int nextEntityIdx(const int i, const int x) { return (i + x) % 4; }

  inline auto toCurveIdx = [](const size_t z) { return static_cast<size_t>(std::floor((z/100)-1)); };

  struct ClippingResult {


    explicit ClippingResult(const std::vector<Clipper2Lib::PointD>& oldV) : originalVertices_(oldV) {}

    struct HostVertex {
      size_t hostIdx{};
      Clipper2Lib::PointD pt{};
    };

    struct NewVertex {
      int onEdgeIdx{};
      Clipper2Lib::PointD pt{};
      size_t trimmingCurveZ{};
    };

    using VertexVariant = std::variant<HostVertex, NewVertex>;

    auto isAlreadyThere(const auto& pt) {
      const auto it = std::ranges::find_if(vertices_, [&](const VertexVariant& vertexVariant) {
        return std::visit([&](const auto& vertex) { return FloatCmp::eq(vertex.pt.x, pt.x) and FloatCmp::eq(vertex.pt.y, pt.y); },
                   vertexVariant);
      });
      return it != vertices_.end();
    }
    auto isAlreadyThereHV(const size_t hostIdx) {
      const auto it = std::ranges::find_if(vertices_, [&](const VertexVariant& vV) {
        if (std::holds_alternative<HostVertex>(vV))
          return std::get<HostVertex>(vV).hostIdx == hostIdx;
        return false;
      });
      return it != vertices_.end();
    }

    void addOriginalVertex(const Clipper2Lib::PointD& pt) {
      assert(pt.z < 4);
      if (not isAlreadyThereHV(pt.z)) vertices_.emplace_back(HostVertex(pt.z, originalVertices_[pt.z]));
    }

    void addNewVertex(const int edgeIdx, const Clipper2Lib::PointD& pt, const size_t trimmingCurveZ) {
      vertices_.emplace_back(NewVertex(edgeIdx, pt, trimmingCurveZ));
    }

    void finish() {
      // Sort the points in counter clockwise manner such that the first point is in the lower left
      auto getPoint = [](auto&& vertexVariant) { return vertexVariant.pt; };

      const auto minVertex
          = *std::ranges::min_element(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
              const auto a = std::visit(getPoint, aV);
              const auto b = std::visit(getPoint, bV);
              return (a.y < b.y) || (a.y == b.y && a.x < b.x);
            });

      auto polarAngle
          = [](const auto& p, const auto& reference) -> double { return atan2(p.y - reference.y, p.x - reference.x); };

      // Sort the points based on polar angle with respect to the reference point
      auto comparePoints
          = [&](const VertexVariant& aV, const VertexVariant& bV, const VertexVariant& referenceV) -> bool {
        const auto a         = std::visit(getPoint, aV);
        const auto b         = std::visit(getPoint, bV);
        const auto reference = std::visit(getPoint, referenceV);
        const double angleA  = polarAngle(a, reference);
        const double angleB  = polarAngle(b, reference);

        if (angleA != angleB) return angleA < angleB;

        // If two points have the same angle, sort based on distance to the reference point
        return (a.x - reference.x) * (a.x - reference.x) + (a.y - reference.y) * (a.y - reference.y)
               < (b.x - reference.x) * (b.x - reference.x) + (b.y - reference.y) * (b.y - reference.y);
      };

      // Sort counter-clockwise
      std::ranges::sort(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
        return comparePoints(aV, bV, minVertex);
      });

      // Remove duplicate entries for NewVertices, begin at the first host index
      const auto it = std::ranges::find_if(vertices_, [](const VertexVariant& vV) {
        return std::holds_alternative<HostVertex>(vV);
      });
      if (it == vertices_.end())
        DUNE_THROW(Dune::NotImplemented, "Algorithm needs at least one HostVertex to work");

      std::ranges::rotate(vertices_, it);
      for (const auto& vV : vertices_) {

      }


    }

    std::vector<VertexVariant> vertices_{};
  private:
    std::vector<Clipper2Lib::PointD> originalVertices_;
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

  auto clipElementRectangle(const auto& element, const auto& patchTrimData)
      -> std::tuple<ElementTrimFlag, ClippingResult> {
    using namespace Clipper2Lib;

    auto geo = element.geometry();
    std::array<FieldVector<double, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    PathD eleRect;
    for (const auto i : std::views::iota(0, 4))
      eleRect.push_back({corners[i][0], corners[i][1], i});

    PathsD trimmingCurves;
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
      trimmingCurves.push_back(tempPath);
    }

    // todo should this be hardcoded here?
    std::array<FieldVector<double, 2>, 4> dirs {FieldVector<double, 2>{1.0, 0.0}, FieldVector<double, 2>{0, 1},
      FieldVector<double, 2>{-1, 0}, FieldVector<double, 2>{0, -1}};


    // @todo check z value if same as eleRect -> full
    auto isFullElement = [&](const auto& clippedEdges) { return FloatCmp::eq(Area(clippedEdges), Area(eleRect)); };

    // First determine if element is trimmed
    // @todo try this with clipper for non-boundary elements with ClipperD
    const auto intersectResult = Intersect({eleRect}, trimmingCurves, FillRule::NonZero, 2);

    if (intersectResult.empty())
      return std::make_tuple(ElementTrimFlag::empty, ClippingResult{eleRect});
    if (isFullElement(intersectResult.front())) {
      matplot::rectangle(eleRect[0].x, eleRect[1].y, eleRect[1].x - eleRect[0].x, eleRect[3].y - eleRect[0].y);
      return std::make_tuple(ElementTrimFlag::full, ClippingResult{eleRect});
    }

    // Element is trimmed
    ClipperD clipper(8);
    ClippingResult result(eleRect);

    auto checkParallel = [&](const auto& curve, const int edgeIndex) -> bool {
      if (curve.degree().front() == 1) {
        FindIntersectionCurveAndLineResult rV;
        std::tie(rV, std::ignore, std::ignore) = findIntersectionLinearCurveAndLine(curve, corners[edgeIndex],
          dirs[edgeIndex], {curve.domainMidPoint()[0], 0.5});
        if (rV == FindIntersectionCurveAndLineResult::linesParallel)
          return true;
      }
      return false;
    };

    clipper.SetZCallback([&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top,
                             const PointD& pt) {
      if (not checkParallel(patchTrimData.loops().front().curves()[toCurveIdx(e2bot.z)], giveEdgeIdx(e1bot.z, e1top.z)))
        result.addNewVertex(giveEdgeIdx(e1bot.z, e1top.z), pt, e2bot.z);
    });

    PathsD resultClosedPaths{};
    PathsD resultOpenPaths{};

    for (const auto i : std::views::iota(0, 4)) {
      // std::cout << "Edge " << i << std::endl;

      clipper.Clear();
      clipper.AddClip(trimmingCurves);

      clipper.AddOpenSubject(PathsD{{eleRect[i], eleRect[nextEntityIdx(i, 1)]}});
      clipper.Execute(ClipType::Intersection, FillRule::NonZero, resultClosedPaths, resultOpenPaths);
      // We receive here all intersection points but not the sampled curve points
      // resultClosedPaths should be empty
      assert(resultClosedPaths.empty());
      assert(resultOpenPaths.size()<=1);

      if (not resultOpenPaths.empty())
        for (const auto& p : resultOpenPaths.front()) {
          if (p.z < 4)
            result.addOriginalVertex(p);
          else {
            // Check if this is neeeded
            std::cout << "maybe we need thhis point\n";
          }
        }
    }

    auto isPointInTrimmingCurves = [&](const auto& pt1) -> bool {
      return std::ranges::any_of(trimmingCurves, [&](const auto& curve) {
        auto it = std::ranges::find_if(
            curve, [&](const auto& pt2) { return FloatCmp::eq(pt1.x, pt2.x) and FloatCmp::eq(pt1.y, pt2.y); });
        return it != curve.end();
      });
    };

    // check if any element vertex is a starting point of one of the trimming curves as the algorihtm cannot find them
    // @todo only for boundary elements
    for (const auto& cP : eleRect)
      if (isPointInTrimmingCurves(cP)) result.addOriginalVertex(cP);

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
}  // namespace Dune::IGANEW::DefaultTrim::Impl
