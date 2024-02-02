// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

namespace Dune::IGANEW::DefaultTrim::Util {

  inline int giveEdgeIdx(const std::size_t e1, const std::size_t e2) {
    if ((e1 == 0 and e2 == 1) or (e1 == 1 and e2 == 0)) return 0;
    if ((e1 == 1 and e2 == 2) or (e1 == 2 and e2 == 1)) return 1;
    if ((e1 == 2 and e2 == 3) or (e1 == 3 and e2 == 2)) return 2;
    if ((e1 == 3 and e2 == 0) or (e1 == 0 and e2 == 3)) return 3;
    DUNE_THROW(Dune::GridError, "No corresponding edge Index for Input");
  }

  constexpr std::array edgeDirections{FieldVector<double, 2>{1.0, 0.0}, FieldVector<double, 2>{0, 1},
                                      FieldVector<double, 2>{-1, 0}, FieldVector<double, 2>{0, -1}};
  constexpr std::array<std::array<int, 2>, 4> edgeLookUp{std::array{0, 1}, {1, 3}, {3, 2}, {2, 0}};
  constexpr std::array vertexIndexMapping = {0u, 1u, 3u, 2u};

  auto isCornerVertex(const auto& pt, const auto& eleRect) -> std::pair<bool, ptrdiff_t> {
    auto it = std::ranges::find_if(eleRect, [&](const auto& vertex) {
      return FloatCmp::eq(vertex.x, pt.x, 1e-8) and FloatCmp::eq(vertex.y, pt.y, 1e-8);
    });
    return std::make_pair(it != eleRect.end(), std::ranges::distance(eleRect.begin(), it));
  }

  auto isPointOnLine(const auto& pt, const auto& p1, const auto& p2) -> bool {
    // Check if the points are collinear using the slope formula
    auto x  = pt[0];
    auto y  = pt[1];
    auto x1 = p1[0];
    auto y1 = p1[1];
    auto x2 = p2[0];
    auto y2 = p2[1];

    // Check if the line is horizontal
    if (std::abs(y1 - y2) < 1e-6) {
      return std::abs(y - y1) < 1e-6 && ((x >= x1 && x <= x2) || (x >= x2 && x <= x1));
    }
    // Check if the line is vertical
    if (std::abs(x1 - x2) < 1e-6) {
      return std::abs(x - x1) < 1e-6 && ((y >= y1 && y <= y2) || (y >= y2 && y <= y1));
    }

    // Check if the given point lies on the line
    const double slope = (y2 - y1) / (x2 - x1);
    return std::abs((y - y1) - slope * (x - x1)) < 1e-6;
  }

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
    struct InsideVertex {
      Clipper2Lib::PointD pt{};
      size_t curveIdxI{};
      size_t curveIdxJ{};
      size_t loopIdx{};
    };

    using VertexVariant = std::variant<HostVertex, NewVertex, InsideVertex>;

    void addOriginalVertex(const Clipper2Lib::PointD& pt) {
      assert(pt.z < 4);
      if (FloatCmp::eq(pt.x, originalVertices_[pt.z].x) and FloatCmp::eq(pt.y, originalVertices_[pt.z].y)
          and not isAlreadyThere(pt.z))
        vertices_.emplace_back(HostVertex(pt.z, originalVertices_[pt.z]));
    }
    void addOriginalVertex(const size_t hostIdx) {
      if (not isAlreadyThere(hostIdx)) vertices_.emplace_back(HostVertex(hostIdx, originalVertices_[hostIdx]));
    }

    void addNewVertex(const int edgeIdx, const Clipper2Lib::PointD& pt, const size_t trimmingCurveZ) {
      if (not isAlreadyThere(pt)) vertices_.emplace_back(NewVertex(edgeIdx, pt, trimmingCurveZ));
    }
    void addInsideVertex(const Clipper2Lib::PointD& pt, const size_t curveIndexI, const size_t curveIndexJ,
                         size_t loopIdx) {
      vertices_.emplace_back(InsideVertex(pt, curveIndexI, curveIndexJ, loopIdx));
    }

    void finish(const auto& patchTrimData) {
      // Sort the points in counter clockwise manner such that the first point is in the lower left
      const auto minVertex = lowerLeftVertex();
      std::ranges::sort(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
        return comparePoints(aV, bV, minVertex);
      });

      // Resort that we always start on a hostVertex (if there is one)
      if (const auto it = std::ranges::find_if(
              vertices_, [](const VertexVariant& vV) { return std::holds_alternative<HostVertex>(vV); });
          it != vertices_.end())
        std::ranges::rotate(vertices_, it);

      // Now resort InsideVertices (sometimes they end up on the wrong position)
      for (const auto i : std::views::iota(0ul, vertices_.size())) {
        if (const auto insideVertex = std::get_if<InsideVertex>(&vertices_[i])) {
          auto it = vertices_.begin() + static_cast<long>(i);
          const size_t indexBeforeIt = it == vertices_.begin() ? vertices_.size() - 1 : std::ranges::distance(vertices_.begin(), it) - 1;
          const size_t indexAfterIt = it == std::prev(vertices_.end()) ? 0 : std::ranges::distance(vertices_.begin(), it) + 1;

          // Check vertex before if they line up correctly
          const auto vertexBefore = std::get_if<NewVertex>(&vertices_[indexBeforeIt]);
          const auto vertexAfter = std::get_if<NewVertex>(&vertices_[indexAfterIt]);
          if (vertexBefore && vertexAfter && patchTrimData.getIndices(vertexBefore->trimmingCurveZ)
                     == std::make_pair(insideVertex->loopIdx, insideVertex->curveIdxI) ) {

            continue;
          }
          // If we end up here sort the vertex to a correct place
          auto correctVertexIt = findCorrespondingVertexIt(*insideVertex, patchTrimData);
          assert(correctVertexIt != vertices_.end() && "No corresponding newVertex found");

          // Add the insideVertex at correctVertexIt+1
          const auto correctItForInsideVertex = correctVertexIt + 1;
          if (std::distance(vertices_.begin(), correctItForInsideVertex) < i) {
            vertices_.insert(correctItForInsideVertex, *insideVertex);
            vertices_.erase(vertices_.begin() + static_cast<long>(i+1));
          } else {
            std::ranges::rotate(it, it + 1, correctItForInsideVertex);
          }
        }
      }
    }

    void report() const {
      std::cout << "Vertices found\n";
      struct Visitor {
        void operator()(const HostVertex& v) const {
          std::cout << "Pt: " << v.pt << " Host Idx: " << v.hostIdx << std::endl;
        }
        void operator()(const NewVertex& v) const {
          std::cout << "Edge: " << v.onEdgeIdx << " Pt: " << v.pt << " On TC: " << v.trimmingCurveZ << std::endl;
        }
        void operator()(const InsideVertex& v) const {
          std::cout << " Pt: " << v.pt << " On TC: " << v.curveIdxI << " and " << v.curveIdxJ << std::endl;
        }
      };
      for (auto& vV : vertices_)
        std::visit(Visitor{}, vV);
    }

    std::vector<VertexVariant> vertices_{};

  private:
    std::vector<Clipper2Lib::PointD> originalVertices_;

    auto isAlreadyThere(const Clipper2Lib::PointD& pt) -> bool {
      const auto it = std::ranges::find_if(vertices_, [&](const VertexVariant& vertexVariant) {
        return std::visit(
            [&](const auto& vertex) {
              return FloatCmp::eq(vertex.pt.x, pt.x, 1e-8) and FloatCmp::eq(vertex.pt.y, pt.y, 1e-8);
            },
            vertexVariant);
      });
      return it != vertices_.end();
    }
    auto isAlreadyThere(const size_t hostIdx) -> bool {
      const auto it = std::ranges::find_if(vertices_, [&](const VertexVariant& vV) {
        if (std::holds_alternative<HostVertex>(vV)) return std::get<HostVertex>(vV).hostIdx == hostIdx;
        return false;
      });
      return it != vertices_.end();
    }

    inline static auto getPoint = [](auto&& vertexVariant) { return vertexVariant.pt; };

    // Sort the points based on polar angle with respect to the reference point
    inline static auto comparePoints
        = [](const VertexVariant& aV, const VertexVariant& bV, const VertexVariant& referenceV) -> bool {
      auto polarAngle
          = [](const auto& p, const auto& reference) -> double { return atan2(p.y - reference.y, p.x - reference.x); };

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

    auto findCorrespondingVertexIt(const InsideVertex& insideVertex, const auto& patchTrimData) {
      return std::ranges::find_if(vertices_, [&](const VertexVariant& vV) {
        return std::holds_alternative<NewVertex>(vV)
               && patchTrimData.getIndices(std::get<NewVertex>(vV).trimmingCurveZ)
                      == std::make_pair(insideVertex.loopIdx, insideVertex.curveIdxI);
      });
    }

    [[nodiscard]] const VertexVariant& lowerLeftVertex() const {
      return *std::ranges::min_element(vertices_, [&](const VertexVariant& aV, const VertexVariant& bV) {
        const auto a = std::visit(getPoint, aV);
        const auto b = std::visit(getPoint, bV);
        return (a.y < b.y) || (a.y == b.y && a.x < b.x);
      });
    }
  };

}  // namespace Dune::IGANEW::DefaultTrim::Util
