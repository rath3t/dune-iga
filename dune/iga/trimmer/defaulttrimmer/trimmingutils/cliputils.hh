// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <array>
#include <clipper2/clipper.core.h>
#include <ranges>
#include <variant>

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/exceptions.hh>

namespace Dune::IGANEW::DefaultTrim::Util {

inline int giveEdgeIdx(const std::size_t e1, const std::size_t e2) {
  if ((e1 == 0 and e2 == 1) or (e1 == 1 and e2 == 0))
    return 0;
  if ((e1 == 1 and e2 == 2) or (e1 == 2 and e2 == 1))
    return 1;
  if ((e1 == 2 and e2 == 3) or (e1 == 3 and e2 == 2))
    return 2;
  if ((e1 == 3 and e2 == 0) or (e1 == 0 and e2 == 3))
    return 3;
  DUNE_THROW(Dune::GridError, "No corresponding edge Index for Input");
}

constexpr std::array edgeDirections{
    FieldVector<double, 2>{1.0, 0.0},
     FieldVector<double, 2>{  0,   1},
     FieldVector<double, 2>{ -1,   0},
    FieldVector<double, 2>{  0,  -1}
};
constexpr std::array<std::array<int, 2>, 4> edgeLookUp{
    std::array{0, 1},
    {1, 3},
    {3, 2},
    {2, 0}
};
// @todo use transformations
constexpr std::array vertexIndexMapping = {0u, 1u, 3u, 2u};
constexpr std::array edgeIndexMapping   = {2u, 1u, 3u, 0u};
auto isCornerVertex(const auto& pt, const auto& eleRect) -> std::pair<bool, ptrdiff_t> {
  auto it = std::ranges::find_if(eleRect, [&](const auto& vertex) {
    return FloatCmp::eq(vertex.x, pt.x, 1e-8) and FloatCmp::eq(vertex.y, pt.y, 1e-8);
  });
  return std::make_pair(it != eleRect.end(), std::ranges::distance(eleRect.begin(), it));
}

struct ClippingResult
{
  explicit ClippingResult(const std::vector<Clipper2Lib::PointD>& oldV)
      : originalVertices_(oldV) {}

  struct Vertex
  {
    struct HostVertexImpl
    {
      size_t hostIdx{};
    };
    struct NewVertexImpl
    {
      int onEdgeIdx{};
      u_int64_t trimmingCurveZ{};
    };
    struct InsideVertexImpl
    {
      size_t curveIdxI{};
      size_t curveIdxJ{};
      size_t loopIdx{};
    };
    using VertexVariant = std::variant<HostVertexImpl, NewVertexImpl, InsideVertexImpl>;
    Clipper2Lib::PointD pt{};
    VertexVariant vertexData;
    // std::optional<size_t> hostIdx;
    // std::optional<size_t> onEdgeIdx;

    bool isHost() const {
      return std::holds_alternative<HostVertexImpl>(vertexData);
    }
    bool isInside() const {
      return std::holds_alternative<InsideVertexImpl>(vertexData);
    }
    bool isNew() const {
      return std::holds_alternative<NewVertexImpl>(vertexData);
    }
    auto loop() const {
      assert(isInside());
      return std::get<InsideVertexImpl>(vertexData).loopIdx;
    }
    auto formerCurve() const {
      assert(isInside());
      return std::get<InsideVertexImpl>(vertexData).curveIdxI;
    }
    auto subsequentCurve() const {
      assert(isInside());
      return std::get<InsideVertexImpl>(vertexData).curveIdxJ;
    }
    auto zValue() const {
      assert(isNew());
      return std::get<NewVertexImpl>(vertexData).trimmingCurveZ;
    }
    auto hostId() const {
      assert(isHost());
      return std::get<HostVertexImpl>(vertexData).hostIdx;
    }
    auto edgeId() const {
      assert(isNew());
      return std::get<NewVertexImpl>(vertexData).onEdgeIdx;
    }
    static Vertex HostVertex(const Clipper2Lib::PointD& pt, size_t hostIdx) {
      return Vertex{.pt = pt, .vertexData = HostVertexImpl(hostIdx)};
    }
    static Vertex NewVertex(int onEdgeIdx, const Clipper2Lib::PointD& pt, u_int64_t trimmingCurveZ) {
      return Vertex{.pt = pt, .vertexData = NewVertexImpl(onEdgeIdx, trimmingCurveZ)};
    }
    static Vertex InsideVertex(const Clipper2Lib::PointD& pt, size_t curveIdxI, size_t curveIdxJ, size_t loopIdx) {
      return Vertex{.pt = pt, .vertexData = InsideVertexImpl(curveIdxI, curveIdxJ, loopIdx)};
    }
    friend std::ostream& operator<<(std::ostream& os, const Vertex& vertex) {
      os << "Vertex: " << vertex.pt << " ";
      if (vertex.isHost())
        os << "HostVertex: " << std::get<HostVertexImpl>(vertex.vertexData).hostIdx;
      if (vertex.isNew())
        os << "NewVertex: " << std::get<NewVertexImpl>(vertex.vertexData).onEdgeIdx << " "
           << std::get<NewVertexImpl>(vertex.vertexData).trimmingCurveZ;
      if (vertex.isInside())
        os << "InsideVertex: " << std::get<InsideVertexImpl>(vertex.vertexData).curveIdxI << " "
           << std::get<InsideVertexImpl>(vertex.vertexData).curveIdxJ << " "
           << std::get<InsideVertexImpl>(vertex.vertexData).loopIdx;
      return os;
    }
  };

  void addOriginalVertex(const Clipper2Lib::PointD& pt) {
    assert(pt.z < 5);
    if (FloatCmp::eq(pt.x, originalVertices_[pt.z].x) and FloatCmp::eq(pt.y, originalVertices_[pt.z].y) and
        not isAlreadyThere(pt.z))
      vertices_.emplace_back(Vertex::HostVertex(originalVertices_[pt.z], pt.z));
  }
  void addOriginalVertex(const size_t hostIdx) {
    if (not isAlreadyThere(hostIdx))
      vertices_.emplace_back(Vertex::HostVertex(originalVertices_[hostIdx], hostIdx));
  }

  void addNewVertex(const int edgeIdx, const Clipper2Lib::PointD& pt, const size_t trimmingCurveZ) {
    if (not isAlreadyThere(pt))
      vertices_.emplace_back(Vertex::NewVertex(edgeIdx, pt, trimmingCurveZ));
  }
  void addInsideVertex(const Clipper2Lib::PointD& pt, const size_t curveIndexI, const size_t curveIndexJ,
                       size_t loopIdx) {
    vertices_.emplace_back(Vertex::InsideVertex(pt, curveIndexI, curveIndexJ, loopIdx));
  }

  void finish(const auto& patchTrimData) {
    // Sort the points in counter clockwise manner such that the first point is in the lower left
    const auto minVertex = lowerLeftVertex();
    std::ranges::sort(vertices_, [&](const Vertex& aV, const Vertex& bV) { return comparePoints(aV, bV, minVertex); });

    // Resort that we always start on a hostVertex (if there is one)
    if (const auto it = std::ranges::find_if(vertices_, [](const Vertex& vV) { return vV.isHost(); });
        it != vertices_.end())
      std::ranges::rotate(vertices_, it);

    // Now resort InsideVertices (sometimes they end up on the wrong position)
    for (const auto i : std::views::iota(0ul, vertices_.size())) {
      if (vertices_[i].isInside()) {
        const auto insideVertex = vertices_[i];
        auto it                 = vertices_.begin() + static_cast<long>(i);
        const size_t indexBeforeIt =
            it == vertices_.begin() ? vertices_.size() - 1 : std::ranges::distance(vertices_.begin(), it) - 1;
        const size_t indexAfterIt =
            it == std::prev(vertices_.end()) ? 0 : std::ranges::distance(vertices_.begin(), it) + 1;

        // Check vertex before if they line up correctly
        const auto vertexBefore = vertices_[indexBeforeIt];
        const auto vertexAfter  = vertices_[indexAfterIt];

        if (vertexBefore.isNew() and vertexAfter.isNew()) {
          auto indices = patchTrimData.getIndices(vertexBefore.zValue());
          if (indices.loop == insideVertex.loop() and indices.curve == insideVertex.formerCurve())
            continue;
        }
        // If we end up here sort the vertex to a correct place
        auto correctVertexIt = findCorrespondingVertexIt(insideVertex, patchTrimData);
        assert(correctVertexIt != vertices_.end() && "No corresponding newVertex found");

        // Add the insideVertex at correctVertexIt+1
        const auto correctItForInsideVertex = correctVertexIt + 1;
        if (std::distance(vertices_.begin(), correctItForInsideVertex) < i) {
          vertices_.insert(correctItForInsideVertex, insideVertex);
          vertices_.erase(vertices_.begin() + static_cast<long>(i + 1));
        } else {
          std::ranges::rotate(it, it + 1, correctItForInsideVertex);
        }
      }
    }
  }

  void report() const {
    std::cout << "Vertices found\n";

    for (auto& vV : vertices_)
      std::cout << vV << std::endl;
  }

  std::vector<Vertex> vertices_{};

private:
  std::vector<Clipper2Lib::PointD> originalVertices_;

  auto isAlreadyThere(const Clipper2Lib::PointD& pt) -> bool {
    const auto it = std::ranges::find_if(vertices_, [&](const Vertex& vertex) {
      return FloatCmp::eq(vertex.pt.x, pt.x, 1e-8) and FloatCmp::eq(vertex.pt.y, pt.y, 1e-8);
    });
    return it != vertices_.end();
  }
  auto isAlreadyThere(const size_t hostIdx) -> bool {
    const auto it = std::ranges::find_if(vertices_, [&](const Vertex& vV) {
      if (vV.isHost())
        return vV.hostId() == hostIdx;
      return false;
    });
    return it != vertices_.end();
  }

  // Sort the points based on polar angle with respect to the reference point
  inline static auto comparePoints = [](const Vertex& aV, const Vertex& bV, const Vertex& referenceV) -> bool {
    auto polarAngle = [](const auto& p, const auto& reference) -> double {
      return atan2(p.y - reference.y, p.x - reference.x);
    };

    const auto a         = aV.pt;
    const auto b         = bV.pt;
    const auto reference = referenceV.pt;
    const double angleA  = polarAngle(a, reference);
    const double angleB  = polarAngle(b, reference);

    if (angleA != angleB)
      return angleA < angleB;

    // If two points have the same angle, sort based on distance to the reference point
    return (a.x - reference.x) * (a.x - reference.x) + (a.y - reference.y) * (a.y - reference.y) <
           (b.x - reference.x) * (b.x - reference.x) + (b.y - reference.y) * (b.y - reference.y);
  };

  auto findCorrespondingVertexIt(const Vertex& insideVertex, const auto& patchTrimData) {
    return std::ranges::find_if(vertices_, [&](const Vertex& vV) {
      if (vV.isNew()) {
        const auto indices = patchTrimData.getIndices(vV.zValue());
        return indices.loop == insideVertex.loop() and indices.curve == insideVertex.formerCurve();
      } else
        return false;
    });
  }

  [[nodiscard]] const Vertex& lowerLeftVertex() const {
    return *std::ranges::min_element(vertices_, [&](const Vertex& aV, const Vertex& bV) {
      const auto a = aV.pt;
      const auto b = bV.pt;
      return (a.y < b.y) || (a.y == b.y && a.x < b.x);
    });
  }
};

} // namespace Dune::IGANEW::DefaultTrim::Util
