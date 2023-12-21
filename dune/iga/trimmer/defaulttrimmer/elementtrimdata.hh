// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <variant>

namespace Dune::IGANEW::DefaultTrim {

  enum class ElementTrimFlag;

  template <typename GridFamily>
  struct ElementTrimDataImpl {
    static constexpr int dim      = GridFamily::Trimmer::mydimension;
    static constexpr int dimworld = GridFamily::Trimmer::dimensionworld;
    using ctype                   = typename GridFamily::Trimmer::ctype;

    using HostElement = typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity;

    using EdgeTrimmedParameterSpaceGeometry =
        typename GridFamily::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
    using EdgePatchGeometry = typename EdgeTrimmedParameterSpaceGeometry::PatchGeometry;

    using VertexGeometry = FieldVector<ctype, dim>;

    using EdgeVariant = std::variant<EdgeTrimmedParameterSpaceGeometry, int>;
    using EdgeData    = std::pair<int, EdgeVariant>;

    explicit ElementTrimDataImpl(auto flag) : flag_(flag) {}

    void addEdge(int hostEdgeIdx) { edgeData.emplace_back(std::make_pair(hostEdgeIdx, hostEdgeIdx)); }

    void addEdge(EdgePatchGeometry& edge, int hostEdgeIdx) {
      edgeData.emplace_back(std::make_pair(hostEdgeIdx, edge));
      vertexGeometries_.emplace_back(edge.corners(0));
      vertexGeometries_.emplace_back(edge.corners(1));
    }

    [[nodiscard]] size_t numEdges() const { return edgeData.size(); }
    [[nodiscard]] EdgeData edgeInfo(int i) const {
      assert(edgeData.size() > i);
      return edgeData[i];
    }
    [[nodiscard]] ElementTrimFlag flag() const { return flag_; }

   private:
    ElementTrimFlag flag_;
    std::vector<EdgeData> edgeData{};
    std::vector<VertexGeometry> vertexGeometries_{};
  };

}  // namespace Dune::IGANEW::DefaultTrim
