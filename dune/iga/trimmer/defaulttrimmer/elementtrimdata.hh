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

    using Vertex = FieldVector<ctype, dim>;

    struct VertexInfo {
      bool isHost;
      int idx;

      std::optional<Vertex> geometry;
    };

    struct EdgeInfo {
      bool isHost;
      bool isTrimmed;
      int idx;

      std::optional<EdgePatchGeometry> geometry;
    };

    explicit ElementTrimDataImpl(auto flag) : flag_(flag) {}

    void addEdge(int idx) {
      edges_.emplace_back(true, false, idx, std::nullopt);
      vertices_.emplace_back(true, idx + 1, std::nullopt);
    }

    void addEdgeHostNew(int idx, EdgePatchGeometry& geometry, Vertex& v2) {
      edges_.emplace_back(true, true, idx, geometry);
      vertices_.emplace_back(false, newVertexCounter_++, v2);
    }
    void addEdgeNewNew(EdgePatchGeometry& geometry, Vertex& v2) {
      edges_.emplace_back(false, true, newEdgeCounter_++, geometry);
      vertices_.emplace_back(false, newVertexCounter_++, v2);
    }
    void addEdgeNewHost(int idx, EdgePatchGeometry& geometry, int v2Idx) {
      edges_.emplace_back(true, true, idx, geometry);
      vertices_.emplace_back(true, v2Idx, std::nullopt);
    }

    // Getter
    [[nodiscard]] ElementTrimFlag flag() const { return flag_; }

   private:
    ElementTrimFlag flag_;
    std::vector<VertexInfo> vertices_{};
    std::vector<EdgeInfo> edges_;

    int newVertexCounter_ = 4;
    int newEdgeCounter_ = 4;
  };

}  // namespace Dune::IGANEW::DefaultTrim
