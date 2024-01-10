// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <matplot/matplot.h>

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

    /* this is for testing purposes only */
    void drawResult(const std::string& filename, auto eleGeometry) {
      if (static_cast<int>(flag_) != 2) return;

      auto figure = matplot::figure(true);
      matplot::hold("on");

      constexpr std::array<std::array<int, 2>, 4> edgeLookUp{std::array{0, 1}, {1, 3}, {3, 2}, {2, 0}};
      constexpr std::array idxLookUp = {0, 1, 3, 2};

      auto plotEllipse = [](const Vertex& v) {
        constexpr auto w = 0.025;
        const auto c     = matplot::ellipse(v[0] - (w / 2), v[1] - (w / 2), w, w);
        c->color("blue");
      };

      auto plotLine = [](std::vector<Vertex>& vs, bool thin = false) {
        std::vector<double> x;
        std::ranges::transform(vs, std::back_inserter(x), [](auto& v) { return v[0]; });

        std::vector<double> y;
        std::ranges::transform(vs, std::back_inserter(y), [](auto& v) { return v[1]; });
        if (thin)
          matplot::plot(x, y)->line_width(0.5).color("grey");
        else
          matplot::plot(x, y)->line_width(2).color("red");
      };

      for (auto c : edgeLookUp) {
        std::vector<Vertex> vs{eleGeometry.corner(c[0]), eleGeometry.corner(c[1])};
        plotLine(vs, true);
      }

      for (auto& v : vertices_) {
        if (v.geometry.has_value())
          plotEllipse(v.geometry.value());
        else
          plotEllipse(eleGeometry.corner(idxLookUp[v.idx]));
      }

      for (auto& curve : edges_) {
        if (curve.geometry.has_value()) {
          std::vector<Vertex> vs;
          for (double u : Utilities::linspace(curve.geometry.value().domain().front(), 20))
            vs.push_back(curve.geometry.value().global(u));
          plotLine(vs);
        } else {
          auto cornerIdx = edgeLookUp[curve.idx];
          std::vector<Vertex> vs{eleGeometry.corner(cornerIdx.front()), eleGeometry.corner(cornerIdx.back())};
          plotLine(vs);
        }
      }

      matplot::axis(matplot::equal);
      matplot::save(filename, "gif");
    }

    // Getter
    [[nodiscard]] ElementTrimFlag flag() const { return flag_; }

  private:
    ElementTrimFlag flag_;
    std::vector<VertexInfo> vertices_{};
    std::vector<EdgeInfo> edges_;

    int newVertexCounter_ = 4;
    int newEdgeCounter_   = 4;
  };

}  // namespace Dune::IGANEW::DefaultTrim
