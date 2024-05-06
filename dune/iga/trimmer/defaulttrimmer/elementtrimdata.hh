// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <clipper2/clipper.core.h>
#include <matplot/matplot.h>

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>

namespace Dune::IGANEW::DefaultTrim {

enum class ElementTrimFlag
{
  full,
  empty,
  trimmed
};

template <typename Grid>
struct ElementTrimDataImpl
{
  using GridFamily              = typename Grid::GridFamily;
  static constexpr int dim      = GridFamily::Trimmer::mydimension;
  static constexpr int dimworld = GridFamily::Trimmer::dimensionworld;
  using ctype                   = typename GridFamily::Trimmer::ctype;

  using HostEntity = typename GridFamily::Trimmer::TrimmerTraits::YASPGridType::Traits::template Codim<0>::Entity;

  using EdgeTrimmedParameterSpaceGeometry =
      typename GridFamily::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
  using EdgePatchGeometry = typename EdgeTrimmedParameterSpaceGeometry::PatchGeometry;

  using Vertex = FieldVector<ctype, dim>;

  enum class TrimmedHostEdgeDirection
  {
    HostNew,
    NewHost,
    NewNew
  };

  struct VertexInfo
  {
    bool isHost{};
    int idx{};

    std::optional<Vertex> geometry{};
  };

  struct EdgeInfo
  {
    bool isHost{};
    bool isTrimmed{};
    int idx{};

    std::optional<EdgePatchGeometry> geometry{};
    TrimmedHostEdgeDirection direction{};
  };

  explicit ElementTrimDataImpl(auto flag, const HostEntity& hostEntity)
      : flag_(flag),
        hostEntity_(hostEntity) {}

  // Delete default constructor, I think its implicitly delted anyway or not viable at least
  ElementTrimDataImpl() = delete;

  bool operator==(const ElementTrimDataImpl& other) const {
    return hostEntity_ == other.hostEntity_;
  }

  void addEdge(int idx) {
    edges_.emplace_back(EdgeInfo{.isHost = true, .isTrimmed = false, .idx = idx});
    vertices_.emplace_back(true, idx == 3 ? 0 : idx + 1, std::nullopt);
  }

  void addEdgeHostNew(int idx, EdgePatchGeometry& geometry, Vertex& v2) {
    edges_.emplace_back(EdgeInfo{.isHost    = true,
                                 .isTrimmed = true,
                                 .idx       = idx,
                                 .geometry  = transform(geometry),
                                 .direction = TrimmedHostEdgeDirection::HostNew});
    vertices_.emplace_back(false, newVertexCounter_++, v2);
  }

  void addEdgeNewNew(EdgePatchGeometry& geometry, Vertex& v2) {
    edges_.emplace_back(false, true, newEdgeCounter_++, transform(geometry));
    vertices_.emplace_back(false, newVertexCounter_++, v2);
  }

  void addEdgeNewHost(int idx, EdgePatchGeometry& geometry, int v2Idx) {
    edges_.emplace_back(EdgeInfo{.isHost    = true,
                                 .isTrimmed = true,
                                 .idx       = idx,
                                 .geometry  = transform(geometry),
                                 .direction = TrimmedHostEdgeDirection::NewHost});
    vertices_.emplace_back(true, v2Idx, std::nullopt);
  }
  void addEdgeNewNewOnHost(int idx, EdgePatchGeometry& geometry, Vertex& v2) {
    edges_.emplace_back(EdgeInfo{.isHost    = true,
                                 .isTrimmed = true,
                                 .idx       = idx,
                                 .geometry  = transform(geometry),
                                 .direction = TrimmedHostEdgeDirection::NewNew});
    vertices_.emplace_back(false, newVertexCounter_++, v2);
  }

  /* this is for testing purposes only */
  void drawResult(const std::string& filename, bool inParameterSpace, bool newFig = true) {
    if (flag_ == ElementTrimFlag::empty)
      return;
    auto eleGeometry      = hostEntity_.geometry();
    auto lowerLeftCorner  = inParameterSpace ? eleGeometry.corner(0) : Dune::FieldVector<double, 2>({0, 0});
    auto upperRightCorner = eleGeometry.corner(3);
    Dune::DiagonalMatrix<double, 2> scaling({inParameterSpace ? 1.0 / (upperRightCorner[0] - lowerLeftCorner[0]) : 1,
                                             inParameterSpace ? 1.0 / (upperRightCorner[1] - lowerLeftCorner[1]) : 1});

    if (newFig)
      auto figure = matplot::figure(true);
    matplot::hold("on");

    matplot::text(eleGeometry.center()[0], eleGeometry.center()[1], std::to_string(edges_.size()))
        ->font("Arial")
        .font_size(18);

    constexpr std::array idxLookUp = {0, 1, 3, 2};

    const double scale = eleGeometry.volume() / 5;

    auto plotEllipse = [&](const Vertex& v) {
      const auto w = 0.025 * scale;
      const auto c = matplot::ellipse(v[0] - (w / 2), v[1] - (w / 2), w, w);
      c->color("blue");
    };

    auto plotLine = [&](std::vector<Vertex>& vs, bool thin = false) {
      std::vector<double> x;
      std::ranges::transform(vs, std::back_inserter(x), [](auto& v) { return v[0]; });

      std::vector<double> y;
      std::ranges::transform(vs, std::back_inserter(y), [](auto& v) { return v[1]; });
      if (thin)
        matplot::plot(x, y)->line_width(0.5).color("grey");
      else
        matplot::plot(x, y)->line_width(2).color("red");

      // Annotate
      if (thin)
        return;
      matplot::text(x.front(), y.front(), "S")->font("Arial").color("blue").font_size(18);
    };

    for (auto c : edgeLookUp) {
      Dune::FieldVector<double, 2> firstVertexPos  = eleGeometry.corner(c[0]) - lowerLeftCorner;
      Dune::FieldVector<double, 2> secondVertexPos = eleGeometry.corner(c[1]) - lowerLeftCorner;
      scaling.mv(firstVertexPos, firstVertexPos);
      scaling.mv(secondVertexPos, secondVertexPos);
      std::vector<Vertex> vs{firstVertexPos, secondVertexPos};
      plotLine(vs, true);
    }

    for (auto& v : vertices_) {
      if (v.geometry.has_value()) {
        Dune::FieldVector<double, 2> pos = (v.geometry.value() - lowerLeftCorner);
        scaling.mv(pos, pos);
        plotEllipse(pos);
      } else {
        Dune::FieldVector<double, 2> pos = eleGeometry.corner(idxLookUp[v.idx]) - lowerLeftCorner;
        scaling.mv(pos, pos);
        plotEllipse(pos);
      }
    }

    for (auto& curve : edges_) {
      if (curve.geometry.has_value()) {
        std::vector<Vertex> vs;
        for (double u : Utilities::linspace(curve.geometry.value().domain().front(), 40)) {
          Dune::FieldVector<double, 2> pos = (curve.geometry.value().global(u) - lowerLeftCorner);
          scaling.mv(pos, pos);
          vs.push_back(pos);
        }
        plotLine(vs);
      } else {
        auto cornerIdx                               = edgeLookUp[curve.idx];
        Dune::FieldVector<double, 2> firstVertexPos  = (eleGeometry.corner(cornerIdx.front()) - lowerLeftCorner);
        Dune::FieldVector<double, 2> secondVertexPos = (eleGeometry.corner(cornerIdx.back()) - lowerLeftCorner);
        scaling.mv(firstVertexPos, firstVertexPos);
        scaling.mv(secondVertexPos, secondVertexPos);
        std::vector<Vertex> vs{firstVertexPos, secondVertexPos};
        plotLine(vs);
      }
    }

    matplot::axis(matplot::equal);
    static int elleC = 0;
    if (newFig)
      matplot::save(filename + std::to_string(elleC++) + (inParameterSpace ? "inParameterSpace" : ""), "gif");
  }
  bool checkInside(const Dune::FieldVector<double, 2>& local) const {
    Clipper2Lib::PointD p(local[0], local[1]);
    const auto result = Clipper2Lib::PointInPolygon<double>(p, path_);
    return result == Clipper2Lib::PointInPolygonResult::IsInside or result == Clipper2Lib::PointInPolygonResult::IsOn;
  }

  double volume() const {
    return Clipper2Lib::Area(path_);
  }

  // Getter
  [[nodiscard]] ElementTrimFlag flag() const {
    return flag_;
  }
  [[nodiscard]] const HostEntity& hostEntity() const {
    return hostEntity_;
  }
  [[nodiscard]] const std::vector<VertexInfo>& vertices() const {
    return vertices_;
  }
  [[nodiscard]] const std::vector<EdgeInfo>& edges() const {
    return edges_;
  }

  VertexInfo vertex(int i) const {
    return vertices_[i];
  }
  const EdgeInfo& edge(int i) const {
    return edges_[i];
  }

  [[nodiscard]] unsigned int size(unsigned int codim) const {
    if (codim == 2)
      return vertices_.size();
    if (codim == 1)
      return edges_.size();
    DUNE_THROW(NotImplemented, "Size Function is only implemented for edges and vertices");
  }

private:
  static constexpr std::array<std::array<int, 2>, 4> edgeLookUp{
      std::array{0, 1},
      {1, 3},
      {3, 2},
      {2, 0}
  };

public:
  auto& finalize() {
    sampleCurves();
    return *this;
  }

private:
  auto sampleCurves() {
    auto eleGeo           = hostEntity_.geometry();
    auto lowerLeftCorner  = eleGeo.corner(0);
    auto upperRightCorner = eleGeo.corner(3);
    Dune::DiagonalMatrix<double, 2> scaling(
        {1.0 / (upperRightCorner[0] - lowerLeftCorner[0]), 1.0 / (upperRightCorner[1] - lowerLeftCorner[1])});

    for (auto& edge : edges_) {
      if (edge.geometry.has_value()) {
        for (const auto v : Utilities::linspace(edge.geometry->domain()[0], 100)) {
          auto fV = edge.geometry->global({v}) - lowerLeftCorner;
          scaling.mv(fV, fV);
          path_.emplace_back(fV[0], fV[1]);
        }
      } else {
        auto cornerIdx = edgeLookUp[edge.idx];

        auto firstVertex  = eleGeo.corner(cornerIdx.front()) - lowerLeftCorner;
        auto secondVertex = eleGeo.corner(cornerIdx.back()) - lowerLeftCorner;
        scaling.mv(firstVertex, firstVertex);
        scaling.mv(secondVertex, secondVertex);

        path_.emplace_back(firstVertex[0], firstVertex[1]);
        path_.emplace_back(secondVertex[0], secondVertex[1]);
      }
    }
  }
  ElementTrimFlag flag_;
  HostEntity hostEntity_;
  std::vector<VertexInfo> vertices_{};
  std::vector<EdgeInfo> edges_{};
  Clipper2Lib::PathD path_;

  int newVertexCounter_ = 4;
  int newEdgeCounter_   = 4;
};

} // namespace Dune::IGANEW::DefaultTrim
