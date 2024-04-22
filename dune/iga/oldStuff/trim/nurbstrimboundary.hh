// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <clipper2/clipper.core.h>
#include <ranges>
#include <utility>

#include "dune/iga/io/ibra/ibrageometry.hh"
#include "dune/iga/nurbspatchgeometry.hh"

namespace Dune::IGA {

class Boundary
{
public:
  static constexpr int worldDim = 2;
  static constexpr int dim      = 1;

  using Geometry            = NURBSPatchGeometry<dim, worldDim>;
  using ControlPointType    = typename Geometry::ControlPointType;
  using ControlPointNetType = Dune::IGA::MultiDimensionalNet<dim, ControlPointType>;
  using PatchData           = Dune::IGA::NURBSPatchData<dim, worldDim>;

  using Point = Dune::FieldVector<double, worldDim>;

  Geometry nurbsGeometry;
  using DomainType = Utilities::Domain<double>;
  DomainType domain{};
  std::array<Point, 2> endPoints;

  Boundary() = default;

  Boundary(Geometry _nurbsGeometry, const Utilities::Domain<double>& _domain)
      : nurbsGeometry(std::move(_nurbsGeometry)),
        domain(_domain),
        endPoints({nurbsGeometry(domain.left()), nurbsGeometry(domain.right())}) {
    assert(domain.right() > domain.left());
  }

  explicit Boundary(Ibra::BrepTrim& _trim)
      : nurbsGeometry(geometryFromTrim(_trim.geometry)),
        domain(_trim.domain),
        endPoints({nurbsGeometry(domain.left()), nurbsGeometry(domain.right())}) {
  }

  explicit Boundary(const Point& a, const Point& b)
      : nurbsGeometry(lineGeometryFromPoints(a, b)),
        domain(nurbsGeometry.domain()[0]),
        endPoints({a, b}) {
  }

  // Helper classes for construction of nurbsGeometry
private:
  static Geometry geometryFromTrim(Ibra::Curve2D& _curve) {
    const auto _cp = _curve.transformControlPoints()[0];
    std::array<int, dim> dimSize{static_cast<int>(_cp.size())};

    // Construct patch Data
    std::array<std::vector<double>, dim> knotSpans{_curve.compileKnotVectors()};
    auto controlNet{ControlPointNetType(dimSize, _cp)};

    return Geometry(PatchData(knotSpans, controlNet, std::array<int, 1>{_curve.degree}));
  }

  /** @brief creates a line geometry from a to b */
  static Geometry lineGeometryFromPoints(const Point& a, const Point b) {
    std::vector<Point> _controlPoints{a, b};
    std::array<ControlPointType, worldDim> _cp{
        ControlPointType{.p{_controlPoints.front()}, .w = 1},
        { .p{_controlPoints.back()}, .w = 1}
    };
    std::array<int, 1> dimSize{2};

    // Construct patch Data
    std::array<std::vector<double>, dim> knotSpans{
        std::vector<double>{0.0, 0.0, 1.0, 1.0}
    };
    auto controlNet{ControlPointNetType(dimSize, _cp)};

    return Geometry(PatchData(knotSpans, controlNet, std::array<int, 1>{1}));
  }

public:
  [[nodiscard]] int degree() const {
    return nurbsGeometry.degree()[0];
  };

  enum class EdgeOrientation
  {
    u,
    v,
    Unknown
  };
  [[nodiscard]] EdgeOrientation getEdgeOrientation(const double tolerance = 1e-6) const {
    if (Dune::FloatCmp::eq(endPoints[0][0], endPoints[1][0], tolerance))
      return EdgeOrientation::v;
    else if (Dune::FloatCmp::eq(endPoints[0][1], endPoints[1][1], tolerance))
      return EdgeOrientation::u;
    else
      return EdgeOrientation::Unknown;
  }

  template <typename ctype = double>
  [[nodiscard]] Clipper2Lib::Path<ctype> path(unsigned int samples = 200, const int intScale = 0) const {
    Clipper2Lib::Path<ctype> path;

    auto scaler = [intScale](const auto x) -> ctype { return x * static_cast<ctype>(std::pow(10, intScale)); };

    if (degree() == 1 && endPoints.size() == 2) {
      path.emplace_back(scaler(endPoints.front()[0]), scaler(endPoints.front()[1]));
      path.emplace_back(scaler(endPoints.back()[0]), scaler(endPoints.back()[1]));
      return path;
    }

    path.reserve(samples);
    auto linS = Utilities::linspace(domain, samples);
    for (auto u : linS) {
      auto vertex = nurbsGeometry(u);
      path.emplace_back(scaler(vertex[0]), scaler(vertex[1]));
    }
    return path;
  }
};

using BoundaryLoop = std::vector<Boundary>;

class TrimData
{
public:
  std::vector<BoundaryLoop> boundaryLoops;

  TrimData() = default;

  void addLoop(const std::vector<Boundary>& _boundaries) {
    boundaryLoops.push_back({_boundaries});
  }

  [[nodiscard]] size_t numBoundaries() const {
    const auto joinedView = std::ranges::join_view(boundaryLoops);
    return std::ranges::distance(joinedView.begin(), joinedView.end());
  }
};

} // namespace Dune::IGA
