// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace Dune::IGA::GeometryKernel {

/**
 *  @brief This returns the center of mass for a polygon
 *
 *  @details See https://en.wikipedia.org/wiki/Centroid#Of_a_polygon for details
 */
template <GeometryConcept Geometry>
auto centerOfMass(const Geometry& geometry) -> FieldVector<typename Geometry::ctype, Geometry::coorddimension> {
  int n = geometry.corners();

  using Point = FieldVector<typename Geometry::ctype, Geometry::coorddimension>;
  std::vector<Point> corners;
  corners.reserve(n);

  for (const auto i : Dune::range(n))
    corners.push_back(geometry.corner(i));

  double a  = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  for (const auto i : Dune::range(n - 1)) {
    cx += (corners[i][0] + corners[i + 1][0]) * (corners[i][0] * corners[i + 1][1] - corners[i + 1][0] * corners[i][1]);
    cy += (corners[i][1] + corners[i + 1][1]) * (corners[i][0] * corners[i + 1][1] - corners[i + 1][0] * corners[i][1]);
    a += corners[i][0] * corners[i + 1][1] - corners[i + 1][0] * corners[i][1];
  }
  // Last index
  cx += (corners[n - 1][0] + corners[0][0]) * (corners[n - 1][0] * corners[0][1] - corners[0][0] * corners[n - 1][1]);
  cy += (corners[n - 1][1] + corners[0][1]) * (corners[n - 1][0] * corners[0][1] - corners[0][0] * corners[n - 1][1]);
  a += corners[n - 1][0] * corners[0][1] - corners[0][0] * corners[n - 1][1];

  a *= 0.5;
  cx *= 1 / (6 * a);
  cy *= 1 / (6 * a);
  return {cx, cy};
}

template <int dim, int dimworld, typename ctype = double>
auto centerOfMass(const std::vector<MultiLinearGeometry<ctype, dim, dimworld>>& elements)
    -> FieldVector<ctype, dimworld> {
  double xs        = 0.0;
  double ys        = 0.0;
  double totalArea = 0.0;

  for (const auto& triangle : elements) {
    auto center = triangle.center();
    auto area   = triangle.volume();

    totalArea += area;
    xs += area * center[0];
    ys += area * center[1];
  }
  return {xs / totalArea, ys / totalArea};
}

} // namespace Dune::IGA::GeometryKernel