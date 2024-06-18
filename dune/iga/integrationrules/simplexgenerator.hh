// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <mapbox/earcut.hpp>

#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>

namespace Dune::IGA {

template <typename GridImp>
struct SimplexGenerator
{
  using ParameterSpaceElement = typename GridImp::ParameterSpace::template Codim<0>::ParameterSpaceGridEntity;

  static constexpr int dim = GridImp::dimension;

  using Point   = FieldVector<double, dim>;
  using Element = MultiLinearGeometry<double, dim, dim>;
  using Index   = std::uint64_t;

  struct Parameters
  {
    int boundaryDivisions{5};
    double targetAccuracy{std::numeric_limits<double>::max()};

    bool useAdaptiveDivisions() const {
      return targetAccuracy < 1;
    }
  };

  static auto createSimplicies(const ParameterSpaceElement& element, const Parameters& parameters = Parameters{})
      -> std::tuple<std::vector<Element>, std::vector<Point>, std::vector<Index>> {
    std::vector<Point> vertices{};
    std::vector<Element> elements{};

    auto& trimData  = element.trimData();
    auto hostEntity = element.getHostEntity();

    for (const auto& edgeInfo : trimData.edges()) {
      if (edgeInfo.isTrimmed) {
        auto localGeometry = GeometryKernel::transformToSpan(edgeInfo.geometry.value(), hostEntity.geometry());
        splitBoundary(localGeometry, parameters, vertices);
      } else {
        auto refElement = referenceElement(hostEntity);
        switch (edgeInfo.idx) {
          case 0:
            vertices.push_back(refElement.template geometry<2>(0).center());
            break;
          case 1:
            vertices.push_back(refElement.template geometry<2>(1).center());
            break;
          case 2:
            vertices.push_back(refElement.template geometry<2>(3).center());
            break;
          case 3:
            vertices.push_back(refElement.template geometry<2>(2).center());
            break;
          default:
            assert(edgeInfo.idx < 4);
            __builtin_unreachable();
        }
      }
    }

    auto indices = triangulate(vertices);
    for (auto it = indices.begin(); it < indices.end(); it += 3)
      elements.emplace_back(GeometryTypes::triangle, std::vector<FieldVector<double, dim>>{
                                                         vertices[*it], vertices[*(it + 1)], vertices[*(it + 2)]});

    return std::make_tuple(elements, vertices, indices);
  }

private:
private:
  static void splitBoundary(const auto& localGeometry, const Parameters& parameters, std::vector<Point>& points) {
    if (localGeometry.affine()) {
      points.push_back(localGeometry.global({0.0}));
      return;
    }

    if (not parameters.useAdaptiveDivisions()) {
      points.reserve(parameters.boundaryDivisions);
      for (auto local : Utilities::linspace(0.0, 1.0, parameters.boundaryDivisions))
        points.push_back(localGeometry.global({local}));
    } else {
      // Guess initial amount of segments (a linear segment of length 1 should be divided by 1 divions)
      double curveLength = computeCurveLength(localGeometry);
      auto segments      = static_cast<unsigned int>(std::ceil(1 * curveLength * localGeometry.degree()[0]));

      bool converged = false;
      std::vector<Point> newPoints;
      while (not converged) {
        std::tie(converged, newPoints) =
            checkForConvergence(localGeometry, segments, curveLength, parameters.targetAccuracy);
        ++segments;
        if (converged)
          points.insert(points.end(), newPoints.begin(), newPoints.end());
      }
    }
    // points.push_back(localGeometry.global({1.0}));

    // Remove last element to avoid duplication
    points.pop_back();
  }

  static auto triangulate(const std::vector<Point>& points) -> std::vector<Index> {
    std::vector<std::vector<Point>> polygonInput;
    polygonInput.push_back(points);

    return mapbox::earcut<Index>(polygonInput);
  }

  // Function to calculate the distance between two points
  static auto checkForConvergence(const auto& localGeometry, unsigned int numSegments, double curveLength,
                                  double targetAccuracy) -> std::pair<bool, std::vector<Point>> {
    auto [totalDistance, points] = getPointsAndDistance(localGeometry, numSegments);
    auto converged               = std::abs((totalDistance - curveLength) / curveLength) < targetAccuracy;
    return std::make_pair(converged, points);
  }
  static auto computeCurveLength(const auto& localGeometry) -> double {
    unsigned int numSegments = localGeometry.degree()[0] * 25;

    return std::get<0>(getPointsAndDistance(localGeometry, numSegments));
  }

  static auto getPointsAndDistance(const auto& localGeometry,
                                   unsigned int numSegments) -> std::pair<double, std::vector<Point>> {
    auto distance = [](const Point& p1, const Point& p2) {
      return std::sqrt(Dune::power(p2[0] - p1[0], 2) + Dune::power(p2[1] - p1[1], 2));
    };

    std::vector<Point> points;
    for (auto local : Utilities::linspace(0.0, 1.0, numSegments + 1))
      points.push_back(localGeometry.global({local}));

    double totalDistance = 0.0;
    for (auto i : Dune::range(1ul, points.size()))
      totalDistance += distance(points[i - 1], points[i]);

    return std::make_pair(totalDistance, points);
  }
};

} // namespace Dune::IGA

// Add support for Dune::FieldVector in Earcut
namespace mapbox::util {

template <typename T>
struct nth<0, Dune::FieldVector<T, 2>>
{
  inline static auto get(const Dune::FieldVector<T, 2>& t) {
    return t[0];
  };
};

template <typename T>
struct nth<1, Dune::FieldVector<T, 2>>
{
  inline static auto get(const Dune::FieldVector<T, 2>& t) {
    return t[1];
  };
};
} // namespace mapbox::util