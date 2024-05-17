// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <mapbox/earcut.hpp>

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmerpreferences.hh>

namespace Dune::IGA::DefaultTrim {

template <typename GridImp>
struct SimplexIntegrationRuleGenerator
{
  using PatchElement = typename GridImp::Traits::template Codim<0>::Entity;

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

  static auto createIntegrationRule(const PatchElement& element, int quadratureOrder,
                                    const Parameters& parameters  = Parameters{},
                                    const QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
    if (not element.impl().getLocalEntity().isTrimmed()) {
      return QuadratureRules<double, dim>::rule(element.type(), quadratureOrder, qt);
    }

    auto [elements, vertices, indices] = createSimplicies(element, parameters);
    return makeQuadratureRule(elements, quadratureOrder, qt);
  }

  static auto createSimplicies(const PatchElement& element, const Parameters& parameters = Parameters{})
      -> std::tuple<std::vector<Element>, std::vector<Point>, std::vector<Index>> {
    std::vector<Point> vertices{};
    std::vector<Element> elements{};

    auto& trimData  = element.impl().getLocalEntity().trimData();
    auto hostEntity = element.impl().getLocalEntity().getHostEntity();

    for (const auto& edgeInfo : trimData.edges()) {
      if (edgeInfo.isTrimmed) {
        auto localGeometry         = GeometryKernel::transformToSpan(edgeInfo.geometry.value(), hostEntity.geometry());
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
      // Guess initial amount of segments (a quadratic segment of lenght 1 should be devided by 5 divions)
      double curveLength   = localGeometry.curveLength();
      auto initialSegments = static_cast<unsigned int>(std::ceil(5 * curveLength * localGeometry.degree()[0]));

      std::vector<Point> initialPoints;
      for (const auto i : Dune::range(initialSegments)) {
        double t = static_cast<double>(i) / initialSegments;
        initialPoints.push_back(localGeometry.global({t}));
      }

      points.push_back(initialPoints.front());

      for (size_t i_ = 0; i_ < initialPoints.size() - 1; ++i_) {
        const auto i = static_cast<double>(i_);
        refineSegments(points, localGeometry, i * (1.0 / initialSegments), (i + 1) * (1.0 / initialSegments),
                       parameters.targetAccuracy);
      }
    }
    points.push_back(localGeometry.global({1.0}));

    // Remove last element to avoid duplication
    points.pop_back();
  }

  static auto triangulate(const std::vector<Point>& points) -> std::vector<Index> {
    std::vector<std::vector<Point>> polygonInput;
    polygonInput.push_back(points);

    return mapbox::earcut<Index>(polygonInput);
  }
  static auto makeQuadratureRule(const std::vector<Element>& elements, int quadratureOrder,
                                 const QuadratureType::Enum qt) {
    assert(qt == QuadratureType::GaussLegendre);

    Dune::QuadratureRule<double, dim> vector{};

    for (auto& subElement : elements) {
      const auto& rule = Dune::QuadratureRules<double, dim>::rule(subElement.type(), quadratureOrder, qt);
      for (auto ip : rule) {
        auto globalInSpan = subElement.global(ip.position());
        vector.emplace_back(globalInSpan, ip.weight() * subElement.integrationElement(ip.position()));
      }
    }
    return vector;
  }

  // Function to calculate the distance between two points
  static double distance(const Point& p1, const Point& p2) {
    return std::sqrt(Dune::power(p2[0] - p1[0], 2) + Dune::power(p2[1] - p1[1], 2));
  }

  static void refineSegments(std::vector<Point>& points, const auto& localGeometry, double t_start, double t_end,
                             double tolerance) {
    Point p1 = localGeometry.global({t_start});
    Point p2 = localGeometry.global({t_end});

    // Calculate the midpoint of the parameter
    double t_mid = (t_start + t_end) / 2;
    Point p_mid  = localGeometry.global({t_mid});

    // Calculate the perpendicular distance from the midpoint to the line segment
    double numerator =
        std::abs((p2[1] - p1[1]) * p_mid[0] - (p2[0] - p1[0]) * p_mid[1] + p2[0] * p1[1] - p2[1] * p1[0]);
    double denominator     = distance(p1, p2);
    double dist_to_segment = numerator / denominator;

    if (dist_to_segment > tolerance) {
      // If the distance is greater than the tolerance, subdivide recursively
      refineSegments(points, localGeometry, t_start, t_mid, tolerance);
      refineSegments(points, localGeometry, t_mid, t_end, tolerance);
    } else {
      points.push_back(p2); // Include the endpoint of the segment
    }
  }
};

} // namespace Dune::IGA::DefaultTrim

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
