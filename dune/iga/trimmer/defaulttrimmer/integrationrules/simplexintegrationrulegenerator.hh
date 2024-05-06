// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <mapbox/earcut.hpp>

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmerpreferences.hh>

namespace Dune::IGANEW::DefaultTrim {

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
    int maxBoundaryDivisions{5};
    double targetTolerance{1e-10}; // not used atm
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
        std::vector<Point> pointsT = splitBoundary(localGeometry, parameters);
        vertices.insert(vertices.end(), pointsT.begin(), pointsT.end());
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
  static std::vector<Point> splitBoundary(const auto& localGeometry, const Parameters& parameters) {
    std::vector<Point> points;

    if (localGeometry.affine()) {
      points.push_back(localGeometry.global({0.0}));
      return points;
    }

    points.reserve(parameters.maxBoundaryDivisions);

    // todo at the moment only boundaryDivisions
    for (auto local : Utilities::linspace(0.0, 1.0, parameters.maxBoundaryDivisions))
      points.push_back(localGeometry.global({local}));

    // Remove last element to avoid duplication
    points.pop_back();
    return points;
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
};

} // namespace Dune::IGANEW::DefaultTrim

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
