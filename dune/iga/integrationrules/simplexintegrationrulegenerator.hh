// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/integrationrules/simplexgenerator.hh>

namespace Dune::IGA {

template <typename GridImp>
struct SimplexIntegrationRuleGenerator
{
  using ParameterSpaceElement = typename GridImp::ParameterSpace::template Codim<0>::ParameterSpaceGridEntity;

  static constexpr int dim = GridImp::dimension;

  using Point   = FieldVector<double, dim>;
  using Element = MultiLinearGeometry<double, dim, dim>;
  using Index   = std::uint64_t;

  using SimplexGeneratorImpl = SimplexGenerator<GridImp>;
  using Parameters           = typename SimplexGeneratorImpl::Parameters;

  static auto createIntegrationRule(const ParameterSpaceElement& element, int quadratureOrder,
                                    const Parameters& parameters  = Parameters{},
                                    const QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
    if constexpr (GridImp::ParameterSpace::isAlwaysTrivial)
      return QuadratureRules<double, dim>::rule(element.type(), quadratureOrder, qt);
    else {
      if (not element.isTrimmed()) {
        return QuadratureRules<double, dim>::rule(element.type(), quadratureOrder, qt);
      }
      auto [elements, vertices, indices] = SimplexGeneratorImpl::createSimplicies(element, parameters);
      return makeQuadratureRule(elements, quadratureOrder, qt);
    }
  }

private:
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

} // namespace Dune::IGA
