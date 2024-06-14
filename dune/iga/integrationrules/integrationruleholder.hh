// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/integrationrules/simplexintegrationrulegenerator.hh>

namespace Dune::IGA {
namespace DefaultTrim {
  template <typename GridImp>
  struct DefaultIntegrationRuleGenerator
  {
    using Generator = SimplexIntegrationRuleGenerator<GridImp>;
    auto operator()(const auto& element, int order, QuadratureType::Enum qt = QuadratureType::GaussLegendre) const {
      return Generator::createIntegrationRule(element, order, parameters, qt);
    }

  private:
    typename Generator::Parameters parameters{.boundaryDivisions = Preferences::getInstance().boundaryDivisions(),
                                              .targetAccuracy    = Preferences::getInstance().targetAccuracy()};
  };
}; // namespace DefaultTrim

template <typename GridImp>
struct IntegrationRuleHolder
{
  using ParameterSpaceElement = typename GridImp::ParameterSpace::template Codim<0>::ParameterSpaceGridEntity;

  static constexpr int dim = GridImp::dimension;

  using FunctionType =
      std::function<QuadratureRule<double, dim>(const ParameterSpaceElement&, int, QuadratureType::Enum)>;

  IntegrationRuleHolder()
      : generator_(DefaultTrim::DefaultIntegrationRuleGenerator<GridImp>{}) {}

  void integrationRule(FunctionType generator) {
    generator_ = generator;
  }

  FunctionType integrationRule() const {
    return generator_;
  }

private:
  FunctionType generator_;
};
} // namespace Dune::IGA