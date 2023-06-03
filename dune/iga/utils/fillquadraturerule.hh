// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/geometry/quadraturerules.hh>

namespace Dune::IGA {
  template <typename TrimmedElementRepresentationType, int dim>
  void fillQuadratureRuleImpl(Dune::QuadratureRule<double, dim>& vector,
                              const TrimmedElementRepresentationType& elementRepr, int order,
                              const QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
    vector.clear();

    auto gridView = elementRepr.gridView();

    for (auto subElement : elements(gridView)) {
      auto subElementGeo = subElement.geometry();

      const auto& rule = Dune::QuadratureRules<double, dim>::rule(subElement.type(), order, qt);
      for (auto ip : rule) {
        auto globalInSpan = subElementGeo.global(ip.position());
        vector.emplace_back(globalInSpan, ip.weight() * subElementGeo.integrationElement(ip.position()));
      }
    }

    if (qt != QuadratureType::Enum::GaussLobatto) return;
    struct Comparator {
      bool operator()(const Dune::FieldVector<double, 2>& p1, const Dune::FieldVector<double, 2>& p2) const {
        if (p1[0] < p2[0]) return true;
        if (p1[0] > p2[0]) return false;
        return p1[1] < p2[1];
      };
    };
    std::map<Dune::FieldVector<double, 2>, double, Comparator> uniquePoints;

    for (const auto& gp : vector) {
      auto [it, success] = uniquePoints.try_emplace(gp.position(), gp.weight());
      if (not success) it->second += gp.weight();
    }

    vector.clear();
    for (const auto& [pos, weight] : uniquePoints) {
      vector.emplace_back(pos, weight);
    }
  }
}  // namespace Dune::IGA
