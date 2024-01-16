// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>
namespace Dune::IGANEW::GeometryKernel {
  template <int dim_, int dimworld_, typename ScalarType = double>
  auto transform(const NURBSPatch<dim_, dimworld_, ScalarType>& patch,
                 const std::array<Utilities::Domain<double>, dim_>& newDomain
                 = std::array<Utilities::Domain<double>, dim_>{}) {
    auto patchData = patch.patchData();

    for (auto i = 0; i < dim_; ++i) {
      for (int j = 0; j < patchData.degree[i] + 1; ++j) {
        patchData.knotSpans[i][j]               = newDomain[i].left();
        *(patchData.knotSpans[i].end() - 1 - j) = newDomain[i].right();
      }

      std::for_each(patchData.knotSpans[i].begin() + patchData.degree[i] + 1,
                    patchData.knotSpans[i].end() - patchData.degree[i] - 1,
                    [&](auto& knotSpan) { knotSpan = mapToRange(knotSpan, patch.domain()[i], newDomain[i]); });
    }
    return NURBSPatch<dim_, dimworld_, ScalarType>(patchData);
  }
}  // namespace Dune::IGANEW::GeometryKernel
