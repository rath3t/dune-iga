// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>

// todo is this tested?
namespace Dune::IGA::GeometryKernel {
template <int dim_, int dimworld_, typename ScalarType = double>
auto transform(
    const NURBSPatch<dim_, dimworld_, ScalarType>& patch,
    const std::array<Utilities::Domain<double>, dim_>& newDomain = std::array<Utilities::Domain<double>, dim_>{}) {
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

template <int dim, int dimworld, typename ParameterSpaceGeometry, typename ScalarType = double>
auto transformToSpan(const NURBSPatch<dim, dimworld, ScalarType>& patch, const ParameterSpaceGeometry& span) {
  assert(span.corners() == 4);

  Dune::FieldVector<ScalarType, dimworld> offset = span.corner(0);
  std::array<ScalarType, dimworld> scaling       = {span.corner(1)[0] - span.corner(0)[0],
                                                    span.corner(2)[1] - span.corner(0)[1]};

  using PatchData     = typename NURBSPatch<dim, dimworld, ScalarType>::PatchData;
  PatchData patchData = patch.patchData();

  std::ranges::for_each(patchData.controlPoints.directGetAll(), [&]<typename CP>(CP& cp) {
    typename CP::VectorType local{};
    for (const auto i : Dune::range(dimworld)) {
      local[i] = (cp.p[i] - offset[i]) / scaling[i];
      local[i] = std::clamp(local[i], 0.0, 1.0);
    }
    cp.p = local;
  });

  return NURBSPatch<dim, dimworld, ScalarType>(patchData);
}

} // namespace Dune::IGA::GeometryKernel
