// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/hierarchicpatch/geometrykernel/controlpoint.hh>
#include <dune/iga/hierarchicpatch/utils/mdnet.hh>

namespace Dune::IGANEW {
  /** \brief class that holds all data regarding the NURBS structure
   *
   * @tparam dim Dimension of the patch
   * @tparam dimworld Dimension of the control point coordinates , i.e. where the patch lives in
   * @tparam NurbsGridLinearAlgebraTraits Traits where FixedVectorType is derived
   */
  template <std::size_t dim, std::size_t dimworld_, typename ScalarType = double>
  struct NURBSPatchData {
    static constexpr int patchDim = dim;
    static constexpr int dimworld = dimworld_;
    using GlobalCoordinateType    = FieldVector<ScalarType, dimworld>;
    using ControlPointType        = ControlPoint<GlobalCoordinateType>;
    using ControlPointNetType     = MultiDimensionNet<dim, ControlPointType>;

    NURBSPatchData() = default;
    NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const ControlPointNetType& controlPointsI,
                   const std::array<int, dim>& degreeInput)
        : knotSpans(knotSpansI), controlPoints(controlPointsI), degree(degreeInput) {}

    std::array<std::vector<double>, dim> knotSpans{};
    ControlPointNetType controlPoints{};
    std::array<int, dim> degree{};
  };
}  // namespace Dune::IGANEW
