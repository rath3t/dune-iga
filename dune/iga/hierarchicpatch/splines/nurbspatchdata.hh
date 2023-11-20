// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/hierarchicpatch/geometrykernel/controlpoint.hh>
#include <dune/iga/hierarchicpatch/utils/mdnet.hh>

namespace Dune::IGANEW {
  /** \brief struct that holds all data regarding the NURBS geometric structure
   *
   * @tparam dim Dimension of the patch
   * @tparam dimworld Dimension of the control point coordinates , i.e. where the patch lives in
   * @tparam NurbsGridLinearAlgebraTraits Traits where FixedVectorType is derived
   */
  template <int dim, int dimworld_, typename ScalarType = double>
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

  //Deduction guide, since std::array, takes size_t as second template argument
  template <std::size_t dim, std::size_t dimworld_, typename ScalarType = double>
  NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const MultiDimensionNet<dim, ControlPoint<FieldVector<ScalarType, dimworld_>>>& controlPointsI,
                 const std::array<int, dim>& degreeInput) -> NURBSPatchData<dim,dimworld_,ScalarType>;


}  // namespace Dune::IGANEW
