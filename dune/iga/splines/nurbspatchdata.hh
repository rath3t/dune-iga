// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/geometrykernel/controlpoint.hh>
#include <dune/iga/utils/mdnet.hh>

namespace Dune::IGANEW {

/**
 * @brief Struct that holds all data regarding the NURBS geometric structure
 *
 * @tparam dim Dimension of the patch
 * @tparam dimworld_ Dimension of the control point coordinates, i.e., where the patch lives in
 * @tparam ScalarType The type for the functions values and arguments, defaults to double
 */
template <int dim, int dimworld_, typename ScalarType = double>
struct NURBSPatchData
{
  static constexpr int patchDim = dim;                                ///< Dimension of the patch
  static constexpr int dimworld = dimworld_;                          ///< Dimension of the control point coordinates
  using GlobalCoordinateType    = FieldVector<ScalarType, dimworld>;  ///< Type for global coordinates
  using ControlPointType        = ControlPoint<GlobalCoordinateType>; ///< Type for control points
  using ControlPointNetType     = MultiDimensionalNet<dim, ControlPointType>; ///< Type for the net of control points

  /* @brief Default constructor */
  NURBSPatchData() = default;

  /**
   * @brief Constructor with explicit initialization
   *
   * @param knotSpansI Knot spans for each dimension
   * @param controlPointsI Net of control points
   * @param degreeInput Degree of the NURBS patch in each dimension
   */
  NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const ControlPointNetType& controlPointsI,
                 const std::array<int, dim>& degreeInput)
      : knotSpans(knotSpansI),
        controlPoints(controlPointsI),
        degree(degreeInput) {}

  std::array<std::vector<double>, dim> knotSpans{}; ///< Knot spans for each dimension
  ControlPointNetType controlPoints{};              ///< Net of control points
  std::array<int, dim> degree{};                    ///< Degree of the NURBS patch in each dimension
};

/**
 * @brief Deduction guide for creating NURBSPatchData instances without explicitly specifying template arguments
 */
template <std::size_t dim, std::size_t dimworld_, typename ScalarType = double>
NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI,
               const MultiDimensionalNet<dim, ControlPoint<FieldVector<ScalarType, dimworld_>>>& controlPointsI,
               const std::array<int, dim>& degreeInput) -> NURBSPatchData<dim, dimworld_, ScalarType>;

} // namespace Dune::IGANEW
