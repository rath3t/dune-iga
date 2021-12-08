//
// Created by lex on 15.11.21.
//

#pragma once

#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/multidimensionNet.hh>

namespace Dune::IGA {
  /** \brief class that holds all data regarding the NURBS structure
   *
   * @tparam dim Dimension of the patch
   * @tparam dimworld Dimension of the control point coordinates , i.e. where the patch lives in
   * @tparam NurbsGridLinearAlgebraTraits Traits where FixedVectorType is derived
   */
  template <std::size_t dim, std::size_t dimworld, LinearAlgebra NurbsGridLinearAlgebraTraits = DuneLinearAlgebraTraits<double>>
  struct NURBSPatchData {
    using GlobalCoordinateType = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<dimworld>;
    using ControlPointType     = ControlPoint<GlobalCoordinateType>;
    using ControlPointNetType  = MultiDimensionNet<dim, ControlPointType>;

    NURBSPatchData() = default;
    NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const ControlPointNetType& controlPointsI,
                   const std::array<int, dim>& degreeInput)
        : knotSpans(knotSpansI), controlPoints(controlPointsI), degree(degreeInput) {}

    std::array<std::vector<double>, dim> knotSpans;
    ControlPointNetType controlPoints;
    std::array<int, dim> degree;
  };
}  // namespace Dune::IGA