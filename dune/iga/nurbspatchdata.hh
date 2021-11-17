//
// Created by lex on 15.11.21.
//

#pragma once

#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/traits.hh>
#include <dune/iga/multidimensionNet.hh>



namespace Dune::IGA{
/** \brief class that holds all data regarding the NURBS structure */
template <std::integral auto  dim, std::integral auto  dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double, dim, dimworld>>
class NURBSPatchData {
public:
  using GlobalCoordinateType      = typename NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
  using ControlPointType = ControlPoint<GlobalCoordinateType>;

  using ControlPointNetType = MultiDimensionNet<dim, ControlPointType>;
  NURBSPatchData()          = default;
  /** \brief constructor for a NURBSPatchData from knots, control points, weights and order
     *
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order order of the NURBS structure for each dimension
   */
  NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const ControlPointNetType& controlPointsI,
                 const std::array<int, dim>& orderI)
      : knotSpans(knotSpansI), controlPoints(controlPointsI), order(orderI) {}

  std::array<std::vector<double>, dim> knotSpans;
  ControlPointNetType controlPoints;
  std::array<int, dim> order;
};
}