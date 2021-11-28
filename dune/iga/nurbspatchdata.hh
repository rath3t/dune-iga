//
// Created by lex on 15.11.21.
//

#pragma once

#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/multidimensionNet.hh>

namespace Dune::IGA {
  /** \brief class that holds all data regarding the NURBS structure */
  template <std::size_t dim, std::size_t dimworld,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double>>
  class NURBSPatchData {
  public:
    using GlobalCoordinateType = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<dimworld>;
    using ControlPointType     = ControlPoint<GlobalCoordinateType>;
    using ControlPointNetType  = MultiDimensionNet<dim, ControlPointType>;

    NURBSPatchData() = default;
    NURBSPatchData(const std::array<std::vector<double>, dim>& knotSpansI, const ControlPointNetType& controlPointsI,
                   const std::array<int, dim>& orderI)
        : knotSpans(knotSpansI), controlPoints(controlPointsI), order(orderI) {}

    std::array<std::vector<double>, dim> knotSpans;
    ControlPointNetType controlPoints;
    std::array<int, dim> order;
  };
}  // namespace Dune::IGA