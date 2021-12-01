//
// Created by lex on 16.11.21.
//

#pragma once

#include <dune/grid/common/capabilities.hh>
#include <dune/iga/concepts.hh>
//#include <dune/grid/test/checkidset.hh>
//#include <dune/grid/test/gridcheck.hh>

namespace Dune::IGA
{
  template<int dim,int dimworld,NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
  class NURBSGrid;
}
namespace Dune::Capabilities {
  template <std::integral auto dim, std::integral auto dimworld, int codim,Dune::IGA::NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
  requires(codim == 0 || (codim == 1 && dim < 3) || (codim == 2 && dim == 3) || (codim == 1 && dim == 3)
           || codim == dim) struct hasEntity<Dune::IGA::NURBSGrid<dim, dimworld,NurbsGridLinearAlgebraTraits>, codim> {
    static const bool v = true;
  };
}  // namespace Dune::Capabilities

template <std::integral auto dim, std::integral auto dimworld,Dune::IGA::NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGA::NURBSGrid<dim, dimworld,NurbsGridLinearAlgebraTraits>> : public std::true_type {};
template <std::integral auto dim, std::integral auto dimworld,Dune::IGA::NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
struct EnableLevelIntersectionIteratorCheck<Dune::IGA::NURBSGrid<dim, dimworld,NurbsGridLinearAlgebraTraits>> {
  static const bool v = true;
};
