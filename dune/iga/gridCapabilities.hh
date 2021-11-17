//
// Created by lex on 16.11.21.
//

#pragma once

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/checkindexset.hh>
namespace Dune::Capabilities
{

  template< std::integral auto  dim, std::integral auto dimworld >
  struct hasEntity<Dune::IGA::NURBSGrid<dim,dimworld>,0>
  {
    static const bool v = true;
  };


  template< std::integral auto dim, std::integral auto dimworld >
  struct hasEntityIterator<Dune::IGA::NURBSGrid<dim,dimworld>,0>
  {
    static const bool v = hasEntity<Dune::IGA::NURBSGrid<dim,dimworld>, 0>::v;
  };




}
template< std::integral auto dim, std::integral auto dimworld >
struct EnableLevelIntersectionIteratorCheck<Dune::IGA::NURBSGrid<dim,dimworld>>
{
  static const bool v = false;
};
