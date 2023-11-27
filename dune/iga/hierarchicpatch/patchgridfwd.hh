
#pragma once
#include <dune/grid/common/grid.hh>
namespace Dune::IGANEW {

  template <int dim, int dimworld, typename TrimmerType>
  class PatchGrid;

  template <int dim, int dimworld, typename TrimmerType>
  struct PatchGridFamily;

  template <int codim, int dim, class GridImp>
  class PatchGridEntity;
}  // namespace Dune::IGANEW
