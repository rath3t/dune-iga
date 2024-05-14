
#pragma once
#include <dune/grid/common/grid.hh>
namespace Dune::IGA {

template <int dim, int dimworld, template <int, int, typename> typename TrimmerType, typename ScalarType>
class PatchGrid;

template <int dim, int dimworld, template <int, int, typename> typename TrimmerType, typename ScalarType>
struct PatchGridFamily;

template <int codim, int dim, class GridImp>
class PatchGridEntity;
} // namespace Dune::IGA
