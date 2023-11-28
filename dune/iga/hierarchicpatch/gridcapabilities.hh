
#pragma once

#include "patchgridfwd.hh"

#include <dune/grid/common/capabilities.hh>
namespace Dune {

  namespace Capabilities {
    /** \brief has entities for some codimensions as host grid
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename TrimmerType, int codim>
    struct hasEntity<const IGANEW::PatchGrid<dim, dimworld, TrimmerType>, codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = hasEntity<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType, int codim>
    struct hasEntity<IGANEW::PatchGrid<dim, dimworld, TrimmerType>, codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = hasEntity<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType, int codim>
    struct hasEntity<const Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>,
                     codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = hasEntity<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType, int codim>
    struct hasEntityIterator<IGANEW::PatchGrid<dim, dimworld, TrimmerType>, codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;

      static const bool v = hasEntityIterator<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType, int codim>
    struct hasEntityIterator<Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>,
                             codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = hasEntityIterator<ParameterSpaceGrid, codim>::v;
    };

    /** \brief PatchGrid can communicate when the host grid can communicate
     *  \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename TrimmerType, int codim>
    struct canCommunicate<IGANEW::PatchGrid<dim, dimworld, TrimmerType>, codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;

      static const bool v = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType, int codim>
    struct canCommunicate<Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>,
                          codim> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    /** \brief has conforming level grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct isLevelwiseConforming<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;

      static const bool v = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct isLevelwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct isLeafwiseConforming<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;

      static const bool v = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct isLeafwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct hasBackupRestoreFacilities<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;

      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct hasBackupRestoreFacilities<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct threadSafe<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct threadSafe<Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct viewThreadSafe<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct viewThreadSafe<Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename TrimmerType>
    struct isCartesian<IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, typename TrimmerType>
    struct isCartesian<Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, TrimmerType>>> {
      using ParameterSpaceGrid = typename TrimmerType::ParameterSpaceGrid;
      static const bool v      = false;
    };

  }  // end namespace Capabilities

}  // namespace Dune
namespace Dune {
  template <class Grid>
  struct EnableBoundarySegmentIndexCheck;

}
template <class Grid>
struct EnableLevelIntersectionIteratorCheck;

template <int dim, int dimworld, typename TrimmerType>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGANEW::PatchGrid<dim, dimworld, TrimmerType>>
    : Dune::EnableBoundarySegmentIndexCheck<
          typename Dune::IGANEW::PatchGrid<dim, dimworld, TrimmerType>::ParameterSpaceGrid> {};

template <int dim, int dimworld, typename TrimmerType>
struct EnableLevelIntersectionIteratorCheck<Dune::IGANEW::PatchGrid<dim, dimworld, TrimmerType>> {
  static const bool v = true;
};
