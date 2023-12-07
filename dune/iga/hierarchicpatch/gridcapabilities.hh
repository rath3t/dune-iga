
#pragma once

#include "patchgridfwd.hh"

#include <dune/grid/common/capabilities.hh>
namespace Dune {

  namespace Capabilities {
    /** @brief has entities for some codimensions as host grid
     * \ingroup PatchGrid
     */
    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<const IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using GridFamilyTraitsT =  GridFamilyTraits<dim, dimworld,ScalarType>;
      static const bool v      = GridFamilyTraitsT::template hasEntity<codim>;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using GridFamilyTraitsT =  GridFamilyTraits<dim, dimworld,ScalarType>;
      static const bool v      = GridFamilyTraitsT::template hasEntity<codim>;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<
        const Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>,
        codim> {
      using GridFamilyTraitsT =  GridFamilyTraits<dim, dimworld,ScalarType>;
      static const bool v      = GridFamilyTraitsT::template hasEntity<codim>;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntityIterator<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      using GridFamilyTraitsT =  GridFamilyTraits<dim, dimworld,ScalarType>;
      static const bool v      = GridFamilyTraitsT::template hasEntityIterator<codim>;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntityIterator<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>, codim> {
      using GridFamilyTraitsT =  GridFamilyTraits<dim, dimworld,ScalarType>;
      static const bool v      = GridFamilyTraitsT::template hasEntityIterator<codim>;
    };

    /** @brief PatchGrid can communicate when the host grid can communicate
     *  \ingroup PatchGrid
     */
    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct canCommunicate<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct canCommunicate<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    /** @brief has conforming level grids when host grid has
     * \ingroup PatchGrid
     */
    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLevelwiseConforming<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLevelwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    /** @brief has conforming leaf grids when host grid has
     * \ingroup PatchGrid
     */
    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLeafwiseConforming<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLeafwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct hasBackupRestoreFacilities<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct hasBackupRestoreFacilities<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct threadSafe<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct threadSafe<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct viewThreadSafe<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct viewThreadSafe<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isCartesian<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isCartesian<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
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

template <int dim, int dimworld, template <int,int, typename> typename GridFamilyTraits, typename ScalarType>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>>
    : Dune::EnableBoundarySegmentIndexCheck<
          typename Dune::IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>::ParameterSpaceGrid> {};

template <int dim, int dimworld, template <int,int, typename> typename GridFamilyTraits, typename ScalarType>
struct EnableLevelIntersectionIteratorCheck<Dune::IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
  static const bool v = true;
};
