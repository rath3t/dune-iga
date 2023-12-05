
#pragma once

#include "patchgridfwd.hh"

#include <dune/grid/common/capabilities.hh>
namespace Dune {

  namespace Capabilities {
    /** @brief has entities for some codimensions as host grid
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<const IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using Trimmer = typename GridFamilyTraits<dim, dimworld,ScalarType>::Trimmer;
      static const bool v      = Trimmer::template hasEntity<codim>;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using Trimmer = typename GridFamilyTraits<dim, dimworld,ScalarType>::Trimmer;
      static const bool v      = Trimmer::template hasEntity<codim>;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntity<
        const Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>,
        codim> {
      using Trimmer = typename GridFamilyTraits<dim, dimworld,ScalarType>::Trimmer;
      static const bool v      = Trimmer::template hasEntity<codim>;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntityIterator<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      using Trimmer = typename GridFamilyTraits<dim, dimworld,ScalarType>::Trimmer;
      static const bool v      = Trimmer::template hasEntityIterator<codim>;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct hasEntityIterator<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>, codim> {
      using Trimmer = typename GridFamilyTraits<dim, dimworld,ScalarType>::Trimmer;
      static const bool v      = Trimmer::template hasEntityIterator<codim>;
    };

    /** @brief PatchGrid can communicate when the host grid can communicate
     *  \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct canCommunicate<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits,
              int codim>
    struct canCommunicate<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>, codim> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = canCommunicate<ParameterSpaceGrid, codim>::v;
    };

    /** @brief has conforming level grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLevelwiseConforming<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLevelwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = isLevelwiseConforming<ParameterSpaceGrid>::v;
    };

    /** @brief has conforming leaf grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLeafwiseConforming<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isLeafwiseConforming<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = isLeafwiseConforming<ParameterSpaceGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct hasBackupRestoreFacilities<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;

      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct hasBackupRestoreFacilities<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct threadSafe<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct threadSafe<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct viewThreadSafe<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct viewThreadSafe<
        Dune::Grid<dim, dimworld, ScalarType, IGANEW::PatchGridFamily<dim, dimworld, GridFamilyTraits, ScalarType>>> {
      using ParameterSpaceGrid = typename GridFamilyTraits<dim, dimworld,ScalarType>::TrimmerTraits::ParameterSpaceGrid;
      static const bool v      = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
    struct isCartesian<IGANEW::PatchGrid<dim, dimworld, GridFamilyTraits, ScalarType>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, typename ScalarType, template <int,int, typename> typename GridFamilyTraits>
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
