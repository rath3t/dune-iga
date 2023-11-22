
#pragma once

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/gridcheck.hh>
namespace Dune {

  namespace IGANEW {
    template <int dim, int dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    class PatchGrid;
  }

  namespace Capabilities {
    /** \brief has entities for some codimensions as host grid
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid, int codim>
    struct hasEntity<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>, codim> {
      static const bool v = hasEntity<HostGrid, codim>::v;
    };

    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid, int codim>
    struct hasEntityIterator<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>, codim> {
      static const bool v = hasEntityIterator<HostGrid, codim>::v;
    };

    /** \brief PatchGrid can communicate when the host grid can communicate
     *  \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid, int codim>
    struct canCommunicate<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>, codim> {
      static const bool v = canCommunicate<HostGrid, codim>::v;
    };

    /** \brief has conforming level grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct isLevelwiseConforming<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = isLevelwiseConforming<HostGrid>::v;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup PatchGrid
     */
    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct isLeafwiseConforming<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = isLeafwiseConforming<HostGrid>::v;
    };

    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct hasBackupRestoreFacilities<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct threadSafe<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct viewThreadSafe<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = false;
    };

    template <std::size_t dim, std::size_t dimworld, Trimming trim, typename ScalarType, typename HostGrid>
    struct isCartesian<IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
      static const bool v = false;
    };

  }  // end namespace Capabilities

}  // namespace Dune
template <int dim, int dimworld, Trimming trim, typename ScalarType, typename HostGrid>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>>
    : public std::true_type {};

template <int dim, int dimworld, Trimming trim, typename ScalarType, typename HostGrid>
struct EnableLevelIntersectionIteratorCheck<Dune::IGANEW::PatchGrid<dim, dimworld, trim, ScalarType, HostGrid>> {
  static const bool v = true;
};
