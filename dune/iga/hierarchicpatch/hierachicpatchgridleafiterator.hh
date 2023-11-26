// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The PatchGridLeafIterator class
 */

namespace Dune::IGANEW {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup PatchGrid
   */
  template <int codim, PartitionIteratorType pitype, class GridImp>
  class PatchGridLeafIterator {
   private:
    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LeafIterator
        HostGridLeafIterator;

   public:
    constexpr static int codimension = codim;

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    PatchGridLeafIterator() = default;
    //! \todo Please doc me !
    explicit PatchGridLeafIterator(const GridImp* identityGrid)
        : identityGrid_(identityGrid),
          hostLeafIterator_(identityGrid->hostgrid_->leafGridView().template begin<codim, pitype>()) {}

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param identityGrid  pointer to grid instance
     */
    explicit PatchGridLeafIterator(const GridImp* identityGrid, [[maybe_unused]] bool endDummy)
        : identityGrid_(identityGrid),
          hostLeafIterator_(identityGrid->hostgrid_->leafGridView().template end<codim, pitype>()) {}

    //! prefix increment
    void increment() { ++hostLeafIterator_; }

    //! dereferencing
    Entity dereference() const { return Entity{{identityGrid_, *hostLeafIterator_}}; }

    //! equality
    bool equals(const PatchGridLeafIterator& i) const { return hostLeafIterator_ == i.hostLeafIterator_; }

   private:
    const GridImp* identityGrid_;

    HostGridLeafIterator hostLeafIterator_;
  };

}  // namespace Dune::IGANEW
