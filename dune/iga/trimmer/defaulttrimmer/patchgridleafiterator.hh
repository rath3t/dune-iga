// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLeafIterator class
 */

namespace Dune::IGANEW::DefaultTrim {

  /** @brief Iterator over all entities of a given codimension and level of a grid.
   *  @ingroup PatchGrid
   */
  template <int codim, PartitionIteratorType pitype, class GridImp>
  class PatchGridLeafIterator {
   private:
    // LevelIterator to the equivalent entity in the host grid
    using IteratorImpl=typename GridImp::Trimmer::template ParameterSpaceLeafIterator<codim,pitype>;
    typedef typename GridImp::Trimmer::template Codim<codim>::ParameterSpaceGridEntity ParameterSpaceGridEntity;

   public:
    constexpr static int codimension = codim;

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    PatchGridLeafIterator() = default;
    //! @todo Please doc me !
    explicit PatchGridLeafIterator(const GridImp* patchGrid)
        : patchGrid_(patchGrid),
          hostLeafIterator_(patchGrid->parameterSpaceGrid().leafGridView().template begin<codim, pitype>()) {}

    /** @brief Constructor which create the end iterator
     *  @param endDummy      Here only to distinguish it from the other constructor
     *  @param patchGrid  pointer to grid instance
     */
    explicit PatchGridLeafIterator(const GridImp* patchGrid, [[maybe_unused]] bool endDummy)
        : patchGrid_(patchGrid),
          hostLeafIterator_(patchGrid->parameterSpaceGrid().leafGridView().template end<codim, pitype>()) {}

    //! prefix increment
    void increment() { ++hostLeafIterator_; }

    //! dereferencing
    Entity dereference() const {

      if constexpr (codim==0) {
        auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_, *hostLeafIterator_,id_};
        auto realEntity= typename Entity::Implementation{patchGrid_,std::move(parameterSpaceEntity)};
        return Entity{std::move(realEntity)};
      }
      else {
        auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_,id_};
        auto realEntity= typename Entity::Implementation{patchGrid_,std::move(parameterSpaceEntity)};
        DUNE_THROW(NotImplemented,"This is doing the wrong thing");

        return Entity{std::move(realEntity)};
      }
    }

    //! equality
    bool equals(const PatchGridLeafIterator& i) const { return hostLeafIterator_ == i.hostLeafIterator_; }

   private:
    const GridImp* patchGrid_;
    typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId id_;
    IteratorImpl hostLeafIterator_;
  };

}  // namespace Dune::IGANEW
