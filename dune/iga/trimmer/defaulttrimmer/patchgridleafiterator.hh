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
class PatchGridLeafIterator
{
private:
  // LevelIterator to the equivalent entity in the host grid
  using IteratorImplR = typename GridImp::Trimmer::template ParameterSpaceLeafIterator<codim, pitype>;
  using OLDIteratorImpl =
      typename GridImp::ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LeafIterator;
  using IteratorImpl = IteratorImplR; // std::conditional_t<codim==0,IteratorImplR,OLDIteratorImpl>;
  typedef typename GridImp::Trimmer::template Codim<codim>::ParameterSpaceGridEntity ParameterSpaceGridEntity;
  using ElementTrimData = typename GridImp::Trimmer::ElementTrimData;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;
  PatchGridLeafIterator() = default;
  //! @todo Please doc me !
  // template<typename =void> requires (codim!=0)
  // explicit PatchGridLeafIterator(const GridImp* patchGrid)
  //     : patchGrid_(patchGrid),
  //       parameterSpaceLeafIterator(patchGrid->parameterSpaceGrid().leafGridView().template begin<codim, pitype>())
  //       {}

  // template<typename =void> requires (codim==0)
  explicit PatchGridLeafIterator(const GridImp* patchGrid)
      : patchGrid_(patchGrid),
        parameterSpaceLeafIterator(
            patchGrid_->trimmer().entityContainer_.template begin<codim>(patchGrid_->maxLevel())) {
  }

  /** @brief Constructor which create the end iterator
   *  @param endDummy      Here only to distinguish it from the other constructor
   *  @param patchGrid  pointer to grid instance
   */
  // template<typename =void> requires (codim!=0)
  // explicit PatchGridLeafIterator(const GridImp* patchGrid, [[maybe_unused]] bool endDummy)
  //     : patchGrid_(patchGrid),
  //       parameterSpaceLeafIterator(patchGrid->parameterSpaceGrid().leafGridView().template end<codim, pitype>()) {}

  // template<typename =void> requires (codim==0)
  explicit PatchGridLeafIterator(const GridImp* patchGrid, [[maybe_unused]] bool endDummy)
      : patchGrid_(patchGrid),
        parameterSpaceLeafIterator(patchGrid_->trimmer().entityContainer_.template end<codim>(patchGrid_->maxLevel())) {
  }

  //! prefix increment
  void increment() {
    ++parameterSpaceLeafIterator;
  }

  using GlobalIdSetId = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
  //! dereferencing
  Entity dereference() const {
    if constexpr (codim == 0) {
      // auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_, *parameterSpaceLeafIterator,id_};
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLeafIterator};
      return Entity{std::move(realEntity)};
    } else if (not parameterSpaceLeafIterator->stemsFromTrim()) { // subentity is untrimmed

      // auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_,*parameterSpaceLeafIterator};
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLeafIterator};

      return Entity{std::move(realEntity)};
    } else {
      DUNE_THROW(NotImplemented, "This is doing the wrong thing");
      // auto parameterSpaceEntity=
      // ParameterSpaceGridEntity{patchGrid_,id_,ElementTrimData(),parameterSpaceLeafIterator->level()}; auto
      // realEntity= typename Entity::Implementation{patchGrid_,std::move(parameterSpaceEntity)};
      //
      // return Entity{std::move(realEntity)};
    }
  }

  //! equality
  bool equals(const PatchGridLeafIterator& i) const {
    return parameterSpaceLeafIterator == i.parameterSpaceLeafIterator;
  }

private:
  const GridImp* patchGrid_;
  GlobalIdSetId id_;
  IteratorImpl parameterSpaceLeafIterator;
};

} // namespace Dune::IGANEW::DefaultTrim
