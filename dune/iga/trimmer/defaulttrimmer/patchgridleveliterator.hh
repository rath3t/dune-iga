// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLevelIterator class
 */

namespace Dune::IGANEW::DefaultTrim {

/** @brief Iterator over all entities of a given codimension and level of a grid.
 * @ingroup PatchGrid
 */
template <int codim, PartitionIteratorType pitype, class GridImp>
class PatchGridLevelIterator
{
  using IteratorImplR = typename GridImp::Trimmer::template ParameterSpaceLevelIterator<codim, pitype>;
  using OLDIteratorImpl =
      typename GridImp::ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LevelIterator;
  using IteratorImpl = IteratorImplR; // std::conditional_t<codim==0,IteratorImplR,OLDIteratorImpl>;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;
  typedef typename GridImp::Trimmer::template Codim<codim>::ParameterSpaceGridEntity ParameterSpaceGridEntity;

  //! Constructor
  PatchGridLevelIterator() = default;
  // template<typename =void> requires (codim!=0)
  // explicit PatchGridLevelIterator(const GridImp* patchGrid, int level)
  //     : patchGrid_(patchGrid),
  //       parameterSpaceLevelIterator(patchGrid->parameterSpaceGrid().levelGridView(level).template begin<codim,
  //       pitype>()) {}

  // template<typename =void> requires (codim==0)
  explicit PatchGridLevelIterator(const GridImp* patchGrid, int level)
      : patchGrid_(patchGrid),
        parameterSpaceLevelIterator(patchGrid_->trimmer().entityContainer_.template begin<codim>(level)) {
  }

  /** @brief Constructor which create the end iterator
      @param endDummy      Here only to distinguish it from the other constructor
      @param patchGrid  pointer to PatchGrid instance
      @param level         grid level on which the iterator shall be created
   */
  // template<typename =void> requires (codim!=0)
  // explicit PatchGridLevelIterator(const GridImp* patchGrid, int level, [[maybe_unused]] bool endDummy)
  //     : patchGrid_(patchGrid),
  //       parameterSpaceLevelIterator(patchGrid->parameterSpaceGrid().levelGridView(level).template end<codim,
  //       pitype>()) {}

  // template<typename =void> requires (codim==0)
  explicit PatchGridLevelIterator(const GridImp* patchGrid, int level, [[maybe_unused]] bool endDummy)
      : patchGrid_(patchGrid),
        parameterSpaceLevelIterator(patchGrid_->trimmer().entityContainer_.template end<codim>(level)) {
  }

  //! prefix increment
  void increment() {
    ++parameterSpaceLevelIterator;
  }
  using GlobalIdSetId   = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
  using ElementTrimData = typename GridImp::Trimmer::ElementTrimData;

  //! dereferencing
  Entity dereference() const {
    if constexpr (codim == 0) {
      // auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_, *parameterSpaceLevelIterator,id_};
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLevelIterator};
      return Entity{std::move(realEntity)};
    } else /* if(id_.elementState==GlobalIdSetId::ElementState::full)// subentity is untrimmed */ {
      // auto parameterSpaceEntity= ParameterSpaceGridEntity{patchGrid_,*parameterSpaceLevelIterator,id_};
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLevelIterator};

      return Entity{std::move(realEntity)};
    }
    // else {
    //   DUNE_THROW(NotImplemented,"This is doing the wrong thing");
    //   // auto parameterSpaceEntity=
    //   ParameterSpaceGridEntity{patchGrid_,id_,ElementTrimData(),parameterSpaceLevelIterator->level()};
    //   // auto realEntity= typename Entity::Implementation{patchGrid_,std::move(parameterSpaceEntity)};
    //   //
    //   // return Entity{std::move(realEntity)};
    // }
  }

  //! equality
  bool equals(const PatchGridLevelIterator& i) const {
    return parameterSpaceLevelIterator == i.parameterSpaceLevelIterator;
  }

private:
  const GridImp* patchGrid_;
  typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId id_;

  IteratorImpl parameterSpaceLevelIterator;
};

} // namespace Dune::IGANEW::DefaultTrim
