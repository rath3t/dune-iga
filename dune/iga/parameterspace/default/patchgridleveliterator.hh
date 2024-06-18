// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLevelIterator class
 */

namespace Dune::IGA::DefaultParameterSpace {

/** @brief Iterator over all entities of a given codimension and level of a grid.
 * @ingroup PatchGrid
 */
template <int codim, PartitionIteratorType pitype, class GridImp>
class PatchGridLevelIterator
{
  using IteratorImplR = typename GridImp::ParameterSpace::template ParameterSpaceLevelIterator<codim, pitype>;
  using OLDIteratorImpl =
      typename GridImp::ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LevelIterator;
  using IteratorImpl = IteratorImplR; // std::conditional_t<codim==0,IteratorImplR,OLDIteratorImpl>;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;
  typedef typename GridImp::ParameterSpace::template Codim<codim>::ParameterSpaceGridEntity ParameterSpaceGridEntity;

  // Constructor
  PatchGridLevelIterator() = default;

  explicit PatchGridLevelIterator(const GridImp* patchGrid, int level)
      : patchGrid_(patchGrid),
        parameterSpaceLevelIterator(patchGrid_->parameterSpace().entityContainer_.template begin<codim>(level)) {}

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
        parameterSpaceLevelIterator(patchGrid_->parameterSpace().entityContainer_.template end<codim>(level)) {}

  // prefix increment
  void increment() {
    ++parameterSpaceLevelIterator;
  }
  using GlobalIdSetId   = typename GridImp::GridFamily::ParameterSpaceTraits::GlobalIdSetId;
  using ElementTrimData = typename GridImp::ParameterSpace::ElementTrimData;

  // dereferencing
  Entity dereference() const {
    auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLevelIterator};
    return Entity{std::move(realEntity)};
  }

  // equality
  bool equals(const PatchGridLevelIterator& i) const {
    return parameterSpaceLevelIterator == i.parameterSpaceLevelIterator;
  }

private:
  const GridImp* patchGrid_;
  typename GridImp::GridFamily::ParameterSpaceTraits::GlobalIdSetId id_;

  IteratorImpl parameterSpaceLevelIterator;
};

} // namespace Dune::IGA::DefaultParameterSpace
