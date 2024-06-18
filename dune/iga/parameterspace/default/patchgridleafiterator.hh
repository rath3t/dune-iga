// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLeafIterator class
 */

namespace Dune::IGA::DefaultParameterSpace {

/** @brief Iterator over all entities of a given codimension and level of a grid.
 *  @ingroup PatchGrid
 */
template <int codim, PartitionIteratorType pitype, class GridImp>
class PatchGridLeafIterator
{
  using IteratorImplR = typename GridImp::ParameterSpace::template ParameterSpaceLeafIterator<codim, pitype>;
  using OLDIteratorImpl =
      typename GridImp::ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LeafIterator;
  using IteratorImpl = IteratorImplR; // std::conditional_t<codim==0,IteratorImplR,OLDIteratorImpl>;
  typedef typename GridImp::ParameterSpace::template Codim<codim>::ParameterSpaceGridEntity ParameterSpaceGridEntity;
  using ElementTrimData = typename GridImp::ParameterSpace::ElementTrimData;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;
  PatchGridLeafIterator() = default;

  explicit PatchGridLeafIterator(const GridImp* patchGrid)
      : patchGrid_(patchGrid),
        parameterSpaceLeafIterator(
            patchGrid_->parameterSpace().entityContainer_.template begin<codim>(patchGrid_->maxLevel())) {}

  explicit PatchGridLeafIterator(const GridImp* patchGrid, [[maybe_unused]] bool endDummy)
      : patchGrid_(patchGrid),
        parameterSpaceLeafIterator(
            patchGrid_->parameterSpace().entityContainer_.template end<codim>(patchGrid_->maxLevel())) {}

  // prefix increment
  void increment() {
    ++parameterSpaceLeafIterator;
  }

  using GlobalIdSetId = typename GridImp::GridFamily::ParameterSpaceTraits::GlobalIdSetId;
  // dereferencing
  Entity dereference() const {
    if constexpr (codim == 0) {
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLeafIterator};
      return Entity{std::move(realEntity)};
    } else {
      auto realEntity = typename Entity::Implementation{patchGrid_, *parameterSpaceLeafIterator};
      return Entity{std::move(realEntity)};
    }
  }

  // equality
  bool equals(const PatchGridLeafIterator& i) const {
    return parameterSpaceLeafIterator == i.parameterSpaceLeafIterator;
  }

private:
  const GridImp* patchGrid_;
  GlobalIdSetId id_;
  IteratorImpl parameterSpaceLeafIterator;
};

} // namespace Dune::IGA::DefaultParameterSpace
