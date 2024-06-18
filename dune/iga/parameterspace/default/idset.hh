// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <variant>

namespace Dune::IGA::DefaultParameterSpace {

template <class GridImp>
class PatchGridGlobalIdSet : public IdSet<GridImp, PatchGridGlobalIdSet<GridImp>,
                                          typename std::remove_const_t<GridImp>::Traits::GlobalIdSet::IdType>
{
  typedef typename std::remove_const_t<GridImp>::ParameterSpaceGrid ParameterSpaceGrid;

  using ParameterSpace = typename GridImp::ParameterSpace;
  // friend class GridImp::ParameterSpaceType;

  // using TrimmingCurve= typename ParameterSpace::TrimmingCurve;
  using UntrimmedParameterSpaceGrid = typename ParameterSpace::UntrimmedParameterSpaceGrid;

public:
  // constructor stores reference to a grid
  PatchGridGlobalIdSet() = default;
  explicit PatchGridGlobalIdSet(const GridImp& g)
      : grid_(&g) {}
  // PatchGridGlobalIdSet(const GridImp& g, const std::vector<TrimmingCurve>& trimmingCurves) : grid_(&g) {}

  // define the type used for persistent indices
  using IdType              = typename ParameterSpace::ParameterSpaceTraits::GlobalIdSetId;
  using PersistentIndexType = typename ParameterSpace::ParameterSpaceTraits::PersistentIndexType;

  // get id of an entity
  /*
     We use the remove_const to extract the Type from the mutable class,
     because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const_t<GridImp>::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return e.impl().getLocalEntity().id();
  }

  // get id of subEntity

  IdType subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i,
               int codim) const {
    return e.impl().getLocalEntity().subId(i, codim);
  }

  /** @todo Should be private */
  void update() {}

  PersistentIndexType newFreeIndex() {
    lastFreeIndex_ += 1;
    return lastFreeIndex_;
  }
  PersistentIndexType lastFreeIndex_{};

  const GridImp* grid_;
};
} // namespace Dune::IGA::DefaultParameterSpace
