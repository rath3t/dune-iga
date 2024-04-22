// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <variant>

namespace Dune::IGANEW::DefaultTrim {

template <class GridImp>
class PatchGridGlobalIdSet : public IdSet<GridImp, PatchGridGlobalIdSet<GridImp>,
                                          typename std::remove_const<GridImp>::type::Traits::GlobalIdSet::IdType>
{
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

  using Trimmer = typename GridImp::Trimmer;
  // friend class GridImp::TrimmerType;

  // using TrimmingCurve= typename Trimmer::TrimmingCurve;
  using UntrimmedParameterSpaceGrid = typename Trimmer::UntrimmedParameterSpaceGrid;

public:
  //! constructor stores reference to a grid
  PatchGridGlobalIdSet() = default;
  PatchGridGlobalIdSet(const GridImp& g)
      : grid_(&g) {
  }
  // PatchGridGlobalIdSet(const GridImp& g, const std::vector<TrimmingCurve>& trimmingCurves) : grid_(&g) {}

  //! define the type used for persistent indices
  using IdType              = typename Trimmer::TrimmerTraits::GlobalIdSetId;
  using PersistentIndexType = typename Trimmer::TrimmerTraits::PersistentIndexType;

  //! get id of an entity
  /*
     We use the remove_const to extract the Type from the mutable class,
     because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return e.impl().getHostEntity().id();
  }

  //! get id of subEntity
  /*

   */
  IdType subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i,
               int codim) const {
    // @todo Trim, the sub indeces are wrong!!!
    //  Return sub id of the host entity
    return e.impl().getHostEntity().subId(i, codim);
  }

  /** @todo Should be private */
  void update() {
  }

  PersistentIndexType newFreeIndex() {
    lastFreeIndex_ += 1;
    return lastFreeIndex_;
  }
  PersistentIndexType lastFreeIndex_{};

  const GridImp* grid_;
};
} // namespace Dune::IGANEW::DefaultTrim
