// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLevelIterator class
 */

namespace Dune::IGANEW::IdentityTrim {

/** @brief Iterator over all entities of a given codimension and level of a grid.
 * @ingroup PatchGrid
 */
template <int codim, PartitionIteratorType pitype, class GridImp>
class PatchGridLevelIterator
{
  typedef typename GridImp::ParameterSpaceGrid::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator
      HostGridLevelIterator;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;

  // Constructor
  PatchGridLevelIterator() = default;
  explicit PatchGridLevelIterator(const GridImp* patchGrid, int level)
      : patchGrid_(patchGrid),
        hostLevelIterator_(patchGrid->parameterSpaceGrid().levelGridView(level).template begin<codim, pitype>()) {}

  /** @brief Constructor which create the end iterator
      @param endDummy      Here only to distinguish it from the other constructor
      @param patchGrid  pointer to PatchGrid instance
      @param level         grid level on which the iterator shall be created
   */
  explicit PatchGridLevelIterator(const GridImp* patchGrid, int level, [[maybe_unused]] bool endDummy)
      : patchGrid_(patchGrid),
        hostLevelIterator_(patchGrid->parameterSpaceGrid().levelGridView(level).template end<codim, pitype>()) {}

  // prefix increment
  void increment() {
    ++hostLevelIterator_;
  }

  // dereferencing
  Entity dereference() const {
    return Entity{
        {patchGrid_, *hostLevelIterator_}
    };
  }

  // equality
  bool equals(const PatchGridLevelIterator& i) const {
    return hostLevelIterator_ == i.hostLevelIterator_;
  }

private:
  const GridImp* patchGrid_;

  HostGridLevelIterator hostLevelIterator_;
};

} // namespace Dune::IGANEW::IdentityTrim
