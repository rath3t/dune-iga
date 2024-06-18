// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/grid/common/gridenums.hh>

/** \file
 * @brief The PatchGridLeafIterator class
 */

namespace Dune::IGA::IdentityParameterSpace {

/** @brief Iterator over all entities of a given codimension and level of a grid.
 *  @ingroup PatchGrid
 */
template <int codim, PartitionIteratorType pitype, class GridImp>
class PatchGridLeafIterator
{
private:
  // LevelIterator to the equivalent entity in the host grid
  using IteratorImpl = typename GridImp::ParameterSpace::template ParameterSpaceLeafIterator<codim, pitype>;

public:
  constexpr static int codimension = codim;

  typedef typename GridImp::template Codim<codim>::Entity Entity;
  PatchGridLeafIterator() = default;
  // @todo Please doc me !
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

  // prefix increment
  void increment() {
    ++hostLeafIterator_;
  }

  // dereferencing
  Entity dereference() const {
    return Entity{
        {patchGrid_, *hostLeafIterator_}
    };
  }

  // equality
  bool equals(const PatchGridLeafIterator& i) const {
    return hostLeafIterator_ == i.hostLeafIterator_;
  }

private:
  const GridImp* patchGrid_;

  IteratorImpl hostLeafIterator_;
};

} // namespace Dune::IGA::IdentityParameterSpace
