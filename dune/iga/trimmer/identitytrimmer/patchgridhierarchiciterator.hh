// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The PatchGridHierarchicIterator class
 */

namespace Dune::IGANEW::IdentityTrim {

//**********************************************************************
//
/** @brief Iterator over the descendants of an entity.
 * @ingroup PatchGrid
   Mesh entities of codimension 0 ("elements") allow to visit all entities of
   codimension 0 obtained through nested, hierarchic refinement of the entity.
   Iteration over this set of entities is provided by the HierarchicIterator,
   starting from a given entity.
 */
template <class GridImp>
class PatchGridHierarchicIterator
{
  // Type of the corresponding HierarchicIterator in the host grid
  typedef
      typename GridImp::ParameterSpaceGrid::template Codim<0>::Entity::HierarchicIterator HostGridHierarchicIterator;

public:
  constexpr static int codimension = 0;

  typedef typename GridImp::template Codim<0>::Entity Entity;

  //! the default Constructor
  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostHierarchicIterator_(startEntity.impl().getHostEntity().hbegin(maxLevel)) {
  }

  //! @todo Please doc me !
  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel,
                                       [[maybe_unused]] bool endDummy)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostHierarchicIterator_(startEntity.impl().getHostEntity().hend(maxLevel)) {
  }

  //! @todo Please doc me !
  void increment() {
    ++hostHierarchicIterator_;
  }

  //! dereferencing
  Entity dereference() const {
    return Entity{
        {parameterSpaceGrid_, *hostHierarchicIterator_}
    };
  }

  //! equality
  bool equals(const PatchGridHierarchicIterator& i) const {
    return hostHierarchicIterator_ == i.hostHierarchicIterator_;
  }

private:
  const GridImp* parameterSpaceGrid_;

  HostGridHierarchicIterator hostHierarchicIterator_;
};

} // namespace Dune::IGANEW::IdentityTrim
