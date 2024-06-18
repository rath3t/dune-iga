// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

/** \file
 * @brief The PatchGridHierarchicIterator class
 */

namespace Dune::IGA::IdentityParameterSpace {

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

  // the default Constructor
  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostHierarchicIterator_(startEntity.impl().getLocalEntity().hbegin(maxLevel)) {}

  // TODO Please doc me !
  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel,
                                       [[maybe_unused]] bool endDummy)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostHierarchicIterator_(startEntity.impl().getLocalEntity().hend(maxLevel)) {}

  // TODO Please doc me !
  void increment() {
    ++hostHierarchicIterator_;
  }

  // dereferencing
  Entity dereference() const {
    return Entity{
        {parameterSpaceGrid_, *hostHierarchicIterator_}
    };
  }

  // equality
  bool equals(const PatchGridHierarchicIterator& i) const {
    return hostHierarchicIterator_ == i.hostHierarchicIterator_;
  }

private:
  const GridImp* parameterSpaceGrid_;

  HostGridHierarchicIterator hostHierarchicIterator_;
};

} // namespace Dune::IGA::IdentityParameterSpace
