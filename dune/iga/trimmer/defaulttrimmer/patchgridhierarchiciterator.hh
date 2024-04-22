// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The PatchGridHierarchicIterator class
 */

namespace Dune::IGANEW::DefaultTrim {

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
  typedef typename GridImp::Trimmer::template Codim<0>::ParameterSpaceGridEntity ParameterSpaceGridEntity;
  using Trimmer = typename GridImp::Trimmer;

public:
  constexpr static int codimension = 0;

  typedef typename GridImp::template Codim<0>::Entity Entity;

  //! the default Constructor
  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel)
      : parameterSpaceGrid_(parameterSpaceGrid),
        maxLevel_{maxLevel} // , hostHierarchicIterator_(startEntity.impl().getHostEntity().hbegin(maxLevel))
  {
    // extract the implementation of the grid entity and the parameter space host entity
    stackChildren(&startEntity.impl().getHostEntity());
    setCurrentEntity();
  }

  explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel,
                                       [[maybe_unused]] bool endDummy)
      : parameterSpaceGrid_(parameterSpaceGrid),
        maxLevel_{maxLevel}
  // ,          hostHierarchicIterator_(startEntity.impl().getHostEntity().hend(maxLevel)),maxLevel_{maxLevel}
  {
    // sets current entity to nullptr
    setCurrentEntity();
  }

  void increment() {
    // exit if no further descendants exist
    if (parameterSpaceElementStack_.empty())
      return;

    auto target = parameterSpaceElementStack_.top();
    parameterSpaceElementStack_.pop(); // remove current son
    stackChildren(target);             // add descendants of current son

    setCurrentEntity(); // since std::stack is LIFO, we set the current entity to the first son of the old son

    // ++hostHierarchicIterator_;
    // ++descendantLocalIndex_;
  }

  //! dereferencing
  Entity dereference() const {
    auto realEntity = typename Entity::Implementation{parameterSpaceGrid_, *currentEntityPtr_};
    return Entity{std::move(realEntity)};
  }

  //! equality
  bool equals(const PatchGridHierarchicIterator& i) const {
    // if the iterators point to the different entities they are not equal,
    // But if we are nullptr by construction (end iterator) and the incremented iterator also becomes a nullptr due to
    // setCurrentEntity(), we return true and this terminates the iteration
    return currentEntityPtr_ == i.currentEntityPtr_;
  }

private:
  void stackChildren(const ParameterSpaceGridEntity* target) {
    // Load sons of target onto the iterator stack
    // if the given entity is leaf or max level we do not add anything to the stack
    if (target->level() < maxLevel_ && !target->isLeaf())
      for (auto descendantId : target->entityInfo_.decendantIds)
        parameterSpaceElementStack_.push(
            &parameterSpaceGrid_->trimmer().entityContainer_.template entity<0>(descendantId, target->level() + 1));
  }

  void setCurrentEntity() {
    currentEntityPtr_ = parameterSpaceElementStack_.empty() ? nullptr : parameterSpaceElementStack_.top();
  }

  std::stack<const ParameterSpaceGridEntity*> parameterSpaceElementStack_;
  const ParameterSpaceGridEntity* currentEntityPtr_;
  const GridImp* parameterSpaceGrid_;
  // HostGridHierarchicIterator hostHierarchicIterator_;
  int maxLevel_;
};

} // namespace Dune::IGANEW::DefaultTrim
