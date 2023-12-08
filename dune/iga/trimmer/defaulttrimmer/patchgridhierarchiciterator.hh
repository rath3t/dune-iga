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
  class PatchGridHierarchicIterator {
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
        : elementInfo_{startEntity.impl().getHostEntity().entityInfo_},
          parameterSpaceGrid_(parameterSpaceGrid),
          hostHierarchicIterator_(startEntity.impl().getHostEntity().hbegin(maxLevel)) {
      assert(startEntity.level() + 1 == maxLevel && "Only direct descendants are implemented");
    }

    //! @todo Please doc me !
    explicit PatchGridHierarchicIterator(const GridImp* parameterSpaceGrid, const Entity& startEntity, int maxLevel,
                                         [[maybe_unused]] bool endDummy)
        : elementInfo_{},
          parameterSpaceGrid_(parameterSpaceGrid),
          hostHierarchicIterator_(startEntity.impl().getHostEntity().hend(maxLevel)) {
      assert(startEntity.level() + 1 == maxLevel && "Only direct descendants are implemented");
    }

    //! @todo Please doc me !
    void increment() {
      ++hostHierarchicIterator_;
      ++descendantLocalIndex_;
    }

    //! dereferencing
    Entity dereference() const {
      auto parameterSpaceEntity
          = ParameterSpaceGridEntity{parameterSpaceGrid_, *hostHierarchicIterator_,
                                     parameterSpaceGrid_->trimmer().entityContainer_.idToElementInfoMap.at(
                                         elementInfo_.decendantIds[descendantLocalIndex_])};
      auto realEntity = typename Entity::Implementation{parameterSpaceGrid_, std::move(parameterSpaceEntity)};
      return Entity{std::move(realEntity)};
    }

    //! equality
    bool equals(const PatchGridHierarchicIterator& i) const {
      return hostHierarchicIterator_ == i.hostHierarchicIterator_;
    }

   private:
    using ElementInfo = typename Trimmer::TrimmerTraits::ElementInfo;
    ElementInfo elementInfo_;
    const GridImp* parameterSpaceGrid_;
    unsigned int descendantLocalIndex_{};
    HostGridHierarchicIterator hostHierarchicIterator_;
  };

}  // namespace Dune::IGANEW::DefaultTrim
