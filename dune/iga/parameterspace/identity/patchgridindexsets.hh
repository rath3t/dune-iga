// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/** \file
 * @brief The index and id sets for the PatchGrid class
 */

#include <vector>

#include <dune/grid/common/indexidset.hh>

namespace Dune::IGA::IdentityTrim {

template <class GridImp>
class PatchGridLevelIndexSet
    : public IndexSet<GridImp, PatchGridLevelIndexSet<GridImp>,
                      typename std::remove_const<GridImp>::type::ParameterSpaceGrid::LevelGridView::IndexSet::IndexType,
                      typename std::remove_const<GridImp>::type::ParameterSpaceGrid::LevelGridView::IndexSet::Types>
{
public:
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid HostGrid;
  typedef typename HostGrid::LevelGridView::IndexSet::Types Types;

  constexpr static int dim = GridImp::dimension;

  // get index of an entity
  template <int codim>
  int index(const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
    return grid_->parameterSpaceGrid().levelIndexSet(level_).template index<codim>(e.impl().getLocalEntity());
  }

  // get index of subEntity of a codim 0 entity
  template <int cc>
  int subIndex(const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const {
    return grid_->parameterSpaceGrid().levelIndexSet(level_).subIndex(e.impl().getLocalEntity(), i, codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(int codim) const {
    return grid_->parameterSpaceGrid().levelIndexSet(level_).size(codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(GeometryType type) const {
    // @todo Trim, count cube and none types i.e. full and trimmed elements

    return grid_->parameterSpaceGrid().levelIndexSet(level_).size(type);
  }

  /** @brief Deliver all geometry types used in this grid */
  Types types(int codim) const {
    // @todo Trim, this should return none and cube for trimmed geometries

    return grid_->parameterSpaceGrid().levelIndexSet(level_).types(codim);
  }

  /** @brief Return true if the given entity is contained in the index set */
  template <class EntityType>
  bool contains(const EntityType& e) const {
    return grid_->parameterSpaceGrid().levelIndexSet(level_).contains(
        grid_->template getHostEntity<EntityType::codimension>(e));
  }

  /** @brief Set up the index set */
  void update(const GridImp& grid, int level) {
    grid_  = &grid;
    level_ = level;
  }

  GridImp* grid_;

  int level_;
};

template <class GridImp>
class PatchGridLeafIndexSet
    : public IndexSet<GridImp, PatchGridLeafIndexSet<GridImp>,
                      typename std::remove_const<GridImp>::type::ParameterSpaceGrid::LeafGridView::IndexSet::IndexType,
                      typename std::remove_const<GridImp>::type::ParameterSpaceGrid::LeafGridView::IndexSet::Types>
{
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

public:
  typedef typename ParameterSpaceGrid::LevelGridView::IndexSet::Types Types;

  /*
   * We use the remove_const to extract the Type from the mutable class,
   * because the const class is not instantiated yet.
   */
  constexpr static int dim = std::remove_const<GridImp>::type::dimension;

  // constructor stores reference to a grid and level
  explicit PatchGridLeafIndexSet(const GridImp& grid)
      : grid_(&grid) {}

  // get index of an entity
  /*
      We use the RemoveConst to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int codim>
  int index(const typename std::remove_const<GridImp>::type::template Codim<codim>::Entity& e) const {
    return grid_->parameterSpaceGrid().leafIndexSet().template index<codim>(e.impl().getLocalEntity());
  }

  // get index of subEntity of a codim 0 entity
  /*
      We use the RemoveConst to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int cc>
  int subIndex(const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i,
               int codim) const {
    return grid_->parameterSpaceGrid().leafIndexSet().subIndex(e.impl().getLocalEntity(), i, codim);
  }

  // get number of entities of given type
  std::size_t size(GeometryType type) const {
    return grid_->parameterSpaceGrid().leafIndexSet().size(type);
  }

  // get number of entities of given codim
  std::size_t size(int codim) const {
    return grid_->parameterSpaceGrid().leafIndexSet().size(codim);
  }

  /** @brief Deliver all geometry types used in this grid */
  Types types(int codim) const {
    return grid_->parameterSpaceGrid().leafIndexSet().types(codim);
  }

  /** @brief Return true if the given entity is contained in the index set */
  template <class EntityType>
  bool contains(const EntityType& e) const {
    return grid_->parameterSpaceGrid().leafIndexSet().contains(
        grid_->template getHostEntity<EntityType::codimension>(e));
  }

  void update(const GridImp& grid) {
    grid_ = &grid;
  }

  const GridImp* grid_;
};

} // namespace Dune::IGA::IdentityTrim
