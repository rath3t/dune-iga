// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The index and id sets for the PatchGrid class
 */

#include <vector>

#include <dune/grid/common/indexidset.hh>

namespace Dune::IGANEW::IdentityTrim {

/** @todo Take the index types from the host grid */
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

  //! get index of an entity
  template <int codim>
  int index(const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
    return grid_->parameterSpaceGrid().levelIndexSet(level_).template index<codim>(
        grid_->template getHostEntity<codim>(e));
  }

  //! get index of subEntity of a codim 0 entity
  template <int cc>
  int subIndex(const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const {
    // @todo Trim, the subindeces are wrong!
    return grid_->parameterSpaceGrid().levelIndexSet(level_).subIndex(grid_->template getHostEntity<cc>(e), i, codim);
  }

  //! get number of entities of given codim, type and on this level
  std::size_t size(int codim) const {
    // @todo Trim,coun trimmed elements!

    return grid_->parameterSpaceGrid().levelIndexSet(level_).size(codim);
  }

  //! get number of entities of given codim, type and on this level
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

  //! constructor stores reference to a grid and level
  explicit PatchGridLeafIndexSet(const GridImp& grid)
      : grid_(&grid) {
  }

  //! get index of an entity
  /*
      We use the RemoveConst to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int codim>
  int index(const typename std::remove_const<GridImp>::type::template Codim<codim>::Entity& e) const {
    return grid_->parameterSpaceGrid().leafIndexSet().template index<codim>(grid_->template getHostEntity<codim>(e));
  }

  //! get index of subEntity of a codim 0 entity
  /*
      We use the RemoveConst to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int cc>
  int subIndex(const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i,
               int codim) const {
    // @todo Trim, this is wrong for the trimmed case
    return grid_->parameterSpaceGrid().leafIndexSet().subIndex(grid_->template getHostEntity<cc>(e), i, codim);
  }

  //! get number of entities of given type
  std::size_t size(GeometryType type) const {
    // @todo Trim, count cube and none types i.e. full and trimmed elements
    return grid_->parameterSpaceGrid().leafIndexSet().size(type);
  }

  //! get number of entities of given codim
  std::size_t size(int codim) const {
    return grid_->parameterSpaceGrid().leafIndexSet().size(codim);
  }

  /** @brief Deliver all geometry types used in this grid */
  Types types(int codim) const {
    // @todo Trim, this should provide cube and none!
    return grid_->parameterSpaceGrid().leafIndexSet().types(codim);
  }

  /** @brief Return true if the given entity is contained in the index set */
  template <class EntityType>
  bool contains(const EntityType& e) const {
    return grid_->parameterSpaceGrid().leafIndexSet().contains(
        grid_->template getHostEntity<EntityType::codimension>(e));
  }

  /** @todo Currently we support only vertex and element indices */
  void update(const GridImp& grid) {
    grid_ = &grid;
  }

  const GridImp* grid_;
};

template <class GridImp>
class PatchGridGlobalIdSet
    : public IdSet<GridImp, PatchGridGlobalIdSet<GridImp>,
                   typename std::remove_const<GridImp>::type::ParameterSpaceGrid::Traits::GlobalIdSet::IdType>
{
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

public:
  //! constructor stores reference to a grid
  explicit PatchGridGlobalIdSet(const GridImp& g)
      : grid_(&g) {
  }

  //! define the type used for persistent indices
  typedef typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType IdType;

  //! get id of an entity
  /*
     We use the remove_const to extract the Type from the mutable class,
     because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return grid_->parameterSpaceGrid().globalIdSet().id(e.impl().getHostEntity());
  }

  //! get id of subEntity
  /*

   */
  IdType subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i,
               int codim) const {
    // @todo Trim, the sub indeces are wrong!!!
    //  Return sub id of the host entity
    return grid_->parameterSpaceGrid().globalIdSet().subId(e.impl().getHostEntity(), i, codim);
  }

  /** @todo Should be private */
  void update() {
  }

  const GridImp* grid_;
};

template <class GridImp>
class PatchGridLocalIdSet
    : public IdSet<GridImp, PatchGridLocalIdSet<GridImp>,
                   typename std::remove_const<GridImp>::type::ParameterSpaceGrid::Traits::LocalIdSet::IdType>
{
private:
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

public:
  //! define the type used for persistent local ids
  typedef typename ParameterSpaceGrid::Traits::LocalIdSet::IdType IdType;

  //! constructor stores reference to a grid
  PatchGridLocalIdSet(const GridImp& g)
      : grid_(&g) {
  }

  //! get id of an entity
  /*
      We use the remove_const to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return grid_->parameterSpaceGrid().localIdSet().id(e.impl().getHostEntity());
  }

  //! get id of subEntity
  /*
   * We use the remove_const to extract the Type from the mutable class,
   * because the const class is not instantiated yet.
   */
  IdType subId(const typename std::remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const {
    // Return sub id of the host entity
    /* @todo Trim, the sub indices are wrong!!! */

    return grid_->parameterSpaceGrid().localIdSet().subId(e.impl().getHostEntity(), i, codim);
  }

  /** @todo Should be private */
  void update() {
  }

  const GridImp* grid_;
};

} // namespace Dune::IGANEW::IdentityTrim
