// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/** \file
 * @brief The index and id sets for the PatchGrid class
 */

#include <dune/grid/common/indexidset.hh>

namespace Dune::IGANEW::DefaultTrim {

namespace Impl {
  enum class IndexSetType
  {
    Level,
    Leaf
  };

  template <class GridImp, IndexSetType type_>
  struct PatchGridIndexSetImpl
  {
    using GridFamily = typename GridImp::GridFamily;
    using GeoTypes   = typename GridFamily::GeometryTypes;
    using HostGrid   = typename std::remove_const_t<GridImp>::ParameterSpaceGrid;

    using HostIndexSet =
        std::conditional_t<type_ == IndexSetType::Leaf,
                           typename std::remove_const_t<GridImp>::ParameterSpaceGrid::LeafGridView::IndexSet,
                           typename std::remove_const_t<GridImp>::ParameterSpaceGrid::LevelGridView::IndexSet>;
    constexpr static int dim = GridImp::dimension;

    PatchGridIndexSetImpl() = default;

    template <typename = void>
    requires(type_ == IndexSetType::Level)
    PatchGridIndexSetImpl(const GridImp& grid, int level)
        : grid_(&grid),
          level_(level) {}

    template <typename = void>
    requires(type_ == IndexSetType::Leaf)
    explicit PatchGridIndexSetImpl(const GridImp& grid)
        : grid_(&grid) {}

    template <int codim>
    int index(const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
      return e.impl().getLocalEntity().index();
    }

    template <int cc>
    int subIndex(const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const {
      return e.impl().getLocalEntity().subIndex(i, codim);
    }

    std::size_t size(int codim) const {
      return grid_->trimmer().entityContainer_.size(codim, level());
    }

    std::size_t size(GeometryType type) const {
      return grid_->trimmer().entityContainer_.size(type, level());
    }

    GeoTypes types(int codim) const {
      return grid_->trimmer().entityContainer_.types(codim, level());
    }

    template <class EntityType>
    bool contains(const EntityType& e) const {
      if constexpr (EntityType::codimension == 0)
        return hostIndexSet().contains(grid_->template getHostEntity<EntityType::codimension>(e).getHostEntity());
      else {
        if (not e.impl().isTrimmed())
          return hostIndexSet().contains(grid_->template getHostEntity<EntityType::codimension>(e).getHostEntity());
        return grid_->trimmer().entityContainer_.contains(e, level());
      }
    }

    /** @brief Set up the index set */
    template <typename = void>
    requires(type_ == IndexSetType::Level)
    void update(const GridImp& grid, int level) {
      grid_  = &grid;
      level_ = level;
    }

    template <typename = void>
    requires(type_ == IndexSetType::Leaf)
    void update(const GridImp& grid) {
      grid_ = &grid;
    }

  private:
    const HostIndexSet& hostIndexSet() const {
      if constexpr (type_ == IndexSetType::Leaf)
        return grid_->parameterSpaceGrid().leafIndexSet();
      else
        return grid_->parameterSpaceGrid().levelIndexSet(level_);
    }

    int level() const {
      if constexpr (type_ == IndexSetType::Leaf)
        return grid_->maxLevel();
      else
        return level_;
    }

    GridImp* grid_{};
    int level_{};
  };

} // namespace Impl

template <class GridImp>
class PatchGridLevelIndexSet
    : public IndexSet<GridImp, PatchGridLevelIndexSet<GridImp>,
                      typename std::remove_const_t<GridImp>::ParameterSpaceGrid::LevelGridView::IndexSet::IndexType,
                      typename std::remove_const_t<GridImp>::GridFamily::GeometryTypes>
{
  using Implementation = Impl::PatchGridIndexSetImpl<GridImp, Impl::IndexSetType::Level>;

public:
  using GeoTypes = typename Implementation::GeoTypes;
  using HostGrid = typename Implementation::HostGrid;

  constexpr static int dim = Implementation::dimension;

  // get index of an entity
  template <int codim>
  int index(const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
    return impl_.template index<codim>(e);
  }

  // get index of subEntity of a codim 0 entity
  template <int cc>
  int subIndex(const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const {
    return impl_.template subIndex<cc>(e, i, codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(int codim) const {
    return impl_.size(codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(GeometryType type) const {
    return impl_.size(type);
  }

  /** @brief Deliver all geometry types used in this grid */
  GeoTypes types(int codim) const {
    return impl_.types(codim);
  }

  /** @brief Return true if the given entity is contained in the index set */
  template <class EntityType>
  bool contains(const EntityType& e) const {
    return impl_.contains(e);
  }

  /** @brief Set up the index set */
  void update(const GridImp& grid, int level) {
    impl_.update(grid, level);
  }

private:
  Implementation impl_{};
};

template <class GridImp>
class PatchGridLeafIndexSet
    : public IndexSet<GridImp, PatchGridLeafIndexSet<GridImp>,
                      typename std::remove_const_t<GridImp>::ParameterSpaceGrid::LeafGridView::IndexSet::IndexType,
                      typename std::remove_const_t<GridImp>::GridFamily::GeometryTypes>
{
  using Implementation = Impl::PatchGridIndexSetImpl<GridImp, Impl::IndexSetType::Leaf>;

public:
  using GeoTypes = typename Implementation::GeoTypes;
  using HostGrid = typename Implementation::HostGrid;

  constexpr static int dim = Implementation::dimension;

  // constructor stores reference to a grid and level
  explicit PatchGridLeafIndexSet(const GridImp& grid)
      : impl_(grid) {}

  // get index of an entity
  /*
      We use the RemoveConst to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  // get index of an entity
  template <int codim>
  int index(const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
    return impl_.template index<codim>(e);
  }

  // get index of subEntity of a codim 0 entity
  template <int cc>
  int subIndex(const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const {
    return impl_.template subIndex<cc>(e, i, codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(int codim) const {
    return impl_.size(codim);
  }

  // get number of entities of given codim, type and on this level
  std::size_t size(GeometryType type) const {
    return impl_.size(type);
  }

  /** @brief Deliver all geometry types used in this grid */
  GeoTypes types(int codim) const {
    return impl_.types(codim);
  }

  /** @brief Return true if the given entity is contained in the index set */
  template <class EntityType>
  bool contains(const EntityType& e) const {
    return impl_.contains(e);
  }

  void update(const GridImp& grid) {
    impl_.update(grid);
  }

private:
  Implementation impl_{};
};

template <class GridImp>
class PatchGridLocalIdSet
    : public IdSet<GridImp, PatchGridLocalIdSet<GridImp>,
                   typename std::remove_const_t<GridImp>::ParameterSpaceGrid::Traits::LocalIdSet::IdType>
{
  using ParameterSpaceGrid = typename std::remove_const_t<GridImp>::ParameterSpaceGrid;

public:
  // define the type used for persistent local ids
  using IdType = typename ParameterSpaceGrid::Traits::LocalIdSet::IdType;

  // constructor stores reference to a grid
  PatchGridLocalIdSet(const GridImp& g)
      : grid_(&g) {}

  // get id of an entity
  /*
      We use the remove_const to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const_t<GridImp>::Traits::template Codim<cd>::Entity& e) const {
    return grid_->parameterSpaceGrid().globalIdSet().id(e.impl().getHostEntity().getHostEntity());
  }

  // get id of subEntity
  /*
   * We use the remove_const to extract the Type from the mutable class,
   * because the const class is not instantiated yet.
   */
  IdType subId(const typename std::remove_const_t<GridImp>::template Codim<0>::Entity& e, int i, int codim) const {
    return grid_->parameterSpaceGrid().localIdSet().subId(e.impl().getHostEntity().getHostEntity(), i, codim);
  }

private:
  void update() {}

  const GridImp* grid_;
};

} // namespace Dune::IGANEW::DefaultTrim
