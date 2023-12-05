// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The PatchGridEntity class
 */

#include "patchgridgeometry.hh"

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"

namespace Dune::IGANEW {

  // Forward declarations

  template <int codim, int dim, class GridImp>
  class PatchGridEntity;

  template <int codim, PartitionIteratorType pitype, class GridImp>
  class PatchGridLevelIterator;

  template <class GridImp>
  class PatchGridLevelIntersectionIterator;

  template <class GridImp>
  class PatchGridLeafIntersectionIterator;

  template <class GridImp>
  class PatchGridHierarchicIterator;

  // External forward declarations
  template <class Grid>
  struct HostGridAccess;

  //**********************************************************************
  //
  // --PatchGridEntity
  // --Entity
  //
  /** @brief The implementation of entities in a PatchGrid
   *   @ingroup PatchGrid
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template <int codim, int dim, class GridImp>
  class PatchGridEntity : public EntityDefaultImplementation<codim, dim, GridImp, PatchGridEntity> {
    template <class GridImp_>
    friend class PatchGridLevelIndexSet;

    template <class GridImp_>
    friend class PatchGridLeafIndexSet;

    template <class GridImp_>
    friend class PatchGridLocalIdSet;

    template <class GridImp_>
    friend class PatchGridGlobalIdSet;

    friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;

   private:
    // The codimension of this entity wrt the host grid
    constexpr static int CodimInHostGrid = GridImp::ParameterSpaceGrid::dimension - GridImp::dimension + codim;

    // equivalent entity in the host grid
    //@todo Trimmer should also provide a


   public:
    typedef typename GridImp::ctype ctype;
    using TrimmerType = typename GridImp::Trimmer;

    using ParameterSpaceGridEntity = typename TrimmerType::template ParameterSpaceGridEntity<codim>;
    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

    PatchGridEntity() : patchGrid_(nullptr) {}

    PatchGridEntity(const GridImp* patchGrid, const ParameterSpaceGridEntity& hostEntity)
        : hostEntity_(hostEntity), patchGrid_(patchGrid) {}

    PatchGridEntity(const GridImp* patchGrid, ParameterSpaceGridEntity&& hostEntity)
        : hostEntity_(std::move(hostEntity)), patchGrid_(patchGrid) {}

    //! @todo Please doc me !
    PatchGridEntity(const PatchGridEntity& original)
        : hostEntity_(original.hostEntity_), patchGrid_(original.patchGrid_) {}

    PatchGridEntity(PatchGridEntity&& original) noexcept
        : hostEntity_(std::move(original.hostEntity_)), patchGrid_(original.patchGrid_) {}

    //! @todo Please doc me !
    PatchGridEntity& operator=(const PatchGridEntity& original) {
      if (this != &original) {
        patchGrid_  = original.patchGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    //! @todo Please doc me !
    PatchGridEntity& operator=(PatchGridEntity&& original) noexcept {
      if (this != &original) {
        patchGrid_  = original.patchGrid_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    bool equals(const PatchGridEntity& other) const { return untrimmedHostEntity() == other.untrimmedHostEntity(); }

    //! returns true if father entity exists
    bool hasFather() const { return hostEntity_.hasFather(); }

    //! Create EntitySeed
    EntitySeed seed() const { return EntitySeed(untrimmedHostEntity()); }

    //! level of this element
    int level() const { return hostEntity_.level(); }

    /** @brief The partition type for parallel computing
     */
    PartitionType partitionType() const { return hostEntity_.partitionType(); }

    /** @brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities(unsigned int cc) const { return hostEntity_.subEntities(cc); }

      using ParameterSpaceGeometry = typename TrimmerType::template LocalParameterSpaceGeometry<codim>;


    //! geometry of this entity
    Geometry geometry() const {

      auto geo = typename Geometry::Implementation(
          hostEntity_.geometry(), patchGrid_->patchGeometries_[this->level()].template localView<codim, TrimmerType>());
      return Geometry(geo);
    }




    const auto& untrimmedHostEntity()const {
      if constexpr (requires {hostEntity_.untrimmedHostEntity();})
        return hostEntity_.untrimmedHostEntity();
      else
      return hostEntity_;
    }

    const auto& hostEntity()const {
        return hostEntity_;
    }
  private:

    ParameterSpaceGridEntity hostEntity_;
    const GridImp* patchGrid_;
  };

  //***********************
  //
  //  --PatchGridEntity
  //
  //***********************
  /** @brief Specialization for codim-0-entities.
   * @ingroup PatchGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   */
  template <int dim, class GridImp>
  class PatchGridEntity<0, dim, GridImp> : public EntityDefaultImplementation<0, dim, GridImp, PatchGridEntity> {
    friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;


   public:
    typedef typename GridImp::ctype ctype;
    using TrimmerType = typename GridImp::Trimmer;
    // The codimension of this entitypointer wrt the host grid
    constexpr static int dimension       = GridImp::dimension;
    constexpr static int CodimInHostGrid = GridImp::ParameterSpaceGrid::dimension - dimension;
    constexpr static int dimworld        = GridImp::dimensionworld;

    // equivalent entity in the host grid
    using ParameterSpaceGridEntity = typename TrimmerType::template ParameterSpaceGridEntity<0>;


    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on this level
    typedef PatchGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

    //! The Iterator over intersections on the leaf level
    typedef PatchGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef PatchGridHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    PatchGridEntity() : patchGrid_(nullptr) {}

    PatchGridEntity(const GridImp* patchGrid, const ParameterSpaceGridEntity& hostEntity)
        : hostEntity_(hostEntity), patchGrid_(patchGrid) {}

    PatchGridEntity(const GridImp* patchGrid, ParameterSpaceGridEntity&& hostEntity)
        : hostEntity_(std::move(hostEntity)), patchGrid_(patchGrid) {}

    //! @todo Please doc me !
    PatchGridEntity(const PatchGridEntity& original)
        : hostEntity_(original.hostEntity_), patchGrid_(original.patchGrid_) {}

    PatchGridEntity(PatchGridEntity&& original) noexcept
        : hostEntity_(std::move(original.hostEntity_)), patchGrid_(original.patchGrid_) {}

    //! @todo Please doc me !
    PatchGridEntity& operator=(const PatchGridEntity& original) {
      if (this != &original) {
        patchGrid_  = original.patchGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    //! @todo Please doc me !
    PatchGridEntity& operator=(PatchGridEntity&& original) noexcept {
      if (this != &original) {
        patchGrid_  = original.patchGrid_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    [[nodiscard]] bool equals(const PatchGridEntity& other) const { return untrimmedHostEntity() == other.untrimmedHostEntity(); }

    //! returns true if father entity exists
    [[nodiscard]] bool hasFather() const { return hostEntity_.hasFather(); }

    //! Create EntitySeed
    [[nodiscard]] EntitySeed seed() const { return EntitySeed(untrimmedHostEntity()); }

    //! Level of this element
    [[nodiscard]] int level() const { return untrimmedHostEntity().level(); }

    /** @brief The partition type for parallel computing */
    [[nodiscard]] PartitionType partitionType() const { return untrimmedHostEntity().partitionType(); }

    //! Geometry of this entity
    [[nodiscard]] Geometry geometry() const {
      //@todo Trim not hostEntity_
      static_assert(
          std::is_same_v<
              decltype(patchGrid_->patchGeometries_[this->level()].template localView<0, TrimmerType>()),
              typename GeometryKernel::NURBSPatch<dim, dimworld, ctype>::template GeometryLocalView<0, TrimmerType>>);
      // auto referenceEle= referenceElement(*this);
      auto geo = typename Geometry::Implementation(
          hostEntity_.geometry(), patchGrid_->patchGeometries_[this->level()].template localView<0, TrimmerType>());
      return Geometry(geo);
    }

    /** @brief Return the number of subEntities of codimension codim.
     */
    [[nodiscard]] unsigned int subEntities(unsigned int codim) const { return hostEntity_.subEntities(codim); }

    /** @brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    [[nodiscard]] typename GridImp::template Codim<cc>::Entity subEntity(int i) const {
      //@todo how to handout vertices and edges?
      //trimData().subEntity(int i)
      return PatchGridEntity<cc, dim, GridImp>(patchGrid_, hostEntity_.template subEntity<cc>(i));
    }

    //! First level intersection
    [[nodiscard]] PatchGridLevelIntersectionIterator<GridImp> ilevelbegin() const {
      return PatchGridLevelIntersectionIterator<GridImp>(
          patchGrid_, patchGrid_->parameterSpaceGrid().levelGridView(level()).ibegin(untrimmedHostEntity()));
    }

    //! Reference to one past the last neighbor
    PatchGridLevelIntersectionIterator<GridImp> ilevelend() const {
      return PatchGridLevelIntersectionIterator<GridImp>(
          patchGrid_, patchGrid_->parameterSpaceGrid().levelGridView(level()).iend(untrimmedHostEntity()));
    }

    //! First leaf intersection
    PatchGridLeafIntersectionIterator<GridImp> ileafbegin() const {
      return PatchGridLeafIntersectionIterator<GridImp>(
          patchGrid_, patchGrid_->parameterSpaceGrid().leafGridView().ibegin(untrimmedHostEntity()));
    }

    //! Reference to one past the last leaf intersection
    PatchGridLeafIntersectionIterator<GridImp> ileafend() const {
      return PatchGridLeafIntersectionIterator<GridImp>(
          patchGrid_, patchGrid_->parameterSpaceGrid().leafGridView().iend(untrimmedHostEntity()));
    }

    //! returns true if Entity has NO children
    bool isLeaf() const { return untrimmedHostEntity().isLeaf(); }

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::Entity father() const {
      return PatchGridEntity(patchGrid_, untrimmedHostEntity().father());
    }

    /** @brief Location of this element relative to the reference element element of the father.
     * This is sufficient to interpolate all dofs in conforming case.
     * Nonconforming may require access to neighbors of father and
     * computations with local coordinates.
     * On the fly case is somewhat inefficient since dofs  are visited several times.
     * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
     * implementation of numerical algorithms is only done for simple discretizations.
     * Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather() const {
      return LocalGeometry(typename LocalGeometry::Implementation(hostEntity_.geometryInFather()));
    }

    /** @brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    PatchGridHierarchicIterator<GridImp> hbegin(int maxLevel) const {
      return PatchGridHierarchicIterator<const GridImp>(patchGrid_, *this, maxLevel);
    }

    //! Returns iterator to one past the last son
    PatchGridHierarchicIterator<GridImp> hend(int maxLevel) const {
      return PatchGridHierarchicIterator<const GridImp>(patchGrid_, *this, maxLevel, true);
    }

    //! @todo Please doc me !
    bool wasRefined() const {
      if (patchGrid_->adaptationStep != GridImp::adaptDone) return false;

      int level = this->level();
      int index = patchGrid_->levelIndexSet(level).index(*this);
      return patchGrid_->refinementMark_[level][index];
    }

    //! @todo Please doc me !
    bool mightBeCoarsened() const { return true; }

    auto trimData() const { return patchGrid_->trimData(*this); }

    // /////////////////////////////////////////
    //   Internal stuff
    // /////////////////////////////////////////

    const auto& untrimmedHostEntity()const {
      if constexpr (requires {hostEntity_.untrimmedHostEntity();})
        return hostEntity_.untrimmedHostEntity();
      else
      return hostEntity_;
    }
    const auto& hostEntity()const {
      return hostEntity_;
    }

    auto& subId(int i, int codim)const  {
      if (codim==0)
        return hostEntity_.id();
      else if (codim==1)
      return patchGrid_->entityContainer_.globalEdgesIdOfElementsMap_.at(hostEntity_.id())[i];
      else if (codim==2)
        return patchGrid_->entityContainer_.globalVerticesIdOfElementsMap.at(hostEntity_.id())[i];
assert(codim>=0 and codim<3);
    }
  private:
    ParameterSpaceGridEntity hostEntity_;
    const GridImp* patchGrid_;

  };  // end of PatchGridEntity codim = 0

  template <int cd, int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamily,
            template <int, int, class> class PatchGridEntity>
  auto referenceElement(
      const PatchGridEntity<cd, dim, const PatchGrid<dim, dimworld, GridFamily, ScalarType>>& entity) {
    return GridFamily<dim,dimworld, ScalarType>::Trimmer::referenceElement(entity);
  }

  template <int cd, int dim, int dimworld, typename ScalarType, template <int,int, typename> typename GridFamily,
            template <int, int, class> class PatchGridEntity>
  auto referenceElement(
      const Entity<cd, dim, const PatchGrid<dim, dimworld, GridFamily, ScalarType>, PatchGridEntity>& entity) {
    return referenceElement(entity.impl());
  }

}  // namespace Dune::IGANEW
