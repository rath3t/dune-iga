// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file parameterspace.hh
 * @brief Definition of the identity parameterspace class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/hierarchicpatch/patchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridview.hh>
#include <dune/iga/parameterspace/identity/idset.hh>
#include <dune/iga/parameterspace/identity/patchgridentityseed.hh>
#include <dune/iga/parameterspace/identity/patchgridhierarchiciterator.hh>
#include <dune/iga/parameterspace/identity/patchgridindexsets.hh>
#include <dune/iga/parameterspace/identity/patchgridintersectioniterator.hh>
#include <dune/iga/parameterspace/identity/patchgridleafiterator.hh>
#include <dune/iga/parameterspace/identity/patchgridleveliterator.hh>
#include <dune/iga/parameterspace/identity/patchgridlocalgeometry.hh>

namespace Dune::IGA {

namespace GeometryKernel {
  template <int dim_, int dimworld_, typename ScalarType>
  class NURBSPatch;
}

namespace IdentityParameterSpace {

  /**
   * @brief Parameter struct representing parameters for the trimming operation.
   */
  struct Parameter
  {
  };

  /**
   * @brief ElementTrimData struct representing trim data for an element.
   * @tparam mydim_ Dimension of the element.
   * @tparam ScalarType Scalar type for the coordinates.
   */
  template <int mydim_, typename ScalarType>
  struct ElementTrimDataImpl
  {
  };

  /**
   * @brief ElementTrimDataContainer struct representing a container for element trim data.
   * @tparam ParameterSpaceGrid Type of the parameter space grid.
   */
  template <typename ParameterSpaceGrid>
  struct ElementTrimDataContainerImpl
  {
  };

  /**
   * @brief PatchTrimData struct representing trim data for a patch.
   * @tparam dim Dimension of the patch.
   * @tparam ScalarType Scalar type for the coordinates.
   */
  template <int dim, typename ScalarType>
  struct PatchTrimDataImpl
  {
  };
  template <int dim, int dimworld, typename ScalarType>
  class ParameterSpaceImpl;
  template <int dim, int dimworld, typename ScalarType>
  struct PatchGridFamily
  {
    using ctype          = ScalarType;
    using Grid           = PatchGrid<dim, dimworld, PatchGridFamily, ScalarType>;
    using ParameterSpace = ParameterSpaceImpl<dim, dimworld, ScalarType>;

    using PatchTrimData = PatchTrimDataImpl<dim,
                                            ScalarType>; ///< Patch trim data type.

    using GlobalIdSet   = PatchGridGlobalIdSet<const Grid>;
    using LocalIdSet    = PatchGridLocalIdSet<const Grid>;
    using LevelIndexSet = PatchGridLevelIndexSet<const Grid>;
    using LeafIndexSet  = PatchGridLeafIndexSet<const Grid>;
    template <int codim, PartitionIteratorType pitype>
    using LeafIterator = PatchGridLeafIterator<codim, pitype, const Grid>;
    template <int codim, PartitionIteratorType pitype>
    using LevelIterator             = PatchGridLevelIterator<codim, pitype, const Grid>;
    using LeafIntersection          = PatchGridLeafIntersection<const Grid>;
    using LevelIntersection         = PatchGridLevelIntersection<const Grid>;
    using LeafIntersectionIterator  = PatchGridLeafIntersectionIterator<const Grid>;
    using LevelIntersectionIterator = PatchGridLevelIntersectionIterator<const Grid>;
    using HierarchicIterator        = PatchGridHierarchicIterator<const Grid>;

    template <int codim>
    static const bool hasEntity = true;

    template <int codim>
    static const bool hasEntityIterator = true;

    template <int codim>
    static const bool hasHostEntity = true;

    struct ParameterSpaceTraits
    {
      using ParameterType = Parameter; ///< Type for trimming parameters.

      using PatchTrimData = PatchTrimDataImpl<dim, ScalarType>;
      using ParameterSpaceGrid =
          YaspGrid<dim, TensorProductCoordinates<ScalarType, dim>>; ///< Type of the Parametric grid
      template <int codim>
      struct Codim
      {
        // This Geometry maps from the reference Element to knotspans
        using LocalParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
        // This Geometry maps from the reference Element subTypes to 0..1
        using LocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;
        // The entity living in the knotspan space
        using ParameterSpaceGridEntity = typename ParameterSpaceGrid::template Codim<codim>::Entity;

        using ParameterSpaceGridEntitySeed = typename ParameterSpaceGrid::Traits::template Codim<codim>::EntitySeed;
        using EntityImp                    = PatchGridEntity<codim, dim, const Grid>;
        using EntitySeedImp                = PatchGridEntitySeed<codim, const Grid>;
      };

      using ParameterSpaceLeafIntersection  = typename ParameterSpaceGrid::Traits::LeafIntersection;
      using ParameterSpaceLevelIntersection = typename ParameterSpaceGrid::Traits::LevelIntersection;
    };
    // clang-format off
      typedef GridTraits<
        dim, dimworld, Grid,
      PatchGridGeometry,
      PatchGridEntity,
      PatchGridLevelIterator,
      PatchGridLeafIntersection,
      PatchGridLevelIntersection,
      PatchGridLeafIntersectionIterator,
      PatchGridLevelIntersectionIterator,
      PatchGridHierarchicIterator,
      PatchGridLeafIterator,
      LevelIndexSet,
      LeafIndexSet,
      GlobalIdSet,
      typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
      LocalIdSet,
      typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LocalIdSet::IdType,
      Communication<No_Comm>,
      PatchGridLevelGridViewTraits,
      PatchGridLeafGridViewTraits,
      PatchGridEntitySeed,
      PatchGridLocalGeometry,
          typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LevelIndexSet::IndexType,
          typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LevelIndexSet::Types,
          typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LeafIndexSet::IndexType,
          typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LeafIndexSet::Types>
          Traits;
    // clang-format on
    // using GlobalIdSetType =  PatchGridGlobalIdSet<const Grid>;

    // template <int codim, PartitionIteratorType pitype>
    // using LeafIterator = PatchGridLeafIterator<codim,pitype,GridImp>;
  };

  /**
   * @brief ParameterSpace struct representing a parameterspace with identity trimming (no trimming at all).
   * @ingroup ParameterSpace
   * @tparam dim Dimension of the patch.
   * @tparam ScalarType Scalar type for the coordinates.
   */
  template <int dim, int dimworld, typename ScalarType = double>
  class ParameterSpaceImpl
  {
  public:
    using GridFamily           = PatchGridFamily<dim, dimworld, ScalarType>; ///< Scalar type for the coordinates.
    using PatchTrimData        = typename GridFamily::PatchTrimData;
    using GridTraits           = typename GridFamily::Traits;
    using ParameterSpaceTraits = typename GridFamily::ParameterSpaceTraits;

    static constexpr bool isValid = true;

    template <int codim>
    using Codim   = typename ParameterSpaceTraits::template Codim<codim>;
    using GridImp = typename GridTraits::Grid;
    friend GridImp;

    static constexpr int mydimension    = GridImp::dimension;      ///< Dimension of the patch.
    static constexpr int dimensionworld = GridImp::dimensionworld; ///< Dimension of the world.
    using ctype                         = ScalarType;

    template <int codim>
    static constexpr bool isLocalGeometryLinear =
        true; ///< boolean for the linearity of the local geometry, for the untrimmed case this is always true
    static constexpr bool isAlwaysTrivial = true; ///< Boolean indicating if the trimming is always trivial, no
                                                  ///< trimming or simple deletion of element.

    template <int codim>
    using Entity = typename GridFamily::Traits::template Codim<codim>::Entity;
    // First level intersection
    [[nodiscard]] PatchGridLevelIntersectionIterator<const GridImp> ilevelbegin(const Entity<0>& ent) const {
      return PatchGridLevelIntersectionIterator<const GridImp>(
          grid_, parameterSpaceGrid().levelGridView(ent.level()).ibegin(ent.impl().getLocalEntity()));
    }

    // Reference to one past the last neighbor
    PatchGridLevelIntersectionIterator<const GridImp> ilevelend(const Entity<0>& ent) const {
      return PatchGridLevelIntersectionIterator<const GridImp>(
          grid_, parameterSpaceGrid().levelGridView(ent.level()).iend(ent.impl().getLocalEntity()));
    }

    // First leaf intersection
    PatchGridLeafIntersectionIterator<const GridImp> ileafbegin(const Entity<0>& ent) const {
      return PatchGridLeafIntersectionIterator<const GridImp>(
          grid_, parameterSpaceGrid().leafGridView().ibegin(ent.impl().getLocalEntity()));
    }

    // Reference to one past the last leaf intersection
    PatchGridLeafIntersectionIterator<const GridImp> ileafend(const Entity<0>& ent) const {
      return PatchGridLeafIntersectionIterator<const GridImp>(
          grid_, parameterSpaceGrid().leafGridView().iend(ent.impl().getLocalEntity()));
    }

    /**
     * @brief Default constructor for ParameterSpace.
     */
    ParameterSpaceImpl() = default;

    using ParameterSpaceGrid = typename ParameterSpaceTraits::ParameterSpaceGrid; ///< Type of the Parametric grid
    template <int mydim>
    using ReferenceElementType =
        typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement; ///< Reference element type.

    using ElementTrimData = ElementTrimDataImpl<ParameterSpaceGrid::dimension,
                                                typename ParameterSpaceGrid::ctype>; ///< Element trim data type.

    using ElementTrimDataContainer =
        ElementTrimDataContainerImpl<ParameterSpaceGrid>; ///< Container for element trim data.

    template <int codim, PartitionIteratorType pitype>
    using ParameterSpaceLeafIterator =
        typename ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LeafIterator;

    using GlobalIdSet   = typename GridFamily::GlobalIdSet;
    using LocalIdSet    = typename GridFamily::LocalIdSet;
    using LeafIndexSet  = typename GridFamily::LeafIndexSet;
    using LevelIndexSet = typename GridFamily::LevelIndexSet;

    using EntityContainerType = Empty;

    using ParameterType = Parameter; ///< Type for trimming parameters.

    /**
     * @brief Get the reference element for a given entity.
     * @tparam EntityType Type of the entity.
     * @param entity The entity for which the reference element is requested.
     * @return Reference element for the entity.
     */
    template </* Dune::Concept::Entity */ typename EntityType>
    static auto referenceElement(const EntityType& entity) {
      return Dune::referenceElement<ctype, EntityType::mydimension>(entity.type());
    }

    /**
     * @brief Create the parameter space grid based on the patch and trim data.
     * @tparam dimworld Dimension of the world.
     * @param patchData NURBS patch data.
     * @param trimData Optional patch trim data.
     */
    void createParameterSpaceGrid() {
      parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(grid_->patchGeometries_[0].uniqueKnotVector());
    }

    /**
     * @brief Constructor for ParameterSpace with patch and trim data.
     * @tparam dimworld Dimension of the world.
     * @param patch NURBS patch data.
     * @param trimData Optional patch trim data.
     */
    ParameterSpaceImpl(GridImp& grid, const std::optional<typename GridFamily::PatchTrimData>& trimData,
                       const ParameterType& par = {})
        : grid_{&grid},
          leafIndexSet_(std::make_unique<LeafIndexSet>(*grid_)),
          globalIdSet_(std::make_unique<GlobalIdSet>(*grid_)),
          localIdSet_(std::make_unique<LocalIdSet>(*grid_)) {
      createParameterSpaceGrid();
      update(grid_);
    }

    ParameterSpaceImpl& operator=(ParameterSpaceImpl&& other) noexcept {
      this->grid_               = std::move(other.grid_);
      this->parameterSpaceGrid_ = std::move(other.parameterSpaceGrid_);
      leafIndexSet_             = std::make_unique<LeafIndexSet>(this->grid_);
      globalIdSet_              = std::make_unique<GlobalIdSet>(this->grid_);
      localIdSet_               = std::make_unique<LocalIdSet>(this->grid_);

      update(grid_);

      return *this;
    }

    GridImp* grid_;

    /**
     * @brief Pass parameters to the parameterspace.
     * @param par The parameters.
     */
    void setup() {}
    void setParameters(const ParameterType&) {}

    /**
     * @brief Get a const reference to the parameter space grid.
     * @return Const reference to the parameter space grid.
     */
    const ParameterSpaceGrid& parameterSpaceGrid() const {
      return *parameterSpaceGrid_;
    }

    /**
     * @brief Get a reference to the parameter space grid.
     * @return Reference to the parameter space grid.
     */
    ParameterSpaceGrid& parameterSpaceGrid() {
      return *parameterSpaceGrid_;
    }

    template <class EntitySeed>
    typename GridFamily::Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const {
      using EntityImp = typename ParameterSpaceTraits::template Codim<EntitySeed::codimension>::EntityImp;

      return EntityImp(grid_, parameterSpaceGrid().entity(seed.impl().hostEntitySeed()));
    }

    template <class Entity>
    typename GridFamily::Traits::template Codim<Entity::codimension>::EntitySeed seed(const Entity& ent) const {
      using EntitySeedImp = typename ParameterSpaceTraits::template Codim<Entity::codimension>::EntitySeedImp;

      return EntitySeedImp(ent.getLocalEntity());
    }

    /**
     * @brief Refine the grid globally.
     * @param ref Number of refinement levels.
     */
    auto globalRefine(int refCount) {
      parameterSpaceGrid().globalRefine(refCount);
      update(grid_);
    }

  protected:
  protected:
    // compute the grid indices and ids
    void update(GridImp* grid) {
      grid_ = grid;
      localIdSet_->update();

      globalIdSet_->update();

      // //////////////////////////////////////////
      //   Create the index sets
      // //////////////////////////////////////////
      for (int i = levelIndexSets_.size(); i <= maxLevel(); i++) {
        auto p = std::make_unique<LevelIndexSet>();
        levelIndexSets_.emplace_back(std::move(p));
      }

      for (int i = 0; i <= maxLevel(); i++)
        if (levelIndexSets_[i])
          levelIndexSets_[i]->update(*grid_, i);

      leafIndexSet_->update(*grid_);
    }

    /** @brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    [[nodiscard]] int maxLevel() const {
      return parameterSpaceGrid().maxLevel();
    }

    size_t numBoundarySegments() const {
      return parameterSpaceGrid_->numBoundarySegments();
    }

    // Our set of level indices
    std::vector<std::unique_ptr<LevelIndexSet>> levelIndexSets_;

    // TODO Please doc me !
    std::unique_ptr<LeafIndexSet> leafIndexSet_;

    // TODO Please doc me !
    std::unique_ptr<GlobalIdSet> globalIdSet_;

    // TODO Please doc me !
    std::unique_ptr<LocalIdSet> localIdSet_;

    std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_; ///< The parameter space grid.
  };

} // namespace IdentityParameterSpace
} // namespace Dune::IGA
