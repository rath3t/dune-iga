// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file trimmer.hh
 * @brief Definition of the identity trimmer class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include "entitycontainer.hh"
#include "idset.hh"
#include "patchgridentityseed.hh"
#include "patchgridhierarchiciterator.hh"
#include "patchgridindexsets.hh"
#include "patchgridintersectioniterator.hh"
#include "patchgridintersections.hh"
#include "patchgridleafiterator.hh"
#include "patchgridleveliterator.hh"
#include "trimmedentity.hh"
#include "trimmedlocalgeometry.hh"
#include "patchtrimdata.hh"

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include <dune/iga/hierarchicpatch/patchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridview.hh>
#include <dune/iga/trimmer/identitytrimmer/patchgridlocalgeometry.hh>
#include <dune/iga/trimmer/intersectionvariants.hh>
#include <dune/iga/trimmer/localgeometryvariant.hh>

#include <dune/subgrid/subgrid.hh>

namespace Dune::IGANEW::DefaultTrim {
  template <typename HostIdType>
  struct IdType;
}
template <typename HostIdType>
struct std::hash<Dune::IGANEW::DefaultTrim::IdType<HostIdType>> {
  std::size_t operator()(const Dune::IGANEW::DefaultTrim::IdType<HostIdType>& k) const {
    using std::hash;

    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:

    return ((hash<HostIdType>()(k.id)
             ^ (hash<typename Dune::IGANEW::DefaultTrim::IdType<HostIdType>::ElementState>()(k.entityIdType) << 1))
            >> 1);
  }
};

namespace Dune::IGANEW {
  namespace GeometryKernel {
    template <int dim_, int dimworld_, typename ScalarType>
    class NURBSPatch;
  }

  namespace DefaultTrim {
    template <typename HostIdType>
    struct IdType {
      enum class EntityIdType { host, newId };
      EntityIdType entityIdType{EntityIdType::host};
      HostIdType id{};

      bool operator==(const IdType& other) const { return (entityIdType == other.entityIdType) and (id == other.id); }

      friend std::ostream& operator<<(std::ostream& stream, const IdType& id) {
        stream << "Type: " << std::string(id.entityIdType == EntityIdType::host ? "Host" : "new")
               << ", Key: " << id.id << "\n";
        return stream;
      }
    };

    template <typename HostIdType>
    bool operator<(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      if (lhs.entityIdType == rhs.entityIdType)
        return lhs.id < rhs.id;
      else if (lhs.entityIdType < rhs.entityIdType)
        return true;
      else
        return false;
    }

    template <typename HostIdType>
    bool operator>(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      return rhs < lhs;
    }

    template <typename HostIdType>
    bool operator==(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      if (lhs.entityIdType == rhs.entityIdType)
        return lhs.id == rhs.id;
      else
        return false;
    }

    template <typename HostIdType>
    bool operator<=(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      if (lhs.entityIdType == rhs.entityIdType)
        return lhs.id <= rhs.id;
      else if (lhs.entityIdType <= rhs.entityIdType)
        return true;
      else
        return false;
    }

    template <typename HostIdType>
    bool operator>=(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      return rhs >= lhs;
    }

    template <typename HostIdType>
    bool operator!=(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      if (lhs.entityIdType == rhs.entityIdType)
        return lhs.id != rhs.id;
      else
        return true;
    }

    template <typename HostIdType, int codim>
    struct EntityInfoImpl {
      static constexpr int codimension = codim;
      struct Empty {};
      int indexInLvlStorage;
      int lvl;
      bool stemFromTrim{false};
      IdType<HostIdType> id;


      auto stemsFromTrim() const {
        return stemFromTrim;
      }
    };

    template <typename HostIdType>
struct EntityInfoImpl<HostIdType,0> {
      static constexpr int codimension = 0;
      struct Empty {};
      int indexInLvlStorage{-1};
       int unTrimmedIndexInLvl{-1};
       int trimmedIndexInLvl{-1};
      int lvl;
      IdType<HostIdType> id;

      //if the element id is new we know that we are trimmed
      auto stemsFromTrim() const {
        return id.entityIdType==IdType<HostIdType>::EntityIdType::newId;
      }

       std::optional<IdType<HostIdType>> fatherId;
      ReservedVector<IdType<HostIdType>, 4> decendantIds;
    };

    /**
     * @brief Parameter struct representing parameters for the trimming operation.
     */
    struct Parameter {};

    /**
     * @brief ElementTrimData struct representing trim data for an element.
     * @tparam mydim_ Dimension of the element.
     * @tparam ScalarType Scalar type for the coordinates.
     */
    template <int mydim_, typename ScalarType>
    struct ElementTrimDataImpl {};

    /**
     * @brief ElementTrimDataContainer struct representing a container for element trim data.
     * @tparam ParameterSpaceGrid Type of the parameter space grid.
     */
    template <typename ParameterSpaceGrid>
    struct ElementTrimDataContainerImpl {};

    /**
     * @brief PatchTrimData struct representing trim data for a patch.
     * @tparam dim Dimension of the patch.
     * @tparam ScalarType Scalar type for the coordinates.
     */

    template <int dim, int dimworld, typename ScalarType>
    class TrimmerImpl;
    template <int dim, int dimworld, typename ScalarType>
    struct PatchGridFamily {
      using ctype   = ScalarType;
      using Grid    = PatchGrid<dim, dimworld, PatchGridFamily, ScalarType>;
      using Trimmer = TrimmerImpl<dim, dimworld, ScalarType>;

      using GlobalIdSet = PatchGridGlobalIdSet<const Grid>;

      using LocalIdSet    = PatchGridGlobalIdSet<const Grid>;
      using LevelIndexSet = PatchGridLevelIndexSet<const Grid>;
      using LeafIndexSet  = PatchGridLeafIndexSet<const Grid>;
      template <int codim, PartitionIteratorType pitype>
      using LeafIterator = PatchGridLeafIterator<codim, pitype, const Grid>;

      template <int codim, PartitionIteratorType pitype>
      using LevelIterator     = PatchGridLevelIterator<codim, pitype, const Grid>;
      using LeafIntersection  = PatchGridLeafIntersection<const Grid>;
      using LevelIntersection = PatchGridLevelIntersection<const Grid>;

      using LeafIntersectionIterator  = PatchGridLeafIntersectionIterator<const Grid>;
      using LevelIntersectionIterator = PatchGridLevelIntersectionIterator<const Grid>;
      using HierarchicIterator        = PatchGridHierarchicIterator<const Grid>;

      template <int codim>
      static const bool hasEntity = true;//codim == 0;

      template <int codim>
      static const bool hasEntityIterator =true;// codim == 0;

      template <int codim>
      static const bool hasHostEntity = true;

      struct TrimmerTraits {
        using ParameterSpaceGrid
            = Dune::SubGrid<dim,
                            YaspGrid<dim, TensorProductCoordinates<ScalarType, dim>>>;  ///< Type of the Parametric grid
        using HostIdType    = typename ParameterSpaceGrid::GlobalIdSet::IdType;
        using GlobalIdSetId = IdType<HostIdType>;
        using PatchTrimData = PatchTrimDataImpl<const Grid>;  ///< Patch trim data type.
        using TrimmingCurve = GeometryKernel::NURBSPatch<dim-1,dim,ctype>;
        using ElementInfo   = EntityInfoImpl<HostIdType, 0>;

        template <int codim>
        struct Codim {
          using EntityInfo = EntityInfoImpl<HostIdType, codim>;
          // This Geometry maps from the reference Element to knotspans
          using UntrimmedParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
          using TrimmedParameterSpaceGeometry
              = TrimmedLocalGeometryImpl<dim - codim, dim, const Grid, LocalGeometryTag::InParameterSpace>;
          using LocalParameterSpaceGeometry
              = Trim::LocalGeometryVariant<Trimmer, UntrimmedParameterSpaceGeometry, TrimmedParameterSpaceGeometry>;
          // This Geometry maps from the reference Element subTypes to 0..1
          using UntrimmedLocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;
          using TrimmedLocalGeometry
              = TrimmedLocalGeometryImpl<dim - codim, dim, const Grid, LocalGeometryTag::InReferenceElement>;
          using LocalGeometry = Trim::LocalGeometryVariant<Trimmer, UntrimmedLocalGeometry, TrimmedLocalGeometry>;
          // The entity living in the knotspan space
          using ParameterSpaceGridEntity = TrimmedParameterSpaceGridEntity<codim, dim, const Grid>;

          using ParameterSpaceGridEntitySeed = typename ParameterSpaceGrid::Traits::template Codim<codim>::EntitySeed;
          using EntityImp                    = PatchGridEntity<codim, dim, const Grid>;
          using EntitySeedImpl               = PatchGridEntitySeed<codim, const Grid>;
          using HostParameterSpaceGridEntity = typename ParameterSpaceGrid::Traits::template Codim<codim>::Entity;
          using UnTrimmedHostParameterSpaceGridEntity = typename ParameterSpaceGrid::Traits::template Codim<codim>::Entity;
        };

        using HostLeafIntersection                   = typename ParameterSpaceGrid::Traits::LeafIntersection;
        using HostLevelIntersection                  = typename ParameterSpaceGrid::Traits::LevelIntersection;
        using TrimmedParameterSpaceLeafIntersection  = TrimmedLeafIntersection<const Grid>;
        using TrimmedParameterSpaceLevelIntersection = TrimmedLevelIntersection<const Grid>;

        using ParameterSpaceLeafIntersection = TrimmedParameterSpaceLeafIntersection;
        // using ParameterSpaceLeafIntersection =
        // Trim::IntersectionVariant<Trimmer,UntrimmedParameterSpaceLeafIntersection,TrimmedParameterSpaceLeafIntersection>;
        // using ParameterSpaceLevelIntersection =
        // Trim::IntersectionVariant<Trimmer,UntrimmedParameterSpaceLevelIntersection,TrimmedParameterSpaceLevelIntersection>;
        using ParameterSpaceLevelIntersection = TrimmedParameterSpaceLevelIntersection;
      };

      using GeometryTypes = ReservedVector<GeometryType,2>;
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
      typename TrimmerTraits::GlobalIdSetId,
      LocalIdSet,
      typename TrimmerTraits::GlobalIdSetId,
      Communication<No_Comm>,
      PatchGridLevelGridViewTraits,
      PatchGridLeafGridViewTraits,
      PatchGridEntitySeed,
      PatchGridLocalGeometry,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LevelIndexSet::IndexType,
          GeometryTypes,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LeafIndexSet::IndexType,
          GeometryTypes>
          Traits;
      // clang-format on
      // using GlobalIdSetType =  PatchGridGlobalIdSet<const Grid>;

      // template <int codim, PartitionIteratorType pitype>
      // using LeafIterator = PatchGridLeafIterator<codim,pitype,GridImp>;
    };

    /**
     * @brief Trimmer struct representing a trimmer with identity trimming (no trimming at all).
     * @ingroup Trimmer
     * @tparam dim Dimension of the patch.
     * @tparam ScalarType Scalar type for the coordinates.
     */
    template <int dim, int dimworld, typename ScalarType>
    class TrimmerImpl {
     public:
      using GridFamily    = PatchGridFamily<dim, dimworld, ScalarType>;  ///< Scalar type for the coordinates.
      using GridTraits    = typename GridFamily::Traits;
      using TrimmerTraits = typename GridFamily::TrimmerTraits;

      template <int codim>
      using Codim   = typename TrimmerTraits::template Codim<codim>;
      using GridImp = typename GridTraits::Grid;
      friend GridImp;
      static constexpr bool isValid = (dim == 2 and (dimworld == 2 or dimworld == 3));

      template <int cd, int d, typename G>
      friend class TrimmedParameterSpaceGridEntity;
      template <int cd, PartitionIteratorType pi, typename G>
      friend class PatchGridLeafIterator;
      template <int cd, PartitionIteratorType pi, typename G>
      friend class PatchGridLevelIterator;
      friend class PatchGridHierarchicIterator<const GridImp>;
      friend class PatchGridLevelIntersection<const GridImp>;
      friend class PatchGridLeafIntersection<const GridImp>;
      friend class TrimmedLevelIntersection<const GridImp>;
      friend class TrimmedLeafIntersection<const GridImp>;


      template <int cd>
      using EntityInfo = typename GridFamily::TrimmerTraits::template Codim<cd>::EntityInfo;

      static constexpr int mydimension    = GridImp::dimension;       ///< Dimension of the patch.
      static constexpr int dimensionworld = GridImp::dimensionworld;  ///< Dimension of the world.
      using ctype                         = ScalarType;

      template <int codim>
      static constexpr bool isLocalGeometryLinear
          = true;  ///< boolean for the linearity of the local geometry, for the untrimmed case this is always true
      static constexpr bool isAlwaysTrivial = true;  ///< Boolean indicating if the trimming is always trivial, no
      ///< trimming or simple deletion of element.

      template <int codim>
      using Entity = typename GridFamily::Traits::template Codim<codim>::Entity;
      //! First level intersection
      [[nodiscard]] PatchGridLevelIntersectionIterator<const GridImp> ilevelbegin(const Entity<0>& ent) const {
        // DUNE_THROW(NotImplemented, "ilevelbegin");

        return PatchGridLevelIntersectionIterator<const GridImp>(
            grid_, parameterSpaceGrid().levelGridView(ent.level()).ibegin(ent.impl().getHostEntity().getHostEntity()));
      }

      //! Reference to one past the last neighbor
      PatchGridLevelIntersectionIterator<const GridImp> ilevelend(const Entity<0>& ent) const {
        // DUNE_THROW(NotImplemented, "ilevelend");

        return PatchGridLevelIntersectionIterator<const GridImp>(
            grid_, parameterSpaceGrid().levelGridView(ent.level()).iend(ent.impl().getHostEntity().getHostEntity()));
      }

      //! First leaf intersection
      PatchGridLeafIntersectionIterator<const GridImp> ileafbegin(const Entity<0>& ent) const {
        // DUNE_THROW(NotImplemented, "ileafbeginileafbegin");

        return PatchGridLeafIntersectionIterator<const GridImp>(
            grid_, parameterSpaceGrid().leafGridView().ibegin(ent.impl().getHostEntity().getHostEntity()));
      }

      //! Reference to one past the last leaf intersection
      PatchGridLeafIntersectionIterator<const GridImp> ileafend(const Entity<0>& ent) const {
        // DUNE_THROW(NotImplemented, "ileafendileafend");

        return PatchGridLeafIntersectionIterator<const GridImp>(
            grid_, parameterSpaceGrid().leafGridView().iend(ent.impl().getHostEntity().getHostEntity()));
      }

      template <class EntitySeed>
      typename GridFamily::Traits::template Codim<EntitySeed::codimension>::Entity entity(
          const EntitySeed& seed) const {
        using EntityImp               = typename TrimmerTraits::template Codim<EntitySeed::codimension>::EntityImp;
        auto [lvl, indexInLvlStorage] = seed.impl().data();
        return EntityImp{grid_,entityContainer_.template entity<EntitySeed::codimension>(lvl, indexInLvlStorage)};
      }

      template <int codim>
      using EntityImp = typename TrimmerTraits::template Codim<codim>::EntityImp;

      template <typename EntityImpl>
      typename GridFamily::Traits::template Codim<EntityImpl::codimension>::EntitySeed seed(
          const EntityImpl& ent) const {
        using EntitySeedImp = typename TrimmerTraits::template Codim<EntityImpl::codimension>::EntitySeedImpl;

        return EntitySeedImp(ent);
      }

      /**
       * @brief Default constructor for Trimmer.
       */
      TrimmerImpl() = default;

      using ParameterSpaceGrid = typename TrimmerTraits::ParameterSpaceGrid;  ///< Type of the Parametric grid
      template <int mydim>
      using ReferenceElementType =
          typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement;  ///< Reference element type.

      using ElementTrimData = ElementTrimDataImpl<ParameterSpaceGrid::dimension,
                                                  typename ParameterSpaceGrid::ctype>;  ///< Element trim data type.
      using PatchTrimData   =typename GridFamily::TrimmerTraits::PatchTrimData;  ///< Patch trim data type.
      using TrimmingCurve   =typename GridFamily::TrimmerTraits::TrimmingCurve;  ///< Patch trim data type.
      using ElementTrimDataContainer
          = ElementTrimDataContainerImpl<ParameterSpaceGrid>;  ///< Container for element trim data.

      using EntityContainer = VectorEntityContainer<GridImp>;
      template <int codim, PartitionIteratorType pitype>
      using ParameterSpaceLeafIterator = typename EntityContainer::template EntityConstInteratorImpl<codim>;

      template <int codim, PartitionIteratorType pitype>
      using ParameterSpaceLevelIterator = ParameterSpaceLeafIterator<codim, pitype>;
      using GlobalIdSet                 = typename GridFamily::GlobalIdSet;
      using LocalIdSet                  = typename GridFamily::LocalIdSet;
      using LeafIndexSet                = typename GridFamily::LeafIndexSet;
      using LevelIndexSet               = typename GridFamily::LevelIndexSet;
      friend LeafIndexSet;
      friend LevelIndexSet;

      using ParameterType = Parameter;  ///< Type for trimming parameters.

      static auto trimElement(const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element, const PatchTrimData& trimmingCurves);

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

      using UntrimmedParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

      /**
       * @brief Create the parameter space grid based on the patch and trim data.
       * @tparam dimworld Dimension of the world.
       * @param patchData NURBS patch data.
       * @param trimData Optional patch trim data.
       */
      void createParameterSpaceGrid() {
        untrimmedParameterSpaceGrid_
            = std::make_unique<UntrimmedParameterSpaceGrid>(grid_->patchGeometries_[0].uniqueKnotVector());

        parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);

        refineParameterSpaceGrid(0, true);
      }

      /**
       * @brief Constructor for Trimmer with patch and trim data.
       * @tparam dimworld Dimension of the world.
       * @param patch NURBS patch data.
       * @param trimData Optional patch trim data.
       */
      TrimmerImpl(GridImp& grid, const std::optional<PatchTrimData>& trimData)
          : grid_{&grid},
            leafIndexSet_(std::make_unique<LeafIndexSet>(*grid_)),
            globalIdSet_(std::make_unique<GlobalIdSet>(*grid_)),
            localIdSet_(std::make_unique<LocalIdSet>(*grid_)) ,trimData_{trimData}{
        createParameterSpaceGrid();
        setIndices();
      }



      TrimmerImpl& operator=(TrimmerImpl&& other) noexcept {
        this->grid_               = other.grid_;
        this->parameterSpaceGrid_ = other.parameterSpaceGrid_;
        leafIndexSet_             = std::make_unique<LeafIndexSet>(this->grid_);
        globalIdSet_              = std::make_unique<GlobalIdSet>(this->grid_);
        localIdSet_               = std::make_unique<LocalIdSet>(this->grid_);

        setIndices();

        return *this;
      }

      GridImp* grid_;

      /**
       * @brief Pass parameters to the trimmer.
       * @param par The parameters.
       */
      void setup(const ParameterType&) {}

      /**
       * @brief Get a const reference to the parameter space grid.
       * @return Const reference to the parameter space grid.
       */
      const ParameterSpaceGrid& parameterSpaceGrid() const { return *parameterSpaceGrid_; }

      /**
       * @brief Get a reference to the parameter space grid.
       * @return Reference to the parameter space grid.
       */
      ParameterSpaceGrid& parameterSpaceGrid() { return *parameterSpaceGrid_; }

      /**
       * @brief Refine the grid globally.
       * @param ref Number of refinement levels.
       */
      void globalRefine(int refCount) {
        if (refCount == 0) return;
        refineParameterSpaceGrid(refCount);
        setIndices();

        // @todo Trim move the refine here from the grid
        ;
      }

     protected:
     protected:
      //! compute the grid indices and ids
      void setIndices() {
        localIdSet_->update();

        globalIdSet_->update();

        // //////////////////////////////////////////
        //   Create the index sets
        // //////////////////////////////////////////
        for (int i = levelIndexSets_.size(); i <= maxLevel(); i++) {
          auto p = std::make_unique<typename GridFamily::LevelIndexSet>();
          levelIndexSets_.emplace_back(std::move(p));
        }

        for (int i = 0; i <= maxLevel(); i++)
          if (levelIndexSets_[i]) levelIndexSets_[i]->update(*grid_, i);

        leafIndexSet_->update(*grid_);
      }

      void createLevel(GridImp& grid, int lvl);
      void refineParameterSpaceGrid(int refCount, bool initFlag = false);

      /** @brief Return maximum level defined in this grid.
       *
       * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
       */
      [[nodiscard]] int maxLevel() const { return parameterSpaceGrid().maxLevel(); }

      EntityContainer entityContainer_;

      //! Our set of level indices
      std::vector<std::unique_ptr<typename GridFamily::LevelIndexSet>> levelIndexSets_;

      //! @todo Please doc me !
      std::unique_ptr<LeafIndexSet> leafIndexSet_;

      //! @todo Please doc me !
      std::unique_ptr<GlobalIdSet> globalIdSet_;

      //! @todo Please doc me !
      std::unique_ptr<LocalIdSet> localIdSet_;

      std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;  ///< The parameter space grid.

      std::unique_ptr<UntrimmedParameterSpaceGrid> untrimmedParameterSpaceGrid_;

      std::optional<PatchTrimData> trimData_;
    };

  }  // namespace DefaultTrim
}  // namespace Dune::IGANEW

#include "createlevel.hh"
#include "trimelement.hh"
