// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file parameterspace.hh
 * @brief Definition of the identity parameterspace class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/hierarchicpatch/patchgridfwd.hh>
#include <dune/iga/hierarchicpatch/patchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridview.hh>
#include <dune/iga/parameterspace/default/elementtrimdata.hh>
#include <dune/iga/parameterspace/default/entitycontainer.hh>
#include <dune/iga/parameterspace/default/entityinfo.hh>
#include <dune/iga/parameterspace/default/idset.hh>
#include <dune/iga/parameterspace/default/patchgridentityseed.hh>
#include <dune/iga/parameterspace/default/patchgridhierarchiciterator.hh>
#include <dune/iga/parameterspace/default/patchgridindexsets.hh>
#include <dune/iga/parameterspace/default/patchgridintersectioniterator.hh>
#include <dune/iga/parameterspace/default/patchgridintersections.hh>
#include <dune/iga/parameterspace/default/patchgridleafiterator.hh>
#include <dune/iga/parameterspace/default/patchgridleveliterator.hh>
#include <dune/iga/parameterspace/default/patchtrimdata.hh>
#include <dune/iga/parameterspace/default/referenceelement.hh>
#include <dune/iga/parameterspace/default/parameterspaceentity.hh>
#include <dune/iga/parameterspace/default/trimmedlocalgeometry.hh>
#include <dune/iga/parameterspace/identity/patchgridlocalgeometry.hh>
#include <dune/iga/parameterspace/localgeometryvariant.hh>
#include <dune/subgrid/subgrid.hh>

namespace Dune::IGA::DefaultTrim {
template <typename HostIdType>
struct IdType;
}
template <typename HostIdType>
struct std::hash<Dune::IGA::DefaultTrim::IdType<HostIdType>>
{
  std::size_t operator()(const Dune::IGA::DefaultTrim::IdType<HostIdType>& k) const {
    using std::hash;

    // Compute individual hash values for first, second and third and combine them using XOR
    // and bit shifting:
    // todo
    // return (hash<HostIdType>()(k.id) ^
    //         hash<typename Dune::IGA::DefaultTrim::IdType<HostIdType>::ElementState>()(k.entityIdType) << 1) >>
    //        1;
    return (hash<HostIdType>()(k.id));
  }
};

namespace Dune::IGA {
namespace GeometryKernel {
  template <int dim_, int dimworld_, typename ScalarType>
  class NURBSPatch;
}

namespace DefaultTrim {

  template <typename HostIdType>
  struct IdType
  {
    enum class EntityIdType
    {
      host,
      newId
    };
    EntityIdType entityIdType{EntityIdType::host};
    HostIdType id{};
    std::optional<HostIdType> hostId{};

    bool operator==(const IdType& other) const {
      return (entityIdType == other.entityIdType) and (id == other.id);
    }

    friend std::ostream& operator<<(std::ostream& stream, const IdType& id) {
      stream << "Type: " << std::string(id.entityIdType == EntityIdType::host ? "Host" : "new") << ", Key: " << id.id
             << "\n";
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

  /**
   * @brief Parameter struct representing parameters for the trimming operation.
   */
  struct Parameter
  {
    size_t splitter      = 100;
    int dummy            = 7;     ///< Dummy variable.
    double trimPrecision = 1e-10; ///< Precision for trimming.
  };

  /**
   * @brief PatchTrimData struct representing trim data for a patch.
   * @tparam dim Dimension of the patch.
   * @tparam ScalarType Scalar type for the coordinates.
   */
  template <int dim, int dimworld, typename ScalarType>
  class ParameterSpaceImpl;
  template <int dim, int dimworld, typename ScalarType>
  struct PatchGridFamily
  {
    using ctype                   = ScalarType;
    static constexpr int patchDim = dim;
    static constexpr int worldDim = dimworld;
    using Grid                    = PatchGrid<dim, dimworld, PatchGridFamily, ScalarType>;
    using ParameterSpace          = ParameterSpaceImpl<dim, dimworld, ScalarType>;

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
    static const bool hasEntity = true; // codim == 0;

    template <int codim>
    static const bool hasEntityIterator = true; // codim == 0;

    template <int codim>
    static const bool hasHostEntity = true;

    struct ParameterSpaceTraits
    {
      using ParameterType       = Parameter; ///< Type for trimming parameters.
      using YASPGridType        = YaspGrid<dim, TensorProductCoordinates<ScalarType, dim>>;
      using ParameterSpaceGrid  = Dune::SubGrid<dim, YASPGridType>; ///< Type of the Parametric grid
      using HostIdType          = typename ParameterSpaceGrid::GlobalIdSet::IdType;
      using PersistentIndexType = typename YASPGridType::PersistentIndexType;

      using GlobalIdSetId        = IdType<HostIdType>;
      using PatchTrimData        = PatchTrimDataImpl<const Grid>; ///< Patch trim data type.
      using TrimmingCurve        = GeometryKernel::NURBSPatch<dim - 1, dim, ctype>;
      using ElementInfo          = EntityInfoImpl<ParameterSpaceTraits, 0>;
      using ReferenceElementType = TrimmedReferenceElement<2, const Grid>;

      using ElementTrimData = ElementTrimDataImpl<const Grid>; ///< Element trim data type.

      template <int codim>
      struct Codim
      {
        using EntityInfo = EntityInfoImpl<ParameterSpaceTraits, codim>;
        // This Geometry maps from the reference Element to knotspans
        using UntrimmedParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
        using TrimmedParameterSpaceGeometry =
            TrimmedLocalGeometryImpl<dim - codim, dim, const Grid, LocalGeometryTag::InParameterSpace>;
        using LocalParameterSpaceGeometry =
            Trim::LocalGeometryVariant<ParameterSpace, UntrimmedParameterSpaceGeometry, TrimmedParameterSpaceGeometry>;
        // This Geometry maps from the reference Element subTypes to 0..1
        using UntrimmedLocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;
        using TrimmedLocalGeometry =
            TrimmedLocalGeometryImpl<dim - codim, dim, const Grid, LocalGeometryTag::InReferenceElement>;
        using LocalGeometry = Trim::LocalGeometryVariant<ParameterSpace, UntrimmedLocalGeometry, TrimmedLocalGeometry>;
        // The entity living in the knotspan space
        using ParameterSpaceGridEntity = TrimmedParameterSpaceGridEntity<codim, dim, const Grid>;

        using ParameterSpaceGridEntitySeed = typename ParameterSpaceGrid::Traits::template Codim<codim>::EntitySeed;
        using EntityImp                    = PatchGridEntity<codim, dim, const Grid>;
        using EntitySeedImpl               = PatchGridEntitySeed<codim, const Grid>;
        using HostParameterSpaceGridEntity = typename ParameterSpaceGrid::Traits::template Codim<codim>::Entity;
      };

      using HostLeafIntersection                   = typename ParameterSpaceGrid::Traits::LeafIntersection;
      using HostLevelIntersection                  = typename ParameterSpaceGrid::Traits::LevelIntersection;
      using TrimmedParameterSpaceLeafIntersection  = TrimmedLeafIntersection<const Grid>;
      using TrimmedParameterSpaceLevelIntersection = TrimmedLevelIntersection<const Grid>;

      using ParameterSpaceLeafIntersection  = TrimmedParameterSpaceLeafIntersection;
      using ParameterSpaceLevelIntersection = TrimmedParameterSpaceLevelIntersection;
    };

    using GeometryTypes = ReservedVector<GeometryType, 2>;

    typedef GridTraits<
        dim, dimworld, Grid, PatchGridGeometry, PatchGridEntity, PatchGridLevelIterator, PatchGridLeafIntersection,
        PatchGridLevelIntersection, PatchGridLeafIntersectionIterator, PatchGridLevelIntersectionIterator,
        PatchGridHierarchicIterator, PatchGridLeafIterator, LevelIndexSet, LeafIndexSet, GlobalIdSet,
        typename ParameterSpaceTraits::GlobalIdSetId, LocalIdSet, typename ParameterSpaceTraits::GlobalIdSetId,
        Communication<No_Comm>, PatchGridLevelGridViewTraits, PatchGridLeafGridViewTraits, PatchGridEntitySeed,
        PatchGridLocalGeometry, typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LevelIndexSet::IndexType,
        GeometryTypes, typename ParameterSpaceTraits::ParameterSpaceGrid::Traits::LeafIndexSet::IndexType,
        GeometryTypes>
        Traits;
  };

  /**
   * @brief ParameterSpace struct representing a parameterspace with identity trimming (no trimming at all).
   * @ingroup ParameterSpace
   * @tparam dim Dimension of the patch.
   * @tparam ScalarType Scalar type for the coordinates.
   */
  template <int dim, int dimworld, typename ScalarType>
  class ParameterSpaceImpl
  {
  public:
    using GridFamily           = PatchGridFamily<dim, dimworld, ScalarType>; ///< Scalar type for the coordinates.
    using GridTraits           = typename GridFamily::Traits;
    using ParameterSpaceTraits = typename GridFamily::ParameterSpaceTraits;

    template <int codim>
    using Codim   = typename ParameterSpaceTraits::template Codim<codim>;
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
    using EntityInfo = typename GridFamily::ParameterSpaceTraits::template Codim<cd>::EntityInfo;

    template <int cd>
    using HostEntity = typename GridFamily::ParameterSpaceTraits::template Codim<cd>::HostParameterSpaceGridEntity;

    static constexpr int mydimension    = GridImp::dimension;      ///< Dimension of the patch.
    static constexpr int dimensionworld = GridImp::dimensionworld; ///< Dimension of the world.
    using ctype                         = ScalarType;

    template <int codim>
    static constexpr bool isLocalGeometryLinear =
        true; ///< boolean for the linearity of the local geometry, for the untrimmed case this is always true
    static constexpr bool isAlwaysTrivial = false; ///< Boolean indicating if the trimming is always trivial, no
    ///< trimming or simple deletion of element.

    template <int codim>
    using Entity = typename GridFamily::Traits::template Codim<codim>::Entity;

    template <int codim>
    using YASPEntity =
        typename GridFamily::ParameterSpace::ParameterSpaceTraits::YASPGridType::Traits::template Codim<codim>::Entity;

    // First level intersection
    [[nodiscard]] PatchGridLevelIntersectionIterator<const GridImp> ilevelbegin(const Entity<0>& ent) const {
      using IntersectionIterator = PatchGridLevelIntersectionIterator<const GridImp>;
      auto& localEntity          = ent.impl().getLocalEntity();

      if (not localEntity.isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
        return IntersectionIterator(
            grid_, parameterSpaceGrid().levelGridView(ent.level()).ibegin(localEntity.getHostEntity()));
      return IntersectionIterator(grid_, localEntity, IntersectionIterator::PositionToken::Begin,
                                  parameterSpaceGrid().levelGridView(ent.level()).ibegin(localEntity.getHostEntity()));
    }

    // Reference to one past the last neighbor
    PatchGridLevelIntersectionIterator<const GridImp> ilevelend(const Entity<0>& ent) const {
      using IntersectionIterator = PatchGridLevelIntersectionIterator<const GridImp>;
      auto& localEntity          = ent.impl().getLocalEntity();

      if (not localEntity.isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
        return IntersectionIterator(grid_,
                                    parameterSpaceGrid().levelGridView(ent.level()).iend(localEntity.getHostEntity()));
      return IntersectionIterator(grid_, localEntity, IntersectionIterator::PositionToken::End,
                                  parameterSpaceGrid().levelGridView(ent.level()).ibegin(localEntity.getHostEntity()));
    }

    // First leaf intersection
    PatchGridLeafIntersectionIterator<const GridImp> ileafbegin(const Entity<0>& ent) const {
      using IntersectionIterator = PatchGridLeafIntersectionIterator<const GridImp>;
      auto& localEntity          = ent.impl().getLocalEntity();

      if (not localEntity.isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
        return IntersectionIterator(grid_, parameterSpaceGrid().leafGridView().ibegin(localEntity.getHostEntity()));
      return IntersectionIterator(grid_, localEntity, IntersectionIterator::PositionToken::Begin,
                                  parameterSpaceGrid().leafGridView().ibegin(localEntity.getHostEntity()));
    }

    // Reference to one past the last leaf intersection
    PatchGridLeafIntersectionIterator<const GridImp> ileafend(const Entity<0>& ent) const {
      using IntersectionIterator = PatchGridLeafIntersectionIterator<const GridImp>;
      auto& localEntity          = ent.impl().getLocalEntity();

      if (not localEntity.isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
        return IntersectionIterator(grid_, parameterSpaceGrid().leafGridView().iend(localEntity.getHostEntity()));
      return IntersectionIterator(grid_, localEntity, IntersectionIterator::PositionToken::End,
                                  parameterSpaceGrid().leafGridView().ibegin(localEntity.getHostEntity()));
    }

    template <class EntitySeed>
    typename GridFamily::Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const {
      assert(seed.isValid());
      using EntityImp               = typename ParameterSpaceTraits::template Codim<EntitySeed::codimension>::EntityImp;
      auto [lvl, indexInLvlStorage] = seed.impl().data();
      return EntityImp{grid_, entityContainer_.template entity<EntitySeed::codimension>(lvl, indexInLvlStorage)};
    }

    template <int codim>
    using EntityImp = typename ParameterSpaceTraits::template Codim<codim>::EntityImp;

    template <typename EntityImpl>
    typename GridFamily::Traits::template Codim<EntityImpl::codimension>::EntitySeed seed(const EntityImpl& ent) const {
      using EntitySeedImp = typename ParameterSpaceTraits::template Codim<EntityImpl::codimension>::EntitySeedImpl;

      return EntitySeedImp(ent);
    }

    /**
     * @brief Default constructor for ParameterSpace.
     */
    ParameterSpaceImpl() = default;

    using ParameterSpaceGrid = typename ParameterSpaceTraits::ParameterSpaceGrid; ///< Type of the Parametric grid
    template <int mydim>
    using ReferenceElementType =
        typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement; ///< Reference element type.

    using ElementTrimData = typename GridFamily::ParameterSpaceTraits::ElementTrimData; ///< Element trim data type.
    using PatchTrimData   = typename GridFamily::ParameterSpaceTraits::PatchTrimData;   ///< Patch trim data type.
    using TrimmingCurve   = typename GridFamily::ParameterSpaceTraits::TrimmingCurve;   ///< Patch trim data type.

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

    using ParameterType = Parameter; ///< Type for trimming parameters.

    using GlobalIdType = typename GridFamily::ParameterSpaceTraits::GlobalIdSetId;

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
     */
    void createParameterSpaceGrid() {
      untrimmedParameterSpaceGrid_ =
          std::make_unique<UntrimmedParameterSpaceGrid>(grid_->patchGeometries_[0].uniqueKnotVector());

      parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);

      refineParameterSpaceGrid(0, true);
    }

    /**
     * @brief Constructor for ParameterSpace with patch and trim data.
     * @param grid patch grid
     * @param trimData Optional patch trim data.
     */
    ParameterSpaceImpl(GridImp& grid, const std::optional<PatchTrimData>& trimData, const ParameterType& par = {})
        : grid_{&grid},
          leafIndexSet_(std::make_unique<LeafIndexSet>(*grid_)),
          globalIdSet_(std::make_unique<GlobalIdSet>(*grid_)),
          localIdSet_(std::make_unique<LocalIdSet>(*grid_)),
          trimData_{trimData},
          parameters_(par) {
      setup();
      createParameterSpaceGrid();
      update(grid_);
    }

    ParameterSpaceImpl& operator=(ParameterSpaceImpl&& other) noexcept {
      this->grid_               = other.grid_;
      this->parameterSpaceGrid_ = other.parameterSpaceGrid_;
      leafIndexSet_             = std::make_unique<LeafIndexSet>(this->grid_);
      globalIdSet_              = std::make_unique<GlobalIdSet>(this->grid_);
      localIdSet_               = std::make_unique<LocalIdSet>(this->grid_);

      update(grid_);

      return *this;
    }

    GridImp* grid_;

    void setup() {
      if (trimData_.has_value())
        trimData_->prepare(parameters_, grid_->tensorProductCoordinates(0));
    }
    /**
     * @brief Change the parameters to the parameterspace.
     * @param par The parameters.
     */
    void setParameters(const ParameterType& par) {
      parameters_ = par;
    }

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

    /**
     * @brief Refine the grid globally.
     * @param refCount Number of refinement levels.
     */
    void globalRefine(int refCount) {
      if (refCount == 0)
        return;
      refineParameterSpaceGrid(refCount);
      update(grid_);

      // @todo Trim move the refine here from the grid
    }
    auto& patchTrimData() const {
      return *trimData_;
    }

    /**
     * \brief Creates trimInfos for each element at the requested level (or max level if not specified)
     * @todo Add execution policy
     * @param level_
     * @return
     */
    std::vector<ElementTrimData> trimElements(std::optional<int> level_ = std::nullopt) {
      int level     = level_.value_or(maxLevel());
      bool initFlag = level == 0;
      if (initFlag)
        numBoundarySegments_ = untrimmedParameterSpaceGrid_->numBoundarySegments();

      std::vector<ElementTrimData> elementTrimDatas;
      auto gv = untrimmedParameterSpaceGrid_->levelGridView(level);
      for (const auto& ele : elements(gv)) {
        if (trimData_.has_value())
          elementTrimDatas.emplace_back(trimElement(ele, gv, trimData_.value(), initFlag));
        else
          elementTrimDatas.emplace_back(ElementTrimFlag::full, ele);
      }
      return elementTrimDatas;
    }

    // compute the grid indices and ids
    void update(GridImp* grid) {
      grid_ = grid;
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
        if (levelIndexSets_[i])
          levelIndexSets_[i]->update(*grid_, i);

      leafIndexSet_->update(*grid_);
    }

    void refineParameterSpaceGrid(int refCount, bool initFlag = false);

    // The following are helper methods for `refineParameterSpaceGrid`
    ElementTrimData trimElement(const YASPEntity<0>& element, const auto& gv, const PatchTrimData& patchTrimData,
                                bool initFlag);

    GlobalIdType makeElementID(const HostEntity<0>& ele);
    void createAndSaveElementInfo(const std::tuple<unsigned int, unsigned int, int>& indices, const HostEntity<0>& ele,
                                  bool trimmed);

    void collectElementEdges(int level, const HostEntity<0>& ele, const ElementTrimData& eleTrimData);
    void collectElementVertices(int level, const HostEntity<0>& ele, const ElementTrimData& eleTrimData);
    void createSubEntities(int level);
    void createElements(int level, const std::vector<ElementTrimData>& trimDatas);
    GlobalIdType idForTrimmedHostEdge(typename ParameterSpaceTraits::PersistentIndexType hostEdgeId,
                                      const typename ElementTrimData::EdgeInfo& trimmedEdge);
    GlobalIdType idForTrimmedVertex(const FieldVector<double, 2>& vertex);

    /** @brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    [[nodiscard]] int maxLevel() const {
      return parameterSpaceGrid().maxLevel();
    }

    size_t numBoundarySegments() const {
      assert(numBoundarySegments_ != std::numeric_limits<size_t>::max());
      return numBoundarySegments_;
    }

    EntityContainer entityContainer_;

    // Our set of level indices
    std::vector<std::unique_ptr<typename GridFamily::LevelIndexSet>> levelIndexSets_;

    std::unique_ptr<LeafIndexSet> leafIndexSet_;

    std::unique_ptr<GlobalIdSet> globalIdSet_;

    std::unique_ptr<LocalIdSet> localIdSet_;

    std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_; ///< The parameter space grid.

    std::unique_ptr<UntrimmedParameterSpaceGrid> untrimmedParameterSpaceGrid_;

    std::optional<PatchTrimData> trimData_;

    ParameterType parameters_;

    size_t numBoundarySegments_{std::numeric_limits<size_t>::max()};
    std::map<typename ParameterSpaceTraits::YASPGridType::GlobalIdSetType::IdType,
             std::vector<std::tuple<Impl::CurveLoopIndexEncoder::IndexResult, double, double, size_t>>>
        boundarySegmentsArchive_;
  };

} // namespace DefaultTrim
} // namespace Dune::IGA

#include "createentities.hh"
#include "createlevel.hh"
#include "trimelement.hh"
