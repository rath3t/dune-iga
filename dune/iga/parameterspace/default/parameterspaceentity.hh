// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/parameterspace/default/patchgridentityseed.hh>

namespace Dune::IGA::DefaultParameterSpace {

template <class GridImp>
class PatchGridHierarchicIterator;
// : TrimPatchEntitiy
template <int codim_, int dim, class GridImp>
class TrimmedParameterSpaceGridEntity
{
  using ctype = typename GridImp::ctype;

  static constexpr int mydimension = dim;

  using LocalCoordinate = FieldVector<ctype, mydimension>;

  friend PatchGridHierarchicIterator<const GridImp>;
  friend PatchGridEntitySeed<codim_, const GridImp>;

  using ParameterSpace = typename GridImp::ParameterSpace;
  friend ParameterSpace;
  using GlobalIdSetIdType = typename ParameterSpace::ParameterSpaceTraits::GlobalIdSetId;
  using EntityInfo        = typename ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::EntityInfo;

  using TrimInfo =
      std::conditional_t<codim_ == 0, typename ParameterSpace::ElementTrimData,
                         std::conditional_t<codim_ == 1, typename ParameterSpace::ElementTrimData::EdgeInfo,
                                            typename ParameterSpace::ElementTrimData::VertexInfo>>;

  using HostParameterSpaceGridEntity =
      typename ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::HostParameterSpaceGridEntity;
  using UntrimmedParameterSpaceGeometry =
      typename GridImp::ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::UntrimmedParameterSpaceGeometry;
  using TrimmedParameterSpaceGeometry =
      typename GridImp::ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::TrimmedParameterSpaceGeometry;
  using LocalParameterSpaceGeometry =
      typename GridImp::ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::LocalParameterSpaceGeometry;
  using ParameterSpaceGridEntitySeed =
      typename GridImp::ParameterSpace::ParameterSpaceTraits::template Codim<codim_>::ParameterSpaceGridEntitySeed;

public:
  TrimmedParameterSpaceGridEntity()                                                            = default;
  TrimmedParameterSpaceGridEntity(const TrimmedParameterSpaceGridEntity& other) noexcept       = default;
  TrimmedParameterSpaceGridEntity(TrimmedParameterSpaceGridEntity&& other) noexcept            = default;
  TrimmedParameterSpaceGridEntity& operator=(const TrimmedParameterSpaceGridEntity& other)     = default;
  TrimmedParameterSpaceGridEntity& operator=(TrimmedParameterSpaceGridEntity&& other) noexcept = default;

  // @todo tidy up these constructors
  // Entity with codim 0 but trimmed thus needs trimdata
  template <typename = void>
  requires(codim_ == 0)
  TrimmedParameterSpaceGridEntity(const GridImp* grid, const HostParameterSpaceGridEntity& untrimmedElement,
                                  const EntityInfo& entInfo, const TrimInfo& trimData)
      : grid_{grid},
        entityInfo_{entInfo},
        hostEntity_{untrimmedElement},
        trimData_{trimData} {
    assert(entityInfo_.lvl == untrimmedElement.level());
  }

  // Entity with codim 0
  template <typename = void>
  requires(codim_ == 0)
  TrimmedParameterSpaceGridEntity(const GridImp* grid, const HostParameterSpaceGridEntity& untrimmedElement,
                                  const EntityInfo& entInfo)
      : grid_{grid},
        entityInfo_{entInfo},
        hostEntity_{untrimmedElement},
        trimData_{} {
    assert(entityInfo_.lvl == untrimmedElement.level());
  }

  template <typename = void>
  requires(codim_ != 0)
  TrimmedParameterSpaceGridEntity(const GridImp* grid, const HostParameterSpaceGridEntity& untrimmedElement,
                                  const EntityInfo& entInfo)
      : grid_{grid},
        entityInfo_{entInfo},
        hostEntity_{untrimmedElement},
        trimData_{entInfo.trimInfo} {
    assert(entityInfo_.lvl == untrimmedElement.level());
  }

  template <typename = void>
  requires(codim_ != 0)
  TrimmedParameterSpaceGridEntity(const GridImp* grid, EntityInfo entInfo)
      : grid_{grid},
        entityInfo_{entInfo},
        trimData_{entInfo.trimInfo} {}

  auto& id() const {
    return entityInfo_.id;
  }

  auto isTrimmed() const {
    return entityInfo_.isTrimmed();
  }

  auto index() const {
    // TODO
    // if constexpr (codim_ == 0)
    //   return isTrimmed() ? entityInfo_.trimmedIndexInLvl : entityInfo_.unTrimmedIndexInLvl;
    // else
    return entityInfo_.indexInLvlStorage;
  }

  auto subIndex(int i, int codim) const {
    if (codim == 0)
      return index();

    return grid_->parameterSpace().entityContainer_.template subIndexFromId<codim_>(entityInfo_.id, i, codim,
                                                                                    this->level());
  }

  template <typename = void>
  requires(codim_ == 0)
  auto& subId(int i, int codim) const {
    return grid_->parameterSpace().entityContainer_.subId(entityInfo_.id, i, codim);
  }

  // TODO isn't this the same as hostEntity()
  HostParameterSpaceGridEntity getHostEntity() const {
    if (codim_ == 0 or not isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
      return hostEntity_;
    DUNE_THROW(NotImplemented, "getHostEntity");
  }

  const EntityInfo& entityInfo() const {
    return entityInfo_;
  }

  const TrimInfo& trimData() const {
    assert(trimData_.has_value());
    return trimData_.value();
  }

public:
  [[nodiscard]] bool operator==(const TrimmedParameterSpaceGridEntity& other) const {
    if constexpr (codim_ == 0)
      return hostEntity_ == other.hostEntity_;
    else
      return entityInfo_.id == other.entityInfo_.id;
  }

  // returns true if father entity exists
  template <typename T = void>
  requires(codim_ == 0)
  [[nodiscard]] bool hasFather() const {
    return hostEntity_.hasFather();
  }

  // Create EntitySeed
  [[nodiscard]] ParameterSpaceGridEntitySeed seed() const {
    DUNE_THROW(NotImplemented, " seed");
    if constexpr (codim_ == 0)
      return hostEntity_.seed();
    return {};
  }

  // Level of this element
  [[nodiscard]] int level() const {
    return entityInfo_.lvl;
  }

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    if (isTrimmed() or Preferences::getInstance().reportTrimmedElementGeometryTypeAsNone())
      return GeometryTypes::none(mydimension);
    return GeometryTypes::cube(mydimension);
  }

  /** @brief The partition type for parallel computing */
  [[nodiscard]] PartitionType partitionType() const {
    if constexpr (codim_ == 0)
      return hostEntity_.partitionType();
    DUNE_THROW(NotImplemented, "partitionType not implemented for codim!=0 objects");
  }

  // Geometry of this entity
  [[nodiscard]] LocalParameterSpaceGeometry geometry() const {
    if (not isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
      return hostEntity_.geometry();
    if constexpr (codim_ == 1 or codim_ == 2) /* edge, vertex */ {
      return TrimmedParameterSpaceGeometry(trimData_->geometry.value());
    } else if constexpr (codim_ == 0) /* element */ {
      return TrimmedParameterSpaceGeometry(hostEntity_.geometry(), this);
    }
    __builtin_unreachable();
  }

  /** @brief Return the number of subEntities of codimension codim.
   */
  [[nodiscard]] unsigned int subEntities(unsigned int codim) const {
    if (codim_ == codim)
      return Dune::referenceElement<double, mydimension>(Dune::GeometryTypes::cube(mydimension)).size(codim);
    if constexpr (codim_ == 0) {
      if (not isTrimmed() or not Preferences::getInstance().reconstructTrimmedLocalGeometry())
        return hostEntity_.subEntities(codim);

      return trimData_.value().size(codim);
    }
    if constexpr (codim_ == 1) {
      return 2;
    }
    assert(codim < 3);
    DUNE_THROW(NotImplemented, "subEntities(): Codim has to be 2 or less");
  }

  /** @brief Provide access to sub entity i of given codimension. Entities
   *  are numbered 0 ... subEntities(cc)-1
   */
  template <int cc>
  requires(codim_ == 0)
  [[nodiscard]] TrimmedParameterSpaceGridEntity<cc, mydimension, GridImp> subEntity(int i) const {
    if constexpr (cc == 0)
      return *this;
    auto entity = grid_->parameterSpace().entityContainer_.template entity<cc>(subId(i, cc), this->level());
    return entity;
  }

  // First level intersection
  template <typename = void>
  requires(codim_ == 0)
  [[nodiscard]] decltype(auto) ilevelbegin() const {
    // if(trimData_)
    //   return trimData_.template ilevelbegin<codim_>(localId_);
    return hostEntity_.ilevelbegin();
  }

  // Reference to one past the last neighbor
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ilevelend() const {
    return hostEntity_.ilevelend();
  }

  // First leaf intersection
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ileafbegin() const {
    return hostEntity_.ileafbegin();
  }

  // Reference to one past the last leaf intersection
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ileafend() const {
    return hostEntity_.ileafend();
  }

  // returns true if Entity has NO children
  template <typename = void>
  requires(codim_ == 0)
  bool isLeaf() const {
    return hostEntity_.isLeaf();
  }

  // Inter-level access to father element on coarser grid.
  // Assumes that meshes are nested.
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) father() const {
    assert(entityInfo_.fatherId.has_value());
    return grid_->parameterSpace().entityContainer_.template entity<0>(entityInfo_.fatherId.value(), this->level());
    // return TrimmedParameterSpaceGridEntity(grid_, hostEntity_.father(),
    // grid_->parameterspace().entityContainer_.idToElementInfoMap.at( entityInfo_.fatherId.value()));
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
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) geometryInFather() const {
    return hostEntity_.geometryInFather();
  }

  // /** @brief Inter-level access to son elements on higher levels<=maxlevel.
  //  * This is provided for sparsely stored nested unstructured meshes.
  //  * Returns iterator to first son.
  //  */
  // template <typename = void>
  // requires(codim_ == 0) decltype(auto) hbegin(int maxLevel) const { return hostEntity_.hbegin(maxLevel); }
  //
  // //  Returns iterator to one past the last son
  // template <typename = void>
  // requires(codim_ == 0) decltype(auto) hend(int maxLevel) const { return hostEntity_.hend(maxLevel); }

  template <typename = void>
  requires(codim_ == 0)
  bool wasRefined() const {
    return hostEntity_.wasRefined();
  }

  template <typename = void>
  requires(codim_ == 0)

  bool mightBeCoarsened() const {
    return hostEntity_.mightBeCoarsened();
  }

  const auto& hostEntity() const {
    return hostEntity_;
  }

  template <typename = void>
  requires(codim_ == 0)
  unsigned int hostIndexInLvl() const {
    return entityInfo_.hostIndexInLvl;
  }

  QuadratureRule<double, dim> getQuadratureRule(const std::optional<int>& p_order = std::nullopt,
                                                const QuadratureType::Enum qt = QuadratureType::GaussLegendre) const {
    // The parameterSpace has linear geometry, so if no order is specified, we assume it to be myDimension
    int order = p_order.value_or(mydimension);
    if (not isTrimmed())
      return QuadratureRules<double, mydimension>::rule(this->type(), order, qt);

    return grid_->integrationRule()(*this, order, qt);
  }

  const GridImp* grid_;

private:
  EntityInfo entityInfo_;
  HostParameterSpaceGridEntity hostEntity_;
  std::optional<TrimInfo> trimData_{};
};

} // namespace Dune::IGA::DefaultParameterSpace
