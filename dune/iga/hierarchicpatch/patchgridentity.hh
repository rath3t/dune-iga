// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

/** \file
 * @brief The PatchGridEntity class
 */

#include "patchgridgeometry.hh"

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include <dune/iga/integrationrules/simplexintegrationrulegenerator.hh>
#include <dune/iga/utils/igahelpers.hh>

namespace Dune::IGA {

// Forward declarations

template <int codim, int dim, class GridImp>
class PatchGridEntity;

template <class GridImp>
class PatchGridLeafIntersectionIterator;

// External forward declarations
template <class Grid>
struct HostGridAccess;

/** @brief The implementation of entities in a PatchGrid
 *   @ingroup PatchGrid
 *
 *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
 *  An entity of codimension c in dimension d is a d-c dimensional object.
 *
 */
template <int codim, int dim, class GridImp>
class PatchGridEntity : public EntityDefaultImplementation<codim, dim, GridImp, PatchGridEntity>
{
  template <class GridImp_>
  friend class PatchGridLevelIndexSet;

  template <class GridImp_>
  friend class PatchGridLeafIndexSet;

  template <class GridImp_>
  friend class PatchGridLocalIdSet;

  template <class GridImp_>
  friend class PatchGridGlobalIdSet;

  friend struct HostGridAccess<typename std::remove_const_t<GridImp>>;

  using ParameterSpace           = typename GridImp::ParameterSpace;
  using ParameterSpaceGridEntity = typename ParameterSpace::template Codim<codim>::ParameterSpaceGridEntity;

public:
  typedef typename GridImp::ctype ctype;

  typedef typename GridImp::template Codim<codim>::Geometry Geometry;

  // The type of the EntitySeed interface class
  typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

  PatchGridEntity()
      : patchGrid_(nullptr) {}

  PatchGridEntity(const GridImp* patchGrid, const ParameterSpaceGridEntity& hostEntity)
      : localEntity_(hostEntity),
        patchGrid_(patchGrid) {}

  PatchGridEntity(const GridImp* patchGrid, ParameterSpaceGridEntity&& hostEntity)
      : localEntity_(std::move(hostEntity)),
        patchGrid_(patchGrid) {}

  PatchGridEntity(const PatchGridEntity& original)
      : localEntity_(original.localEntity_),
        patchGrid_(original.patchGrid_) {}

  PatchGridEntity(PatchGridEntity&& original) noexcept
      : localEntity_(std::move(original.localEntity_)),
        patchGrid_(original.patchGrid_) {}

  PatchGridEntity& operator=(const PatchGridEntity& original) {
    if (this != &original) {
      patchGrid_   = original.patchGrid_;
      localEntity_ = original.localEntity_;
    }
    return *this;
  }

  PatchGridEntity& operator=(PatchGridEntity&& original) noexcept {
    if (this != &original) {
      patchGrid_   = original.patchGrid_;
      localEntity_ = std::move(original.localEntity_);
    }
    return *this;
  }

  bool equals(const PatchGridEntity& other) const {
    return getLocalEntity() == other.getLocalEntity();
  }

  // returns true if father entity exists
  bool hasFather() const {
    return localEntity_.hasFather();
  }

  /** @brief returns  EntitySeed */
  EntitySeed seed() const {
    return patchGrid_->parameterSpace().seed(*this);
  }

  /** @brief  level of this element */
  int level() const {
    return localEntity_.level();
  }

  /** @brief The partition type for parallel computing
   */
  PartitionType partitionType() const {
    return localEntity_.partitionType();
  }

  /** @brief Return the number of subEntities of codimension codim.
   */
  unsigned int subEntities(unsigned int cc) const {
    return localEntity_.subEntities(cc);
  }

  // geometry of this entity
  Geometry geometry() const {
    auto geo = typename Geometry::Implementation(
        localEntity_.geometry(),
        patchGrid_->patchGeometries_[this->level()].template localView<codim, ParameterSpace>());
    return Geometry(geo);
  }

  const auto& getLocalEntity() const {
    return localEntity_;
  }

  bool isTrimmed() const {
    return localEntity_.entityInfo().trimmed;
  }

private:
  ParameterSpaceGridEntity localEntity_;
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
class PatchGridEntity<0, dim, GridImp> : public EntityDefaultImplementation<0, dim, GridImp, PatchGridEntity>
{
  friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;

public:
  typedef typename GridImp::ctype ctype;
  using ParameterSpace = typename GridImp::ParameterSpace;
  // The codimension of this entitypointer wrt the host grid
  constexpr static int dimension       = GridImp::dimension;
  constexpr static int CodimInHostGrid = GridImp::ParameterSpaceGrid::dimension - dimension;
  constexpr static int dimworld        = GridImp::dimensionworld;

  // equivalent entity in the host grid
  using ParameterSpaceGridEntity = typename ParameterSpace::template Codim<0>::ParameterSpaceGridEntity;

  typedef typename GridImp::template Codim<0>::Geometry Geometry;

  typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

  // The Iterator over intersections on this level
  typedef typename GridImp::GridFamily::LevelIntersectionIterator LevelIntersectionIterator;

  // The Iterator over intersections on the leaf level
  typedef typename GridImp::GridFamily::LeafIntersectionIterator LeafIntersectionIterator;

  // Iterator over descendants of the entity
  typedef typename GridImp::GridFamily::HierarchicIterator HierarchicIterator;

  // The type of the EntitySeed interface class
  typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;
  typedef typename GridImp::ParameterSpace::ParameterSpaceTraits::template Codim<0>::ParameterSpaceGridEntitySeed
      ParameterSpaceGridEntitySeed;

  PatchGridEntity()
      : patchGrid_(nullptr) {}

  PatchGridEntity(const GridImp* patchGrid, const ParameterSpaceGridEntity& hostEntity)
      : localEntity_(hostEntity),
        patchGrid_(patchGrid) {}

  PatchGridEntity(const GridImp* patchGrid, ParameterSpaceGridEntity&& hostEntity)
      : localEntity_(std::move(hostEntity)),
        patchGrid_(patchGrid) {}

  PatchGridEntity(const PatchGridEntity& original)
      : localEntity_(original.localEntity_),
        patchGrid_(original.patchGrid_) {}

  PatchGridEntity(PatchGridEntity&& original) noexcept
      : localEntity_(std::move(original.localEntity_)),
        patchGrid_(original.patchGrid_) {}

  PatchGridEntity& operator=(const PatchGridEntity& original) {
    if (this != &original) {
      patchGrid_   = original.patchGrid_;
      localEntity_ = original.localEntity_;
    }
    return *this;
  }

  PatchGridEntity& operator=(PatchGridEntity&& original) noexcept {
    if (this != &original) {
      patchGrid_   = original.patchGrid_;
      localEntity_ = std::move(original.localEntity_);
    }
    return *this;
  }

  [[nodiscard]] bool equals(const PatchGridEntity& other) const {
    return localEntity_ == other.localEntity_;
  }

  /** @brief  returns true if father entity exists */
  [[nodiscard]] bool hasFather() const {
    return localEntity_.hasFather();
  }

  /** @brief Return EntitySeed */
  [[nodiscard]] EntitySeed seed() const {
    return patchGrid_->parameterSpace().seed(*this);
  }

  /** @brief  Level of this element */
  [[nodiscard]] int level() const {
    return getLocalEntity().level();
  }

  /** @brief The partition type for parallel computing */
  [[nodiscard]] PartitionType partitionType() const {
    return getLocalEntity().partitionType();
  }

  /** @brief  Geometry of this entity */
  [[nodiscard]] Geometry geometry() const {
    static_assert(
        std::is_same_v<
            decltype(patchGrid_->patchGeometries_[this->level()].template localView<0, ParameterSpace>()),
            typename GeometryKernel::NURBSPatch<dim, dimworld, ctype>::template GeometryLocalView<0, ParameterSpace>>);
    auto geo = typename Geometry::Implementation(
        localEntity_.geometry(), patchGrid_->patchGeometries_[this->level()].template localView<0, ParameterSpace>());
    return Geometry(geo);
  }

  /** @brief Return the number of subEntities of codimension codim. */
  [[nodiscard]] unsigned int subEntities(unsigned int codim) const {
    return localEntity_.subEntities(codim);
  }

  /** @brief Provide access to sub entity i of given codimension. Entities
   *  are numbered 0 ... subEntities(cc)-1
   */
  template <int cc>
  [[nodiscard]] typename GridImp::template Codim<cc>::Entity subEntity(int i) const {
    return PatchGridEntity<cc, dim, GridImp>(patchGrid_, localEntity_.template subEntity<cc>(i));
  }

  // First level intersection
  [[nodiscard]] LevelIntersectionIterator ilevelbegin() const {
    return patchGrid_->parameterSpace().ilevelbegin(*this);
  }

  // Reference to one past the last neighbor
  LevelIntersectionIterator ilevelend() const {
    return patchGrid_->parameterSpace().ilevelend(*this);
  }

  // First leaf intersection
  LeafIntersectionIterator ileafbegin() const {
    return patchGrid_->parameterSpace().ileafbegin(*this);
  }

  // Reference to one past the last leaf intersection
  LeafIntersectionIterator ileafend() const {
    return patchGrid_->parameterSpace().ileafend(*this);
  }

  // returns true if Entity has NO children
  bool isLeaf() const {
    return localEntity_.isLeaf();
  }

  // Inter-level access to father element on coarser grid.
  // Assumes that meshes are nested.
  typename GridImp::template Codim<0>::Entity father() const {
    return PatchGridEntity(patchGrid_, localEntity_.father());
  }

  /** @brief Location of this element relative to the reference element element of the father.
   * This is sufficient to interpolate all dofs in conforming case.
   * Nonconforming may require access to neighbors of father and
   * computations with local coordinates.
   * On the fly case is somewhat inefficient since dofs are visited several times.
   * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
   * implementation of numerical algorithms is only done for simple discretizations.
   * Assumes that meshes are nested.
   */
  LocalGeometry geometryInFather() const {
    return LocalGeometry(typename LocalGeometry::Implementation(localEntity_.geometryInFather()));
  }

  /** @brief Inter-level access to son elements on higher levels<=maxlevel.
   * This is provided for sparsely stored nested unstructured meshes.
   * Returns iterator to first son.
   */
  HierarchicIterator hbegin(int maxLevel) const {
    return HierarchicIterator(patchGrid_, *this, maxLevel);
  }

  // Returns iterator to one past the last son
  HierarchicIterator hend(int maxLevel) const {
    return HierarchicIterator(patchGrid_, *this, maxLevel, true);
  }

  bool wasRefined() const {
    if (patchGrid_->adaptationStep != GridImp::adaptDone)
      return false;

    int level = this->level();
    int index = patchGrid_->levelIndexSet(level).index(*this);
    return patchGrid_->refinementMark_[level][index];
  }

  bool mightBeCoarsened() const {
    return true;
  }

  bool isTrimmed() const {
    if constexpr (requires { localEntity_.isTrimmed(); })
      return localEntity_.isTrimmed();
    else
      return false;
  }

  auto trimData() const {
    return patchGrid_->trimData(*this);
  }

  const auto& getLocalEntity() const {
    return localEntity_;
  }

  QuadratureRule<double, dim> getQuadratureRule(const std::optional<int>& p_order = std::nullopt,
                                                const QuadratureType::Enum qt     = QuadratureType::GaussLegendre) const
  requires(not ParameterSpace::isAlwaysTrivial and dimension == 2)
  {
    auto degree = patchGrid_->patchGeometry(this->level()).degree();
    int order   = p_order.value_or(dimension * *std::ranges::max_element(degree));
    if (not isTrimmed())
      return QuadratureRules<double, dimension>::rule(this->type(), order, qt);

    return patchGrid_->integrationRule()(localEntity_, order, qt);
  }

private:
  ParameterSpaceGridEntity localEntity_;
  const GridImp* patchGrid_;
};

template <int cd, int dim, int dimworld, typename ScalarType, template <int, int, typename> typename GridFamily,
          template <int, int, class> class PatchGridEntity>
auto referenceElement(const PatchGridEntity<cd, dim, const PatchGrid<dim, dimworld, GridFamily, ScalarType>>& entity) {
  return GridFamily<dim, dimworld, ScalarType>::ParameterSpace::referenceElement(entity);
}

template <int cd, int dim, int dimworld, typename ScalarType, template <int, int, typename> typename GridFamily,
          template <int, int, class> class PatchGridEntity>
auto referenceElement(
    const Entity<cd, dim, const PatchGrid<dim, dimworld, GridFamily, ScalarType>, PatchGridEntity>& entity) {
  return referenceElement(entity.impl());
}

} // namespace Dune::IGA
