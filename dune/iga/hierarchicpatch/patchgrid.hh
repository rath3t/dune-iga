// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The PatchGrid class
 */

#include "concepts.hh"
#include "enums.hh"
#include "gridcapabilities.hh"
#include "patchgridentity.hh"
#include "patchgridfactory.hh"
#include "patchgridfwd.hh"
#include "patchgridgeometry.hh"
#include "patchgridview.hh"
#include "traits.hh"

#include <map>
#include <string>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/common/grid.hh>

namespace Dune::Functions {
template <typename GV, typename ScalarType>
class NurbsPreBasis;
}
namespace Dune::IGANEW {

namespace IdentityTrim {
  template <int dim, int dimworld, typename ScalarType = double>
  struct PatchGridFamily;
}
namespace Impl {
  template <int dim>
  class NurbsPreBasisFactoryFromDegreeElevation;
}

// External forward declarations
template <class Grid>
struct HostGridAccess;

/**
 * @brief Provides a NURBS grid based on a single NURBS patch
 * @ingroup PatchGrid
 * @tparam dim The dimension of the grid
 * @tparam dimworld The dimension of the embedding space
 * @tparam TrimmerType_ The trimmer of the trimmer
 * @tparam ScalarType The type for the coordinates
 * Example Create surface in 3D space:
 * @code
 * using namespace Dune::IGANEW;
 *
 * // Define a NURBS patch data
 * const int dim = 2;
 * const int dimworld = 3;
 * const std::array<int, dim> order = {2, 2};
 * const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
 *
 * using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;
 * const std::vector<std::vector<ControlPoint>> controlPoints = {
 *     {{.p = {0, 0}, .w = 1}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
 *     {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 1}, {.p = {1, 0.5}, .w = 1}},
 *     {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}
 * };
 * std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
 * auto controlNet = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);
 *
 * // Create a PatchGrid from the NURBS patch data
 * PatchGrid<dim, dimworld> grid({knotSpans, controlNet, order});
 *
 * // Perform global refinement
 * grid.globalRefine(2);
 *
 * @endcode
 */
template <int dim, int dimworld, template <int, int, typename> typename GridFamily_ = IdentityTrim::PatchGridFamily,
          typename ScalarType = double>
class PatchGrid : public GridDefaultImplementation<dim, dimworld, ScalarType, GridFamily_<dim, dimworld, ScalarType>>
{
  friend class PatchGridLeafGridView<const PatchGrid>;
  friend class PatchGridLevelGridView<const PatchGrid>;
  friend GridFamily_<dim, dimworld, ScalarType>::LevelIndexSet;
  friend GridFamily_<dim, dimworld, ScalarType>::LeafIndexSet;
  friend GridFamily_<dim, dimworld, ScalarType>::GlobalIdSet;
  friend GridFamily_<dim, dimworld, ScalarType>::LocalIdSet;

  //! type of the used GridFamily for this grid
public:
  using GridFamily = GridFamily_<dim, dimworld, ScalarType>;

private:
  using LevelIndexSetImpl = GridFamily::LevelIndexSet;
  using LeafIndexSetImpl  = GridFamily::LeafIndexSet;
  using GlobalIdSetImpl   = GridFamily::GlobalIdSet;
  using LocalIdSetImpl    = GridFamily::LocalIdSet;
  template <int codim, PartitionIteratorType pitype>
  using LeafIteratorImpl = GridFamily::template LeafIterator<codim, pitype>;
  template <int codim, PartitionIteratorType pitype>
  using LevelIteratorImpl = GridFamily::template LevelIterator<codim, pitype>;
  friend GridFamily::HierarchicIterator;
  friend GridFamily::LevelIntersection;
  friend GridFamily::LeafIntersection;
  friend GridFamily::LeafIntersectionIterator;
  friend GridFamily::LevelIntersectionIterator;

  // friend class PatchGridLevelIntersection<const PatchGrid>;
  // friend class PatchGridLeafIntersectionIterator<const PatchGrid>;
  // friend class PatchGridLeafIntersection<const PatchGrid>;
  friend class PatchGridLevelGridView<PatchGrid>;
  friend class PatchGridLeafGridView<PatchGrid>;
  friend struct HostGridAccess<PatchGrid>;
  friend class GridFactory<PatchGrid>;

public:
  using Trimmer = typename GridFamily::Trimmer;
  //! The type used to store coordinates, inherited from the Trimmer
  using ctype = typename Trimmer::ctype;

  friend Trimmer;

private:
  friend class Impl::NurbsPreBasisFactoryFromDegreeElevation<dim>;
  friend class Functions::NurbsPreBasis<typename PatchGrid::LeafGridView, ctype>;
  friend class Functions::NurbsPreBasis<typename PatchGrid::LevelGridView, ctype>;

  template <int codim, PartitionIteratorType pitype, class GridImp_>
  friend class PatchGridLevelIterator;

  template <int codim, PartitionIteratorType pitype, class GridImp_>
  friend class PatchGridLeafIterator;

  template <int codim_, int dim_, class GridImp_>
  friend class PatchGridEntity;

public:
  using PatchTensorProductCoordinatesType =
      typename GeometryKernel::NURBSPatch<dim, dimworld, ctype>::TensorProductCoordinatesType;

  //**********************************************************
  // The Interface Methods
  //**********************************************************

  //! the Traits
  using Traits = typename GridFamily::Traits;

  using ParameterSpaceGrid = typename Trimmer::ParameterSpaceGrid;

  using PatchTrimData = typename GridFamily::TrimmerTraits::PatchTrimData;

  /** @brief Constructor
   *
   * @param hostgrid The host grid wrapped by the PatchGrid
   */
  explicit PatchGrid(const NURBSPatchData<dim, dimworld, ctype>& patchData,
                     const std::optional<PatchTrimData>& patchTrimData            = std::nullopt,
                     const typename GridFamily::TrimmerTraits::ParameterType& par = {})
      : patchGeometries_(1, GeometryKernel::NURBSPatch<dim, dimworld, ctype>(patchData)),
        trimmer_(std::make_unique<Trimmer>(*this, patchTrimData, par)) {
    patchGeometriesUnElevated = patchGeometries_;
  }

  PatchGrid& operator=(PatchGrid&& other) noexcept {
    patchGeometries_          = std::move(other.patchGeometries_);
    patchGeometriesUnElevated = std::move(other.patchGeometriesUnElevated);
    trimmer_                  = std::move(other.trimmer_);
    trimmer_->update(this);
    return *this;
  }

  /** @brief Return maximum level defined in this grid.
   *
   * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
   */
  [[nodiscard]] int maxLevel() const {
    return trimmer_->maxLevel();
  }

  //! Iterator to first entity of given codim on level
  template <int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin(int level) const {
    return LevelIteratorImpl<codim, All_Partition>(this, level);
  }

  //! one past the end on this level
  template <int codim>
  typename Traits::template Codim<codim>::LevelIterator lend(int level) const {
    return LevelIteratorImpl<codim, All_Partition>(this, level, true);
  }

  //! Iterator to first entity of given codim on level
  template <int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin(int level) const {
    return LevelIteratorImpl<codim, PiType>(this, level);
  }

  //! one past the end on this level
  template <int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend(int level) const {
    return LevelIteratorImpl<codim, PiType>(this, level, true);
  }

  //! Iterator to first leaf entity of given codim
  template <int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
    return LeafIteratorImpl<codim, All_Partition>(this);
  }

  //! one past the end of the sequence of leaf entities
  template <int codim>
  typename Traits::template Codim<codim>::LeafIterator leafend() const {
    return LeafIteratorImpl<codim, All_Partition>(this, true);
  }

  //! Iterator to first leaf entity of given codim
  template <int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
    return LeafIteratorImpl<codim, PiType>(this);
  }

  //! one past the end of the sequence of leaf entities
  template <int codim, PartitionIteratorType PiType>
  typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
    return LeafIteratorImpl<codim, PiType>(this, true);
  }

  /** @brief Number of grid entities per level and codim
   */
  [[nodiscard]] int size(int level, int codim) const {
    return levelIndexSet(level).size(codim);
  }

  /** @brief returns the number of boundary segments within the macro grid
   */
  [[nodiscard]] size_t numBoundarySegments() const {
    // @todo Trim this is wrong another trimmer functionality should care about this
    return trimmer_->parameterSpaceGrid().numBoundarySegments();
  }

  //! number of leaf entities per codim in this process
  [[nodiscard]] int size(int codim) const {
    return leafIndexSet().size(codim);
  }

  //! number of entities per level, codim and geometry type in this process
  int size(int level, GeometryType type) const {
    return levelIndexSet(level).size(type);
  }

  //! number of leaf entities per codim and geometry type in this process
  int size(GeometryType type) const {
    return leafIndexSet().size(type);
  }

  /** @brief Access to the GlobalIdSet */
  const typename Traits::GlobalIdSet& globalIdSet() const {
    return *trimmer_->globalIdSet_;
  }

  /** @brief Access to the LocalIdSet */
  const typename Traits::LocalIdSet& localIdSet() const {
    return *trimmer_->localIdSet_;
  }

  /** @brief Access to the LevelIndexSets */
  const typename Traits::LevelIndexSet& levelIndexSet(int level) const {
    if (level < 0 || level > maxLevel()) {
      DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    }
    return *trimmer_->levelIndexSets_[level];
  }

  /** @brief Access to the LeafIndexSet */
  const typename Traits::LeafIndexSet& leafIndexSet() const {
    return *trimmer_->leafIndexSet_;
  }

  /** @brief Create Entity from EntitySeed */
  template <class EntitySeed>
  typename Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const {
    return trimmer_->entity(seed);
  }

  /** @name Grid Refinement Methods */
  /*@{*/

  /** global refinement
   */
  void globalRefine(int refCount) {
    for (int i = 0; i < refCount; ++i) {
      const auto& finestPatchData = patchGeometries_.back().patchData();
      auto newfinestPatchData     = finestPatchData;

      auto& olduniqueKnotVector = patchGeometries_.back().uniqueKnotVector();
      std::array<std::vector<ctype>, dim> newUniqueKnotVecs;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots = Splines::generateRefinedKnots(finestPatchData.knotSpans, refDirection, 1);
        newfinestPatchData   = Splines::knotRefinement<dim>(newfinestPatchData, additionalKnots, refDirection);
        auto& newKnotVec     = newUniqueKnotVecs[refDirection];

        // compute new unique Knot vectors from refinement
        newKnotVec.resize(olduniqueKnotVector[refDirection].size() + additionalKnots.size());
        std::ranges::merge(olduniqueKnotVector[refDirection], additionalKnots, std::begin(newKnotVec));
      }

      patchGeometries_.emplace_back(newfinestPatchData, newUniqueKnotVecs);
      patchGeometriesUnElevated.emplace_back(patchGeometries_.back());
    }
    // Here it is exploited that the knot refinem above is in sync with the globaRefine of the trimmer
    trimmer_->globalRefine(refCount);
  }

  /**
   * @brief Returns the coordinates of the untrimmed tensor product grid
   * @param lvl The grid level of the requested coordinates
   * \return the array of the coordinates
   */
  const PatchTensorProductCoordinatesType& tensorProductCoordinates(int lvl) const {
    return patchGeometries_[lvl].uniqueKnotVector();
  }

  /**
   * @brief This refines the grid in one specific direction, this resets all multilevel structure, since YaspGrid does
   * not support this!
   * @param refines how often should the element be
   * refined in the given directions. This splits the element in half, quarters ,... in the given direction
   */
  void globalRefineInDirection(const std::array<int, dim>& refines) {
    const auto& finestPatchData = patchGeometries_.back().patchData();
    auto newfinestPatchData     = finestPatchData;
    for (int dir = 0; auto refinesInDirection : refines) {
      if (refinesInDirection == 0) {
        ++dir;
        continue;
      }
      auto additionalKnots = Splines::generateRefinedKnots(finestPatchData.knotSpans, dir, refinesInDirection);
      newfinestPatchData   = Splines::knotRefinement<dim>(newfinestPatchData, additionalKnots, dir);
      ++dir;
    }
    auto newGrid = PatchGrid(newfinestPatchData);
    *this        = std::move(newGrid);
  }

  /**
   * @brief Increases the polynomial degree of the NURBS geometry of the given level
   * @param elevationFactors the increase in polynomial degree of the underlying NURBS
   * @param lvl the level where the degree elevation should be performed
   */
  void degreeElevate(const std::array<int, dim>& elevationFactors, int lvl) {
    if (lvl > maxLevel() and lvl >= 0)
      DUNE_THROW(Dune : RangeError, "This level does not exist");
    auto& patchData                = patchGeometries_[lvl].patchData();
    patchGeometriesUnElevated[lvl] = GeometryKernel::NURBSPatch<dim, dimworld, ctype>(patchData);
    for (int dir = 0; auto elevatesInDirection : elevationFactors) {
      if (elevatesInDirection == 0) {
        ++dir;
        continue;
      }
      patchData = Splines::degreeElevate(patchData, dir, elevatesInDirection);
      ++dir;
    }
    patchGeometries_[lvl] = GeometryKernel::NURBSPatch<dim, dimworld, ctype>(patchData);
  }

  /**
   * @brief Increases the polynomial degree of the NURBS geometry of the given level
   * @param elevationFactors the increase in polynomial degree of the underlying NURBS
   */
  void degreeElevateOnAllLevels(const std::array<int, dim>& elevationFactors) {
    for (int lvl = 0; lvl < maxLevel() + 1; ++lvl)
      degreeElevate(elevationFactors, lvl);
  }

  /** @brief Mark entity for refinement
   *
   * This only works for entities of codim 0.
   * The parameter is currently ignored
   *
   * \return <ul>
   * <li> true, if marking was successful </li>
   * <li> false, if marking was not possible </li>
   * </ul>
   */
  bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
    // @todo trim this does not do the right thing! the knotspans should also be aware of this change
    return false; // trimmer_->parameterSpaceGrid().mark(refCount, getHostEntity<0>(e));
  }

  /** @brief Return refinement mark for entity
   *
   * \return refinement mark (1,0,-1)
   */
  int getMark(const typename Traits::template Codim<0>::Entity& e) const {
    return 0; // trimmer_->parameterSpaceGrid().getMark(getHostEntity<0>(e));
  }

  /** @brief returns true, if at least one entity is marked for adaption */
  bool preAdapt() {
    return trimmer_->paramterSpaceGrid().preAdapt();
  }

  //! Triggers the grid refinement process
  bool adapt() {
    return trimmer_->paramterSpaceGrid().adapt();
  }

  /** @brief Clean up refinement markers */
  void postAdapt() {
    return trimmer_->paramterSpaceGrid().postAdapt();
  }

  /*@}*/

  /** @brief Size of the overlap on the leaf level */
  unsigned int overlapSize(int codim) const {
    return trimmer_->parameterSpaceGrid().leafGridView().overlapSize(codim);
  }

  /** @brief Size of the ghost cell layer on the leaf level */
  unsigned int ghostSize(int codim) const {
    return trimmer_->parameterSpaceGrid().leafGridView().ghostSize(codim);
  }

  /** @brief Size of the overlap on a given level */
  unsigned int overlapSize(int level, int codim) const {
    return trimmer_->parameterSpaceGrid().levelGridView(level).overlapSize(codim);
  }

  /** @brief Size of the ghost cell layer on a given level */
  unsigned int ghostSize(int level, int codim) const {
    return trimmer_->parameterSpaceGrid().levelGridView(level).ghostSize(codim);
  }

#if 0
    /** @brief Distributes this grid over the available nodes in a distributed machine
     *
     * @param minlevel The coarsest grid level that gets distributed
     * @param maxlevel does currently get ignored
     */
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "PatchGrid::loadBalance()");
    }
#endif

  /** @brief dummy communication */
  const Communication<No_Comm>& comm() const {
    return ccobj;
  }

  /** @brief Communicate data of level gridView */
  template <class DataHandle>
  void communicate(DataHandle& handle, InterfaceType iftype, CommunicationDirection dir, int level) const {
    trimmer_->parameterSpaceGrid().levelGridView(level).communicate(handle, iftype, dir);
  }

  /** @brief Communicate data of leaf gridView */
  template <class DataHandle>
  void communicate(DataHandle& handle, InterfaceType iftype, CommunicationDirection dir) const {
    trimmer_->parameterSpaceGrid().leafGridView().communicate(handle, iftype, dir);
  }

  // **********************************************************
  // End of Interface Methods
  // **********************************************************

  //! Returns the hostgrid this PatchGrid lives in
  const ParameterSpaceGrid& parameterSpaceGrid() const {
    return trimmer_->parameterSpaceGrid();
  }
  ParameterSpaceGrid& parameterSpaceGrid() {
    return trimmer_->parameterSpaceGrid();
  }

  //! Returns the hostgrid entity encapsulated in given PatchGrid entity
  template <int codim>
  requires(GridFamily::template hasHostEntity<codim>)
  const typename GridFamily::TrimmerTraits::template Codim<codim>::ParameterSpaceGridEntity& getHostEntity(
      const typename Traits::template Codim<codim>::Entity& e) const {
    return e.impl().getHostEntity();
  }

  auto untrimmedElementNumbers(int lvl) const {
    return patchGeometries_[lvl].numberOfSpans();
  }

  const auto& trimmer() const {
    return *trimmer_;
  }
  const auto& patchGeometry(int i) const {
    return patchGeometries_.at(i);
  }

private:
  PatchGrid() = default;
  std::vector<GeometryKernel::NURBSPatch<dim, dimworld, ctype>> patchGeometries_;
  std::vector<GeometryKernel::NURBSPatch<dim, dimworld, ctype>> patchGeometriesUnElevated;

  std::unique_ptr<Trimmer> trimmer_;

private:
  //! @todo Please doc me !
  Communication<No_Comm> ccobj;

}; // end Class PatchGrid

} // namespace Dune::IGANEW
