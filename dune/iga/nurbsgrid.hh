// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <cstdlib>
#include <mutex>
#include <optional>
#include <utility>

#include "dune/iga/gridcapabilities.hh"
#include "dune/iga/nurbsalgorithms.hh"
#include "dune/iga/nurbsgridindexsets.hh"
#include "dune/iga/nurbsgridleafiterator.hh"
#include "dune/iga/nurbsgridtraits.hh"
#include "dune/iga/nurbsidset.hh"
#include "dune/iga/nurbsintersection.hh"
#include "dune/iga/nurbsleafgridview.hh"
#include "dune/iga/nurbslocalgeometry.hh"
#include "dune/iga/nurbspatch.hh"
#include "dune/iga/utils/concepts.hh"

namespace Dune::IGA {

  template <typename GridImpl, int griddim, std::integral auto... codim>
  std::tuple<std::vector<typename NurbsLeafGridViewTraits<GridImpl>::template Codim<codim>::Entity>...>
      gridEntityTupleGenerator(std::integer_sequence<std::common_type_t<decltype(codim)...>, codim...>);

  template <int dim, int dimworld, typename ScalarType>
  class NURBSGrid;

  template <int cd, typename GridImpl>
  class EntitySeedStruct {
   public:
    [[nodiscard]] bool isValid() const { return valid_; }
    static constexpr int codimension = cd;

   private:
    bool valid_{false};
    template <int codim, int dim, typename GridImpl1>
    friend class NURBSGridEntity;
    template <int dim, int dimworld, typename ScalarType>
    friend class NURBSGrid;

    int index_{-1};
  };

  template <int dim, int dimworld, typename ScalarType>
  struct NurbsGridFamily;

  /** \brief NURBS grid manager */
  template <int dim, int dimworld, typename ScalarType = double>
  class NURBSGrid : public Dune::Grid<dim, dimworld, ScalarType, NurbsGridFamily<dim, dimworld, ScalarType>> {
   public:
    static constexpr std::integral auto dimension      = dim;
    static constexpr std::integral auto dimensionworld = dimworld;
    using ctype                                        = ScalarType;

    using NURBSPatchDataType  = NURBSPatchData<dim, dimworld, ScalarType>;
    using ControlPointNetType = typename NURBSPatchDataType::ControlPointNetType;

    using GridFamily = NurbsGridFamily<dim, dimworld, ScalarType>;

    using Traits = typename GridFamily::Traits;
    template <int cd>
    using Codim = typename GridFamily::Traits::template Codim<cd>;

    using LeafIndexSet  = typename Traits::LeafIndexSet;
    using LevelIndexSet = typename Traits::LevelIndexSet;
    using GridView      = typename Traits::LeafGridView;
    using ElementEntity = typename Traits::template Codim<0>::Entity;
    using GlobalIdSet   = typename Traits::GlobalIdSet;
    NURBSGrid()         = default;

    /** \brief  constructor
     *
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order degree of the B-Spline structure for each dimension
     */
    NURBSGrid(const std::array<std::vector<double>, dim>& knotSpans, const ControlPointNetType& controlPoints,
              const std::array<int, dim>& order)
        : NURBSGrid(NURBSPatchData<dim, dimworld, ScalarType>(knotSpans, controlPoints, order)) {}

    explicit NURBSGrid(const NURBSPatchData<dim, dimworld, ScalarType>& nurbsPatchData,
                       std::optional<std::shared_ptr<TrimData>> _trimData = std::nullopt)
        : entityVector{std::make_unique<decltype(gridEntityTupleGenerator<NURBSGrid, dimension>(
            std::make_integer_sequence<int, dimension + 1>()))>()},
          coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          leafPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, ScalarType>>>(
              1, NURBSPatch<dim, dimworld, ScalarType>(currentPatchRepresentation_, _trimData))},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<const NURBSGrid>(*this, 0))},
          indexSet_{std::make_unique<NURBSGridLeafIndexSet<const NURBSGrid>>((this->leafGridView().impl()))},
          idSet_{std::make_unique<IgaIdSet<const NURBSGrid>>(this->leafGridView())},
          trimData_(_trimData) {
      static_assert(dim <= 3, "Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size() - nurbsPatchData.degree[0] - 1
                 == nurbsPatchData.controlPoints.strideSizes()[0]
             && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim > 1)
        assert(nurbsPatchData.knotSpans[1].size() - nurbsPatchData.degree[1] - 1
                   == nurbsPatchData.controlPoints.strideSizes()[1]
               && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim > 2)
        assert(nurbsPatchData.knotSpans[2].size() - nurbsPatchData.degree[2] - 1
                   == nurbsPatchData.controlPoints.strideSizes()[2]
               && "The size of the controlpoints and the knotvector strideSizes do not match in the third direction");
      // FIXME check sanity of knotvector and degree
      createEntities();
    }

    bool loadBalance() { return false; }

    void globalRefine(int refinementLevel) { globalRefine(refinementLevel, false); }
    void globalRefine(int refinementLevel, bool omitTrim) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots
            = generateRefinedKnots(currentPatchRepresentation_.knotSpans, refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
      updateStateAfterRefinement(omitTrim);
    }

    void globalRefineInDirection(const int dir, const int refinementLevel, bool omitTrim = false) {
      if (refinementLevel == 0) return;
      auto additionalKnots        = generateRefinedKnots(currentPatchRepresentation_.knotSpans, dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
      updateStateAfterRefinement(omitTrim);
    }

    void degreeElevateInDirection(const int dir, const int elevationFactor, bool omitTrim = false) {
      if (elevationFactor == 0) return;
      currentPatchRepresentation_ = degreeElevate(currentPatchRepresentation_, dir, elevationFactor);
      updateStateAfterRefinement(omitTrim);
    }

    void globalDegreeElevate(const int elevationFactor, bool omitTrim = false) {
      if (elevationFactor == 0) return;
      oldPatchRepresentation_.push_back(currentPatchRepresentation_);
      for (int refDirection = 0; refDirection < dim; ++refDirection)
        currentPatchRepresentation_ = degreeElevate(currentPatchRepresentation_, refDirection, elevationFactor);
      updateStateAfterRefinement(omitTrim);
    }

    void globalMultiRefine(const int global, const int uDir, const int vDir) {
      this->globalRefine(global, uDir > 0 or vDir > 0);
      this->globalRefineInDirection(0, uDir, vDir > 0);
      this->globalRefineInDirection(1, vDir);
    }

    int size(int codim) const { return leafPatches_.get()->front().size(codim); }

    const auto& patchData(int i = 0) const { return currentPatchRepresentation_; }
    const auto& lowerOrderPatchData(int i = 0) const { return oldPatchRepresentation_.at(i); }

    bool reportTrimError() const {
      for (const auto& patch : *leafPatches_)
        if (patch.reportTrimError()) return true;
      return false;
    }

    /** \brief returns the number of boundary segments within the macro grid */
    int numBoundarySegments() const { return getPatch().numBoundarySegments(); }
    int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const { return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark(const ElementEntity& element) const { return 0; }
    bool mark(int refCount, const ElementEntity& element) { return false; }

    template <typename Seed>
    typename Codim<Seed::codimension>::Entity entity(const Seed& seed) const {
      return leafGridView_->impl().template getEntity<Seed::codimension>(seed.impl().index_);
    }

    int size(const GeometryType& type) const { return this->leafGridView().size(type); }

    int size(int lvl, const GeometryType& type) const { return this->size(type); }

    const GlobalIdSet& globalIdSet() const { return *idSet_; }
    const LevelIndexSet& levelIndexSet(int lvl) const { return *indexSet_; }
    const LeafIndexSet& leafIndexSet() const { return *indexSet_; }

    int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

    // TODO Fix publicness and privateness
   private:
    void updateStateAfterRefinement(bool omitTrim = false) {
      leafPatches_->clear();
      leafPatches_->emplace_back(currentPatchRepresentation_, omitTrim ? std::nullopt : trimData_);

      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<const NURBSGrid>(*this, 0));
      updateEntities();
      indexSet_ = std::make_unique<NURBSGridLeafIndexSet<const NURBSGrid>>((this->leafGridView().impl()));
      idSet_    = std::make_unique<IgaIdSet<const NURBSGrid>>(this->leafGridView());
    }

    void createEntities() {
      for (int currentPatchId = 0; auto&& patch : *leafPatches_.get()) {
        Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dimension + 1>()), [&](const auto i) {
          std::get<i>(*entityVector.get()).reserve(patch.size(i));
          for (unsigned int j = 0; j < patch.size(i); ++j) {
            std::get<i>(*entityVector.get())
                .emplace_back(NURBSGridEntity<i, dimension, const NURBSGrid>(*leafGridView_, j, currentPatchId));
          }
        });
        ++currentPatchId;
      }
    }
    void updateEntities() {
      entityVector = std::make_unique<decltype(gridEntityTupleGenerator<NURBSGrid, dimension>(
          std::make_integer_sequence<int, dimension + 1>()))>();

      createEntities();
    }

   public:
    using EntityVectorType
        = decltype(gridEntityTupleGenerator<NURBSGrid, dimension>(std::make_integer_sequence<int, dimension + 1>()));
    std::unique_ptr<EntityVectorType> entityVector{};
    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, ScalarType> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, ScalarType> currentPatchRepresentation_;
    std::vector<NURBSPatchData<(size_t)dim, (size_t)dimworld, ScalarType>> oldPatchRepresentation_;
    std::shared_ptr<std::vector<NURBSPatch<dim, dimworld, ScalarType>>> leafPatches_;
    std::shared_ptr<GridView> leafGridView_;
    std::unique_ptr<NURBSGridLeafIndexSet<const NURBSGrid>> indexSet_;
    std::unique_ptr<IgaIdSet<const NURBSGrid>> idSet_;
    std::optional<std::shared_ptr<TrimData>> trimData_ = std::nullopt;

    auto& getPatch() const { return leafPatches_->front(); }
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  auto levelGridView(const NURBSGrid<dim, dimworld, ScalarType>& grid, int level) {
    return grid.levelGridView(level);
  }

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  auto leafGridView(const NURBSGrid<dim, dimworld, ScalarType>& grid) {
    return grid.leafGridView();
  }

  template <int dim, int dimworld, typename ScalarType>
  struct NurbsGridFamily {
    using GridImpl = Dune::IGA::NURBSGrid<dim, dimworld, ScalarType>;
    using Traits
        = NurbsGridTraits<dim, dimworld, GridImpl, NURBSGeometry, NURBSGridEntity, NURBSGridLeafIterator,
                          NURBSintersection, NURBSintersection, NURBSGridInterSectionIterator,
                          NURBSGridInterSectionIterator, NurbsHierarchicIterator, NURBSGridLeafIterator,
                          NURBSGridLeafIndexSet<const GridImpl>, NURBSGridLeafIndexSet<const GridImpl>,
                          IgaIdSet<const GridImpl>, int, IgaIdSet<const GridImpl>, int, Communication<No_Comm>,
                          NurbsLeafGridViewTraits, NurbsLeafGridViewTraits, EntitySeedStruct, NURBSLocalGeometry>;
  };

}  // namespace Dune::IGA
