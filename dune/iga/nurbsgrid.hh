// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once
#include <optional>
#include <utility>

#include <dune/iga/concepts.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsgridtraits.hh>
#include <dune/iga/nurbsidset.hh>
#include <dune/iga/nurbsintersection.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbslocalgeometry.hh>
#include <dune/iga/nurbspatch.hh>
#include <stdlib.h>
#include <mutex>

std::once_flag onceFlag;

namespace Dune::IGA {
  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
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
    template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
    friend class NURBSGrid;

    int index_{-1};
  };

  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily;

  /** \brief NURBS grid manager */
  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  class NURBSGrid : public Dune::Grid<dim, dimworld, typename NurbsGridLinearAlgebraTraitsImpl::value_type,
                                      NurbsGridFamily<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>> {
   public:
    using LinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;

    static constexpr std::integral auto dimension      = dim;
    static constexpr std::integral auto dimensionworld = dimworld;
    using ctype                                        = typename LinearAlgebraTraits::value_type;

    using ControlPointNetType =
        typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType;

    using GridFamily = NurbsGridFamily<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;

    using Traits = typename GridFamily::Traits;
    template <int cd>
    using Codim = typename GridFamily::Traits::template Codim<cd>;

    using LeafIndexSet  = typename Traits::LeafIndexSet;
    using LevelIndexSet = typename Traits::LevelIndexSet;
    using GridView      = typename Traits::LeafGridView;
    using ElementEntity = typename Traits::template Codim<0>::Entity;
    using GlobalIdSet   = typename Traits::GlobalIdSet;
    NURBSGrid()         = default;

    explicit NURBSGrid(const NURBSPatchData<dim, dimworld, LinearAlgebraTraits>& nurbsPatchData,
                       std::optional<std::shared_ptr<TrimData>> _trimData = std::nullopt)
        : coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()},
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this))},
          idSet_{std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView())},
          trimData_(_trimData) {
      static_assert(dim <= 3, "Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size() - nurbsPatchData.degree[0] - 1 == nurbsPatchData.controlPoints.size()[0]
             && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim > 1)
        assert(nurbsPatchData.knotSpans[1].size() - nurbsPatchData.degree[1] - 1
                   == nurbsPatchData.controlPoints.size()[1]
               && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim > 2)
        assert(nurbsPatchData.knotSpans[2].size() - nurbsPatchData.degree[2] - 1
                   == nurbsPatchData.controlPoints.size()[2]
               && "The size of the controlpoints and the knotvector size do not match in the third direction");
      // FIXME check sanity of knotvector and degree
      finestPatches_->emplace_back(currentPatchRepresentation_, std::move(_trimData));
      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      idSet_        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
      indexdSet_    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());
    }

    /** \brief  constructor
     *
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order degree of the B-Spline structure for each dimension
     */
    NURBSGrid(const std::array<std::vector<double>, dim>& knotSpans, const ControlPointNetType& controlPoints,
              const std::array<int, dim>& order)
        : coarsestPatchRepresentation_{NURBSPatchData<dim, dimworld, LinearAlgebraTraits>(knotSpans, controlPoints,
                                                                                          order)},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(currentPatchRepresentation_, *this))},
          idSet_{std::make_unique<Dune::IGA::IgaIdSet<NURBSGrid>>(*this)},
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())} {}

    void globalRefine(int refinementLevel) { globalRefine(refinementLevel, false); }

    void globalRefine(int refinementLevel, bool omitTrim) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots
            = generateRefinedKnots(currentPatchRepresentation_.knotSpans, refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
      updateStateAfterRefinement(omitTrim);
      std::call_once ( onceFlag, [ ]{ putenv(strdup("ALUGRID_VERBOSITY_LEVEL=0")); } );
    }

    void globalRefineInDirection(const int dir, const int refinementLevel, bool omitTrim = false) {
      if (refinementLevel == 0) return;
      auto additionalKnots        = generateRefinedKnots(currentPatchRepresentation_.knotSpans, dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
      updateStateAfterRefinement(omitTrim);
    }

    void globalMultiRefine(const int global, const int uDir, const int vDir) {
      this->globalRefine(global, uDir > 0 or vDir > 0);
      this->globalRefineInDirection(0, uDir, vDir > 0);
      this->globalRefineInDirection(1, vDir);
    }

    [[nodiscard]] int size(int codim) const { return finestPatches_.get()->front().size(codim); }

    [[nodiscard]] bool reportTrimError() const {
      for (const auto& patch : *finestPatches_)
        if (patch.reportTrimError()) return true;
      return false;
    }

    /** \brief returns the number of boundary segments within the macro grid */
    [[nodiscard]] int numBoundarySegments() const { return getPatch().numBoundarySegments(); }
    [[nodiscard]] int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const { return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark(const ElementEntity& element) const { return 0; }
    bool mark(int refCount, const ElementEntity& element) { return false; }

    template <typename Seed>
    typename Codim<Seed::codimension>::Entity entity(const Seed& seed) const {
      return leafGridView_->impl().template getEntity<Seed::codimension>(seed.impl().index_);
    }

    [[nodiscard]] int size(const GeometryType& type) const { return this->leafGridView().size(type); }

    int size(int lvl, const GeometryType& type) const { return this->size(type); }

    const GlobalIdSet& globalIdSet() const { return *idSet_; }
    const LevelIndexSet& levelIndexSet(int lvl) const { return *indexdSet_; }
    const LeafIndexSet& leafIndexSet() const { return *indexdSet_; }

    [[nodiscard]] int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    [[nodiscard]] const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

    // TODO Fix publicness and privateness
   private:

    void updateStateAfterRefinement(bool omitTrim = false) {
      finestPatches_->clear();
      finestPatches_->emplace_back(currentPatchRepresentation_, omitTrim ? std::nullopt : trimData_);

      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      indexdSet_    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());
      idSet_        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
    }

   public:
    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> currentPatchRepresentation_;
    std::shared_ptr<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>> finestPatches_;
    std::shared_ptr<GridView> leafGridView_;
    std::unique_ptr<LeafIndexSet> indexdSet_;
    std::unique_ptr<IgaIdSet<NURBSGrid>> idSet_;
    std::optional<std::shared_ptr<TrimData>> trimData_ = std::nullopt;

    auto& getPatch() const { return finestPatches_->front(); }
  };

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto levelGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid, int level) {
    return grid.levelGridView(level);
  }

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto leafGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid) {
    return grid.leafGridView();
  }

  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily {
    using GridImpl = Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using Traits   = NurbsGridTraits<dim, dimworld, GridImpl, NURBSGeometry, NURBSGridEntity, NURBSGridLeafIterator,
                                   NURBSintersection, NURBSintersection, NURBSGridInterSectionIterator,
                                   NURBSGridInterSectionIterator, NurbsHierarchicIterator, NURBSGridLeafIterator,
                                   NURBSGridLeafIndexSet<GridImpl>, NURBSGridLeafIndexSet<GridImpl>, IgaIdSet<GridImpl>,
                                   int, IgaIdSet<GridImpl>, int, Communication<No_Comm>, NurbsLeafGridViewTraits,
                                   NurbsLeafGridViewTraits, EntitySeedStruct, NURBSLocalGeometry>;
  };

}  // namespace Dune::IGA
