// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once
#include <clipper2/clipper.h>
#include <utility>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/grid/yaspgrid.hh>
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
#include <dune/iga/nurbspatchgeometry.h>
#include <dune/iga/nurbstrimboundary.hh>
#include <dune/iga/nurbstrimfunctionality.hh>
#include <dune/iga/nurbstrimmedpatch.hh>
#include <dune/iga/reconstructedgridhandler.hh>

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

    // Boundaries can be given as an optional
    using BoundariesOpt = std::optional<std::vector<Boundary>>;

    explicit NURBSGrid(const NURBSPatchData<dim, dimworld, LinearAlgebraTraits>& nurbsPatchData,
                       BoundariesOpt _boundaries = std::nullopt)
        : coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()},
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this))},
          idSet_{std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView())},
          boundaries(std::move(_boundaries)),
          nurbsSurface{Dune::IGA::Nurbs<dim>(currentPatchRepresentation_)} {
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
      finestPatches_->emplace_back(currentPatchRepresentation_);
      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      idSet_        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
      indexdSet_    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());

      trimFlags = std::vector<ElementTrimFlag>(size(0));
      std::fill(trimFlags.begin(), trimFlags.end(), ElementTrimFlag::full);

      if (boundaries.has_value()) trimElement();
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
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())},
          nurbsSurface{Dune::IGA::Nurbs<dim>(currentPatchRepresentation_)} {}

    void globalRefine(int refinementLevel) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots
            = generateRefinedKnots(currentPatchRepresentation_.knotSpans, refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
      updateStateAfterRefinement();
    }

    void globalRefineInDirection(const int dir, const int refinementLevel) {
      if (refinementLevel == 0) return;
      auto additionalKnots        = generateRefinedKnots(currentPatchRepresentation_.knotSpans, dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
      updateStateAfterRefinement();
    }

    [[nodiscard]] int size(int codim) const { return finestPatches_.get()->front().size(codim); }

    /** \brief returns the number of boundary segments within the macro grid */
    [[nodiscard]] int numBoundarySegments() const {
      if constexpr (dimension == 1)
        return 2;
      else if constexpr (dimension == 2)
        return (finestPatches_.get()->front().validKnotSize()[0] + finestPatches_.get()->front().validKnotSize()[1])
               * 2;
      else if constexpr (dimension == 3)
        return 2
               * (finestPatches_.get()->front().validKnotSize()[0] * finestPatches_.get()->front().validKnotSize()[1]
                  + finestPatches_.get()->front().validKnotSize()[1] * finestPatches_.get()->front().validKnotSize()[2]
                  + finestPatches_.get()->front().validKnotSize()[0]
                        * finestPatches_.get()->front().validKnotSize()[2]);
      __builtin_unreachable();
    }
    [[nodiscard]] int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const { return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark(const ElementEntity& element) const { return 0; }
    bool mark(int refCount, const ElementEntity& element) { return false; }

    template <typename Seed>
    typename Codim<Seed::codimension>::Entity entity(const Seed& seed) const {
      return leafGridView_->impl().template getEntity<Seed::codimension>(seed.impl().index_);
    }

    [[nodiscard]] int size(const GeometryType& type) const {
      if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1)
          || type == Dune::GeometryTypes::cube(2) || type == Dune::GeometryTypes::cube(3))
        return this->leafGridView().size(dimension - type.dim());
      else
        return 0;
    }
    int size(int lvl, const GeometryType& type) const { return this->size(type); }

    const GlobalIdSet& globalIdSet() const { return *idSet_; }
    const LevelIndexSet& levelIndexSet(int lvl) const { return *indexdSet_; }
    const LeafIndexSet& leafIndexSet() const { return *indexdSet_; }

    [[nodiscard]] int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    [[nodiscard]] const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

    // Get if available the reconstructedGridView for a particular Element (only valid for dim == 2 (and also only for
    // dimwold == 2))
    std::optional<typename ReconstructedGridHandler<dimworld>::GridView> getReconstructedGridViewForTrimmedElement(
        int index) {
      if constexpr (dimworld == 2 && dim == 2) {
        if ((boundaries.has_value()) && (trimFlags[index] == ElementTrimFlag::trimmed))
          return std::make_optional<typename ReconstructedGridHandler<dimworld>::GridView>(
              trimResultMap[index]->grid->leafGridView());
        else
          return std::nullopt;
      } else
        return std::nullopt;
    }

   private:
    void updateStateAfterRefinement() {
      finestPatches_.get()->front() = NURBSPatch<dim, dimworld, LinearAlgebraTraits>(currentPatchRepresentation_);
      leafGridView_                 = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      indexdSet_                    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());
      idSet_                        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());

      if (boundaries.has_value()) {
        nurbsSurface = Dune::IGA::Nurbs<dim>(currentPatchRepresentation_);
        trimResultMap.clear();
        trimElement();
      }
    }

    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> currentPatchRepresentation_;
    std::shared_ptr<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>> finestPatches_;
    std::shared_ptr<GridView> leafGridView_;
    std::unique_ptr<LeafIndexSet> indexdSet_;
    std::unique_ptr<IgaIdSet<NURBSGrid>> idSet_;

    Dune::IGA::Nurbs<dim> nurbsSurface;
    std::vector<ElementTrimFlag> trimFlags;
    BoundariesOpt boundaries;
    std::map<int, std::unique_ptr<ReconstructedGridHandler<dimworld>>> trimResultMap;

    auto parameterSpaceGrid() {
      auto gV        = leafGridView();
      auto patchData = gV.impl().getPatchData();

      using TensorSpaceCoordinates = TensorProductCoordinates<double, dimension>;
      using ParametricGrid         = YaspGrid<dimension, TensorSpaceCoordinates>;

      std::array<std::vector<double>, dimension> tensorProductCoordinates;
      for (int i = 0; i < dimension; ++i)
        std::ranges::unique_copy(patchData.knotSpans[i], std::back_inserter(tensorProductCoordinates[i]));

      return std::make_unique<ParametricGrid>(tensorProductCoordinates);
    }

    // For dim != 2
    void trimElement() {}

    void trimElement()
      requires(dim == 2 && dimworld == 2)
    {
      // Set up trimFlags
      trimFlags = std::vector<ElementTrimFlag>(size(0));
      std::ranges::fill(trimFlags, ElementTrimFlag::full);

      if (!boundaries.has_value()) {
        leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this, trimFlags));
        return;
      }

      // Trimming is done in th 2d parameterSpace of the Grid
      auto paraGrid               = parameterSpaceGrid();
      auto parameterSpaceGridView = paraGrid->leafGridView();

      // Use Trim namespace for more concise function names
      using namespace Impl::Trim;


      // Get Clip as ClipperPath
      Clipper2Lib::PathsD clip = getClip(*boundaries);

      const auto& indexSet = parameterSpaceGridView.indexSet();
      for (auto& element : elements(parameterSpaceGridView)) {
        auto index{indexSet.index(element)};
        auto corners = getElementCorners(element);

        using namespace Dune::IGA::Impl::Trim;

        auto [trimFlag, clippingResultOpt] = clipElement(element, clip);

        trimFlags[index] = trimFlag;

        if (clippingResultOpt.has_value()) {
          auto elementBoundariesOpt = constructElementBoundaries(*clippingResultOpt, corners, *boundaries);

          if (elementBoundariesOpt.has_value()) {
            // ReconstructGrid and save in a GridHandler
            trimResultMap[index] = std::move(std::make_unique<ReconstructedGridHandler<dimworld>>(
                *elementBoundariesOpt,
                NURBSPatchGeometry<dim, dimworld>{
                    std::make_shared<NURBSPatchData<dim, dimworld>>(currentPatchRepresentation_)},
                currentPatchRepresentation_.degree));
          } else
            trimFlags[index] = ElementTrimFlag::empty;
        }

      }  // Element Loop End

      // Update the GridView with the set trimFlags
      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this, trimFlags));
    }
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
