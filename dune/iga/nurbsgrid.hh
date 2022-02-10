// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsgridtraits.hh>
#include <dune/iga/nurbsidset.hh>
#include <dune/iga/nurbsintersection.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbslocalgeometry.hh>
#include <dune/iga/nurbspatch.hh>
#include <dune/common/parallel/communication.hh>

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

    using ControlPointNetType = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType;

    using GridFamily = NurbsGridFamily<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;

    using Traits = typename GridFamily::Traits;
    template <int cd>
    using Codim = typename GridFamily::Traits::template Codim<cd>;

    using LeafIndexSet  = typename Traits::LeafIndexSet;
    using GridView      = typename Traits::LeafGridView;
    using ElementEntity = typename Traits::template Codim<0>::Entity;
    NURBSGrid()         = default;

    explicit NURBSGrid(const NURBSPatchData<dim, dimworld, LinearAlgebraTraits>& nurbsPatchData)
        : coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()} {
      static_assert(dim <= 3, "Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size() - nurbsPatchData.degree[0] - 1 == nurbsPatchData.controlPoints.size()[0]
             && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim > 1)
        assert(nurbsPatchData.knotSpans[1].size() - nurbsPatchData.degree[1] - 1 == nurbsPatchData.controlPoints.size()[1]
               && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim > 2)
        assert(nurbsPatchData.knotSpans[2].size() - nurbsPatchData.degree[2] - 1 == nurbsPatchData.controlPoints.size()[2]
               && "The size of the controlpoints and the knotvector size do not match in the third direction");

      finestPatches_->emplace_back(currentPatchRepresentation_);
      leafGridView_ = std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(finestPatches_, *this);
      idSet_        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
      indexdSet_    = std::make_unique<LeafIndexSet>(this->leafGridView());
    }

    void globalRefine(int refinementLevel) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots        = generateRefinedKnots(currentPatchRepresentation_.knotSpans, refDirection, refinementLevel);
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
        return (finestPatches_.get()->front().validKnotSize()[0] + finestPatches_.get()->front().validKnotSize()[1]) * 2;
      else if constexpr (dimension == 3)
        return 2
               * (finestPatches_.get()->front().validKnotSize()[0] * finestPatches_.get()->front().validKnotSize()[1]
                  + finestPatches_.get()->front().validKnotSize()[1] * finestPatches_.get()->front().validKnotSize()[2]
                  + finestPatches_.get()->front().validKnotSize()[0] * finestPatches_.get()->front().validKnotSize()[2]);
      __builtin_unreachable();
    }
    [[nodiscard]] int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const { return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark(const ElementEntity& element) const { return 0; }
    bool mark(int refCount, const ElementEntity& element) { return false; }

    template <int cd>
    typename Codim<cd>::Entity entity(EntitySeedStruct<cd, NURBSGrid>& seed) const {
      return leafGridView_->template getEntity<cd>(seed.index_);
    }

    [[nodiscard]] int size(const GeometryType& type) const {
      if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1) || type == Dune::GeometryTypes::cube(2)
          || type == Dune::GeometryTypes::cube(3))
        return this->leafGridView().size(dimension - type.dim());
      else
        return 0;
    }
    int size(int lvl, const GeometryType& type) const { return this->size(type); }

    const auto& globalIdSet() const { return *idSet_; }
    const auto& levelIndexSet(int lvl) const { return *indexdSet_; }
    const auto& leafIndexSet() const { return *indexdSet_; }

    [[nodiscard]] int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    [[nodiscard]] const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

  private:
    void updateStateAfterRefinement() {
      finestPatches_.get()->front() = NURBSPatch<dim, dimworld, LinearAlgebraTraits>(currentPatchRepresentation_);
      leafGridView_                 = std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(finestPatches_, *this);
      indexdSet_                    = std::make_unique<LeafIndexSet>(this->leafGridView());
      idSet_                        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
    }

    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> currentPatchRepresentation_;
    std::shared_ptr<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>> finestPatches_;
    std::shared_ptr<NURBSLeafGridView<NURBSGrid>> leafGridView_;
    std::unique_ptr<IgaIdSet<NURBSGrid>> idSet_;
    std::unique_ptr<LeafIndexSet> indexdSet_;
  };

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  NURBSLeafGridView<NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>> levelGridView(
      const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid, int level) {
    return grid.levelGridView(level);
  }

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  NURBSLeafGridView<NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>> leafGridView(
      const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid) {
    return grid.leafGridView();
  }

  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily {
    using GridImpl = Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    typedef NurbsGridTraits<dim, dimworld, GridImpl, NURBSGeometry, NURBSGridEntity, NURBSGridLeafIterator, NURBSintersection,
                            NURBSintersection, NURBSGridInterSectionIterator, NURBSGridInterSectionIterator, NurbsHierarchicIterator,
                            NURBSGridLeafIterator, NURBSGridLeafIndexSet<GridImpl>, NURBSGridLeafIndexSet<GridImpl>, IgaIdSet<GridImpl>,
                            int, IgaIdSet<GridImpl>, int, CollectiveCommunication<No_Comm>, NurbsLeafGridViewTraits,
                            NurbsLeafGridViewTraits, EntitySeedStruct, NURBSLocalGeometry>
        Traits;
  };

}  // namespace Dune::IGA
