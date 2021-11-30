// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/iga/concepts.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbspatch.hh>
//#include <dune/iga/gridcapabilities.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/igaidset.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsintersection.hh>
#include <dune/iga/nurbslocalgeometry.hh>
#include <dune/iga/nurbsgridtraits.hh>

namespace Dune::IGA {
  template<int  cd, typename GridImpl>
  struct EntitySeedStruct {
    static constexpr int codimension = cd;
    [[nodiscard]] bool isValid() const { return true;}
    int index; };


  template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily;


  /** \brief NURBS grid manager */
  template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
      class NURBSGrid : public Dune::Grid<dim,dimworld,typename NurbsGridLinearAlgebraTraitsImpl::value_type,NurbsGridFamily<dim,dimworld,NurbsGridLinearAlgebraTraitsImpl>>{
  public:
    using NurbsGridLinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;
    using GlobalCoordinateType         = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<dimworld>;
    using LocalCoordinateType          = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<dim>;
    using JacobianTransposedType       = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<dim,dimworld>;
    using JacobianInverseTransposed    = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<dimworld,dim>;

    static constexpr std::integral auto dimension      = dim;
    static constexpr std::integral auto dimensionworld = dimworld;
    using ctype                                        = typename NurbsGridLinearAlgebraTraits::value_type;

    using ControlPointNetType = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType;

//    using Comm       = Communication<No_Comm>;
//    using CollectiveCommunication       = Communication<No_Comm>;
    using GridFamily = NurbsGridFamily<dim,dimworld,NurbsGridLinearAlgebraTraitsImpl>;

    using Traits =  typename GridFamily::Traits;
    template <int cd>
    using Codim  =  typename GridFamily::Traits::template Codim<cd>;

    //    struct Traits {
//      using Grid = NURBSGrid;
//      using GridView                 = NURBSLeafGridView<NURBSGrid<dim, dimworld>>;
//      using LeafGridView                 = GridView;
//      using LevelGridView                 = GridView;
//      using IndexSet                 = NURBSGridLeafIndexSet<GridView>;
//      using LevelIndexSet = IndexSet;
//      using LeafIndexSet = IndexSet;
//      using LocalIdSet               = IgaIdSet<NURBSGrid>;
//      using GlobalIdSet               = LocalIdSet;
//      using LeafIntersection = NURBSintersection<NURBSGrid>;
//      using LevelIntersection = LeafIntersection;
//      using LeafIntersectionIterator = NURBSGridInterSectionIterator<NURBSGrid>;
//      using LevelIntersectionIterator     = LeafIntersectionIterator;
//      using HierarchicIterator     = NurbsHierarchicIterator<NURBSGrid>;
//      using CollectiveCommunication = Communication<No_Comm>;
//      template <int cd>
//      struct Codim {
//        using Entity   = NURBSGridEntity<cd, NURBSLeafGridView<NURBSGrid>>;
//        using Geometry = NURBSGeometry<dim - cd, dimworld, dim, NurbsGridLinearAlgebraTraitsImpl>;
//        using LevelIterator = NURBSGridLeafIterator<NURBSGridEntity<cd, NURBSLeafGridView<NURBSGrid>>>;
//        using LeafIterator  = NURBSGridLeafIterator<NURBSGridEntity<cd, NURBSLeafGridView<NURBSGrid>>>;
//        using LocalGeometry =  NURBSLocalGeometry<dim-cd, dim, dim, NurbsGridLinearAlgebraTraitsImpl>;
//        using EntitySeed = EntitySeedStruct<cd>;
//        template <PartitionIteratorType pitype>
//        struct Partition {
//          /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
//          using LeafIterator = NURBSGridLeafIterator<Entity>;
//
//          /** \brief The type of the iterator over the level entities of this codim on this partition. */
//          using LevelIterator = LeafIterator;
//        };
//      };
//
//    };
//    template <int cd>
//    using Codim         = typename GridFamily::Traits::template Codim<cd>;
//    using LeafGridView  = typename GridFamily::Traits::LeafGridView;
//    using LevelGridView = typename GridFamily::Traits::LevelGridView;
//    using LocalIdSet    = typename GridFamily::Traits::LocalIdSet;
//    using GlobalIdSet    = typename GridFamily::Traits::LocalIdSet;
//    using LevelIndexSet = typename GridFamily::Traits::LevelIndexSet;
    using LeafIndexSet = typename Traits::LeafIndexSet;
    using GridView = typename Traits::LeafGridView;
    using ElementEntity = typename Traits::template Codim<0>::Entity;
//    using HierarchicIterator = typename GridFamily::Traits::HierarchicIterator;
//    using Intersection = typename GridFamily::Traits::LeafIntersection;
//    using LevelIntersectionIterator = typename GridFamily::Traits::LeafIntersectionIterator;

    NURBSGrid() = default;


    NURBSGrid(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>& nurbsPatchData)
        : coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatch_{currentPatchRepresentation_},
          idSet_{std::make_unique<IgaIdSet<NURBSGrid>>(*this)},
          indexdSet_{std::make_unique<LeafIndexSet>(this->leafGridView())},
          leafGridView_{std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(currentPatchRepresentation_, *this)}{
      static_assert(dim <= 3, "Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size() - nurbsPatchData.order[0] - 1 == nurbsPatchData.controlPoints.size()[0]
             && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim > 1)
        assert(nurbsPatchData.knotSpans[1].size() - nurbsPatchData.order[1] - 1 == nurbsPatchData.controlPoints.size()[1]
               && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim > 2)
        assert(nurbsPatchData.knotSpans[2].size() - nurbsPatchData.order[2] - 1 == nurbsPatchData.controlPoints.size()[2]
               && "The size of the controlpoints and the knotvector size do not match in the third direction");
    }

    /** \brief  constructor
     *
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order order of the B-Spline structure for each dimension
     */
    NURBSGrid(const std::array<std::vector<double>, dim>& knotSpans, const ControlPointNetType& controlPoints,
              const std::array<int, dim>& order)
        : coarsestPatchRepresentation_{NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>(knotSpans, controlPoints, order)},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatch_{currentPatchRepresentation_},
          leafGridView_{std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(currentPatchRepresentation_, *this)},
          idSet_{std::make_unique<IgaIdSet<NURBSGrid>>(*this)},
          indexdSet_{std::make_unique<LeafIndexSet>(this->leafGridView())}{}

    void globalRefine(int refinementLevel) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots        = generateRefinedKnots(refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
      leafGridView_ =  std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(currentPatchRepresentation_, *this);
      idSet_       = std::make_unique<IgaIdSet<NURBSGrid>>(*this);
      indexdSet_       = std::make_unique<LeafIndexSet>(this->leafGridView());
      finestPatch_ = NURBSPatch<dim, dimworld, NurbsGridLinearAlgebraTraits>(currentPatchRepresentation_);
    }

    void globalRefineInDirection(const int dir, const int refinementLevel) {
      auto additionalKnots        = generateRefinedKnots(dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
      leafGridView_ =  std::make_shared<NURBSLeafGridView<NURBSGrid<dim, dimworld>>>(currentPatchRepresentation_, *this);
      idSet_                      = std::make_unique<IgaIdSet<NURBSGrid>>(*this);
      indexdSet_       = std::make_unique<LeafIndexSet>(this->leafGridView());
      finestPatch_                = NURBSPatch<dim, dimworld, NurbsGridLinearAlgebraTraits>(currentPatchRepresentation_);
    }

    [[nodiscard]] int size(int codim) const { return finestPatch_.size(codim); }
    [[nodiscard]] int numBoundarySegments() const { return 4;}
    [[nodiscard]] int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const {
//      leafGridView_ = std::make_shared<NURBSLeafGridView<NURBSGrid>>(currentPatchRepresentation_, *this);
      return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark (const ElementEntity& element) const {return 0;}
    bool mark(int refCount, const ElementEntity& element){ return false;}

    template<int  cd>
    typename Codim<cd>::Entity entity( EntitySeedStruct<cd,NURBSGrid>& seed ) const  {
      return leafGridView_->template getEntity<cd>(seed.index);
    }

    int size(const GeometryType& type) const{
      if(type==Dune::GeometryTypes::vertex || type==Dune::GeometryTypes::cube(1) || type==Dune::GeometryTypes::cube(2) || type==Dune::GeometryTypes::cube(3))
        return this->leafGridView().size(dimension-type.dim());
      else
        return 0;
    }
    int size(int lvl,const GeometryType& type) const{
      return this->size(type);
    }

    const auto& globalIdSet() const { return *idSet_; }
    const auto& levelIndexSet(int lvl) const { return *indexdSet_; }
    const auto& leafIndexSet() const { return *indexdSet_; }

    [[nodiscard]] int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    [[nodiscard]] const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

  private:
    auto generateRefinedKnots(const int dir, const int refinementLevel) {
      const int newKnotsSizeForEachSpan = Dune::power(2, refinementLevel);
      const auto& knotSpans             = currentPatchRepresentation_.knotSpans;
      auto unique_Knots                 = knotSpans;
      std::vector<double> additionalKnots;
      auto& unique_KnotPerDim = unique_Knots[dir];
      unique_KnotPerDim.erase(std::ranges::begin(std::ranges::unique(unique_KnotPerDim)), std::end(unique_KnotPerDim));

      for (int i = 0; i < unique_KnotPerDim.size() - 1; ++i) {
        const double spanLength = unique_KnotPerDim[i + 1] - unique_KnotPerDim[i];
        const double increment  = spanLength / newKnotsSizeForEachSpan;
        for (int j = 1; j < newKnotsSizeForEachSpan; ++j) {
          additionalKnots.emplace_back(unique_KnotPerDim[i] + increment * j);
        }
      }
      return additionalKnots;
    }

    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, NurbsGridLinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, NurbsGridLinearAlgebraTraits> currentPatchRepresentation_;
    NURBSPatch<dim, dimworld, NurbsGridLinearAlgebraTraits> finestPatch_;
    std::shared_ptr<NURBSLeafGridView<NURBSGrid>> leafGridView_;
    std::unique_ptr<IgaIdSet<NURBSGrid>> idSet_;
    std::unique_ptr<LeafIndexSet> indexdSet_;
  };


  template<std::integral auto  dim, std::integral auto dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
   NURBSLeafGridView<NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>>
  levelGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid, int level)
  {
    return grid.levelGridView(level);
  }

  template<std::integral auto  dim, std::integral auto dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  NURBSLeafGridView<NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>>
  leafGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid)
  {
    return grid.leafGridView();
  }


  template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily
  {
using GridImpl = Dune::IGA::NURBSGrid<dim,dimworld,NurbsGridLinearAlgebraTraitsImpl>;
    typedef NurbsGridTraits<dim,dimworld,GridImpl,
                   NURBSGeometry,
                   NURBSGridEntity,
                   NURBSGridLeafIterator,
                       NURBSintersection,
                       NURBSintersection,
                       NURBSGridInterSectionIterator,
                       NURBSGridInterSectionIterator,
                       NurbsHierarchicIterator,
                       NURBSGridLeafIterator,
                       NURBSGridLeafIndexSet< GridImpl>,
                       NURBSGridLeafIndexSet< GridImpl>,
                       IgaIdSet< GridImpl>,
                       int,
                       IgaIdSet< GridImpl>,
                       int ,
                       CollectiveCommunication<No_Comm>,
                       NurbsLeafGridViewTraits,
                       NurbsLeafGridViewTraits,
                       EntitySeedStruct,
                       NURBSLocalGeometry>
        Traits;
  };

//  template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
//class NURBSGrid : public  NURBSGrid<(std::size_t)dim,(std::size_t)dimworld,NurbsGridLinearAlgebraTraitsImpl>;
}  // namespace Dune::IGA
