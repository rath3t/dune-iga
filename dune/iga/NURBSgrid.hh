// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/iga/NURBSleafgridview.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/igaidset.hh>
#include <dune/iga/gridCapabilities.hh>

namespace Dune::IGA {

  /** \brief NURBS grid manager */
  template <std::integral auto dim, std::integral auto dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  class NURBSGrid {
  public:
    using NurbsGridLinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;
    using GlobalCoordinateType         = typename NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
    using LocalCoordinateType          = typename NurbsGridLinearAlgebraTraits::LocalCoordinateType;
    using JacobianTransposedType       = typename NurbsGridLinearAlgebraTraits::JacobianTransposedType;
    using JacobianInverseTransposed    = typename NurbsGridLinearAlgebraTraits::JacobianInverseTransposed;



        static constexpr std::integral auto dimension      = dim;
    static constexpr std::integral auto dimensionworld = dimworld;
    using ctype                                        = typename GlobalCoordinateType::value_type;

    using ControlPointNetType = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType;

    using Comm = Communication<No_Comm>;

    struct Traits {
      using GridView = NURBSLeafGridView<NURBSGrid<dim, dimworld>>;
      using IndexSet = NURBSGridLeafIndexSet<GridView>;
      using LocalIdSet = IgaIdSet<NURBSGrid>;
      using LeafIntersectionIterator = NURBSGridLeafIterator<NURBSGridEntity<dim-1UL, NURBSLeafGridView<NURBSGrid<dim, dimworld>>>>;
      using IntersectionIterator = NURBSGridLeafIterator<NURBSGridEntity<dim-1UL, NURBSLeafGridView<NURBSGrid<dim, dimworld>>>>;
      template <std::integral auto cd>
      struct Codim {
        using Entity   = NURBSGridEntity<cd, NURBSLeafGridView<NURBSGrid<dim, dimworld>>>;
        using Geometry = NURBSGeometry<dim-cd,  dimworld,dim, NurbsGridLinearAlgebraTraitsImpl>;
        template <PartitionIteratorType pitype>
        struct Partition {
          /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
          using LeafIterator = NURBSGridLeafIterator<Entity>;

          /** \brief The type of the iterator over the level entities of this codim on this partition. */
          using LevelIterator = LeafIterator;
        };
      };
    };
    template <std::integral auto cd>
    using Codim = typename Traits::template Codim<cd>;

    NURBSGrid(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>& nurbsPatchData)
        : coarsestPatchRepresentation_{nurbsPatchData}, currentPatchRepresentation_{coarsestPatchRepresentation_} {
      static_assert(dim<=3,"Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size()-nurbsPatchData.order[0]-1==nurbsPatchData.controlPoints.size()[0] && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim>1)
      assert(nurbsPatchData.knotSpans[1].size()-nurbsPatchData.order[1]-1==nurbsPatchData.controlPoints.size()[1] && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim>2)
      assert(nurbsPatchData.knotSpans[2].size()-nurbsPatchData.order[2]-1==nurbsPatchData.controlPoints.size()[2] && "The size of the controlpoints and the knotvector size do not match in the third direction");
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
          currentPatchRepresentation_{coarsestPatchRepresentation_} {}

    void globalRefine(int refinementLevel) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots        = generateRefinedKnots(refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
    }


    void globalRefineInDirection(const int dir, const int refinementLevel) {
      auto additionalKnots        = generateRefinedKnots(dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
    }

    auto leafGridView() const { return NURBSLeafGridView<NURBSGrid<dim, dimworld>>(currentPatchRepresentation_, *this); }

    const auto& globalIdSet() const
    {
      return Dune::IGA::IgaIdSet(*this);
    }

    const auto& localIdSet() const
    {
      return this->globalIdSet();
    }

    [[nodiscard]] const Comm& comm() const { return ccobj; }

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

    Comm ccobj;
    NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits> currentPatchRepresentation_;
  };
}  // namespace Dune::IGA
