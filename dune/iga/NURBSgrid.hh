// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <dune/iga/igaalgorithms.hh>

namespace Dune::IGA
  {
    /** \brief NURBS grid manager */
    template<int dim, int dimworld>
    class NURBSGrid
    {
    public:

        static constexpr int dimension = dim;
        static constexpr int dimensionworld = dimworld;
        using ctype = double;

        using Comm = Communication<No_Comm>;
        struct Traits{
            template< int cd >
            struct Codim {

                using Entity =  NURBSGridEntity<cd, NURBSLeafGridView<NURBSGrid<dim,dimworld>>>;
                using Geometry = typename Entity :: Geometry;
                template <PartitionIteratorType pitype>
                struct Partition
                {
                    /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
                    using  LeafIterator = NURBSGridLeafIterator<Entity>;
                    /** \brief The type of the iterator over the level entities of this codim on this partition. */
                    using LevelIterator = LeafIterator;

                };
            };
        };


      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      NURBSGrid(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNetFVd<dim,dimworld>& controlPoints,
                   const MultiDimensionNetFVd<dim,1>& weights,
                   const std::array<int,dim>& order)
      :        coarsestPatchRepresentation_{NURBSPatchData<dim,dimworld>(knotSpans, controlPoints, weights, order)},
              currentPatchRepresentation_{coarsestPatchRepresentation_}
      {
      }

      void globalRefine(int refinementLevel)
      {
        if (refinementLevel==0)
          return;
        const int newKnotsSizeForEachSpan = Dune::power(2,refinementLevel);
        const auto& knotSpans = coarsestPatchRepresentation_.getKnots();
        auto unique_Knots = knotSpans;
        std::array<std::vector<double>,dim> additionalKnots;
        for(int curdim=0; curdim<dim; ++curdim)
        {
          auto& unique_KnotPerDim = unique_Knots[curdim];
          unique_KnotPerDim.erase(std::unique(unique_KnotPerDim.begin(), unique_KnotPerDim.end()),
                                  unique_KnotPerDim.end());
          for (int i= 0; i<unique_KnotPerDim.size()-1;++i)
          {
            const double spanLength = unique_KnotPerDim[i+1] - unique_KnotPerDim[i];
            const double increment = spanLength/newKnotsSizeForEachSpan;
            for (int j=1; j<newKnotsSizeForEachSpan; ++j)
            {
              additionalKnots[curdim].emplace_back(unique_KnotPerDim[i]+increment*j);
            }
          }
        }
        for(int curdim=0; curdim<dim; ++curdim)
        { std::cerr<<std::endl;
        for(auto& addKnot : additionalKnots[curdim])
            std::cerr << addKnot << " ";
          }
          std::cerr<<"======================="<<std::endl;
        //        std::cerr<<additionalKnots[1]<<std::endl;
        if constexpr (dim==1)
        {
//
currentPatchRepresentation_ = curveKnotRefinement(coarsestPatchRepresentation_, additionalKnots);
////          std::vector<double> additionalKnots;
////          currentPatchRepresentation_ = curveKnotInsertion(coarsestPatchRepresentation_, additionalKnots);
        }
        else if constexpr (dim==2)
        {
//currentPatchRepresentation_ = curveKnotRefinement(coarsestPatchRepresentation_, additionalKnots);
        }
      }

      auto leafGridView()
      {
        return NURBSLeafGridView<NURBSGrid<dim,dimworld>>(currentPatchRepresentation_,*this);
      }

        [[nodiscard]] const Comm& comm() const
        {
            return ccobj;
        }
    private:
      Comm ccobj;
//      std::shared_ptr <NURBSLeafGridView<NURBSGrid<dim,dimworld>>> leafGridView_;
      NURBSPatchData<dim,dimworld> coarsestPatchRepresentation_;
      NURBSPatchData<dim,dimworld> currentPatchRepresentation_;
    };
  }

