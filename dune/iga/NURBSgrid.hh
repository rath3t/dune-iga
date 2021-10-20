// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>

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
                    using  LeafIterator = NURBSGridLeafIterator<cd, NURBSLeafGridView<NURBSGrid<dim,dimworld>>,Entity>;
                    /** \brief The type of the iterator over the level entities of this codim on this partition. */
                    using LevelIterator = LeafIterator;

                };
            };
        };



        void globalRefine(int refinementevel)
        {

        }
      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      NURBSGrid(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNet<dim,dimworld> controlPoints,
                   const MultiDimensionNet<dim,1> weights,
                   const std::array<int,dim> order)
      : leafGridView_(std::make_shared<NURBSLeafGridView<NURBSGrid<dim,dimworld>>>(knotSpans, controlPoints, weights, order))
      {
      }

      NURBSLeafGridView<NURBSGrid<dim,dimworld>>& leafGridView()
      {
        return *(this->leafGridView_);
      }

        [[nodiscard]] const Comm& comm() const
        {
            return ccobj;
        }
    private:
      Comm ccobj;
      std::shared_ptr <NURBSLeafGridView<NURBSGrid<dim,dimworld>>> leafGridView_;
    };
  }

