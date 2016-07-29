// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRID_HH
#define DUNE_IGA_NURBSGRID_HH

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>

namespace Dune
{
  namespace IGA
  {
    /** \brief NURBS grid manager */
    template<int dim, int dimworld>
    class NURBSGrid
    {
    public:
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
      : leafGridView_(std::make_shared<NURBSLeafGridView<dim,dimworld>>(knotSpans, controlPoints, weights, order))
      {
      }

      NURBSLeafGridView<dim,dimworld>& leafGridView()
      {
        return *(this->leafGridView_);
      }

    private:

      std::shared_ptr <NURBSLeafGridView<dim,dimworld>> leafGridView_;
    };
  }
}

#endif  // DUNE_IGA_BSPLINEGRID_HH
