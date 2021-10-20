// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRID_HH
#define DUNE_IGA_BSPLINEGRID_HH

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplineleafgridview.hh>

namespace Dune::IGA
  {
    /** \brief b-spline grid manager */
    template<int dim, int dimworld>
    class BSplineGrid
    {
    public:

      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      BSplineGrid(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNet<dim,dimworld> controlPoints,
                   const std::array<int,dim> order)
      : leafGridView_(std::make_shared<BSplineLeafGridView<dim,dimworld>>(knotSpans, controlPoints, order))
      {

      }

      BSplineLeafGridView<dim,dimworld> leafGridView()
      {
        return *(this->leafGridView_);
      }



    private:

      std::shared_ptr <BSplineLeafGridView<dim,dimworld>> leafGridView_;
    };
  }

#endif  // DUNE_IGA_BSPLINEGRID_HH
