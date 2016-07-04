// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRIDVIEW_HH
#define DUNE_IGA_BSPLINEGRIDVIEW_HH

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplinegridleafiterator.hh>

namespace Dune
{
  namespace IGA
  {
    /** \brief b-spline grid manager */
    template<int dim, int dimworld>
    class BSplineLeafGridView
    {
    public:
      typedef BSplineLeafGridView<dim, dimworld> BSplineGridViewImp;

      template<int codim>
      struct Codim
      {
        typedef BSplineGridLeafIterator<codim,BSplineGridViewImp> Iterator;
      };

      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      BSplineLeafGridView(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNet<dim,dimworld> controlPoints,
                   const std::array<int,dim> order)
      : BSplinepatch_(std::make_shared<BSplinePatch<dim,dimworld>>(knotSpans, controlPoints, order))
      {
      }

    private:
       std::shared_ptr <BSplinePatch<dim,dimworld>> BSplinepatch_;
    };
  }
}

#endif  // DUNE_IGA_BSPLINEGRIDVIEW_HH
