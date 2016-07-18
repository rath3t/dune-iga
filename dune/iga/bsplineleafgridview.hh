// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRIDVIEW_HH
#define DUNE_IGA_BSPLINEGRIDVIEW_HH

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplinegridleafiterator.hh>
#include <dune/iga/bsplinegridentity.hh>

namespace Dune
{
  namespace IGA
  {
    /** \brief b-spline grid manager */
    template<int dim, int dimworld>
    class BSplineLeafGridView
    {
    public:
      //static const int dim = dim;
      //static const int dimWorld = dimWorld;

      template<int codim, class GridViewImp>
      friend class BSplineGridEntity;


      template<int codim, typename BSplineGridView, typename BSplineEntity>
      friend class BSplineGridLeafIterator;


      typedef BSplineLeafGridView<dim, dimworld> BSplineGridView;
      typedef BSplineGeometry<dim, dimworld> Geometry;

      template<int codim>
      struct Codim
      {
        typedef BSplineGridEntity<codim,BSplineGridView> Entity;
        typedef BSplineGridLeafIterator<codim,BSplineGridView, Entity> Iterator;

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
        int elementSize = BSplinepatch_ ->knotElementNet_->directSize();
        entityVector_.reserve(elementSize);

        //Fill the element vector of codim 0
        for(unsigned int i=0; i<elementSize; ++i)
        {
          entityVector_.push_back(std::make_shared<BSplineGridEntity<0, BSplineGridView>>(*this, i));
        }
      }

      template<int codim>
      typename Codim<codim>::Entity& getEntity(unsigned int directIndex) const
      {
        //need to be rewrite for other codims
        if (codim==0)
        {
          return *(entityVector_.at(directIndex));
        }

      }

//      template<>
//      typename Codim<0>::Entity getEntity<0>(unsigned int directIndex)
//      {
//         return *(entityVector_.at(directIndex));
//      }

      typename Codim<0>::Iterator begin () const
      {
        return BSplineGridLeafIterator<0,BSplineGridView, typename Codim<0>::Entity>(*this, 0);
      }

      typename Codim<0>::Iterator end () const
      {
        int elementSize = entityVector_.size();
        return BSplineGridLeafIterator<0,BSplineGridView, typename Codim<0>::Entity>(*this, elementSize);
      }

    private:
       std::shared_ptr <BSplinePatch<dim,dimworld>> BSplinepatch_;
       std::vector<std::shared_ptr<BSplineGridEntity<0, BSplineGridView>>> entityVector_;
    };
  }
}

#endif  // DUNE_IGA_BSPLINEGRIDVIEW_HH
