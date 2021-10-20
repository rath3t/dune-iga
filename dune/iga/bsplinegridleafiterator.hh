// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_IGA_BSPLINEGRID_LEAFITERATOR_HH
#define DUNE_IGA_BSPLINEGRID_LEAFITERATOR_HH

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplineleafgridview.hh>
#include <dune/common/iteratorfacades.hh>

/** \file
 * \brief The BSplineIterator class
 */

namespace Dune::IGA
  {
    /** \brief b-spline leaf iterator */
    template<int codim, typename BSplineGridView, typename BSplineEntity>
    class BSplineGridLeafIterator:
      public ForwardIteratorFacade<BSplineGridLeafIterator<codim, BSplineGridView, BSplineEntity>,BSplineEntity,BSplineEntity&, int>
    {
    public:
      //typedef typename BSplineGridView::template Codim<codim>::Entity Entity;

      BSplineGridLeafIterator(const BSplineGridView& GridView, int index)
        : BSplineGridView_(&GridView), directIndex_(index)
      {
      }

      BSplineEntity& dereference() const
      {
        return BSplineGridView_-> template getEntity<codim>(directIndex_);
      }

      bool equals(const BSplineGridLeafIterator<codim,typename std::remove_const<BSplineGridView>::type, typename std::remove_const<BSplineEntity>::type>& other) const
      {
        return directIndex_ == other.directIndex_ && BSplineGridView_->BSplinepatch_ == other.BSplineGridView_->BSplinepatch_;
      }

      void increment()
      {
        ++directIndex_;
      }

      void advance(int n)
      {
        directIndex_=directIndex_+n;
      }





    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////

    // B-Spline patch
    private:
      const BSplineGridView* BSplineGridView_;
      //Global index of the knot span, can be transferred into multidimensional indices
      unsigned int directIndex_;
      //std::array<unsigned int,dim> multiIndex;


    };
  }

#endif  // DUNE_IGA_BSPLINEGRID_LEAFITERATOR_HH
