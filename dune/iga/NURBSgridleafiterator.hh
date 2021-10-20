// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_IGA_NURBSGRID_LEAFITERATOR_HH
#define DUNE_IGA_NURBSGRID_LEAFITERATOR_HH

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <dune/common/iteratorfacades.hh>

/** \file
 * \brief The NURBSGridIterator class
 */

 namespace Dune::IGA
  {
    /** \brief NURBS gird leaf iterator */
    template<int codim, typename NURBSGridView, typename NURBSEntity>
    class NURBSGridLeafIterator:
      public ForwardIteratorFacade<NURBSGridLeafIterator<codim, NURBSGridView, NURBSEntity>,NURBSEntity,NURBSEntity&, int>
    {
    public:
      //typedef typename NURBSGridView::template Codim<codim>::Entity Entity;

      NURBSGridLeafIterator(const NURBSGridView& GridView, int index)
        : NURBSGridView_(&GridView), directIndex_(index)
      {
      }

      NURBSEntity& dereference() const
      {
        return NURBSGridView_-> template getEntity<codim>(directIndex_);
      }

      bool equals(const NURBSGridLeafIterator<codim,typename std::remove_const<NURBSGridView>::type, typename std::remove_const<NURBSEntity>::type>& other) const
      {
        return directIndex_ == other.directIndex_ && NURBSGridView_->NURBSpatch_ == other.NURBSGridView_->NURBSpatch_;
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
      const NURBSGridView* NURBSGridView_;
      //Global index of the knot span, can be transferred into multidimensional indices
      unsigned int directIndex_;
      //std::array<unsigned int,dim> multiIndex;


    };
  }

#endif  // DUNE_IGA_NURBSGRID_LEAFITERATOR_HH
