// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_IGA_NURBSGRID_LEAFITERATOR_HH
#define DUNE_IGA_NURBSGRID_LEAFITERATOR_HH

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <dune/common/iteratorfacades.hh>
#include <span>
/** \file
 * \brief The NURBSGridIterator class
 */

 namespace Dune::IGA
  {
    /** \brief NURBS gird leaf iterator */
    template<typename NURBSEntity>
  class NURBSGridLeafIterator : public std::vector<NURBSEntity>::iterator
    {
    public:
      //typedef typename NURBSGridView::template Codim<codim>::Entity Entity;

      using Reference =  NURBSEntity;
      explicit NURBSGridLeafIterator(const typename std::vector<NURBSEntity>::iterator& spanIter)
          :  std::vector<NURBSEntity>::iterator(spanIter)
      {}
//      NURBSGridLeafIterator(const NURBSGridView& GridView, int index)
//        : NURBSGridView_(&GridView), directIndex_(index)
//      {
//      }
//
//      NURBSEntity& dereference() const
//      {
//        return NURBSGridView_-> template getEntity<codim>(directIndex_);
//      }
//
//        NURBSEntity& operator*() const
//        {
//            return  dereference();
//        }
//
//      bool operator==(const NURBSGridLeafIterator<codim,typename std::remove_const<NURBSGridView>::type, typename std::remove_const<NURBSEntity>::type>& other) const
//      {
//        return directIndex_ == other.directIndex_ && NURBSGridView_->NURBSpatch_ == other.NURBSGridView_->NURBSpatch_;
//      }
//
//      void operator++()
//      {
//        ++directIndex_;
//      }
//
//      void advance(int n)
//      {
//        directIndex_=directIndex_+n;
//      }
//
//
//        auto operator->()
//        {
//          return &NURBSGridView_-> template getEntity<codim>(directIndex_);
//        }


    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////

    // B-Spline patch
    private:
//      const NURBSGridView* NURBSGridView_;
      //Global index of the knot span, can be transferred into multidimensional indices
//      unsigned int directIndex_;
      //std::array<unsigned int,dim> multiIndex;


    };
  }

#endif  // DUNE_IGA_NURBSGRID_LEAFITERATOR_HH
