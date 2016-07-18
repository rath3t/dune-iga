// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRIDENTITY_HH
#define DUNE_IGA_BSPLINEGRIDENTITY_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplineleafgridview.hh>

/** \file
 * \brief The BSplineGridEntity class
 */

 namespace Dune
 {
   namespace IGA
   {

   /** \brief
   *.
   */
    template<int codim, class GridViewImp>
    class BSplineGridEntity
    {
//      template <int dim, int dimworld>
//      friend class BSplinePatch;
//
//      template <int dim, int dimworld>
//      friend class BSplineLeafGridView;


      /* !Index set and intersection iterator are not implemented yet

      template <class GridImp_>
      friend class BSplinePatchIntersectionIterator;

      friend class BSplinePatchIndexSet<GridImp>;

      */



  /*public:
    typedef typename <int dim, int dimworld>template BSplinePatchIterator ValidKnotIterator;
    */
    public:

      //! Default Constructor
      BSplineGridEntity ()
        : BSplineGridView_(nullptr)
      {}

      BSplineGridEntity (const GridViewImp& gridView, unsigned int directIndex)
        : BSplineGridView_(&gridView), directIndex_(directIndex)
      {}

      //! Geometry of this entity
      typename GridViewImp::Geometry geometry () const
      {
        auto const &knotElementNet = BSplineGridView_->BSplinepatch_->knotElementNet_;
        auto const &multiIndex = knotElementNet->directToMultiIndex(0);
        return BSplineGridView_->BSplinepatch_->geometry(multiIndex);
      }

      unsigned int  getIndex() const
      {
        return directIndex_;
      }


//    bool equals(const BSplinePatchEntity& other) const
//    {
//      return target_ == other.target_;
//    }


    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    private:
     const GridViewImp* BSplineGridView_;
     unsigned int directIndex_;

    }; // end of OneDGridEntity codim = 0
  }
}

#endif  //DUNE_IGA_BSPLINEGRIDENTITY_HH
