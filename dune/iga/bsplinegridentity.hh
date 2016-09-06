// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEGRIDENTITY_HH
#define DUNE_IGA_BSPLINEGRIDENTITY_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/bsplineleafgridview.hh>
#include <dune/grid/common/gridenums.hh>


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

      //typedef BSplineGeometry<dim, dimworld> Geometry;
      typedef typename GridViewImp::template Codim<codim>::Geometry Geometry;

      //! Default Constructor
      BSplineGridEntity ()
        : BSplineGridView_(nullptr)
      {}

      BSplineGridEntity (const GridViewImp& gridView, unsigned int directIndex)
        : BSplineGridView_(&gridView), directIndex_(directIndex)
      {}

      //! return the element type identifier (segment)
      GeometryType type () const {return this->geometry().type();}

      //! Geometry of this entity
      Geometry geometry () const
      {
        auto const &knotElementNet = BSplineGridView_->BSplinepatch_->knotElementNet_;
        auto const &multiIndex = knotElementNet->directToMultiIndex(0);
        return BSplineGridView_->BSplinepatch_->geometry(multiIndex);
      }

      unsigned int  getIndex() const
      {
        return directIndex_;
      }

      //! only interior entities
      PartitionType partitionType () const { return InteriorEntity; }

      //! Return the number of subentities of codimension codim
      unsigned int subEntities (unsigned int cd) const
      {
        assert(cd==0 || cd==1);
        return (cd==0) ? 1 : 4;
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
