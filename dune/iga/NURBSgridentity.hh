// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRIDENTITY_HH
#define DUNE_IGA_NURBSGRIDENTITY_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>

/** \file
 * \brief The NURBSGridEntity class
 */

namespace Dune::IGA
  {
  /** \brief
  *.
  */
  template<int codim, class GridViewImp>
  class NURBSGridEntity
  {
  public:

      using Geometry = NURBSGeometry<GridViewImp::dimension,GridViewImp::dimensionworld>;
      //! Default Constructor
      NURBSGridEntity ()
        : NURBSGridView_(nullptr)
      {}

      NURBSGridEntity (const GridViewImp& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView), directIndex_(directIndex)
      {}

      //! Geometry of this entity
      typename GridViewImp::Geometry geometry () const
      {
        auto const &knotElementNet = NURBSGridView_->NURBSpatch_->knotElementNet_;
        auto const &multiIndex = knotElementNet->directToMultiIndex(0);
        return NURBSGridView_->NURBSpatch_->geometry(multiIndex);
      }

      [[nodiscard]] unsigned int  getIndex() const
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
     const GridViewImp* NURBSGridView_;
     unsigned int directIndex_;

    }; // end of OneDGridEntity codim = 0
  }

#endif  //DUNE_IGA_NURBSGRIDENTITY_HH
