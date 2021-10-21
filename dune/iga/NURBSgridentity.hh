// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRIDENTITY_HH
#define DUNE_IGA_NURBSGRIDENTITY_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSleafgridview.hh>
#include <numeric>
#include <dune/grid/common/gridenums.hh>






/** \file
 * \brief The NURBSGridEntity class
 */

namespace Dune::IGA
  {
  /** \brief
  *.
  */
  template<int codim, typename GridViewImp>
  class NURBSGridEntity
  {
  public:

      using Geometry = NURBSGeometry<GridViewImp::dimension,GridViewImp::dimensionworld>;
      //! Default Constructor
      NURBSGridEntity ()
        : NURBSGridView_(nullptr)
      {}

      NURBSGridEntity (const GridViewImp& gridView, unsigned int directIndex)
        : NURBSGridView_(&gridView), directIndex_(directIndex),
          parType_{(NURBSGridView_->NURBSpatch_->isBorderElement(directIndex_)? PartitionType::BorderEntity : PartitionType::InteriorEntity)}
      {

      }

      //! Geometry of this entity
      typename GridViewImp::Geometry geometry () const
      {
        auto const &knotElementNet = NURBSGridView_->NURBSpatch_->knotElementNet_;
        auto const &multiIndex = knotElementNet->directToMultiIndex(directIndex_);
        return NURBSGridView_->NURBSpatch_->geometry(multiIndex);
      }

      [[nodiscard]] unsigned int  getIndex() const
      {
        return directIndex_;
      }


      [[nodiscard]] unsigned int  subEntities(unsigned int codim1) const
      {
          if (codim1==0)
              return 0;
          else if (GridViewImp::dimension - codim1 == 0)
              return (1<<GridViewImp::dimension);
          else if (GridViewImp::dimension - codim1 == 1)
              return 3*GridViewImp::dimension+1;
          else if (GridViewImp::dimension - codim1 == 2)
                  return 6;

      }

      template<int codimSub>
      requires(codim==0)  typename GridViewImp::template Codim<codimSub>::Entity subEntity(int i) const
      {
//            if constexpr(codimSub==GridViewImp::dimension) //vertices


      }

      [[nodiscard]] auto type() const
      {
          return GeometryTypes::cube(GridViewImp::dimension-codim);
      }

    [[nodiscard]] PartitionType partitionType() const
    {
          return parType_;
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
      PartitionType parType_;

    }; // end of OneDGridEntity codim = 0
  }

#endif  //DUNE_IGA_NURBSGRIDENTITY_HH
