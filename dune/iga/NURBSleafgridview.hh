// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_NURBSGRID_LEAFGRIDVIEW_HH
#define DUNE_IGA_NURBSGRID_LEAFGRIDVIEW_HH

#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSgridentity.hh>
#include <dune/iga/NURBSgridleafiterator.hh>
#include <dune/iga/NURBSgridindexsets.hh>

namespace Dune
{
  namespace IGA
  {
    /** \brief NURBS grid manager */
    template<int dim, int dimworld>
    class NURBSLeafGridView
    {
    public:

      template<int codim, class GridViewImp>
      friend class NURBSGridEntity;


      template<int codim, typename NURBSGridView, typename NURBSEntity>
      friend class NURBSGridLeafIterator;

      typedef NURBSLeafGridView<dim, dimworld> NURBSGridView;
      typedef NURBSGeometry<dim, dimworld> Geometry;
      typedef NURBSGridLeafIndexSet<NURBSGridView> IndexSet;

      template<int codim>
      struct Codim
      {
        typedef NURBSGridEntity<codim,NURBSGridView> Entity;
        typedef NURBSGridLeafIterator<codim,NURBSGridView, Entity> Iterator;

      };

      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      NURBSLeafGridView(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNet<dim,dimworld> controlPoints,
                   const MultiDimensionNet<dim,1> weights,
                   const std::array<int,dim> order)
      : NURBSpatch_(std::make_shared<NURBSPatch<dim,dimworld>>(knotSpans, controlPoints, weights, order))
      {
        int elementSize = NURBSpatch_ ->knotElementNet_->directSize();
        entityVector_.reserve(elementSize);

        //Fill the element vector of codim 0
        for(unsigned int i=0; i<elementSize; ++i)
        {
          entityVector_.push_back(std::make_shared<NURBSGridEntity<0, NURBSGridView>>(*this, i));
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



      typename Codim<0>::Iterator begin () const
      {
        return NURBSGridLeafIterator<0,NURBSGridView, typename Codim<0>::Entity>(*this, 0);
      }

      typename Codim<0>::Iterator end () const
      {
        int elementSize = entityVector_.size();
        return NURBSGridLeafIterator<0,NURBSGridView, typename Codim<0>::Entity>(*this, elementSize);
      }

      IndexSet indexSet() const
      {
        return NURBSGridLeafIndexSet<NURBSGridView>(*this);
      }

    private:
       std::shared_ptr <NURBSPatch<dim,dimworld>> NURBSpatch_;
       std::vector<std::shared_ptr<NURBSGridEntity<0, NURBSGridView>>> entityVector_;
    };
  }
}


#endif  // DUNE_IGA_NURBSGRIDVIEW_HH
