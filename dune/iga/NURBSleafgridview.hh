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


      /** \brief Collect several types associated to OneDGrid LeafGridViews */
      template< class GridImp>
      struct NurbsLeafGridViewTraits
      {
          typedef NurbsLeafGridViewTraits< GridImp > GridViewImp;

          /** \brief type of the grid */
          typedef typename std::remove_const<GridImp>::type Grid;

//          /** \brief type of the index set */
//          typedef typename GridImpl :: Traits :: LeafIndexSet IndexSet;
//
//          /** \brief type of the intersection */
//          typedef typename GridImpl :: Traits :: LeafIntersection Intersection;
//
//          /** \brief type of the intersection iterator */
//          typedef typename GridImpl :: Traits :: LeafIntersectionIterator IntersectionIterator;
//
//          /** \brief type of the collective communication */
//          typedef typename GridImpl :: Traits :: CollectiveCommunication CollectiveCommunication;

          template< int cd >
          struct Codim
          {
              typedef typename Grid :: Traits
              :: template Codim< cd > :: template Partition< All_Partition > :: LeafIterator
                      Iterator;

              typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;

              typedef typename Grid :: Traits:: template Codim< cd > :: Geometry Geometry;
              typedef typename Grid :: Traits :: template Codim< cd > :: Geometry LocalGeometry;

              /** \brief Define types needed to iterate over entities of a given partition type */
              template <PartitionIteratorType pit >
              struct Partition
              {
                  /** \brief iterator over a given codim and partition type */
                  typedef typename Grid :: Traits :: template Codim< cd >
                  :: template Partition< pit > :: LeafIterator
                          Iterator;
              };
          };

          enum { conforming = true };
      };




    /** \brief NURBS grid manager */
    template<typename GridImpl>
    class NURBSLeafGridView
    {
    public:

      template<int codim, class GridViewImp>
      friend class NURBSGridEntity;

      using Traits = NurbsLeafGridViewTraits<GridImpl>;

      using ctype = double;
      static constexpr int dimension = GridImpl::dimension;
      static constexpr int dimensionworld = GridImpl::dimensionworld;

      template<int codim, typename NURBSGridView, typename NURBSEntity>
      friend class NURBSGridLeafIterator;

      using Grid = typename Traits::Grid;
      typedef NURBSLeafGridView<GridImpl> NURBSGridView;
      typedef NURBSGeometry<dimension, dimensionworld> Geometry;
      typedef NURBSGridLeafIndexSet<NURBSGridView> IndexSet;

        template< int cd >
        struct Codim : public Traits :: template Codim<cd> {};

      /** \brief  constructor
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      NURBSLeafGridView(const std::array<std::vector<double>,dimension>& knotSpans,
                   const MultiDimensionNet<dimensionworld,dimensionworld> controlPoints,
                   const MultiDimensionNet<dimensionworld,1> weights,
                   const std::array<int,dimension> order)
      : NURBSpatch_(std::make_shared<NURBSPatch<dimension,dimensionworld>>(knotSpans, controlPoints, weights, order))
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

        /** \brief obtain collective communication object */
        const CollectiveCommunication &comm () const
        {
            return grid().comm();
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
       std::shared_ptr <NURBSPatch<dimension,dimensionworld>> NURBSpatch_;
       std::vector<std::shared_ptr<NURBSGridEntity<0, NURBSGridView>>> entityVector_;
    };
  }
}


#endif  // DUNE_IGA_NURBSGRID_LEAFGRIDVIEW_HH
