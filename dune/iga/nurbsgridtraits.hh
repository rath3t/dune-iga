//
// Created by lex on 29/11/2021.
//

#pragma

#include <dune/grid/common/gridenums.hh>
#include <dune/iga/nurbsleafgridview.hh>

/** \brief This class is a copy of GridTraits dune-grid/grid/common/grid.hh.
 * This class removes all fascade class from GridTraits since some exported types do not work with c++20
 * ranges ,eg.  Dune::IteratorRange<Dune::EntityIterator<0, con*/
template <int dim, int dimw, class GridImp,
          template<int,int,class> class GeometryImp,
          template<int,int,class> class EntityImp,
          template<int,Dune::PartitionIteratorType,class> class LevelIteratorImp,
          template<class> class LeafIntersectionImp,
          template<class> class LevelIntersectionImp,
          template<class> class LeafIntersectionIteratorImp,
          template<class> class LevelIntersectionIteratorImp,
          template<class> class HierarchicIteratorImp,
          template<int,Dune::PartitionIteratorType,class> class LeafIteratorImp,
          class LevelIndexSetImp, class LeafIndexSetImp,
          class GlobalIdSetImp, class GIDType, class LocalIdSetImp, class LIDType, class CCType,
          template<class> class LevelGridViewTraits,
          template<class> class LeafGridViewTraits,
          template<int,class> class EntitySeedImp,
          template<int,int,class> class LocalGeometryImp = GeometryImp
          >
struct NurbsGridTraits
{
  /** \brief The type that implements the grid. */
  typedef GridImp Grid;

  /** \brief The type of the intersection at the leafs of the grid. */
  typedef  LeafIntersectionImp<  GridImp >  LeafIntersection;
  /** \brief The type of the intersection at the levels of the grid. */
  typedef LevelIntersectionImp<  GridImp >  LevelIntersection;
  /** \brief The type of the intersection iterator at the leafs of the grid. */
  typedef LeafIntersectionIteratorImp<  GridImp >  LeafIntersectionIterator;
  /** \brief The type of the intersection iterator at the levels of the grid. */
  typedef  LeafIntersectionIteratorImp<  GridImp >  LevelIntersectionIterator;

  /** \brief The type of the  hierarchic iterator. */
  typedef  HierarchicIteratorImp<  GridImp > HierarchicIterator;

  /**
     * \brief Traits associated with a specific codim.
     * \tparam cd The codimension.
   */
  template <int cd>
  struct Codim
  {
  public:
    typedef GridImp GeometryImpl;
    typedef GridImp LocalGeometryImpl;
    //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
    /** \brief The type of the geometry associated with the entity.*/
    typedef  GeometryImp<dim-cd,dimw,GridImp> Geometry;
    /** \brief The type of the local geometry associated with the entity.*/
    typedef LocalGeometryImp<dim-cd, dim,  GridImp> LocalGeometry;
    /** \brief The type of the entity. */
    // we could - if needed - introduce another struct for dimglobal of Geometry
    typedef EntityImp<cd, dim,  GridImp> Entity;

    /** \brief The type of the entity seed of this codim.*/
    typedef EntitySeedImp<cd,  GridImp > EntitySeed;

    /**
       * \brief Traits associated with a specific grid partition type.
       * \tparam pitype The type of the grid partition.
     */
    template <Dune::PartitionIteratorType pitype>
    struct Partition
    {
      /** \brief The type of the iterator over the level entities of this codim on this partition. */
      typedef LevelIteratorImp< cd, pitype,  GridImp >  LevelIterator;
      /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
      typedef  LeafIteratorImp< cd, pitype,  GridImp  > LeafIterator;
    };

    /** \brief The type of the iterator over all leaf entities of this codim. */
    typedef typename Partition< Dune::All_Partition >::LeafIterator LeafIterator;

    /** \brief The type of the entity pointer for entities of this codim.*/
    typedef typename Partition< Dune::All_Partition >::LevelIterator LevelIterator;

  private:
    friend class EntityImp<cd, dim,  GridImp>;
  };

  /** \brief type of view for leaf grid */
  typedef typename Dune::IGA::NURBSLeafGridView<GridImp>  LeafGridView;
  /** \brief type of view for level grid */
  typedef typename Dune::IGA::NURBSLeafGridView<GridImp> LevelGridView;

  /** \brief The type of the level index set. */
  typedef LevelIndexSetImp LevelIndexSet;
  /** \brief The type of the leaf index set. */
  typedef LeafIndexSetImp LeafIndexSet;
  /** \brief The type of the global id set. */
  typedef GlobalIdSetImp GlobalIdSet;
  /** \brief The type of the local id set. */
  typedef LocalIdSetImp LocalIdSet;

  /** \brief The type of the collective communication. */
  typedef CCType CollectiveCommunication;
};