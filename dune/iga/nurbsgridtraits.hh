// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "dune/iga/nurbsleafgridview.hh"
#include <dune/grid/common/gridenums.hh>

/** \brief This class is a modified copy of GridTraits dune-grid/grid/common/grid.hh. */
template <int dim, int dimw, class GridImp, template <int, int, class> class GeometryImp,
          template <int, int, class> class EntityImp,
          template <int, Dune::PartitionIteratorType, class> class LevelIteratorImp,
          template <class> class LeafIntersectionImp, template <class> class LevelIntersectionImp,
          template <class> class LeafIntersectionIteratorImp, template <class> class LevelIntersectionIteratorImp,
          template <class> class HierarchicIteratorImp,
          template <int, Dune::PartitionIteratorType, class> class LeafIteratorImp, class LevelIndexSetImp,
          class LeafIndexSetImp, class GlobalIdSetImp, class GIDType, class LocalIdSetImp, class LIDType, class CCType,
          template <class> class LevelGridViewTraits, template <class> class LeafGridViewTraits,
          template <int, class> class EntitySeedImp, template <int, int, class> class LocalGeometryImp = GeometryImp>
struct NurbsGridTraits {
  /** \brief The type that implements the grid. */
  typedef GridImp Grid;

  /** \brief The type of the intersection at the leafs of the grid. */
  using LeafIntersection = Dune::Intersection<const GridImp, LeafIntersectionImp<const GridImp>>;
  /** \brief The type of the intersection at the levels of the grid. */
  using LevelIntersection = Dune::Intersection<const GridImp, LeafIntersectionImp<const GridImp>>;
  /** \brief The type of the intersection iterator at the leafs of the grid. */
  using LeafIntersectionIterator = Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorImp<const GridImp>,
                                                              LeafIntersectionImp<const GridImp>>;
  /** \brief The type of the intersection iterator at the levels of the grid. */
  using LevelIntersectionIterator
      = Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorImp<const GridImp>,
                                   LeafIntersectionImp<const GridImp>>;

  /** \brief The type of the  hierarchic iterator. */
  using HierarchicIterator = Dune::EntityIterator<0, const GridImp, HierarchicIteratorImp<const GridImp>>;

  /**
   * \brief Traits associated with a specific codim.
   * \tparam cd The codimension.
   */
  template <int cd>
  struct Codim {
   public:
    //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
    /** \brief The type of the geometry associated with the entity.*/
    using Geometry = Dune::Geometry<dim - cd, dimw, const GridImp, GeometryImp>;
    /** \brief The type of the local geometry associated with the entity.*/
    using LocalGeometry = Dune::Geometry<dim - cd, dim, const GridImp, LocalGeometryImp>;
    /** \brief The type of the entity. */
    // we could - if needed - introduce another struct for dimglobal of Geometry
    using Entity = Dune::Entity<cd, dim, const GridImp, EntityImp>;

    /** \brief The type of the entity seed of this codim.*/
    using EntitySeed = Dune::EntitySeed<const GridImp, EntitySeedImp<cd, const GridImp>>;

    /**
     * \brief Traits associated with a specific grid partition type.
     * \tparam pitype The type of the grid partition.
     */
    template <Dune::PartitionIteratorType pitype>
    struct Partition {
      /** \brief The type of the iterator over the level entities of this codim on this partition. */
      using LevelIterator = Dune::EntityIterator<cd, const GridImp, LevelIteratorImp<cd, pitype, const GridImp>>;
      /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
      using LeafIterator = Dune::EntityIterator<cd, const GridImp, LeafIteratorImp<cd, pitype, const GridImp>>;
    };

    /** \brief The type of the iterator over all leaf entities of this codim. */
    typedef typename Partition<Dune::All_Partition>::LeafIterator LeafIterator;

    /** \brief The type of the entity pointer for entities of this codim.*/
    typedef typename Partition<Dune::All_Partition>::LevelIterator LevelIterator;

   private:
    friend class Dune::Entity<cd, dim, GridImp, EntityImp>;
  };

  /** \brief type of view for leaf grid */
  using LeafGridView = Dune::GridView<LeafGridViewTraits<const GridImp>>;
  /** \brief type of view for level grid */
  //  typedef typename Dune::IGA::NURBSLeafGridView<const GridImp> LevelGridView;
  using LevelGridView = Dune::GridView<LevelGridViewTraits<const GridImp>>;
  /** \brief The type of the level index set. */
  //  typedef LevelIndexSetImp LevelIndexSet;
  using LevelIndexSet = Dune::IndexSet<const GridImp, LevelIndexSetImp, GIDType, std::vector<Dune::GeometryType>>;
  /** \brief The type of the leaf index set. */
  using LeafIndexSet = Dune::IndexSet<const GridImp, LevelIndexSetImp, GIDType, std::vector<Dune::GeometryType>>;
  /** \brief The type of the global id set. */
  //  typedef GlobalIdSetImp GlobalIdSet;
  using GlobalIdSet = Dune::IdSet<const GridImp, GlobalIdSetImp, GIDType>;
  /** \brief The type of the local id set. */
  //  typedef LocalIdSetImp LocalIdSet;
  using LocalIdSet = Dune::IdSet<const GridImp, LocalIdSetImp, GIDType>;

  /** \brief The type of the collective communication. */
  typedef CCType CollectiveCommunication;
  typedef CCType Communication;

  using SubGridType = Dune::IGA::TrimmedSubGrid<2>;
};
