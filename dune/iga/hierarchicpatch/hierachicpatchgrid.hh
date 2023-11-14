// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IDENTITYGRID_HH
#define DUNE_GRID_IDENTITYGRID_HH

/** \file
 * \brief The PatchGrid class
 */

#include <string>
#include <map>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/subgrid/subgrid.hh>

// The components of the PatchGrid interface
#include "hierachicpatchgridgeometry.hh"
#include "hierachicpatchgridentity.hh"
#include "hierachicpatchgridentityseed.hh"
#include "hierachicpatchgridintersectioniterator.hh"
#include "hierachicpatchgridleveliterator.hh"
#include "hierachicpatchgridleafiterator.hh"
#include "hierachicpatchgridhierarchiciterator.hh"
#include "hierachicpatchgridindexsets.hh"

namespace Dune::IGA
{
  // Forward declaration
  template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid>
  class PatchGrid;

  // External forward declarations
  template< class Grid >
  struct HostGridAccess;

  template<int dim, int dimworld,  bool trim,typename ScalarType, typename HostGrid=YaspGrid<dim,TensorProductCoordinates<ScalarType,dim>>>
  struct PatchGridFamily
  {

  public:

    typedef GridTraits<
        dim,
        dimworld,
        PatchGrid<dim,dimworld,trim,ScalarType,HostGrid>,
        PatchGridGeometry,
        PatchGridEntity,
        PatchGridLevelIterator,
        PatchGridLeafIntersection,
        PatchGridLevelIntersection,
        PatchGridLeafIntersectionIterator,
        PatchGridLevelIntersectionIterator,
        PatchGridHierarchicIterator,
        PatchGridLeafIterator,
        PatchGridLevelIndexSet< const PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >,
        PatchGridLeafIndexSet< const PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >,
        PatchGridGlobalIdSet< const PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >,
        typename HostGrid::Traits::GlobalIdSet::IdType,
        PatchGridLocalIdSet< const PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >,
        typename HostGrid::Traits::LocalIdSet::IdType,
        Communication<No_Comm>,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        PatchGridEntitySeed,
        PatchGridGeometry,
        typename HostGrid::Traits::LevelIndexSet::IndexType,
        typename HostGrid::Traits::LevelIndexSet::Types,
        typename HostGrid::Traits::LeafIndexSet::IndexType,
        typename HostGrid::Traits::LeafIndexSet::Types
        > Traits;

  };
  template<std::size_t dim,typename ScalarType>
  auto createUniqueKnotSpans(const std::array<std::vector<ScalarType>, dim>& knotSpans) {
    std::array<std::vector<double>, dim>uniqueKnotVector;
    for (int i = 0; i < dim; ++i)  // create unique knotspan vectors
      std::ranges::unique_copy(knotSpans[i], std::back_inserter(uniqueKnotVector[i]),
                               [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

    return uniqueKnotVector;
  }

  template<std::size_t dim,typename ScalarType>
  auto createParameterPatchGridFromKnotVectors(const std::array<std::vector<ScalarType>, dim>& uniqueKnotVector_) {

    auto grid=std::make_unique< YaspGrid<dim,TensorProductCoordinates<ScalarType,dim>>> (uniqueKnotVector_);
    return grid;
  }


  //**********************************************************************
  //
  // --PatchGrid
  //
  //************************************************************************
  /*!
   * \brief Provides a meta grid that is identical to its host
   * \ingroup GridImplementations
   * \ingroup PatchGrid
   *
   * \tparam HostGrid The host grid type wrapped by the PatchGrid
   */
  template<std::size_t dim, std::size_t dimworld,  bool trim=true,typename ScalarType=double, typename HostGrid=YaspGrid<dim,TensorProductCoordinates<ScalarType,dim>>>
  class PatchGrid
  : public GridDefaultImplementation<dim, dimworld,
                                     ScalarType, PatchGridFamily<dim,dimworld,trim,ScalarType, HostGrid> >
  {
    using ParameterSpaceUntrimmedGrid = HostGrid;
    using SubGridParameterSpaceGridView = std::conditional_t<trim,SubGrid<dim,ParameterSpaceUntrimmedGrid>,ParameterSpaceUntrimmedGrid>;

    friend class PatchGridLevelIndexSet<const PatchGrid >;
    friend class PatchGridLeafIndexSet<const PatchGrid >;
    friend class PatchGridGlobalIdSet<const PatchGrid >;
    friend class PatchGridLocalIdSet<const PatchGrid >;
    friend class PatchGridHierarchicIterator<const PatchGrid >;
    friend class PatchGridLevelIntersectionIterator<const PatchGrid >;
    friend class PatchGridLeafIntersectionIterator<const PatchGrid >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class PatchGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class PatchGridLeafIterator;


    template<int codim_, int dim_, class GridImp_>
    friend class PatchGridEntity;

    friend struct HostGridAccess< PatchGrid >;

  public:

    /** \todo Should not be public */
    typedef HostGrid HostGridType;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    using GridFamily=  PatchGridFamily<dim,dimworld,trim,ScalarType,HostGrid>;

    //! the Traits
    using Traits= typename  PatchGridFamily<dim,dimworld,trim,ScalarType,HostGrid>::Traits ;

    //! The type used to store coordinates, inherited from the HostGrid
    typedef typename ParameterSpaceUntrimmedGrid::ctype ctype;


    /** \brief Constructor
     *
     * \param hostgrid The host grid wrapped by the PatchGrid
     */
    explicit PatchGrid(const NURBSPatchData<dim,dimworld,ctype>& patchData) :
      uniqueCoarseKnotSpans(createUniqueKnotSpans(patchData.knotSpans)),
      hostgrid_(createParameterPatchGridFromKnotVectors(uniqueCoarseKnotSpans)),
      leafIndexSet_(*this),
      globalIdSet_(*this),
      localIdSet_(*this)
    {
      setIndices();
      patchGeometries.emplace_back(patchData);
    }

    //! Destructor
    ~PatchGrid()
    {
      // Delete level index sets
      for (size_t i=0; i<levelIndexSets_.size(); i++)
        if (levelIndexSets_[i])
          delete (levelIndexSets_[i]);
    }


    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {
      return hostgrid_->maxLevel();
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return PatchGridLevelIterator<codim,All_Partition, const PatchGrid >(this, level);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return PatchGridLevelIterator<codim,All_Partition, const PatchGrid >(this, level, true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return PatchGridLevelIterator<codim,PiType, const PatchGrid >(this, level);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return PatchGridLevelIterator<codim,PiType, const PatchGrid >(this, level, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return PatchGridLeafIterator<codim,All_Partition, const PatchGrid >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return PatchGridLeafIterator<codim,All_Partition, const PatchGrid >(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return PatchGridLeafIterator<codim,PiType, const PatchGrid >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return PatchGridLeafIterator<codim,PiType, const PatchGrid >(this, true);
    }


    /** \brief Number of grid entities per level and codim
     */
    [[nodiscard]] int size (int level, int codim) const {
      return hostgrid_->size(level,codim);
    }

    /** \brief returns the number of boundary segments within the macro grid
     */
    [[nodiscard]] size_t numBoundarySegments () const {
      return hostgrid_->numBoundarySegments();
    }

    //! number of leaf entities per codim in this process
    [[nodiscard]] int size (int codim) const {
      return leafIndexSet().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return levelIndexSets_[level]->size(type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return globalIdSet_;
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return localIdSet_;
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level < 0 || level > maxLevel())
      {
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      }
      return *levelIndexSets_[level];
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      typedef PatchGridEntity<
        EntitySeed::codimension,
        HostGrid::dimension,
        const typename Traits::Grid
        > EntityImp;

      return EntityImp(this, hostgrid_->entity(seed.impl().hostEntitySeed()));
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount) {
      for (int i = 0; i < refCount; ++i) {
        const auto& finestPatchData = patchGeometries.back().patchData_;
        auto newfinestPatchData = finestPatchData;
        for (int refDirection = 0; refDirection < dim; ++refDirection) {
          auto additionalKnots  = generateRefinedKnots(finestPatchData.knotSpans, refDirection, 1);
          newfinestPatchData = knotRefinement<dim>(newfinestPatchData, additionalKnots, refDirection);
        }
        patchGeometries.push_back(std::move(newfinestPatchData));

      }

      hostgrid_->globalRefine(refCount);
    }

    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     * The parameter is currently ignored
     *
     * \return <ul>
     * <li> true, if marking was successful </li>
     * <li> false, if marking was not possible </li>
     * </ul>
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e)
    {
      return hostgrid_->mark(refCount, getHostEntity<0>(e));
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      return hostgrid_->getMark(getHostEntity<0>(e));
    }

    /** \brief returns true, if at least one entity is marked for adaption */
    bool preAdapt() {
      return hostgrid_->preAdapt();
    }


    //! Triggers the grid refinement process
    bool adapt()
    {
      return hostgrid_->adapt();
    }

    /** \brief Clean up refinement markers */
    void postAdapt() {
      return hostgrid_->postAdapt();
    }

    /*@}*/

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return hostgrid_->leafGridView().overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return hostgrid_->leafGridView().ghostSize(codim);
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return hostgrid_->levelGridView(level).overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return hostgrid_->levelGridView(level).ghostSize(codim);
    }


#if 0
    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "PatchGrid::loadBalance()");
    }
#endif


    /** \brief dummy communication */
    const Communication<No_Comm>& comm () const
    {
      return ccobj;
    }

    /** \brief Communicate data of level gridView */
    template <class DataHandle>
    void communicate (DataHandle& handle, InterfaceType iftype,
                      CommunicationDirection dir, int level) const
    {
      hostgrid_->levelGridView(level).communicate(handle,iftype,dir);
    }

    /** \brief Communicate data of leaf gridView */
    template <class DataHandle>
    void communicate (DataHandle& handle, InterfaceType iftype,
                      CommunicationDirection dir) const
    {
      hostgrid_->leafGridView().communicate(handle,iftype,dir);
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the hostgrid this PatchGrid lives in
    HostGridType& getHostGrid() const
    {
      return *hostgrid_;
    }


    //! Returns the hostgrid entity encapsulated in given PatchGrid entity
    template <int codim>
    const typename HostGrid::Traits::template Codim<codim>::Entity& getHostEntity(const typename Traits::template Codim<codim>::Entity& e) const
    {
      return e.impl().hostEntity_;
    }

  public:
    std::array<std::vector<ScalarType>,dim> uniqueCoarseKnotSpans;
    std::vector<NURBSPatchGeometry<dim,dimworld,ScalarType>> patchGeometries;
  protected:

    //! The host grid which contains the actual grid hierarchy structure
    std::unique_ptr<HostGrid> hostgrid_;

  private:

    //! compute the grid indices and ids
    void setIndices()
    {
      localIdSet_.update();

      globalIdSet_.update();

      // //////////////////////////////////////////
      //   Create the index sets
      // //////////////////////////////////////////
      for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
        auto* p
          = new PatchGridLevelIndexSet<const PatchGrid >();
        levelIndexSets_.push_back(p);
      }

      for (int i=0; i<=maxLevel(); i++)
        if (levelIndexSets_[i])
          levelIndexSets_[i]->update(*this, i);

      leafIndexSet_.update(*this);

    }



    //! \todo Please doc me !
    Communication<No_Comm> ccobj;

    //! Our set of level indices
    std::vector<PatchGridLevelIndexSet<const PatchGrid >*> levelIndexSets_;

    //! \todo Please doc me !
    PatchGridLeafIndexSet<const PatchGrid > leafIndexSet_;

    //! \todo Please doc me !
    PatchGridGlobalIdSet<const PatchGrid > globalIdSet_;

    //! \todo Please doc me !
    PatchGridLocalIdSet<const PatchGrid > localIdSet_;



  }; // end Class PatchGrid


} // namespace Dune::IGA

namespace Dune {


  namespace Capabilities
  {
    /** \brief has entities for some codimensions as host grid
     * \ingroup PatchGrid
     */
    template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid,int codim>
    struct hasEntity<IGA::PatchGrid<dim,dimworld,trim,ScalarType,HostGrid>, codim>
    {
      static const bool v = hasEntity<HostGrid,codim>::v;
    };

    template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid,int codim>
    struct hasEntityIterator<IGA::PatchGrid<dim,dimworld,trim,ScalarType,HostGrid>, codim>
    {
      static const bool v = hasEntityIterator<HostGrid, codim>::v;
    };

    /** \brief PatchGrid can communicate when the host grid can communicate
     *  \ingroup PatchGrid
     */
    template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid,int codim>
    struct canCommunicate<IGA::PatchGrid<dim,dimworld,trim,ScalarType,HostGrid>, codim>
    {
      static const bool v = canCommunicate<HostGrid, codim>::v;
    };

    /** \brief has conforming level grids when host grid has
     * \ingroup PatchGrid
     */
    template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid>
    struct isLevelwiseConforming<IGA::PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >
    {
      static const bool v = isLevelwiseConforming<HostGrid>::v;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup PatchGrid
     */
    template<std::size_t dim, std::size_t dimworld,  bool trim,typename ScalarType, typename HostGrid>
    struct isLeafwiseConforming<IGA::PatchGrid<dim,dimworld,trim,ScalarType,HostGrid> >
    {
      static const bool v = isLeafwiseConforming<HostGrid>::v;
    };
  } // end namespace Capabilities

}

#endif // DUNE_GRID_IDENTITYGRID_HH
