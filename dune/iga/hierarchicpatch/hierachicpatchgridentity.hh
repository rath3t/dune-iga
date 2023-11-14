// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDENTITY_HH
#define DUNE_IDENTITYGRIDENTITY_HH

/** \file
 * \brief The PatchGridEntity class
 */

#include <dune/grid/common/grid.hh>

namespace Dune::IGA {


  // Forward declarations

  template<int codim, int dim, class GridImp>
  class PatchGridEntity;

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class PatchGridLevelIterator;

  template<class GridImp>
  class PatchGridLevelIntersectionIterator;

  template<class GridImp>
  class PatchGridLeafIntersectionIterator;

  template<class GridImp>
  class PatchGridHierarchicIterator;


  // External forward declarations
  template< class Grid >
  struct HostGridAccess;


  //**********************************************************************
  //
  // --PatchGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a PatchGrid
   *   \ingroup PatchGrid
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int codim, int dim, class GridImp>
  class PatchGridEntity :
    public EntityDefaultImplementation <codim,dim,GridImp,PatchGridEntity>
  {

    template <class GridImp_>
    friend class PatchGridLevelIndexSet;

    template <class GridImp_>
    friend class PatchGridLeafIndexSet;

    template <class GridImp_>
    friend class PatchGridLocalIdSet;

    template <class GridImp_>
    friend class PatchGridGlobalIdSet;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;


  private:

    typedef typename GridImp::ctype ctype;

    // The codimension of this entity wrt the host grid
    constexpr static int CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim;

    // equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Entity HostGridEntity;


  public:

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

    PatchGridEntity()
      : patchGrid_(nullptr)
    {}

    PatchGridEntity(const GridImp* patchGrid, const HostGridEntity& hostEntity)
      : hostEntity_(hostEntity)
      , patchGrid_(patchGrid)
    {}

    PatchGridEntity(const GridImp* patchGrid, HostGridEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity))
      , patchGrid_(patchGrid)
    {}

    //! \todo Please doc me !
    PatchGridEntity(const PatchGridEntity& original)
      : hostEntity_(original.hostEntity_)
      , patchGrid_(original.patchGrid_)
    {}

    PatchGridEntity(PatchGridEntity&& original)
 noexcept       : hostEntity_(std::move(original.hostEntity_))
      , patchGrid_(original.patchGrid_)
    {}

    //! \todo Please doc me !
    PatchGridEntity& operator=(const PatchGridEntity& original)
    {
      if (this != &original)
      {
        patchGrid_ = original.patchGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    //! \todo Please doc me !
    PatchGridEntity& operator=(PatchGridEntity&& original) noexcept
    {
      if (this != &original)
      {
        patchGrid_ = original.patchGrid_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    bool equals(const PatchGridEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    //! returns true if father entity exists
    bool hasFather () const {
      return hostEntity_.hasFather();
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return EntitySeed(hostEntity_);
    }

    //! level of this element
    int level () const {
      return hostEntity_.level();
    }


    /** \brief The partition type for parallel computing
     */
    PartitionType partitionType () const {
      return hostEntity_.partitionType();
    }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cc) const
    {
      return hostEntity_.subEntities(cc);
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      auto geo = typename Geometry::Implementation( hostEntity_.geometry() ,patchGrid_->patchGeometries[this->level()].localView());
      return Geometry( geo);
    }


    HostGridEntity hostEntity_;

  private:

    const GridImp* patchGrid_;

  };




  //***********************
  //
  //  --PatchGridEntity
  //
  //***********************
  /** \brief Specialization for codim-0-entities.
   * \ingroup PatchGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   */
  template<int dim, class GridImp>
  class PatchGridEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, PatchGridEntity>
  {
    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

  public:

    // The codimension of this entitypointer wrt the host grid
    constexpr static int CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;

    // equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Entity HostGridEntity;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on this level
    typedef PatchGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

    //! The Iterator over intersections on the leaf level
    typedef PatchGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef PatchGridHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;



    PatchGridEntity()
      : patchGrid_(nullptr)
    {}

    PatchGridEntity(const GridImp* patchGrid, const HostGridEntity& hostEntity)
      : hostEntity_(hostEntity)
      , patchGrid_(patchGrid)
    {}

    PatchGridEntity(const GridImp* patchGrid, HostGridEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity))
      , patchGrid_(patchGrid)
    {}

    //! \todo Please doc me !
    PatchGridEntity(const PatchGridEntity& original)
      : hostEntity_(original.hostEntity_)
      , patchGrid_(original.patchGrid_)
    {}

    PatchGridEntity(PatchGridEntity&& original)
 noexcept       : hostEntity_(std::move(original.hostEntity_))
      , patchGrid_(original.patchGrid_)
    {}

    //! \todo Please doc me !
    PatchGridEntity& operator=(const PatchGridEntity& original)
    {
      if (this != &original)
      {
        patchGrid_ = original.patchGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    //! \todo Please doc me !
    PatchGridEntity& operator=(PatchGridEntity&& original)
 noexcept     {
      if (this != &original)
      {
        patchGrid_ = original.patchGrid_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    [[nodiscard]] bool equals(const PatchGridEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    //! returns true if father entity exists
    [[nodiscard]] bool hasFather () const {
      return hostEntity_.hasFather();
    }

    //! Create EntitySeed
    [[nodiscard]] EntitySeed seed () const
    {
      return EntitySeed(hostEntity_);
    }

    //! Level of this element
    [[nodiscard]] int level () const
    {
      return hostEntity_.level();
    }


    /** \brief The partition type for parallel computing */
    [[nodiscard]] PartitionType partitionType () const {
      return hostEntity_.partitionType();
    }


    //! Geometry of this entity
    [[nodiscard]] Geometry geometry () const
    {
      static_assert(std::is_same_v<decltype(patchGrid_->patchGeometries[this->level()].localView()),typename NURBSPatchGeometry<dim,dimworld,ctype>::LocalView>);
      auto geo = typename Geometry::Implementation( hostEntity_.geometry() ,patchGrid_->patchGeometries[this->level()].localView());
      return Geometry(geo);
    }


    /** \brief Return the number of subEntities of codimension codim.
     */
    [[nodiscard]] unsigned int subEntities (unsigned int codim) const
    {
      return hostEntity_.subEntities(codim);
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template<int cc>
    [[nodiscard]] typename GridImp::template Codim<cc>::Entity subEntity (int i) const {
      return PatchGridEntity<cc,dim,GridImp>(patchGrid_, hostEntity_.template subEntity<cc>(i));
    }


    //! First level intersection
    [[nodiscard]] PatchGridLevelIntersectionIterator<GridImp> ilevelbegin () const {
      return PatchGridLevelIntersectionIterator<GridImp>(
        patchGrid_,
        patchGrid_->getHostGrid().levelGridView(level()).ibegin(hostEntity_));
    }


    //! Reference to one past the last neighbor
    PatchGridLevelIntersectionIterator<GridImp> ilevelend () const {
      return PatchGridLevelIntersectionIterator<GridImp>(
        patchGrid_,
        patchGrid_->getHostGrid().levelGridView(level()).iend(hostEntity_));
    }


    //! First leaf intersection
    PatchGridLeafIntersectionIterator<GridImp> ileafbegin () const {
      return PatchGridLeafIntersectionIterator<GridImp>(
        patchGrid_,
        patchGrid_->getHostGrid().leafGridView().ibegin(hostEntity_));
    }


    //! Reference to one past the last leaf intersection
    PatchGridLeafIntersectionIterator<GridImp> ileafend () const {
      return PatchGridLeafIntersectionIterator<GridImp>(
        patchGrid_,
        patchGrid_->getHostGrid().leafGridView().iend(hostEntity_));
    }


    //! returns true if Entity has NO children
    bool isLeaf() const {
      return hostEntity_.isLeaf();
    }


    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::Entity father () const {
      return PatchGridEntity(patchGrid_, hostEntity_.father());
    }


    /** \brief Location of this element relative to the reference element element of the father.
     * This is sufficient to interpolate all dofs in conforming case.
     * Nonconforming may require access to neighbors of father and
     * computations with local coordinates.
     * On the fly case is somewhat inefficient since dofs  are visited several times.
     * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
     * implementation of numerical algorithms is only done for simple discretizations.
     * Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather () const
    {
      return LocalGeometry( hostEntity_.geometryInFather() );
    }


    /** \brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    PatchGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
    {
      return PatchGridHierarchicIterator<const GridImp>(patchGrid_, *this, maxLevel);
    }


    //! Returns iterator to one past the last son
    PatchGridHierarchicIterator<GridImp> hend (int maxLevel) const
    {
      return PatchGridHierarchicIterator<const GridImp>(patchGrid_, *this, maxLevel, true);
    }


    //! \todo Please doc me !
    bool wasRefined () const
    {
      if (patchGrid_->adaptationStep!=GridImp::adaptDone)
        return false;

      int level = this->level();
      int index = patchGrid_->levelIndexSet(level).index(*this);
      return patchGrid_->refinementMark_[level][index];
    }


    //! \todo Please doc me !
    bool mightBeCoarsened () const
    {
      return true;
    }


    // /////////////////////////////////////////
    //   Internal stuff
    // /////////////////////////////////////////


    HostGridEntity hostEntity_;
    const GridImp* patchGrid_;

  private:

    typedef typename GridImp::ctype ctype;

  }; // end of PatchGridEntity codim = 0


} // namespace Dune


#endif
