// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include "hierachicpatchgridintersections.hh"
#include "hierachicpatchgridentity.hh"

#include <dune/grid/common/intersection.hh>

/** \file
 * \brief The PatchGridLeafIntersectionIterator and PatchGridLevelIntersectionIterator classes
 */

namespace Dune::IGANEW {

  /** \brief Iterator over all element neighbors
   * \ingroup PatchGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class PatchGridLeafIntersectionIterator
  {

    constexpr static int dim = GridImp::dimension;

    constexpr static int dimworld = GridImp::dimensionworld;

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::LeafGridView::IntersectionIterator HostLeafIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, PatchGridLeafIntersection<GridImp> > Intersection;

    PatchGridLeafIntersectionIterator()
    {}

    PatchGridLeafIntersectionIterator(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid)
      , hostIterator_(hostIterator)
    {}

    //! equality
    bool equals(const PatchGridLeafIntersectionIterator& other) const {
      return hostIterator_ == other.hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++hostIterator_;
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return PatchGridLeafIntersection<GridImp>(identityGrid_,*hostIterator_);
    }

  private:
    //**********************************************************
    //  private data
    //**********************************************************

    const GridImp* identityGrid_ = nullptr;
    HostLeafIntersectionIterator hostIterator_ = {};
  };




  //! \todo Please doc me !
  template<class GridImp>
  class PatchGridLevelIntersectionIterator
  {
    constexpr static int dim = GridImp::dimension;

    constexpr static int dimworld = GridImp::dimensionworld;

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::HostGridType::LevelGridView::IntersectionIterator HostLevelIntersectionIterator;

  public:

    typedef Dune::Intersection<const GridImp, PatchGridLevelIntersection<GridImp> > Intersection;

    PatchGridLevelIntersectionIterator()
    {}

    PatchGridLevelIntersectionIterator(const GridImp* identityGrid,
                                          const HostLevelIntersectionIterator& hostIterator)
      : identityGrid_(identityGrid)
      , hostIterator_(hostIterator)
    {}

    //! equality
    bool equals(const PatchGridLevelIntersectionIterator<GridImp>& other) const {
      return hostIterator_ == other.hostIterator_;
    }


    //! prefix increment
    void increment() {
      ++hostIterator_;
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return PatchGridLevelIntersection<GridImp>(identityGrid_,*hostIterator_);
    }

  private:


    const GridImp* identityGrid_ = nullptr;
    HostLevelIntersectionIterator hostIterator_ = {};

  };


}  // namespace Dune
