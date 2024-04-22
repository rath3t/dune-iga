// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/grid/common/intersection.hh>
#include <dune/iga/hierarchicpatch/patchgridentity.hh>
#include <dune/iga/hierarchicpatch/patchgridintersections.hh>

/** \file
 * @brief The PatchGridLeafIntersectionIterator and PatchGridLevelIntersectionIterator classes
 */

namespace Dune::IGANEW::IdentityTrim {

/** @brief Iterator over all element neighbors
 * @ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <class GridImp>
class PatchGridLeafIntersectionIterator
{
  constexpr static int dim = GridImp::dimension;

  constexpr static int dimworld = GridImp::dimensionworld;

  // The type used to store coordinates
  typedef typename GridImp::ctype ctype;

  typedef typename GridImp::ParameterSpaceGrid::LeafGridView::IntersectionIterator HostLeafIntersectionIterator;

public:
  typedef Dune::Intersection<const GridImp, PatchGridLeafIntersection<GridImp> > Intersection;

  PatchGridLeafIntersectionIterator() {
  }

  PatchGridLeafIntersectionIterator(const GridImp* parameterSpaceGrid, const HostLeafIntersectionIterator& hostIterator)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostIterator_(hostIterator) {
  }

  //! equality
  bool equals(const PatchGridLeafIntersectionIterator& other) const {
    return hostIterator_ == other.hostIterator_;
  }

  //! prefix increment
  void increment() {
    ++hostIterator_;
  }

  //! @brief dereferencing
  Intersection dereference() const {
    return PatchGridLeafIntersection<GridImp>(parameterSpaceGrid_, *hostIterator_);
  }

private:
  //**********************************************************
  //  private data
  //**********************************************************

  const GridImp* parameterSpaceGrid_         = nullptr;
  HostLeafIntersectionIterator hostIterator_ = {};
};

//! @todo Please doc me !
template <class GridImp>
class PatchGridLevelIntersectionIterator
{
  constexpr static int dim = GridImp::dimension;

  constexpr static int dimworld = GridImp::dimensionworld;

  // The type used to store coordinates
  typedef typename GridImp::ctype ctype;

  typedef typename GridImp::ParameterSpaceGrid::LevelGridView::IntersectionIterator HostLevelIntersectionIterator;

public:
  typedef Dune::Intersection<const GridImp, PatchGridLevelIntersection<GridImp> > Intersection;

  PatchGridLevelIntersectionIterator() {
  }

  PatchGridLevelIntersectionIterator(const GridImp* parameterSpaceGrid,
                                     const HostLevelIntersectionIterator& hostIterator)
      : parameterSpaceGrid_(parameterSpaceGrid),
        hostIterator_(hostIterator) {
  }

  //! equality
  bool equals(const PatchGridLevelIntersectionIterator<GridImp>& other) const {
    return hostIterator_ == other.hostIterator_;
  }

  //! prefix increment
  void increment() {
    ++hostIterator_;
  }

  //! @brief dereferencing
  Intersection dereference() const {
    return PatchGridLevelIntersection<GridImp>(parameterSpaceGrid_, *hostIterator_);
  }

private:
  const GridImp* parameterSpaceGrid_          = nullptr;
  HostLevelIntersectionIterator hostIterator_ = {};
};

} // namespace Dune::IGANEW::IdentityTrim
