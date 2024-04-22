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

namespace Dune::IGANEW::DefaultTrim {

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
  using ParameterSpaceLeafIntersection = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection;
  using LeafIntersection               = typename GridImp::Traits::LeafIntersection;

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
  LeafIntersection dereference() const {
    auto parameterspaceIntersection = ParameterSpaceLeafIntersection(parameterSpaceGrid_, *hostIterator_);
    auto realIntersection = typename LeafIntersection::Implementation(parameterSpaceGrid_, parameterspaceIntersection);
    return LeafIntersection(realIntersection);
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
  using ParameterSpaceLevelIntersection = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection;
  using LevelIntersection               = typename GridImp::Traits::LevelIntersection;

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
  LevelIntersection dereference() const {
    auto parameterspaceIntersection = ParameterSpaceLevelIntersection(parameterSpaceGrid_, *hostIterator_);
    auto realIntersection = typename LevelIntersection::Implementation(parameterSpaceGrid_, parameterspaceIntersection);
    return LevelIntersection(realIntersection);
  }

private:
  const GridImp* parameterSpaceGrid_          = nullptr;
  HostLevelIntersectionIterator hostIterator_ = {};
};

} // namespace Dune::IGANEW::DefaultTrim
