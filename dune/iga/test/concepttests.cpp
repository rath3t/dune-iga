// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/geometrykernel/makecirculararc.hh>
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/parameterspace/concepts.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/subgrid/test/common.hh>

template <typename>
struct IsDefaultReferenceElement : std::false_type
{
};

template <typename S>
struct IsDefaultReferenceElement<Dune::Geo::ReferenceElement<S>> : std::true_type
{
};

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity     = typename G::template Codim<0>::Entity;
  using LeafGridView   = typename G::LeafGridView;
  using LevelGridView  = typename G::LevelGridView;
  using GlobalIdSet    = typename G::GlobalIdSet;
  using IndexSet       = typename LeafGridView::IndexSet;
  using ParameterSpace = typename G::ParameterSpace;

  using GridEntityReferenceType = decltype(referenceElement(GridEntity()));
  if constexpr (ParameterSpace::isAlwaysTrivial)
    static_assert(IsDefaultReferenceElement<GridEntityReferenceType>::value);
  // todo
  // else
  //   static_assert(not IsDefaultReferenceElement<GridEntityReferenceType>::value);

  static_assert(Dune::Concept::EntityGeneral<GridEntity>);

  // static_assert(Dune::IGA::Concept::NurbsGeometry<typename GridEntityGeometry::Implementation>);
  static_assert(Dune::Concept::GridView<LeafGridView>);
  static_assert(Dune::Concept::GridView<LevelGridView>);
  static_assert(Dune::Concept::IndexSet<IndexSet>);
  static_assert(Dune::Concept::IdSet<GlobalIdSet>);
}

int main() {
  using namespace Dune::IGA;
  checkConcepts<PatchGrid<1, 1>>();
  checkConcepts<PatchGrid<1, 2>>();
  checkConcepts<PatchGrid<1, 3>>();

  checkConcepts<PatchGrid<2, 2>>();
  checkConcepts<PatchGrid<2, 3>>();

  checkConcepts<PatchGrid<3, 3>>();

  // Check concepts with trim

  checkConcepts<PatchGrid<2, 2, DefaultParameterSpace::PatchGridFamily>>();
  checkConcepts<PatchGrid<2, 3, DefaultParameterSpace::PatchGridFamily>>();
  // checkConcepts<PatchGrid<3, 3, DefaultParameterSpace::PatchGridFamily>>();

  static_assert(
      not Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<1, 1, double>::ParameterSpace>);
  static_assert(
      not Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<1, 2, double>::ParameterSpace>);
  static_assert(
      not Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<1, 3, double>::ParameterSpace>);
  // static_assert(not Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<3, 3,
  // double>::ParameterSpace>);

  static_assert(
      Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<2, 2, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<DefaultParameterSpace::PatchGridFamily<2, 3, double>::ParameterSpace>);

  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<1, 1, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<1, 2, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<1, 3, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<2, 2, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<2, 3, double>::ParameterSpace>);
  static_assert(
      Dune::IGA::Concept::ParameterSpace<IdentityParameterSpace::PatchGridFamily<3, 3, double>::ParameterSpace>);

  return 0;
}
