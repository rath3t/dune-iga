// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/grid/concepts.hh>
// #include <dune/iga/hierarchicpatch/hierachicpatchgridgeometry.hh>
#include <dune/common/classname.hh>

#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/iga/hierarchicpatch/concepts.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/defaulttrimmer.hh>
template <typename>
struct IsDefaultReferenceElement : std::false_type {};

template <typename S>
struct IsDefaultReferenceElement<Dune::Geo::ReferenceElement<S>> : std::true_type {};

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity    = typename G::template Codim<0>::Entity;
  using LeafGridView  = typename G::LeafGridView;
  using LevelGridView = typename G::LevelGridView;
  using GlobalIdSet   = typename G::GlobalIdSet;
  using IndexSet      = typename LeafGridView::IndexSet;
  using TrimmerType   = typename G::TrimmerType;

  using GridEntityReferenceType = decltype(referenceElement(GridEntity()));
  if constexpr (TrimmerType::isAlwaysTrivial)
    static_assert(IsDefaultReferenceElement<GridEntityReferenceType>::value);
  else
    static_assert(not IsDefaultReferenceElement<GridEntityReferenceType>::value);

  static_assert(Dune::Concept::EntityGeneral<GridEntity>);

  // static_assert(Dune::IGANEW::Concept::NurbsGeometry<typename GridEntityGeometry::Implementation>);
  static_assert(Dune::Concept::GridView<LeafGridView>);
  static_assert(Dune::Concept::GridView<LevelGridView>);
  static_assert(Dune::Concept::IndexSet<IndexSet>);
  static_assert(Dune::Concept::IdSet<GlobalIdSet>);
}

int main() {
  checkConcepts<Dune::IGANEW::PatchGrid<1, 1>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 2>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 3>>();

  checkConcepts<Dune::IGANEW::PatchGrid<2, 2>>();
  checkConcepts<Dune::IGANEW::PatchGrid<2, 3>>();

  checkConcepts<Dune::IGANEW::PatchGrid<3, 3>>();

  // Check concepts with trim
  checkConcepts<Dune::IGANEW::PatchGrid<1, 1, Dune::IGANEW::Trim::DefaultTrimmer<1>>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 2, Dune::IGANEW::Trim::DefaultTrimmer<1>>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 3, Dune::IGANEW::Trim::DefaultTrimmer<1>>>();

  checkConcepts<Dune::IGANEW::PatchGrid<2, 2, Dune::IGANEW::Trim::DefaultTrimmer<2>>>();
  checkConcepts<Dune::IGANEW::PatchGrid<2, 3, Dune::IGANEW::Trim::DefaultTrimmer<2>>>();

  checkConcepts<Dune::IGANEW::PatchGrid<3, 3, Dune::IGANEW::Trim::DefaultTrimmer<3>>>();

  using Grid23 = Dune::IGANEW::PatchGrid<2, 3, Dune::IGANEW::Trim::NoOpTrimmer<2>>;
  std::cout << Dune::className<
      typename Grid23::Codim<0>::LocalGeometry::Implementation::GeometryLocalView::ParameterSpaceGeometry>()
            << std::endl;
  std::cout << Dune::className<
      typename Grid23::Codim<1>::LocalGeometry::Implementation::GeometryLocalView::ParameterSpaceGeometry>()
            << std::endl;
  std::cout << Dune::className<
      typename Grid23::Codim<2>::LocalGeometry::Implementation::GeometryLocalView::ParameterSpaceGeometry>()
            << std::endl;

  return 0;
}
