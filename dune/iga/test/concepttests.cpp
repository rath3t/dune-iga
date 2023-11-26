// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/grid/concepts.hh>
// #include <dune/iga/hierarchicpatch/hierachicpatchgridgeometry.hh>
#include <dune/common/classname.hh>

#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/iga/hierarchicpatch/concepts.hh>
#include <dune/iga/hierarchicpatch/hierachicpatchgrid.hh>
template < typename>
struct IsDefaultReferenceElement : std::false_type {};

template < typename S>
struct IsDefaultReferenceElement<Dune::Geo::ReferenceElement<S>> : std::true_type {};

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity         = typename G::template Codim<0>::Entity;
  using LeafGridView       = typename G::LeafGridView;
  using LevelGridView      = typename G::LevelGridView;
  using GlobalIdSet        = typename G::GlobalIdSet;
  using IndexSet           = typename LeafGridView::IndexSet;
  static constexpr Trimming trim = G::trim;


  std::cout << Dune::className<GridEntity>() << std::endl;
  using GridEntityReferenceType= decltype(referenceElement(GridEntity()));
  if constexpr (trim==Trimming::Disabled)
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

  //Check concepts with trim
  checkConcepts<Dune::IGANEW::PatchGrid<1, 1,Trimming::Enabled>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 2,Trimming::Enabled>>();
  checkConcepts<Dune::IGANEW::PatchGrid<1, 3,Trimming::Enabled>>();

  checkConcepts<Dune::IGANEW::PatchGrid<2, 2,Trimming::Enabled>>();
  checkConcepts<Dune::IGANEW::PatchGrid<2, 3,Trimming::Enabled>>();

  checkConcepts<Dune::IGANEW::PatchGrid<3, 3,Trimming::Enabled>>();

  return 0;
}
