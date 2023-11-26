// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/grid/concepts.hh>
// #include <dune/iga/hierarchicpatch/hierachicpatchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/concepts.hh>
#include <dune/iga/hierarchicpatch/hierachicpatchgrid.hh>
#include <dune/common/classname.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
  template <template <typename> class Type, typename>
  struct IsDefaultReferenceElement: std::false_type {};

  template <template <typename> class Type, typename S>
  struct IsDefaultReferenceElement<Type, ReferenceElement<S>> : std::true_type {};

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity         = typename G::template Codim<0>::Entity;
  using GridEntityGeometry = typename GridEntity::Geometry;
  using LeafGridView       = typename G::LeafGridView;
  using LevelGridView      = typename G::LevelGridView;
  using GlobalIdSet        = typename G::GlobalIdSet;
  using IndexSet           = typename LeafGridView::IndexSet;



std::cout<<Dune::className<GridEntity>()<<std::endl;
static_assert(std::IsDefaultReferenceElement<decltype(referenceElement(GridEntity())>::value);



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

  return 0;
}
