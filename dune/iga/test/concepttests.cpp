// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/grid/concepts.hh>
// #include <dune/iga/hierarchicpatch/hierachicpatchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/concepts.hh>
#include <dune/iga/hierarchicpatch/hierachicpatchgrid.hh>

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity         = typename G::template Codim<0>::Entity;
  using GridEntityGeometry = typename GridEntity::Geometry;
  using LeafGridView       = typename G::LeafGridView;
  using LevelGridView      = typename G::LevelGridView;
  using GlobalIdSet        = typename G::GlobalIdSet;
  using IndexSet           = typename LeafGridView::IndexSet;

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
