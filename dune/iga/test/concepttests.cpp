// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/grid/concepts.hh>
#include <dune/iga/nurbsgeometry.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/utils/concepts.hh>

template <typename G>
void checkConcepts() {
  static_assert(Dune::Concept::Grid<G>);

  using GridEntity         = typename G::template Codim<0>::Entity;
  using GridEntityGeometry = typename GridEntity::Geometry;
  using GridView           = typename G::GridView;
  using GlobalIdSet        = typename G::GlobalIdSet;
  using IndexSet           = typename GridView::IndexSet;

  static_assert(Dune::Concept::EntityGeneral<GridEntity>);

  static_assert(Dune::IGA::Concept::NurbsGeometry<typename GridEntityGeometry::Implementation>);
  static_assert(Dune::Concept::GridView<GridView>);
  static_assert(Dune::Concept::IndexSet<IndexSet>);
  static_assert(Dune::Concept::IdSet<GlobalIdSet>);
}

int main() {
  checkConcepts<Dune::IGA::NURBSGrid<1, 1>>();
  checkConcepts<Dune::IGA::NURBSGrid<1, 2>>();
  checkConcepts<Dune::IGA::NURBSGrid<1, 3>>();

  checkConcepts<Dune::IGA::NURBSGrid<2, 2>>();
  checkConcepts<Dune::IGA::NURBSGrid<2, 3>>();

  checkConcepts<Dune::IGA::NURBSGrid<3, 3>>();

  return 0;
}
