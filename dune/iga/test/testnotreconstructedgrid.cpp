// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#define DUNE_CHECK_BOUNDS
#define CHECK_RESERVEDVECTOR
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/subgrid/test/common.hh>

using namespace Dune::IGA;
using namespace Dune;

auto test(int refinement) {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  std::cout << "\n*********************\n";
  std::cout << "Reinement " << refinement << std::endl;
  std::cout << "*********************\n\n";

  using PatchGrid   = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(typename GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/quarter_plate.ibra", true, {refinement, refinement});

  auto grid = gridFactory.createGrid();

  gridcheck(*grid);
  checkIntersectionIterator(*grid);
  checkLeafIntersections(*grid);

  using namespace Functions::BasisFactory;

  auto gridView = grid->leafGridView();

  std::cout << "\nTesting NURBS basis\n";
  auto basis = makeBasis(gridView, nurbs());
  t.subTest(checkBasis(basis, EnableContinuityCheck()));

  std::cout << "\nTesting Lagrange (1) basis\n";
  auto lagrangeBasis1 = makeBasis(gridView, lagrange<1>());
  t.subTest(checkBasis(lagrangeBasis1, EnableContinuityCheck()));

  std::cout << "\nTesting Lagrange (2) basis\n";
  auto lagrangeBasis2 = makeBasis(gridView, lagrange<2>());
  t.subTest(checkBasis(lagrangeBasis2, EnableContinuityCheck()));

  return t;
}

int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  Preferences::getInstance().targetAccuracy(1e-3);
  Preferences::getInstance().reconstructTrimmedLocalGeometry(false);
  Preferences::getInstance().reportTrimmedElementGeometryTypeAsNone(false);

  t.subTest(test(0));
  t.subTest(test(1));
  t.subTest(test(2));
  t.subTest(test(3));

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
