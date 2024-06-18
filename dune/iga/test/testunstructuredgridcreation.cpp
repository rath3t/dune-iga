// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "testhelper.hh"

#include <cfenv>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/uggrid.hh>
#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/io/griddrawer.hh>
#include <dune/iga/io/vtk/igadatacollector.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune::IGA;

auto testUnstructuredGridCreation() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  using PatchGrid   = PatchGrid<2, 2, DefaultParameterSpace::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(typename GridFactory::TrimParameterType{100});

  gridFactory.insertJson("auxiliaryfiles/quarter_plate.ibra", true, {3, 3});

  auto grid = gridFactory.createUnstructedGrid<Dune::UGGrid<2>>();

  Dune::Vtk::VtkWriter writer(grid->leafGridView());
  writer.write("out/unstructuredGrid");

  return t;
}

int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  createOutputFolder("out");
  DefaultParameterSpace::Preferences::getInstance().targetAccuracy(1e-3);

  t.subTest(testUnstructuredGridCreation());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
