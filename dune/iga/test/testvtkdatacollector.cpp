// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/io/vtk/igadatacollector.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::IGA;

auto testIO() {
  TestSuite t;

  using PatchGrid   = IGA::PatchGrid<2, 2, IGA::DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto igaGridFactory = GridFactory();
  igaGridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {2, 2});
  auto igaGrid = igaGridFactory.createGrid();

  auto gridView = igaGrid->leafGridView();

  Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, 0);
  Vtk::UnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.write("element_trim");

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testIO());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}