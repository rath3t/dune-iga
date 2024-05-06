// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
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
#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/io/griddrawer.hh>
#include <dune/iga/io/vtk/igadatacollector.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune::IGANEW;

auto testIbraReader() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  using PatchGrid   = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});

  const std::vector testCases{
      std::tuple<std::string, int, int>{    "auxiliaryfiles/element_trim_xb.ibra", 0, 3},
      {       "auxiliaryfiles/element_trim.ibra", 0, 3},
      {        "auxiliaryfiles/trim_2edges.ibra", 0, 3},
      {         "auxiliaryfiles/trim_multi.ibra", 0, 0},
      {       "auxiliaryfiles/surface-hole.ibra", 1, 3},
      {  "auxiliaryfiles/surface-hole-skew.ibra", 1, 3},
      {"auxiliaryfiles/surface-hole-square.ibra", 1, 3}
  };

  for (auto& [file_name, min, max] : testCases) {
    for (int i = min; i <= max; i++) {
      auto name = file_name.substr(0, file_name.find('.')).substr(file_name.find_last_of('/') + 1);

      std::cout << "Testing now " << name << " (Refinement " << i << ", " << i << ")" << std::endl;
      gridFactory.insertJson(file_name, true, {i, i});
      try {
        auto grid = gridFactory.createGrid();

        auto outputFileName = "out/" + name + +"_" + std::to_string(i) + "_" + std::to_string(i);
        drawGrid(grid.get(), outputFileName + ".gif");

        Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(grid->leafGridView());
        Dune::Vtk::UnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

        vtkWriter.write(outputFileName);

      } catch (Dune::GridError&) {
        t.check(false) << "Grid Creation failed ...\n";
      }
    }
  }

  return t;
}

auto testIbraReader3d() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  using PatchGrid   = PatchGrid<2, 3, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{200});

  const std::vector testCases{
      std::tuple<std::string, int, int>{"auxiliaryfiles/shell-hole.ibra", 0, 2}
  };

  for (auto& [file_name, min, max] : testCases) {
    auto name = file_name.substr(0, file_name.find('.')).substr(file_name.find_last_of('/') + 1);

    for (int i = min; i <= max; i++) {
      gridFactory.insertJson(file_name, true, {i, i});
      auto grid = gridFactory.createGrid();

      auto outputFileName = "out/" + name + +"_" + std::to_string(i) + "_" + std::to_string(i);
      drawGrid(grid.get(), outputFileName + ".gif");

      Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(grid->leafGridView());
      Dune::Vtk::UnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

      vtkWriter.write(outputFileName);
    }
  }

  return t;
}

int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  createOutputFolder("out");

  t.subTest(testIbraReader());
  t.subTest(testIbraReader3d());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
