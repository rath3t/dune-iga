// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <cfenv>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimelement.hh>

using namespace Dune::IGANEW;

auto testExample1(){
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  // Setup
  using GridFactory = Dune::GridFactory<PatchGrid<2, 2, DefaultTrim::PatchGridFamily>>;

  auto gridFactory = GridFactory();
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {0, 0});
  const auto tensorCoordinates = GeometryKernel::NURBSPatch{gridFactory.patchData_}.uniqueKnotVector();
  const Dune::YaspGrid grid{tensorCoordinates};

  const auto& patchTrimData = gridFactory.patchTrimData_.value();
  auto trimmer = DefaultTrim::TrimmerImpl<2, 2, double>{};

  for (auto& ele : elements(grid.leafGridView())) {
    auto elementTrimData = trimmer.trimElement(ele, patchTrimData);
  }

  return t;
}

#include <cfenv>

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  t.subTest(testExample1());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
