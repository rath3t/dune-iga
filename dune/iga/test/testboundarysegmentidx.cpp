// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/tuplevector.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGA;
using namespace Dune;

auto boundaryIndices(size_t refLevel, size_t edgeNr) {
  int nEleEdge = 1 << refLevel;
  return std::views::iota(edgeNr * nEleEdge, (edgeNr + 1) * nEleEdge);
}

bool checkBoundarySegmentIndex(size_t refLevel, size_t edgeNr, size_t indexToCheck) {
  auto indices = boundaryIndices(refLevel, edgeNr);
  return std::ranges::find(indices, indexToCheck) != indices.end();
}

auto runTestHierachic(Dune::TestSuite& t, const std::string& fileName, int initialRefLevel, int maxRefLevel,
                      bool trimmed) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson(fileName, trimmed, {initialRefLevel, initialRefLevel});

  const auto grid = gridFactory.createGrid();
  const auto gv   = grid->leafGridView();

  auto numBoundarySegmentsUntrimmed = grid->trimmer().untrimmedParameterSpaceGrid_->numBoundarySegments();

  for (const auto i : Dune::range(initialRefLevel, maxRefLevel)) {
    if (i > initialRefLevel)
      grid->globalRefine(1);

    // Loop over all intersections, if boundary, check
    std::set<size_t> trimmedBoundarySegmentIndices;
    std::set<size_t> allBoundarySegmentIndices;
    for (const auto& ele : elements(gv)) {
      for (const auto& intersection : intersections(gv, ele)) {
        if (intersection.boundary()) {
          const auto boundarySegmentIdx = intersection.boundarySegmentIndex();
          allBoundarySegmentIndices.insert(boundarySegmentIdx);

          const auto center             = intersection.geometry().center();
          if (FloatCmp::lt(center[0], 1e-8))
            t.check(checkBoundarySegmentIndex(initialRefLevel, 0, boundarySegmentIdx));
          else if (FloatCmp::gt(center[0], 1 - 1e-8))
            t.check(checkBoundarySegmentIndex(initialRefLevel, 1, boundarySegmentIdx));
          else if (FloatCmp::lt(center[1], 1e-8))
            t.check(checkBoundarySegmentIndex(initialRefLevel, 2, boundarySegmentIdx));
          else if (FloatCmp::gt(center[1], 1 - 1e-8))
            t.check(checkBoundarySegmentIndex(initialRefLevel, 3, boundarySegmentIdx));
          else {
            t.check(boundarySegmentIdx >= numBoundarySegmentsUntrimmed);
            trimmedBoundarySegmentIndices.insert(boundarySegmentIdx);
          }
        }
      }
    }
    t.check(trimmedBoundarySegmentIndices.size() + numBoundarySegmentsUntrimmed == grid->numBoundarySegments());

    // TODO Test for consecutive indices (is this enforced by the grid interface?)
    // size_t j = 0;
    // t.check(std::ranges::all_of(allBoundarySegmentIndices, [&j](size_t k) {
    //   return k == j++;
    // }));
  }
}

auto runTestFlat(Dune::TestSuite& t, const std::string& fileName, int refLevel, bool trimmed) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson(fileName, trimmed, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  const auto gv   = grid->leafGridView();

  // Loop over all intersections, if boundary, check
  for (const auto& ele : elements(gv)) {
    for (const auto& intersection : intersections(gv, ele)) {
      if (intersection.boundary()) {
        const auto boundarySegmentIdx = intersection.boundarySegmentIndex();
        const auto center             = intersection.geometry().center();
        if (FloatCmp::lt(center[0], 1e-8))
          t.check(checkBoundarySegmentIndex(refLevel, 0, boundarySegmentIdx));
        else if (FloatCmp::gt(center[0], 1 - 1e-8))
          t.check(checkBoundarySegmentIndex(refLevel, 1, boundarySegmentIdx));
        else if (FloatCmp::lt(center[1], 1e-8))
          t.check(checkBoundarySegmentIndex(refLevel, 2, boundarySegmentIdx));
        else if (FloatCmp::gt(center[1], 1 - 1e-8))
          t.check(checkBoundarySegmentIndex(refLevel, 3, boundarySegmentIdx));
      }
    }
  }
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  std::cout << "Testing element_trim.ibra" << std::endl;
  for (const auto i : Dune::range(2)) {
    runTestFlat(t, "auxiliaryfiles/element_trim.ibra", i, false);
    runTestFlat(t, "auxiliaryfiles/element_trim.ibra", i, true);
  }

  runTestHierachic(t, "auxiliaryfiles/element_trim.ibra", 1, 3, false);
  runTestHierachic(t, "auxiliaryfiles/element_trim.ibra", 1, 3, true);

  std::cout << "Testing surface-hole_normalized.ibra" << std::endl;
  for (const auto i : Dune::range(1, 4)) {
    runTestFlat(t, "auxiliaryfiles/surface-hole_normalized.ibra", i, true);
  }
  runTestHierachic(t, "auxiliaryfiles/surface-hole_normalized.ibra", 1, 3, true);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}