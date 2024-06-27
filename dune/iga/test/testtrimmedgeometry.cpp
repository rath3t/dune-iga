// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
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
using namespace Dune;

template <typename GridView>
double referenceVolume(const GridView& gridView) {
  double volume{};
  for (const auto& ele : elements(gridView)) {
    auto qR  = ele.impl().getQuadratureRule(4);
    auto geo = ele.geometry();
    for (auto [gp, w] : qR) {
      volume += geo.integrationElement(gp) * w;
    }
  }
  return volume;
}

template <typename GridImp>
FieldVector<double, 2> referenceCenter(
    const typename GridImp::ParameterSpace::template Codim<0>::ParameterSpaceGridEntity& element) {
  auto [elements, _, __] =
      IGA::SimplexGenerator<GridImp>::createSimplicies(element, {.boundaryDivisions = 5, .targetAccuracy = 2});

  double xs        = 0.0;
  double ys        = 0.0;
  double totalArea = 0.0;

  for (const auto& triangle : elements) {
    auto center = triangle.center();
    auto area   = triangle.volume();

    totalArea += area;
    xs += area * center[0];
    ys += area * center[1];
  }
  return {xs / totalArea, ys / totalArea};
}

auto testVolumeAndCenter(const std::string& fileName, int ref) {
  TestSuite t;

  auto gridFactory = makePatchGridFactory<2, 2>(withTrimmingCapabilities());
  gridFactory.insertJson(fileName, true, {ref, ref});
  auto grid = gridFactory.createGrid();

  auto gridView = grid->leafGridView();
  double volume{};
  for (const auto& ele : elements(gridView)) {
    volume += ele.geometry().volume();
  }
  double refVolume = referenceVolume(gridView);
  t.check(FloatCmp::eq(refVolume, volume, 1e-8));

  for (const auto& ele : elements(gridView)) {
    if (not ele.impl().isTrimmed())
      continue;

    auto geometry = ele.geometry();
    auto center   = geometry.center();
    auto refCenter =
        geometry.global(referenceCenter<typename decltype(grid)::element_type>(ele.impl().getLocalEntity()));
    t.check(FloatCmp::eq(center, refCenter, 1e-8));
  }

  return t;
}

#include <cfenv>

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testVolumeAndCenter("auxiliaryfiles/element_trim.ibra", 1));
  t.subTest(testVolumeAndCenter("auxiliaryfiles/element_trim_scaled.ibra", 1));

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}