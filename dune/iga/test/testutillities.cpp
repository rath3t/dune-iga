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
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/default/trimmingutils/cliputils.hh>
#include <dune/iga/parameterspace/default/trimmingutils/indextransformations.hh>
#include <dune/iga/utils/expected.hh>
#include <dune/iga/utils/igahelpers.hh>

using namespace Dune;

auto testTransformations() {
  TestSuite t;

  constexpr std::array vertexIndexMapping = {0u, 1u, 3u, 2u};
  constexpr std::array edgeIndexMapping   = {2u, 1u, 3u, 0u};

  using Transformations = Dune::IGA::DefaultParameterSpace::Transformations;

  for (const auto i : Dune::range(4u)) {
    t.check(Transformations::mapToDune(2, i) == vertexIndexMapping[i])
        << "Map to Dune for Vertex not correct, mapping says its " << Transformations::mapToDune(2, i)
        << " but should be " << vertexIndexMapping[i];
    t.check(Transformations::mapToDune(1, i) == edgeIndexMapping[i])
        << "Map to Dune for Edge not correct, mapping says its " << Transformations::mapToDune(1, i)
        << " but should be " << edgeIndexMapping[i];
  }

  constexpr std::array vertexIndexBackMapping = {0u, 1u, 3u, 2u};
  constexpr std::array edgeIndexBackMapping   = {3u, 1u, 0u, 2u};

  for (const auto i : Dune::range(4u)) {
    t.check(Transformations::mapToParameterSpace(2, i) == vertexIndexBackMapping[i])
        << "Map to ParameterSpace for Vertex not correct, mapping says its "
        << Transformations::mapToParameterSpace(2, i) << " but should be " << vertexIndexBackMapping[i];
    ;
    t.check(Transformations::mapToParameterSpace(1, i) == edgeIndexBackMapping[i])
        << "Map to ParameterSpace for Edge not correct, mapping says its " << Transformations::mapToParameterSpace(1, i)
        << " but should be " << edgeIndexBackMapping[i];
    ;
  }

  return t;
}

// todo write a proper test
auto testTransformToSpan() {
  TestSuite t("Test TransformToSpan");

  auto igaGridFactory = makePatchGridFactory<2, 2>(withTrimmingCapabilities());
  igaGridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {1, 1});
  auto igaGrid = igaGridFactory.createGrid();

  for (const auto& ele : Dune::elements(igaGrid->leafGridView()) |
                             std::views::filter([](const auto& ele_) { return ele_.impl().isTrimmed(); })) {
    auto trimEdges  = ele.impl().getLocalEntity().trimData().edges();
    auto hostEntity = ele.impl().getLocalEntity().getHostEntity();

    for (const auto& edgeInfo :
         trimEdges | std::views::filter([](const auto& edgeInfo_) { return edgeInfo_.geometry.has_value(); })) {
      auto edgeGeo            = edgeInfo.geometry.value();
      auto transformedEgdeGeo = transformToSpan(edgeInfo.geometry.value(), hostEntity.geometry());

      for (const auto& cp : transformedEgdeGeo.patchData().controlPoints.directGetAll()) {
        for (auto i : range(2))
          t.check(FloatCmp::ge(cp.p[i], 0.0) and FloatCmp::le(cp.p[i], 1.0))
              << "CP " << cp.p << " not in local reference space [0, 1]^2";
      }
    }
  }

  return t;
}

auto testForEachUntrimmedBoundary(int degreeElevate) {
  TestSuite t("Test testForEachUntrimmedBoundary");

  using PatchGrid   = IGA::PatchGrid<2, 2, IGA::DefaultParameterSpace::PatchGridFamily>;
  using GridView    = PatchGrid::LeafGridView;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto igaGridFactory = GridFactory();
  igaGridFactory.insertJson("auxiliaryfiles/quarter_plate.ibra", true, {0, 0}, {degreeElevate, degreeElevate}, {4, 4});
  auto igaGrid = igaGridFactory.createGrid();

  using namespace Functions::BasisFactory;

  auto gridView = igaGrid->leafGridView();
  Functions::NurbsBasis<GridView> basis2(gridView, nurbs());

  auto size = basis2.size();

  forEachUntrimmedBoundaryDOF(basis2, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    t.check(localView.index(localIndex) < size);
  });

  auto powerB = makeBasis(gridView, Functions::BasisFactory::power<2>(nurbs(), FlatInterleaved()));

  size = powerB.size();
  forEachUntrimmedBoundaryDOF(basis2, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    t.check(localView.index(localIndex) < size);
  });

  return t;
}

auto testClipUtils() {
  TestSuite t("testClipUtils", TestSuite::ThrowPolicy::AlwaysThrow);

  using IGA::DefaultParameterSpace::Util::isConsecutive;

  t.check(isConsecutive(0, 1));
  t.check(isConsecutive(1, 2));
  t.check(isConsecutive(2, 3));
  t.check(isConsecutive(3, 0));

  t.check(not isConsecutive(2, 0));
  t.check(not isConsecutive(0, 3));
  t.check(not isConsecutive(0, 2));
  t.check(not isConsecutive(2, 1));
  t.check(not isConsecutive(3, 1));
  t.check(not isConsecutive(1, 3));

  return t;
}

enum class compute_error
{
  negative_number,
  to_big_number
};

auto testExpected() {
  TestSuite t("", TestSuite::ThrowPolicy::ThrowOnRequired);

  auto computeNumber = [](int i) -> Std::expected<int, compute_error> {
    if (i < 0)
      return Std::unexpected(compute_error::negative_number);
    if (i < 5)
      return i * 5;

    return Std::unexpected(compute_error::to_big_number);
  };

  auto v1 = computeNumber(-1);
  auto v2 = computeNumber(3);
  auto v3 = computeNumber(10);

  t.require(not v1.has_value()) << "v1 reports to have a value, but shouldn't";
  t.check(not v1) << "v1 reports to have a value, but shouldn't";
  t.check(v1.value_or(3) == 3) << "v1 value_or should yield 3, but yields " << v1.value_or(3);
  t.check(v1.error() == compute_error::negative_number) << "v1 should report negative number, but doesn't";

  t.require(v2.has_value()) << "v2 reports to have not have a value, but should";
  // t.check(static_cast<bool>(v2));
  t.check(*v2 == 15) << "v2 should yield a value of 15, but yields " << *v2;
  t.check(v2.value_or(3) == 15) << "v2 value_or should yield 15, but yields " << v2.value_or(3);

  t.require(not v3.has_value()) << "v1 reports to have a value, but shouldn't";
  t.check(not v3) << "v3 reports to have a value, but shouldn't";
  t.check(v3.value_or(3) == 3) << "v3 value_or should yield 3, but yields " << v1.value_or(3);
  t.check(v3.error() == compute_error::to_big_number) << "v3 should report to_big_number, but doesn't";

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testTransformations());
  t.subTest(testTransformToSpan());
  t.subTest(testClipUtils());

  t.subTest(testForEachUntrimmedBoundary(0));
  t.subTest(testForEachUntrimmedBoundary(1));
  t.subTest(testForEachUntrimmedBoundary(2));

  t.subTest(testExpected());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}