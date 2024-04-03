// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#define CHECK_RESERVEDVECTOR
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "dune/iga/gridcapabilities.hh"
#include "dune/iga/io/ibra/ibrareader.hh"
#include "dune/iga/io/igadatacollector.hh"
#include "dune/iga/nurbsbasis.hh"
#include "dune/iga/nurbsgrid.hh"
#include "dune/iga/nurbspatchgeometry.hh"
#include "dune/iga/trim/nurbstrimmer.hh"
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkindexset.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::IGA;

// DEFINITIONS
std::string OUTPUT_FOLDER = "output";

auto testPatchGeometryCurve() {
  TestSuite t;

  const auto dim      = 1;
  const auto dimworld = 2;

  const std::array<int, dim> order                     = {3};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints = {
      { .p = {-4, -4},   .w = 1},
      {.p = {-3, 2.8}, .w = 2.5},
      {  .p = {2, -4},   .w = 1},
      {   .p = {4, 4},   .w = 1}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  NURBSPatchGeometry<dim, dimworld> geometry(std::make_shared<Dune::IGA::NURBSPatchData<dim, dimworld>>(patchData));

  auto p0 = geometry.global({0.0});
  t.check(Dune::FloatCmp::eq(p0, {-4, -4}));

  auto p1 = geometry.global({0.5});
  t.check(Dune::FloatCmp::eq(p1, {-1.32, 0.72}));

  auto p2 = geometry.global({1});
  t.check(Dune::FloatCmp::eq(p2, {4, 4}));

  auto p3 = geometry.global({0.25});
  t.check(Dune::FloatCmp::eq(p3, {-2.7607655502392343, 0.4688995215311005}));

  // Test Operator ()
  auto p4 = geometry({0.4});
  t.check(Dune::FloatCmp::eq(p4, {-1.9854368932038828, 0.7669902912621357}));

  // Check derivative
  auto jc0 = geometry.jacobianTransposed({0});
  t.check(Dune::FloatCmp::eq(jc0[0], {7.5, 51}));

  auto jc1 = geometry.jacobianTransposed({0.5});
  t.check(Dune::FloatCmp::eq(jc1[0], {7.4496, -0.9216}));

  // Check local function
  auto u0 = geometry.local({-4, -4});
  t.check(Dune::FloatCmp::eq(u0, {0}));

  auto u1 = geometry.local({-1.32, 0.72});
  t.check(Dune::FloatCmp::eq(u1, {0.5}));

  // geomdl reports 13.230641820866644 for the length of the curve. The volume function approaches this value,
  // if you use a lot of gauß-points
  auto len = geometry.volume();

  // Check corners
  t.check(geometry.corners() == 2);
  std::array<FieldVector<double, 2>, 2> expectedCorners{
      {FieldVector<double, 2>{-4, -4}, FieldVector<double, 2>{4, 4}}
  };
  for (int i = 0; i < 2; ++i)
    t.check(Dune::FloatCmp::eq(geometry.corner(i), expectedCorners[i]));

  return t;
}

auto testPatchGeometrySurface() {
  TestSuite t("", TestSuite::ThrowPolicy::AlwaysThrow);

  const auto dim                   = 2;
  const auto dimworld              = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0, 1}, .w = 1}, {.p = {1, 0, 1}, .w = 1}, {.p = {2, 0, 2}, .w = 1}},
      {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 1}, {.p = {2, 1, 0}, .w = 1}},
      {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 1}, {.p = {2, 2, 2}, .w = 1}}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  NURBSPatchGeometry<dim, dimworld> geometry(std::make_shared<Dune::IGA::NURBSPatchData<dim, dimworld>>(patchData));

  auto p1 = geometry.global(FieldVector<double, 2>{0.5, 0.5});
  t.check(Dune::FloatCmp::eq(p1, {1.0, 1.0, 0.75})) << "p1 check failed " << p1;

  auto p2 = geometry.global(FieldVector<double, 2>{0, 0});
  t.check(Dune::FloatCmp::eq(p2, {0, 0, 1})) << "p2 check failed " << p2;

  auto p3 = geometry.global(FieldVector<double, 2>{0, 1});
  t.check(Dune::FloatCmp::eq(p3, {2, 0, 2})) << "p3 check failed " << p3;

  // Check derivative
  auto jc1 = geometry.jacobianTransposed({0.5, 0.5});
  t.check(Dune::FloatCmp::eq(jc1[0], {0.0, 2.0, 0.5})) << "jc1[0] check failed " << jc1[0];
  t.check(Dune::FloatCmp::eq(jc1[1], {2.0, 0.0, 0.5})) << "jc1[1] check failed " << jc1[1];

  // Check local function
  auto u1 = geometry.local({1.0, 1.0, 0.75});
  t.check(Dune::FloatCmp::eq(u1, {0.5, 0.5})) << "u1 check failed " << u1;

  auto u2 = geometry.local({2, 0, 2});
  t.check(Dune::FloatCmp::eq(u2, {0, 1})) << "u2 check failed " << u2;

  // Check corners
  t.check(geometry.corners() == 4);

  std::array<FieldVector<double, 3>, 4> expectedCorners{
      {FieldVector<double, 3>{0, 0, 1}, FieldVector<double, 3>{0, 2, 1}, FieldVector<double, 3>{2, 0, 2},
       FieldVector<double, 3>{2, 2, 2}}
  };
  for (int i = 0; i < 4; ++i)
    t.check(Dune::FloatCmp::eq(geometry.corner(i), expectedCorners[i]))
        << "u2 check failed " << geometry.corner(i) << " Expected: " << expectedCorners[i];

  // Check domain
  t.check(Dune::FloatCmp::eq(geometry.domain()[0].left(), 0.0))
      << "Domaincheck failed check failed " << geometry.domain()[0].left();
  t.check(Dune::FloatCmp::eq(geometry.domain()[0].right(), 1.0))
      << "Domaincheck failed check failed " << geometry.domain()[0].right();

  // Check domain midpoint
  t.check(Dune::FloatCmp::eq(geometry.domainMidPoint()[0], 0.5))
      << "Domaincheck failed check failed " << geometry.domainMidPoint()[0];

  return t;
}

auto testIbraReader() {
  TestSuite t;

  auto grid = IbraReader<2, 2>::read("auxiliaryfiles/element.ibra");

  // Check n_ele = 1, n_vert = 4
  t.check(grid->size(0) == 1);
  t.check(grid->size(2) == 4);

  grid->globalRefine(1);

  // Check n_ele = 4, n_vert = 9 after refinement
  t.check(grid->size(0) == 4);
  t.check(grid->size(2) == 9);

  // Test degree (maybe test degree elevate)
  t.check(grid->leafGridView().impl().getPatchData().degree[0] == 1);
  t.check(grid->leafGridView().impl().getPatchData().degree[1] == 1);

  // Enumerate elements and check position of centers
  auto gV = grid->leafGridView();
  std::vector<FieldVector<double, 2>> expectedElementCenters{
      {0.25, 0.25},
      {0.75, 0.25},
      {0.25, 0.75},
      {0.75, 0.75}
  };
  const auto& indexSet = grid->leafGridView().indexSet();
  for (auto& ele : elements(gV))
    t.check(ele.geometry().center() == expectedElementCenters[indexSet.index(ele)]);

  // Test Instantiation
  auto grid3D = IbraReader<2, 3>::read("auxiliaryfiles/shell.ibra", false);
  grid3D->globalRefine(2);

  // Test File not available (error should be caught)

  t.checkThrow<Dune::IOError>([]() { IbraReader<2, 2>::read("fileNotAvailable.ibra"); });

  return t;
}

auto testDataCollectorAndVtkWriter() {
  TestSuite t;

  auto lambdaf = [](auto x) {
    return Dune::FieldVector<double, 2>({std::sin(x[0]), std::cos(3 * x[0]) + std::sin(4 * x[1])});
  };
  std::vector<std::string> geometries{"element_trim_xb", "surface-hole"};

  for (auto& fileName : geometries) {
    auto grid = IbraReader<2, 2>::read("auxiliaryfiles/" + fileName + ".ibra");

    for (auto r : std::views::iota(0, 4)) {
      if (r > 0)
        grid->globalRefine(1);

      const auto gv = grid->leafGridView();
      auto lambaGV  = Dune::Functions::makeAnalyticGridViewFunction(lambdaf, gv);

      for (auto s : std::views::iota(0, 4)) {
        Dune::Vtk::DiscontinuousIgaDataCollector dataCollector1(gv, s);
        Dune::Vtk::UnstructuredGridWriter writer2(dataCollector1, Vtk::FormatTypes::ASCII);

        writer2.addPointData(lambaGV, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
        auto vtkFileName = OUTPUT_FOLDER + "/" + fileName + "_r" + std::to_string(r) + "_s" + std::to_string(s);
        writer2.write(vtkFileName);

        // Check if the file exists
        vtkFileName += ".vtu";
        t.check(std::filesystem::exists(vtkFileName) and std::filesystem::file_size(vtkFileName) > 0);
      }
    }
  }

  return t;
}

auto testExampleSuite() {
  TestSuite t("", TestSuite::AlwaysThrow);

  Dune::GeometryChecker<NURBSGrid<2, 2>> geometryChecker2;
  Dune::GeometryChecker<NURBSGrid<2, 3>> geometryChecker3;

  auto testLoop = [&]<int worldDim = 2>(const auto& _grid, int testRuns, std::string&& name) {
    for (auto i : std::views::iota(0, testRuns)) {
      if (i != 0)
        _grid->globalRefine(1);

      t.check(_grid->reportTrimError(), "Trim Error: " + name);

      if constexpr (worldDim == 2)
        geometryChecker2.checkGeometry(_grid->leafGridView());
      else
        geometryChecker3.checkGeometry(_grid->leafGridView());

      Dune::checkIndexSet(*_grid, _grid->leafGridView(), std::cout);
      gridcheck(*_grid);
    }
  };

  auto grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra");
  testLoop(grid, 4, "surface-hole");

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra");
  testLoop(grid, 4, "plate_quarter");

  grid = IbraReader<2, 2>::read("auxiliaryfiles/pipe_trim.ibra");
  testLoop(grid, 4, "pipe_trim");

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra");
  testLoop(grid, 6, "Element_trim_Xb");

  grid = IbraReader<2, 2>::read("auxiliaryfiles/nurbs_1.ibra");
  testLoop(grid, 4, "nurbs_1");

  grid          = IbraReader<2, 2>::read("auxiliaryfiles/surface-multihole.ibra");
  const auto gv = grid->leafGridView();
  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector1(gv);

  Dune::Vtk::UnstructuredGridWriter writer2(dataCollector1, Vtk::FormatTypes::ASCII);
  writer2.write("TestFileTest");
  testLoop(grid, 2, "surface-multihole");

  auto grid3 = IbraReader<2, 3>::read("auxiliaryfiles/shell-hole.ibra");
  testLoop.template operator()<3>(grid3, 1, "shell-hole");

  return t;
}

auto testMapsInTrimmedPatch() {
  TestSuite t;

  auto getTrimFlagCounts = [](const auto& patch) -> std::tuple<int, int, int> {
    int full    = 0;
    int empty   = patch.emptyElements();
    int trimmed = 0;
    for (auto i : std::views::iota(0, patch.size(0))) {
      ElementTrimFlag flag = patch.getTrimFlag(i);
      if (flag == ElementTrimFlag::full)
        ++full;
      else
        ++trimmed;
    }
    return {full, trimmed, empty};
  };

  // O refinement, 1 trimmed
  auto grid   = IbraReader<2, 2>::read("auxiliaryfiles/element_trim.ibra");
  auto& patch = grid->getPatch();

  t.check(patch.getDirectIndex<0>(0) == 0);
  t.check(patch.getRealIndex<0>(0) == 0);

  // 1 refinement: 3 trimmed, 0 empty, 1 full
  grid->globalRefine(1);

  auto& patch_1_1             = grid->getPatch();
  auto [full, trimmed, empty] = getTrimFlagCounts(patch_1_1);
  t.check(full == 1);
  t.check(empty == 0);
  t.check(trimmed == 3);

  // As n_f + n_t = n, there has to be a 1 to 1 mapping of the indices
  for (int i = 0; i < 4; ++i) {
    t.check(patch_1_1.getDirectIndex<0>(i) == i);
    t.check(patch_1_1.getRealIndex<0>(i) == i);
  }

  // Load next example Grid
  auto grid2 = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra");
  grid2->globalRefine(1);
  auto& patch_2_1 = grid2->getPatch();

  // After one refinement we should have 3 trimmed one empty element
  auto [full2, trimmed2, empty2] = getTrimFlagCounts(patch_2_1);

  t.check(full2 == 0);
  t.check(empty2 == 1);
  t.check(trimmed2 == 3);

  // The second element is the first element with a direct Index and so forth
  for (int i = 1; i < 4; ++i) {
    t.check(patch_2_1.getRealIndex<0>(i) == i - 1);
    t.check(patch_2_1.getDirectIndex<0>(i - 1) == i);
  }

  return t;
}

double calculateArea(const auto& gridView, std::optional<int> order = std::nullopt,
                     QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
  double area = 0;

  Dune::QuadratureRule<double, 2> rule;

  for (auto& ele : elements(gridView)) {
    ele.impl().fillQuadratureRule(rule, order, qt);
    auto geo = ele.geometry();
    for (auto& ip : rule)
      area += geo.integrationElement(ip.position()) * ip.weight();
  }

  return area;
}

/// @brief Test if sum of the weights is the ratio of trimmed to untrimmed surface in parameterspace * A ref (= 1)
void checkSumWeights(const auto& gridView, auto& t) {
  Dune::QuadratureRule<double, 2> rule;

  for (auto& ele : elements(gridView)) {
    if (not ele.impl().isTrimmed())
      continue;

    auto untrimmedArea = 1;
    auto trimmedArea   = ele.impl().elementSubGrid()->calculateArea();
    auto ratio         = trimmedArea / untrimmedArea;

    ele.impl().fillQuadratureRule(rule);
    double weightSum = 0.0;
    for (auto& ip : rule)
      weightSum += ip.weight();

    t.check(Dune::FloatCmp::eq(weightSum, ratio));
  }
}

auto testIntegrationPoints() {
  TestSuite t;

  /// 1. test case A = 10 * 10, r = 3
  auto targetArea = 10 * 10 - (std::numbers::pi * std::pow(3, 2));

  auto grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra");
  grid->globalRefine(1);

  double area = calculateArea(grid->leafGridView(), std::nullopt);
  std::cout << "Area (0): " << area << std::endl;

  area = calculateArea(grid->leafGridView(), std::nullopt, QuadratureType::GaussLobatto);
  std::cout << "Area GaussLobatto (0): " << area << std::endl;

  t.check(Dune::FloatCmp::eq(area, targetArea, 1e-1));

  checkSumWeights(grid->leafGridView(), t);

  grid->globalRefine(1);
  area                = calculateArea(grid->leafGridView());
  auto areaFromVolume = calculateArea(grid->leafGridView());
  t.check(Dune::FloatCmp::eq(areaFromVolume, targetArea, 1e-2))
      << "areaFromVolume " << areaFromVolume << " targetArea " << targetArea;
  std::cout << "Area (1): " << area << std::endl;
  t.check(Dune::FloatCmp::eq(area, targetArea, 1e-2));
  checkSumWeights(grid->leafGridView(), t);

  /// 1. Test case A = 20 * 2, r = 0.3
  targetArea = 20 * 2 - (std::numbers::pi * std::pow(0.3, 2));

  std::cout << "Grid 2 \n";
  auto grid2 = IbraReader<2, 2>::read("auxiliaryfiles/infty_pwh.ibra");
  grid2->globalRefine(1);

  area = calculateArea(grid2->leafGridView());
  std::cout << "Area: (0) " << area << std::endl;

  t.check(Dune::FloatCmp::eq(area, targetArea, 1e-1));
  checkSumWeights(grid->leafGridView(), t);

  grid2->globalMultiRefine(1, 0, 0);
  area = calculateArea(grid2->leafGridView());
  std::cout << "Area: (1) " << area << std::endl;

  t.check(Dune::FloatCmp::eq(area, targetArea, 1e-2));
  checkSumWeights(grid->leafGridView(), t);

  return t;
}

auto testNurbsBasis() {
  TestSuite t;

  auto runBasisChecks = [&](const auto& _grid, int testRuns) {
    for (auto i : std::views::iota(0, testRuns)) {
      if (i != 0)
        _grid->globalRefine(1);

      auto gridView = _grid->leafGridView();
      Dune::Functions::NurbsBasis<decltype(gridView)> basis(gridView);
      t.subTest(checkBasis(basis, EnableContinuityCheck(), EnableContinuityCheck()));
    }
  };

  auto grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra");
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra", true, {1, 1});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra", true, {2, 2});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim_xb.ibra", true, {3, 3});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim.ibra");
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim.ibra", true, {1, 1});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/element_trim.ibra", true, {2, 2});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra", true);
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra", true, {1, 1});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra", true, {2, 2});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/surface-hole.ibra", true, {3, 3});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra", true);
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra", true, {1, 1});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra", true, {2, 2});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra", true, {3, 3});
  runBasisChecks(grid, 4);

  grid = IbraReader<2, 2>::read("auxiliaryfiles/plate_quarter.ibra", true, {4, 4});
  runBasisChecks(grid, 4);

  return t;
}

void createOutputFolder() {
  namespace fs = std::filesystem;
  if (fs::exists(OUTPUT_FOLDER) and fs::is_directory(OUTPUT_FOLDER))
    return;
  else {
    if (not fs::create_directory(OUTPUT_FOLDER))
      DUNE_THROW(Dune::IOError, "Couldn't create output folder!");
  }
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t("", TestSuite::ThrowPolicy::AlwaysThrow);
  createOutputFolder();

  /// Test General stuff
  t.subTest(testPatchGeometryCurve());
  t.subTest(testPatchGeometrySurface());
  t.subTest(testIbraReader());
  t.subTest(testDataCollectorAndVtkWriter());

  /// 1. Test Trimming Functionality
  t.subTest(testExampleSuite());
  t.subTest(testMapsInTrimmedPatch());

  /// 2. Test Integration Points
  t.subTest(testIntegrationPoints());

  /// 3. Test Basis
  t.subTest(testNurbsBasis());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
