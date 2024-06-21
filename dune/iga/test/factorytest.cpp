// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/geometrykernel/makecirculararc.hh>
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>

using namespace Dune;
using namespace Dune::IGA;

auto testFactoryWithTorus() {
  TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  auto gridFactory = makePatchGridFactory<1, 3>();

  const double r = 1.0;
  auto circle    = makeCircularArc(r);
  gridFactory.insertPatch(circle);
  auto grid = gridFactory.createGrid();

  return t;
}

auto diagonalTrimmingLoop() {
  const std::array<std::vector<double>, 1> knotSpan = {
      {
       {0, 0, 1, 1},
       }
  };
  const std::array degree = {1};

  using ControlPoint    = NURBSPatchData<1, 2>::ControlPointType;
  using ControlPointNet = NURBSPatchData<1, 2>::ControlPointNetType;

  const std::vector<ControlPoint> cpsBottom = {
      {{.p = {0, 0}, .w = 1}, {.p = {1, 0}, .w = 1}}
  };
  const std::vector<ControlPoint> cpsRight = {
      {{.p = {1, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}
  };
  const std::vector<ControlPoint> cpsDiagonal = {
      {{.p = {1, 1}, .w = 1}, {.p = {0, 0}, .w = 1}}
  };

  std::vector<NURBSPatchData<1, 2>> loop{};
  loop.emplace_back(knotSpan, ControlPointNet{cpsBottom}, degree);
  loop.emplace_back(knotSpan, ControlPointNet{cpsRight}, degree);
  loop.emplace_back(knotSpan, ControlPointNet{cpsDiagonal}, degree);
  return loop;
}

template <bool manually, int dimworld>
requires(dimworld == 2 or dimworld == 3)
auto testTriangleTrim() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim = 2;
  using Grid            = Dune::IGA::PatchGrid<gridDim, dimworld, DefaultParameterSpace::PatchGridFamily>;

  const double Lx = 1;
  const double Ly = 1;

  auto grid = [&]() {
    if constexpr (manually) {
      const std::array order                                   = {2, 2};
      const std::array<std::vector<double>, gridDim> knotSpans = {
          {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
      };

      using ControlPoint = typename NURBSPatchData<gridDim, dimworld>::ControlPointType;

      auto controlPoints = [&]() {
        if constexpr (dimworld == 2)
          return std::vector<std::vector<ControlPoint>>{
              {     {.p = {0, 0}, .w = 1},      {.p = {Lx / 2, 0}, .w = 1},      {.p = {Lx, 0}, .w = 1}},
              {{.p = {0, Ly / 2}, .w = 1}, {.p = {Lx / 2, Ly / 2}, .w = 1}, {.p = {Lx, Ly / 2}, .w = 1}},
              {    {.p = {0, Ly}, .w = 1},     {.p = {Lx / 2, Ly}, .w = 1},     {.p = {Lx, Ly}, .w = 1}}
          };
        else
          return std::vector<std::vector<ControlPoint>>{
              {     {.p = {0, 0, 0}, .w = 1},      {.p = {Lx / 2, 0, 0}, .w = 1},      {.p = {Lx, 0, 0}, .w = 1}},
              {{.p = {0, Ly / 2, 0}, .w = 1}, {.p = {Lx / 2, Ly / 2, 0}, .w = 1}, {.p = {Lx, Ly / 2, 0}, .w = 1}},
              {    {.p = {0, Ly, 0}, .w = 1},     {.p = {Lx / 2, Ly, 0}, .w = 1},     {.p = {Lx, Ly, 0}, .w = 1}}
          };
      }();

      std::array dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

      auto controlNet = typename NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

      NURBSPatchData<gridDim, dimworld> patchData;
      patchData.knotSpans     = knotSpans;
      patchData.degree        = order;
      patchData.controlPoints = controlNet;

      GridFactory<Grid> gridFactory;
      const double r = 1.0;
      gridFactory.insertPatch(patchData);
      gridFactory.insertTrimLoop(diagonalTrimmingLoop());
      return gridFactory.createGrid();
    } else {
      auto gridFactory = makePatchGridFactory<gridDim, dimworld>(withTrimmingCapabilities());
      gridFactory.insertJson("auxiliaryfiles/triangle.ibra", true, {0, 0});
      return gridFactory.createGrid();
    }
  }();

  auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });

  t.check(grid->size(0, 1) == 3) << "There should only be " << 3 << " edges in this grid, but there are "
                                 << grid->size(0, 1);

  double area = 0;
  for (auto elegeo : elements(grid->leafGridView()) | extractGeo)
    area += elegeo.volume();

  double circumference = 0;
  for (auto edgegeo : edges(grid->leafGridView()) | extractGeo)
    circumference += edgegeo.volume();

  // we created a triangle therefore the area should be a triangle
  double expectedArea = 0.5 * Lx * Ly;
  t.check(FloatCmp::eq(area, expectedArea)) << "The volume is " << area << " but should be " << 0.5 * Lx * Ly;

  double expectedCircumference = Lx + Ly + std::hypot(Lx, Ly);
  t.check(FloatCmp::eq(circumference, expectedCircumference))
      << "The circumference is " << circumference << " but should be " << expectedCircumference;
  return t;
}

auto testFactoryFactories() {
  constexpr int dim      = 2;
  constexpr int worldDim = 2;

  auto identityFactory = makePatchGridFactory<dim, worldDim>();
  auto defaultFactory  = makePatchGridFactory<dim, worldDim>(withTrimmingCapabilities());

  static_assert(std::is_same_v<decltype(identityFactory)::Grid,
                               PatchGrid<dim, worldDim, IdentityParameterSpace::PatchGridFamily>>);
  static_assert(
      std::is_same_v<decltype(defaultFactory)::Grid, PatchGrid<dim, worldDim, DefaultParameterSpace::PatchGridFamily>>);
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t("", TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testFactoryWithTorus());

  t.subTest(testTriangleTrim<true, 2>());
  t.subTest(testTriangleTrim<false, 2>());

  t.subTest(testTriangleTrim<true, 3>());
  t.subTest(testTriangleTrim<false, 3>());

  testFactoryFactories();

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
