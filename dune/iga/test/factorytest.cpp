// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
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
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/iga/trimmer/identitytrimmer/trimmer.hh>

/********
 *TODO This test is currently disabled as the interface in gridfactory cannot handle curve insertion atm (HJ)
 *********/

using namespace Dune;
using namespace Dune::IGANEW;

auto testFactoryWithTorus() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  using Grid = Dune::IGANEW::PatchGrid<1, 3, IdentityTrim::PatchGridFamily>;
  Dune::GridFactory<Grid> gridFactory;
  const double r = 1.0;
  auto circle    = makeCircularArc(r);
  gridFactory.insertPatch(circle);
  auto grid = gridFactory.createGrid();

  return t;
}

auto diagonalTrimmingCurve() {
  const std::array<std::vector<double>, 1> knotSpansCurve = {
      {
       {0, 0, 1, 1},
       }
  };
  using ControlPoint = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointType;

  const std::vector<ControlPoint> controlPointsCurve = {
      {{.p = {0, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}
  };
  const std::array orderCurve = {1};
  auto controlNetCurve        = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<1, 2> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;
  return patchDataCurve;
}

auto testFactoryWithPlateWithTriangularTrim2D() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim   = 2;
  constexpr auto dimworld = 2;
  using Grid              = Dune::IGANEW::PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  const std::array order  = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double Lx                                            = 2;
  const double Ly                                            = 3;
  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {     {.p = {0, 0}, .w = 1},      {.p = {Lx / 2, 0}, .w = 1},      {.p = {Lx, 0}, .w = 1}},
      {{.p = {0, Ly / 2}, .w = 1}, {.p = {Lx / 2, Ly / 2}, .w = 1}, {.p = {Lx, Ly / 2}, .w = 1}},
      {    {.p = {0, Ly}, .w = 1},     {.p = {Lx / 2, Ly}, .w = 1},     {.p = {Lx, Ly}, .w = 1}}
  };

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  gridFactory.insertPatch(patchData);
  gridFactory.insertTrimmingCurve(diagonalTrimmingCurve());
  gridFactory.insertTrimParameters({.dummy = 10, .trimPrecision = 1e-6});
  auto grid       = gridFactory.createGrid();
  auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });

  t.check(grid->size(0, 1) == 3) << "There should only be " << 3 << " edges in this grid, but there are "
                                 << grid->size(0, 1);

  double volume = 0;
  for (auto elegeo : elements(grid->leafGridView()) | extractGeo)
    volume += elegeo.volume();

  double circumference = 0;
  for (auto edgegeo : edges(grid->leafGridView()) | extractGeo)
    circumference += edgegeo.volume();

  // we created a triangle therfore the area should be a triangle
  t.check(FloatCmp::eq(volume, 0.5 * Lx * Ly), "Triangle volume in 2D")
      << "The volume is " << volume << " but should be " << 0.5 * Lx * Ly;
  double expectedCircumference = Lx + Ly + std::hypot(Lx, Ly);
  t.check(FloatCmp::eq(circumference, expectedCircumference), "Triangle circumference in 2D")
      << "The circumference is " << circumference << " but should be " << expectedCircumference;
  return t;
}

auto testFactoryWithPlateWithTriangularTrim3D() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim   = 2;
  constexpr auto dimworld = 3;
  using Grid              = Dune::IGANEW::PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  const std::array order  = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double Lx                                            = 2;
  const double Ly                                            = 3;
  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {     {.p = {0, 0, 0}, .w = 1},      {.p = {Lx / 2, 0, 0}, .w = 1},      {.p = {Lx, 0, 0}, .w = 1}},
      {{.p = {0, Ly / 2, 0}, .w = 1}, {.p = {Lx / 2, Ly / 2, 0}, .w = 1}, {.p = {Lx, Ly / 2, 0}, .w = 1}},
      {    {.p = {0, Ly, 0}, .w = 1},     {.p = {Lx / 2, Ly, 0}, .w = 1},     {.p = {Lx, Ly, 0}, .w = 1}}
  };

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  gridFactory.insertPatch(patchData);
  gridFactory.insertTrimmingCurve(diagonalTrimmingCurve());
  gridFactory.insertTrimParameters({.dummy = 10, .trimPrecision = 1e-6});
  auto grid       = gridFactory.createGrid();
  auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });

  double volume = 0;
  for (auto elegeo : elements(grid->leafGridView()) | extractGeo)
    volume += elegeo.volume();

  double circumference = 0;
  for (auto edgegeo : edges(grid->leafGridView()) | extractGeo)
    circumference += edgegeo.volume();

  t.check(grid->size(0, 1) == 3) << "There should only be " << 3 << " edges in this grid, but there are "
                                 << grid->size(0, 1);

  // we created a triangle therfore the area should be a triangle
  t.check(FloatCmp::eq(volume, 0.5 * Lx * Ly), "Triangle volume in 3D")
      << "The volume is " << volume << " but should be " << 0.5 * Lx * Ly;
  double expectedCircumference = Lx + Ly + std::hypot(Lx, Ly);

  t.check(FloatCmp::eq(circumference, expectedCircumference), "Triangle circumference in 3D")
      << "The circumference is " << circumference << " but should be " << expectedCircumference;

  return t;
}

auto testFactoryWithPlateWithCircularTrim3D() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim   = 2;
  constexpr auto dimworld = 3;
  using Grid              = Dune::IGANEW::PatchGrid<gridDim, dimworld, IdentityTrim::PatchGridFamily>;
  const std::array order  = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double Lx                                            = 2;
  const double Ly                                            = 3;
  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {     {.p = {0, 0, 0}, .w = 1},      {.p = {Lx / 2, 0, 0}, .w = 1},      {.p = {Lx, 0, 0}, .w = 1}},
      {{.p = {0, Ly / 2, 0}, .w = 1}, {.p = {Lx / 2, Ly / 2, 0}, .w = 1}, {.p = {Lx, Ly / 2, 0}, .w = 1}},
      {    {.p = {0, Ly, 0}, .w = 1},     {.p = {Lx / 2, Ly, 0}, .w = 1},     {.p = {Lx, Ly, 0}, .w = 1}}
  };

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  gridFactory.insertPatch(patchData);
  // gridFactory.insertTrimmingCurve(makeCircularArc2D(1.0));

  auto grid       = gridFactory.createGrid();
  auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });

  double volume = 0;
  for (auto elegeo : elements(grid->leafGridView()) | extractGeo)
    volume += elegeo.volume();

  double circumference = 0;
  for (auto edgegeo : edges(grid->leafGridView()) | extractGeo)
    circumference += edgegeo.volume();

  const auto pi = std::numbers::pi_v<double>;

  // we created a rectangle with an ellipse therfore the area should be a rectangle-ellipse
  const double expectedVolume = Lx * Ly - pi * Lx / 2 * Ly / 2;
  t.check(FloatCmp::eq(volume, expectedVolume), "rectangle-ellipse volume in 3D")
      << "The volume is " << volume << " but should be " << expectedVolume;

  double expectedCircumference = 2 * Lx + 2 * Ly;

  // since we ignore inner trims within one element the circumference should be that of the rectangle
  t.check(FloatCmp::eq(circumference, expectedCircumference), "rectangle circumference in 3D")
      << "The circumference is " << circumference << " but should be " << expectedCircumference;

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  t.subTest(testFactoryWithTorus());
  t.subTest(testFactoryWithPlateWithTriangularTrim2D());
  t.subTest(testFactoryWithPlateWithTriangularTrim3D());
  t.subTest(testFactoryWithPlateWithCircularTrim3D());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
