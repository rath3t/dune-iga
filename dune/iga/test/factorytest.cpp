// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
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

using namespace Dune::IGANEW;

auto testFactory() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);
  using Grid = Dune::IGANEW::PatchGrid<1, 3, Dune::IGANEW::DefaultTrim::Trimmer>;
  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  auto circle    = makeCircularArc(r);
  gridFactory.insertPatch(circle);
  gridFactory.setupTrimmer({.dummy = 10, .trimPrecision = 1e-6});
  auto grid = gridFactory.createGrid();

  auto entity       = *grid->leafGridView().begin<0>();
  auto intersection = *grid->leafGridView().ibegin(entity);
  std::cout << intersection.geometryInInside().volume() << std::endl;
  return t;
}

auto testFactoryWithPlateWithCurveInsertion() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 2;
  using Grid = Dune::IGANEW::PatchGrid<gridDim, dimworld, Dune::IGANEW::DefaultTrim::Trimmer>;
  const std::array order = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
         {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 1}, {.p = {1, 0.5}, .w = 1}},
         {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  const std::array<std::vector<double>, gridDim-1> knotSpansCurve = {{{0, 0, 1, 1}, }};
  const std::vector<ControlPoint> controlPointsCurve = {{{.p = {0, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}};
  const std::array orderCurve = {1};
  auto controlNetCurve = Dune::IGANEW::NURBSPatchData<gridDim-1, gridDim>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<gridDim-1, gridDim> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;

  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  gridFactory.insertPatch(patchData);
  gridFactory.insertTrimmingCurve(patchDataCurve);
  gridFactory.setupTrimmer({.dummy = 10, .trimPrecision = 1e-6});
  auto grid = gridFactory.createGrid();

  auto entity       = *grid->leafGridView().begin<0>();
  auto intersection = *grid->leafGridView().ibegin(entity);
  std::cout << intersection.geometryInInside().volume() << std::endl;
  return t;
}

auto testFactoryWithPlateIn3DWithCurveInsertion() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 3;
  using Grid = Dune::IGANEW::PatchGrid<gridDim, dimworld, Dune::IGANEW::DefaultTrim::Trimmer>;
  const std::array order = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0,0}, .w = 1}, {.p = {0.5, 0,0}, .w = 1}, {.p = {1, 0,0}, .w = 1}},
         {{.p = {0, 0.5,0}, .w = 1}, {.p = {0.5, 0.5,0}, .w = 1}, {.p = {1, 0.5,0}, .w = 1}},
         {{.p = {0, 1,0}, .w = 1}, {.p = {0.5, 1,0}, .w = 1}, {.p = {1, 1,0}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  const std::array<std::vector<double>, gridDim-1> knotSpansCurve = {{{0, 0, 1, 1}, }};
  using ControlPointParameterSpace = Dune::IGANEW::NURBSPatchData<gridDim-1, gridDim>::ControlPointType;

  const std::vector<ControlPointParameterSpace> controlPointsCurve = {{{.p = {0, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}};
  const std::array orderCurve = {1};
  auto controlNetCurve = Dune::IGANEW::NURBSPatchData<gridDim-1, gridDim>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<gridDim-1, gridDim> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;

  GridFactory<Grid> gridFactory;
  const double r = 1.0;
  gridFactory.insertPatch(patchData);
  gridFactory.insertTrimmingCurve(patchDataCurve);
  gridFactory.setupTrimmer({.dummy = 10, .trimPrecision = 1e-6});
  auto grid = gridFactory.createGrid();

  auto entity       = *grid->leafGridView().begin<0>();
  auto intersection = *grid->leafGridView().ibegin(entity);
  std::cout << intersection.geometryInInside().volume() << std::endl;
  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);
  t.subTest(testFactory());
  t.subTest(testFactoryWithPlateWithCurveInsertion());
  t.subTest(testFactoryWithPlateIn3DWithCurveInsertion());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
