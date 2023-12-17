// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <clipper2/clipper.h>
#include <unordered_set>

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

#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/idset.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGANEW;

auto diagonalTrimmingCurve(double offset) {
  const std::array<std::vector<double>, 1> knotSpansCurve = {{
      {0, 0, 1, 1},
  }};
  using ControlPoint                                      = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointType;

  const std::vector<ControlPoint> controlPointsCurve
      = {{{.p = {1 - offset, 1 + offset}, .w = 1}, {.p = {-offset, offset}, .w = 1}}};
  const std::array orderCurve = {1};
  auto controlNetCurve        = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<1, 2> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;
  return Dune::IGANEW::GeometryKernel::NURBSPatch(patchDataCurve);
}

auto testFactoryWithPlateWithTriangularTrim2D() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim   = 2;
  constexpr auto dimworld = 2;
  using Grid              = Dune::IGANEW::PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  const std::array order  = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double Lx = 2;
  const double Ly = 3;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {Lx / 2, 0}, .w = 1}, {.p = {Lx, 0}, .w = 1}},
         {{.p = {0, Ly / 2}, .w = 1}, {.p = {Lx / 2, Ly / 2}, .w = 1}, {.p = {Lx, Ly / 2}, .w = 1}},
         {{.p = {0, Ly}, .w = 1}, {.p = {Lx / 2, Ly}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  const auto trimCurve = diagonalTrimmingCurve(0.1);
  using PatchTrimData  = Grid::Trimmer::PatchTrimData;
  PatchTrimData patchTrimData;
  patchTrimData.insertTrimCurve(trimCurve);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  patchData               = Splines::knotRefinement(patchData, {0.5}, 0);
  patchData               = Splines::knotRefinement(patchData, {0.5}, 1);
  Grid grid(patchData, patchTrimData);
  // grid.globalRefine(1);
  auto& parameterGrid = grid.parameterSpaceGrid();
  auto gv             = parameterGrid.leafGridView();

  return t;
}

auto testWithIbraReader() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);



  return t;
}


#include <cfenv>

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  t.subTest(testFactoryWithPlateWithTriangularTrim2D());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
