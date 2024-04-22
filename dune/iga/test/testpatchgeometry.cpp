// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <cfenv>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/patchgrid.hh>

using namespace Dune;
using namespace Dune::IGANEW;

auto testPatchGeometryCurve() {
  TestSuite t;

  constexpr auto dim         = 1;
  constexpr auto dimworld    = 2;
  constexpr std::array order = {3};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 1, 1, 1}}};

  using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints = {
      { .p = {-4, -4},   .w = 1},
      {.p = {-3, 2.8}, .w = 2.5},
      {  .p = {2, -4},   .w = 1},
      {   .p = {4, 4},   .w = 1}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  GeometryKernel::NURBSPatch geometry(patchData);

  auto p0 = geometry.global({0.0});
  t.check(Dune::FloatCmp::eq(p0, {-4, -4}));

  auto p1 = geometry.global({0.5});
  t.check(Dune::FloatCmp::eq(p1, {-1.32, 0.72}));

  auto p2 = geometry.global({1});
  t.check(Dune::FloatCmp::eq(p2, {4, 4}));

  auto p3 = geometry.global({0.25});
  t.check(Dune::FloatCmp::eq(p3, {-2.7607655502392343, 0.4688995215311005}));

  // Test Operator ()
  auto p4 = geometry.global({0.4});
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
  auto len = geometry.curveLength();
  t.check(Dune::FloatCmp::eq(len, 13.230641820, 0.01));

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

  constexpr auto dim         = 2;
  constexpr auto dimworld    = 3;
  constexpr std::array order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0, 1}, .w = 1}, {.p = {1, 0, 1}, .w = 1}, {.p = {2, 0, 2}, .w = 1}},
      {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 1}, {.p = {2, 1, 0}, .w = 1}},
      {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 1}, {.p = {2, 2, 2}, .w = 1}}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  const auto controlNet = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  const GeometryKernel::NURBSPatch geometry(patchData);

  const auto p1 = geometry.global(FieldVector<double, 2>{0.5, 0.5});
  t.check(Dune::FloatCmp::eq(p1, {1.0, 1.0, 0.75})) << "p1 check failed " << p1;

  const auto p2 = geometry.global(FieldVector<double, 2>{0, 0});
  t.check(Dune::FloatCmp::eq(p2, {0, 0, 1})) << "p2 check failed " << p2;

  const auto p3 = geometry.global(FieldVector<double, 2>{0, 1});
  t.check(Dune::FloatCmp::eq(p3, {2, 0, 2})) << "p3 check failed " << p3;

  // Check derivative
  auto jc1 = geometry.jacobianTransposed({0.5, 0.5});
  t.check(Dune::FloatCmp::eq(jc1[0], {0.0, 2.0, 0.5})) << "jc1[0] check failed " << jc1[0];
  t.check(Dune::FloatCmp::eq(jc1[1], {2.0, 0.0, 0.5})) << "jc1[1] check failed " << jc1[1];

  // Check local function
  const auto u1 = geometry.local({1.0, 1.0, 0.75});
  t.check(Dune::FloatCmp::eq(u1, {0.5, 0.5})) << "u1 check failed " << u1;

  const auto u2 = geometry.local({2, 0, 2});
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

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testPatchGeometryCurve());
  t.subTest(testPatchGeometrySurface());

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}