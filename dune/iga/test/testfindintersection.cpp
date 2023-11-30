// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/iga/geometrykernel/findintersection.hh>
#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>
using namespace Dune::IGANEW;


auto diagonalTrimmingCurve() {
  const std::array<std::vector<double>, 1> knotSpansCurve = {{
      {0, 0,0,1, 1, 1},
  }};
  using ControlPoint                                      = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointType;

  const std::vector<ControlPoint> controlPointsCurve = {{{.p = {0, 0}, .w = 1},{.p = {0.6, 0.4}, .w = 1} ,{.p = {1.2, -7}, .w = 1}}};
  const std::array orderCurve                        = {2};
  auto controlNetCurve = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<1, 2> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;
  return GeometryKernel::NURBSPatch(patchDataCurve);
}



auto test1() {
  Dune::TestSuite t;
  auto curve1 = diagonalTrimmingCurve();
  std::cout<<std::setprecision(16)<<std::endl;
  auto [success,tParameter,curvePoint]=Dune::IGANEW::
  findIntersectionCurveAndLine(curve1,Dune::FieldVector<double,2>({0.5,0.5}),{0.0,5.0},{0.5,0.5});
  t.check(success)<<"No intersection found";

  t.check(Dune::FloatCmp::eq(tParameter[1],-0.3041666666666666))<<"tLine is " <<tParameter[1]<< " but should be "<<-0.3041666666666666;
  t.check(Dune::FloatCmp::eq(tParameter[0],0.4166666666666667))<<"tCurve is " <<tParameter[0]<< " but should be "<<0.4166666666666667;
  t.check(Dune::FloatCmp::eq(curvePoint,{0.5,-1.020833333333333}))<<"The obtained point is wrong";
return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  t.subTest(test1());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
