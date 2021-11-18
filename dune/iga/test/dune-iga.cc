// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/iga/NURBSgrid.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/bsplinepatch.hh>
//#include <dune/iga/vtkfile.hh>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/bsplinebasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/nurbsbasis.hh>

using namespace Dune;
using namespace Dune::IGA;

void testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const auto dim      = 1UL;
  const auto dimworld = 3UL;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};
  //  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints
      = {{.p = {1, 3, 4}, .w = 1}, {.p = {2, 2, 2}, .w = 3}, {.p = {3, 4, 5}, .w = 1}, {.p = {5, 1, 7}, .w = 2}, {.p = {4, 7, 2}, .w = 1},
         {.p = {8, 6, 2}, .w = 1}, {.p = {2, 9, 9}, .w = 7}, {.p = {1, 4, 3}, .w = 1}, {.p = {1, 7, 1}, .w = 5}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  //  vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-CurveNewFineResample_knotRefine");
}

void testNURBSGridSurface() {
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2UL;
  const auto dimworld              = 3UL;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  //  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
  //      = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");
}

void testTorusGeometry() {
  const double R = 2.0;
  const double r = 1.0;
  auto circle    = makeCircularArc(r);
  auto nurbsPatchData = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 360.0);

  IGA::NURBSGrid grid(nurbsPatchData);
  grid.globalRefine(0);
  auto gridView = grid.leafGridView();

  const int subSampling = 2;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-SurfaceRevolutionFLAT");

  TestSuite test;
  Dune::GeometryChecker<IGA::NURBSGrid<2UL, 3UL>> geometryChecker;
  geometryChecker.checkGeometry(gridView);
   Dune::checkIndexSet(grid,gridView,std::cout);

  double vol = 0.0;
  for (auto& ele : elements(gridView)) {
    vol += ele.geometry().volume();
  }
  const double pi                   = std::numbers::pi_v<double>;
  const double referenceTorusVolume = 4.0 * pi * pi * r * R;
  test.check(vol - referenceTorusVolume < 1e-4, "The integrated area of the torus surface is wrong!");

  checkEntityLifetime(gridView,gridView.size(0));
//  auto view = std::ranges::transform_view(elements(gridView),[](const auto& ele){return ele.geometry();}; );
  for (auto&& elegeo :  elements(gridView) | std::views::transform([](const auto& ele){return ele.geometry();}))
    checkJacobians(elegeo);

  checkIterators(gridView);

//  gridcheck(grid);

}

void testNURBSSurface() {
  // parameters
  int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim                    = 2;
  const int dimworld               = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  using ControlPoint                                   = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
  //

  //  auto weightNet  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 1);
  testSuite.check(patch.size(1) == 4);
  testSuite.check(patch.size(2) == 9);
}

void testNURBSCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<ControlPoint> controlPoints
      = {{.p = {1, 3, 4}, .w = 2}, {.p = {2, 2, 2}, .w = 2}, {.p = {3, 4, 5}, .w = 1}, {.p = {5, 1, 7}, .w = 1}, {.p = {4, 7, 2}, .w = 4},
         {.p = {8, 6, 2}, .w = 2}, {.p = {2, 9, 9}, .w = 1}, {.p = {1, 4, 3}, .w = 2}, {.p = {1, 7, 1}, .w = 4}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 5);
  testSuite.check(patch.size(1) == controlPoints.size());
}

void testBSplineCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::GlobalCoordinateType;
  const std::vector<ControlPoint> controlPoints
      = {{{1, 3, 4}}, {{2, 2, 2}}, {{3, 4, 5}}, {{5, 1, 7}}, {{4, 7, 2}}, {{8, 6, 2}}, {{2, 9, 9}}, {{1, 4, 3}}, {{1, 7, 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet = MultiDimensionNet<dim, typename Dune::IGA::NURBSPatchData<dim, dimworld>::GlobalCoordinateType>(dimsize, controlPoints);

  IGA::BSplinePatch<dim, dimworld> patch(knotSpans, controlNet, order);
}

void testNurbsGridCylinder() {
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  constexpr double invsqr2         = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  grid.globalRefine(5);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, order);

  testSuite.check(nurbsPatch.size(0) == 1);
  testSuite.check(nurbsPatch.size(1) == 4);
  testSuite.check(nurbsPatch.size(2) == 6);

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");

  ////////////////////////////////////////////////////////////////
}

// template <int dim>
void testNurbsBasis() {
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  constexpr double invsqr2         = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  nurbsPatchData.controlPoints = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
                                  {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
                                  //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
                                  {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};
  nurbsPatchData.order         = order;

  IGA::NURBSGrid<dim, dimworld> grid(nurbsPatchData);
  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite test;

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");
  using GridView = decltype(gridView);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.getPatchData());
  //  Dune::Functions::NurbsPreBasis<GridView,Dune::Functions::FlatMultiIndex<std::size_t>> basis(gridView_, nurbsPatchData);

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView, gridView.getPatchData());
    ////        test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
    test.subTest(checkBasis(basis2));
  }

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView);
    //        test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
    test.subTest(checkBasis(basis2));
  }

  {
    // Check basis created via makeBasis
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, nurbs<dim>(gridView.getPatchData()));
    //          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
    test.subTest(checkBasis(basis2));
  }

  {
    // Check whether a B-Spline basis can be combined with other bases.
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, power<2>(gridView.getPreBasis()));
    //          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
    test.subTest(checkBasis(basis2));
  }
}

void testBsplineBasisFunctions() {
  std::vector<double> knots = {0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3};
  int degree                = 2;
  // std::vector<double> N;

  TestSuite test;

  auto N = Dune::IGA::Bspline<double>::basisFunctions(0.3, knots, degree);
  using Dune::FloatCmp::eq;
  test.check(eq(N[0], 0.16), "P=2,N0,u=0.3");
  test.check(eq(N[1], 0.48), "P=2,N1,u=0.3");
  test.check(eq(N[2], 0.36), "P=2,N2,u=0.3");

  N = Dune::IGA::Bspline<double>::basisFunctions(0.5, knots, degree);  // try knot span boundary

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.5");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.5");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.5");

  N = Dune::IGA::Bspline<double>::basisFunctions(0.0, knots, degree);  // try left end

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.0");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.0");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.0");

  N = Dune::IGA::Bspline<double>::basisFunctions(3.0, knots, degree);  // try right end

  test.check(eq(N[0], 0.0), "P=2,N0,u=3");
  test.check(eq(N[1], 0.0), "P=2,N1,u=3");
  test.check(eq(N[2], 1.0), "P=2,N2,u=3");

  knots  = {0, 0, 0, 0.5, 1, 1, 1};
  degree = 2;

  N = Dune::IGA::Bspline<double>::basisFunctions(0.1, knots, degree);
  test.check(eq(N[0], 0.64), "P=2,N0,u=0.1");
  test.check(eq(N[1], 0.34), "P=2,N1,u=0.1");
  test.check(eq(N[2], 0.02), "P=2,N2,u=0.1");

  N = Dune::IGA::Bspline<double>::basisFunctions(0.01, knots, degree);
  test.check(eq(N[0], 0.9604), "P=2,N0,u=0.01");
  test.check(eq(N[1], 0.0394), "P=2,N1,u=0.01");
  test.check(eq(N[2], 0.0002), "P=2,N2,u=0.01");

  knots  = {0, 0, 0, 0, 0.5, 1, 1, 2, 2, 2, 2};
  degree = 3;

  auto N2 = Dune::IGA::Bspline<double>::basisFunctions(1.45, knots, degree);
  test.check(eq(N2[0], 0.1109166666666667), "P=3,N0,u=1.45");
  test.check(eq(N2[1], 0.4638333333333333), "P=3,N1,u=1.45");
  test.check(eq(N2[2], 0.334125), "P=3,N2,u=1.45");
  test.check(eq(N2[3], 0.091125), "P=3,N3,u=1.45");

  auto dN = Dune::IGA::Bspline<double>::basisFunctionDerivatives(1.45, knots, degree, 3);
  // check values
  test.check(eq(dN[0][0], 0.1109166666666667), "P=3,dN00,u=1.45");
  test.check(eq(dN[0][1], 0.4638333333333333), "P=3,dN01,u=1.45");
  test.check(eq(dN[0][2], 0.334125), "P=3,dN02,u=1.45");
  test.check(eq(dN[0][3], 0.091125), "P=3,dN03,u=1.45");

  // check first derivatives
  test.check(eq(dN[1][0], -0.605), "P=3,dN10,u=1.45");
  test.check(eq(dN[1][1], -0.88), "P=3,dN11,u=1.45");
  test.check(eq(dN[1][2], 0.8775), "P=3,dN12,u=1.45");
  test.check(eq(dN[1][3], 0.6075), "P=3,dN13,u=1.45");

  // check second derivatives
  test.check(eq(dN[2][0], 2.2), "P=3,dN20,u=1.45");
  test.check(eq(dN[2][1], -2.8), "P=3,dN21,u=1.45");
  test.check(eq(dN[2][2], -2.1), "P=3,dN22,u=1.45");
  test.check(eq(dN[2][3], 2.7), "P=3,dN23,u=1.45");

  // check third derivatives
  test.check(eq(dN[3][0], -4.0), "P=3,dN30,u=1.45");
  test.check(eq(dN[3][1], 16.0), "P=3,dN31,u=1.45");
  test.check(eq(dN[3][2], -18.0), "P=3,dN32,u=1.45");
  test.check(eq(dN[3][3], 6.0), "P=3,dN33,u=1.45");
  // https://godbolt.org/z/Ta3fzW553
  auto Nf                           = Dune::IGA::Bspline<double>(knots, degree);
  std::vector<double> NAtEvalPoints = {1,
                                       0.1714677640603567,
                                       0.001371742112482855,
                                       0.0740740740740741,
                                       0.00274348422496571,
                                       0.4682213077274805,
                                       0.1975308641975309,
                                       0.05852766346593506,
                                       0.007315957933241894,
                                       0};
  for (int i = 0; i < NAtEvalPoints.size(); ++i) {
    test.check(eq(Nf(i / (NAtEvalPoints.size() - 1.0) * 2.0)[0], NAtEvalPoints[i]));
  }

  std::array<double, 2> xieta{0.2, 0.25};
  std::array<std::vector<double>, 2> knots2 = {{{0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3}, {0, 0, 0, 2, 2, 2}}};
  std::array<int, 2> degree2{2, 2};
  const std::vector<std::vector<double>> weights2 = {{{1, 2, 3, 4, 5, 6, 7}, {8, 9, 10, 11, 12, 13, 14}, {15, 16, 17, 18, 19, 20, 21}}};
  std::array<int, 2> dimsize                      = {static_cast<int>(weights2.size()), static_cast<int>(weights2[0].size())};
  MultiDimensionNet<2UL, double> weightNet(dimsize, weights2);

  auto N_Nurbs = Dune::IGA::Nurbs<double, 2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();

  test.check(N_Nurbs.size() == (degree2[0] + 1) * (degree2[1] + 1));

  test.check(eq(N_Nurbs[0], 0.04023722627737226), "Nurbs2d P=2,N0");  // check ansatzfunctions in domain
  test.check(eq(N_Nurbs[1], 0.4291970802919708), "Nurbs2d P=2,N1");
  test.check(eq(N_Nurbs[2], 0.2682481751824818), "Nurbs2d P=2,N2");
  test.check(eq(N_Nurbs[3], 0.02299270072992701), "Nurbs2d P=2,N3");
  test.check(eq(N_Nurbs[4], 0.137956204379562), "Nurbs2d P=2,N4");
  test.check(eq(N_Nurbs[5], 0.08175182481751825), "Nurbs2d P=2,N5");
  test.check(eq(N_Nurbs[6], 0.002463503649635036), "Nurbs2d P=2,N6");
  test.check(eq(N_Nurbs[7], 0.01094890510948905), "Nurbs2d P=2,N7");
  test.check(eq(N_Nurbs[8], 0.006204379562043796), "Nurbs2d P=2,N8");
  test.check(eq(std::accumulate(N_Nurbs.begin(), N_Nurbs.end(), 0.0), 1.0), "partition of unity in domain");

  xieta   = {0, 0.1};
  N_Nurbs = Dune::IGA::Nurbs<double, 2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();
  test.check(eq(N_Nurbs[0], 0.8204545454545455), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(N_Nurbs[1], 0.0), "Nurbs P=2,N1");
  test.check(eq(N_Nurbs[2], 0.0), "Nurbs P=2,N2");
  test.check(eq(N_Nurbs[3], 0.1727272727272728), "Nurbs P=2,N3");
  test.check(eq(N_Nurbs[4], 0.0), "Nurbs P=2,N4");
  test.check(eq(N_Nurbs[5], 0.0), "Nurbs P=2,N5");
  test.check(eq(N_Nurbs[6], 0.00681818181818182), "Nurbs P=2,N6");
  test.check(eq(N_Nurbs[7], 0.0), "Nurbs P=2,N7");
  test.check(eq(N_Nurbs[8], 0.0), "Nurbs P=2,N8");

  test.check(eq(std::accumulate(N_Nurbs.begin(), N_Nurbs.end(), 0.0), 1.0), "partition of unity on boundary");
  xieta         = {0, 0.1};
  auto dN_Nurbs = Dune::IGA::Nurbs<double, 2>::basisFunctionDerivatives(xieta, knots2, degree2, weightNet, 5);

  // for(int j= 0; j<dN_Nurbs.directSize() ;++j) {
  //   auto multiIndex = dN_Nurbs.directToMultiIndex(j);
  //   for (int i = 0; i < 2; ++i)
  //     std::cout<<multiIndex[i]<<" ";
  //   std::cout<<std::endl;
  // }

  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 0}), 0.8204545454545455), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 1}), 0.1727272727272728), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 2}), 0.00681818181818182), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // first derivative in u dir
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 0}), -24.166115702479335), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 0}), 26.25454545454545), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 1}), -5.0876033057851231), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 1}), 3.1090909090909089), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 2}), -0.20082644628099172), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 2}), 0.090909090909090925), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in u dir
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 0}), 1236.8386175807659), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 0}), -1441.6132231404956), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 0}), 98.454545454545439), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 1}), 260.38707738542439), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 1}), -170.71735537190079), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 1}), 11.054545454545453), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 2}), 10.278437265214125), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 2}), -4.9917355371900829), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 2}), 0.30909090909090914), "Nurbs P=2,N8");

  // third derivative in u dir
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 0}), -94449.494433440283), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 0}), 110086.82794891056), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 0}), -7518.3471074380141), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 1}), -19884.104091250589), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 1}), 13036.598046581514), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 1}), -844.16528925619821), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 2}), -784.89884570726042), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 2}), 381.18707738542452), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 2}), -23.603305785123968), "Nurbs P=2,N8");

  // fourth derivative in u dir
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 0}), 9616675.7968593743), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 0}), -11208840.663889075), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 0}), 765504.432757325), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 1}), 2024563.3256546052), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 1}), -1327362.7101973905), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 1}), 85951.374906085635), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 2}), 79916.973381102871), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 2}), -38811.775151970498), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 2}), 2403.2456799398947), "Nurbs P=2,N8");

  // fifth derivative in u dir
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 0}), -1223940555.9639204), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 0}), 1426579720.8586094), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 0}), -97427836.896386802), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 1}), -257671695.99240425), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 1}), 168937072.20694059), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 1}), -10939265.897138171), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 2}), -10171251.15759491), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 2}), 4939680.4738871539), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 2}), -305867.63199235022), "Nurbs P=2,N8");

  // first derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 0}), -1.609504132231405), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 1}), 1.479338842975206), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 2}), 0.1301652892561984), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 0}), 3.380916604057099), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 1}), -4.507888805409466), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 2}), 1.126972201352367), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // third derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 0}), -9.220681647428448), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 1}), 12.29424219657127), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 2}), -3.073560549142818), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fourth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 0}), 33.52975144519434), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 1}), -44.70633526025914), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 2}), 11.17658381506479), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fifth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 0}), -152.4079611145197), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 1}), 203.2106148193597), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 2}), -50.80265370483993), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 1,1 mixed derivative
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 0}), 66.39293764087149), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 0}), -51.50413223140495), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 1}), -39.5762584522915), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 1}), 26.62809917355371), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 2}), -3.676183320811419), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 2}), 1.735537190082645), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 2,1 mixed derivative
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 0}), -4511.311932245063), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 0}), 4043.131480090156), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 0}), -193.1404958677686), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 1}), 1791.166723584455), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 1}), -1318.232907588279), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 1}), 94.67768595041321), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 2}), 178.898026091114), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 2}), -91.08940646130728), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 2}), 5.900826446280992), "Nurbs P=2,N8");

  // std::cout<<std::setprecision(16)<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,2})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,2})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,2})<<std::endl;

  // 2,1 mixed derivative
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 0}), 430363.3606745686), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 0}), -408827.1566149851), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 0}), 21583.77160030052), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 1}), -118703.546081676), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 1}), 88813.60562803084), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 1}), -6462.50939143501), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 2}), -12947.75940540574), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 2}), 6609.384604876715), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 2}), -429.1510142749812), "Nurbs P=2,N8");
}

#include <dune/grid/test/checkindexset.hh>
#include <dune/iga/gridCapabilities.hh>
void gridCheck() {
  TestSuite test;

  auto circle = makeCircularArc();
  circle.controlPoints.directGet(0).p[1] += 3;
  const auto patch = makeSurfaceOfRevolution(circle, {2.0, 0, 0}, {0, 1, 0}, 360.0);

  IGA::NURBSGrid<2UL, 3UL> grid(patch);
  grid.globalRefine(2);
  //  grid.globalRefineInDirection(0,2);
  auto gridView = grid.leafGridView();

  //  Dune::checkIndexSet(grid,gridView_,std::cout);
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  //  gridCheck();
  //  std::cout << "test NURBS grid surface" << std::endl;
  //  testNURBSGridSurface();
  //  std::cout << "done with NURBS grid surface" << std::endl;
  //
  //  std::cout << "test NURBS grid Curve" << std::endl;
  //  testNURBSGridCurve();
  //  std::cout << "done with NURBS grid Curve" << std::endl;
  //
  //  std::cout << "test NURBS grid Curve" << std::endl;
  //  testNURBSCurve();
  //  std::cout << "done with NURBS curve" << std::endl;
  //
  //  testNURBSSurface();
  //  std::cout << "done with NURBS surface" << std::endl;
  //
  //  testNurbsGridCylinder();
  //  std::cout << "done with NURBS surface cylinder" << std::endl;
  //
  testTorusGeometry();
  std::cout << "done with NURBS torus " << std::endl;
  //
  testNurbsBasis();
  std::cout << "done with NURBS basis test " << std::endl;

  testBsplineBasisFunctions();
  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
