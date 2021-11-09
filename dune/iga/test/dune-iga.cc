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

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/bsplinebasis.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/iga/nurbsbasis.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/common/float_cmp.hh>



using namespace Dune;
using namespace Dune::IGA;









void testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};
//  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<ControlPoint > controlPoints
      = {{.p = {1, 3, 4}, .w = 1},
         {.p = {2, 2, 2}, .w = 3},
         {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 2},
         {.p = {4, 7, 2}, .w = 1},
         {.p = {8, 6, 2}, .w = 1},
         {.p = {2, 9, 9}, .w = 7},
         {.p = {1, 4, 3}, .w = 1},
         {.p = {1, 7, 1}, .w = 5}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet                       = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

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

  const int dim                    = 2;
  const int dimworld               = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
//  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
//      = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
 using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0, 0, 1}, .w = 2},        {.p = {1, 0, 1}, .w = 2},     {.p = {2, 0, 2}, .w = 1}},
          {{.p = {0, 1, 0}, .w = 1},        {.p = {1, 1, 0}, .w = 4},     {.p = {2, 1, 0}, .w = 1}},
          {{.p = {0, 2, 1}, .w = 1},        {.p = {1, 2, 2}, .w = 2},     {.p = {2, 2, 2}, .w = 4}}};

  std::array< int, dim> dimsize
      = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");
}


void testNURBSGridSurfaceOfRevolution() {
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

//  const int dim                    = 1;
//  const int dimworld               = 3;
//  const std::array<int, dim> order = {2};
//
//  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0,0, 1, 1,1}}};
//  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
//
//  const std::vector<ControlPoint > controlPoints
//      =  {{{.p = {0, 0, 0}, .w = 1},        {.p = {1, -1, 0}, .w = std::sqrt(2.0)/2.0}, {.p = {2, 0, 0}, .w = 1}}};
//
//  std::array<unsigned int, dim> dimsize = { static_cast<unsigned int>(controlPoints.size())};

  auto circle = makeCircularArc({0,0,0},{1,0,0},{0,1,0});

//  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);
//  IGA::NURBSPatchData<dim, dimworld> patchLine(knotSpans, controlNet, order);

  auto patch = makeSurfaceOfRevolution({2.0,0,0},{0,1,0},circle,360.0);

  IGA::NURBSGrid<2, 3> grid(patch.getKnots(), patch.getControlPoints(), patch.getOrder());
  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();


  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-SurfaceRevolution");
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
  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0, 0, 1}, .w = 2},        {.p = {1, 0, 1}, .w = 2},     {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1},        {.p = {1, 1, 0}, .w = 4},     {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1},        {.p = {1, 2, 2}, .w = 2},     {.p = {2, 2, 2}, .w = 4}}};


  std::array< int, dim> dimsize
      = {static_cast< int>(controlPoints.size()), static_cast< int>(controlPoints[0].size())};
  //

//  auto weightNet  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

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

  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
  const std::vector<ControlPoint > controlPoints
      = {{.p = {1, 3, 4}, .w = 2},
         {.p = {2, 2, 2}, .w = 2},
         {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 1},
         {.p = {4, 7, 2}, .w = 4},
         {.p = {8, 6, 2}, .w = 2},
         {.p = {2, 9, 9}, .w = 1},
         {.p = {1, 4, 3}, .w = 2},
         {.p = {1, 7, 1}, .w = 4}};


  std::array< int, dim> dimsize = {static_cast< int>(controlPoints.size())};
  auto controlNet                       = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

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

  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::GlobalCoordinateType;
  const std::vector<ControlPoint > controlPoints
      = {{ {1, 3, 4}},
         { {2, 2, 2}},
         { {3, 4, 5}},
         { {5, 1, 7}},
         { {4, 7, 2}},
         { {8, 6, 2}},
         { {2, 9, 9}},
         { {1, 4, 3}},
         { {1, 7, 1}}};

  std::array< int, dim> dimsize = {static_cast< int>(controlPoints.size())};
  auto controlNet                       = MultiDimensionNet<dim,typename Dune::IGA::NURBSPatchData<dim,dimworld>::GlobalCoordinateType>(dimsize, controlPoints);

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
  constexpr int dim                    = 2;
  constexpr int dimworld               = 3;
  const std::array<int, dim> order = {2, 1};
  constexpr double invsqr2         = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l                                       = 10;
  const double rad                                     = 5;
//  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};


  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0,   0, rad}, .w =       1},  {.p = {0,   l, rad}, .w = 1      }},
          {{.p = {rad, 0, rad}, .w = invsqr2},  {.p = {rad, l, rad}, .w = invsqr2}},
//          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
          {{.p = {rad, 0,   0}, .w =       1},  {.p = {rad, l,   0}, .w = 1      }}};

  std::array< int, dim> dimsize
      = {static_cast< int>(controlPoints.size()), static_cast< int>(controlPoints[0].size())};
  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

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

//template <int dim>
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
  constexpr int dim                    = 2;
  constexpr int dimworld               = 3;
  const std::array<int, dim> order = {2, 1};
  constexpr double invsqr2         = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l                                       = 10;
  const double rad                                     = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};


  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0,   0, rad}, .w =       1},  {.p = {0,   l, rad}, .w = 1      }},
         {{.p = {rad, 0, rad}, .w = invsqr2},  {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0,   0}, .w =       1},  {.p = {rad, l,   0}, .w = 1      }}};

  std::array< int, dim> dimsize
      = {static_cast< int>(controlPoints.size()), static_cast< int>(controlPoints[0].size())};
  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
//  grid.globalRefine(5);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite test;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, order);

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");
  using GridView = decltype(gridView);
  //  Dune::Functions::BasisFactory::nurbs<2>(knotSpans,order);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, knotSpans, order);

    // Test open knot vectors
    std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

      {
        // Check basis created via its constructor
        Functions::NurbsBasis<GridView> basis2(gridView, knotSpans, order);
//        test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
        test.subTest(checkBasis(basis2, AllowZeroBasisFunctions()));

      }

      {
        // Check basis created via makeBasis
        using namespace Functions::BasisFactory;
        auto basis2 = makeBasis(gridView, nurbs<dim>(knotSpans, order));
//          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions()));
      }

      {
        // Check whether a B-Spline basis can be combined with other bases.
        using namespace Functions::BasisFactory;
        auto basis2 = makeBasis(gridView,
                               power<2>(
                                   nurbs<dim>(knotSpans, order)
                                       ));
//          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions(), EnableContinuityCheck()));
          test.subTest(checkBasis(basis2, AllowZeroBasisFunctions()));
      }
}

void testBsplineBasisFunctions()
{
  std::vector<double> knots = {0,0,0,0.5,0.5,2,2,3,3,3};
  int degree = 2;
//std::vector<double> N;

TestSuite test;

auto N = Dune::IGA::Bspline<double>::basisFunctions(0.3,knots,degree);
using Dune::FloatCmp::eq;
test.check(eq(N[0], 0.16) ,"P=2,N0,u=0.3");
test.check(eq(N[1], 0.48) ,"P=2,N1,u=0.3");
test.check(eq(N[2], 0.36) ,"P=2,N2,u=0.3");

N = Dune::IGA::Bspline<double>::basisFunctions(0.5,knots,degree); //try knot span boundary

test.check(eq(N[0], 1.0),"P=2,N0,u=0.5");
test.check(eq(N[1], 0.0),"P=2,N1,u=0.5");
test.check(eq(N[2], 0.0),"P=2,N2,u=0.5");

N = Dune::IGA::Bspline<double>::basisFunctions(0.0,knots,degree); //try left end

test.check(eq(N[0], 1.0),"P=2,N0,u=0.0");
test.check(eq(N[1], 0.0),"P=2,N1,u=0.0");
test.check(eq(N[2], 0.0),"P=2,N2,u=0.0");

N = Dune::IGA::Bspline<double>::basisFunctions(3.0,knots,degree); //try right end

test.check(eq(N[0], 0.0),"P=2,N0,u=3");
test.check(eq(N[1], 0.0),"P=2,N1,u=3");
test.check(eq(N[2], 1.0),"P=2,N2,u=3");

knots = {0,0,0,0.5,1,1,1};
degree = 2;

N = Dune::IGA::Bspline<double>::basisFunctions(0.1,knots,degree);
test.check(eq(N[0], 0.64),"P=2,N0,u=0.1");
test.check(eq(N[1], 0.34),"P=2,N1,u=0.1");
test.check(eq(N[2], 0.02),"P=2,N2,u=0.1");

N = Dune::IGA::Bspline<double>::basisFunctions(0.01,knots,degree);
test.check(eq(N[0], 0.9604),"P=2,N0,u=0.01");
test.check(eq(N[1], 0.0394),"P=2,N1,u=0.01");
test.check(eq(N[2], 0.0002),"P=2,N2,u=0.01");

knots = {0,0,0,0,0.5,1,1,2,2,2,2};
degree = 3;

auto N2 = Dune::IGA::Bspline<double>::basisFunctions(1.45,knots,degree);
test.check(eq(N2[0], 0.1109166666666667),"P=3,N0,u=1.45");
test.check(eq(N2[1], 0.4638333333333333),"P=3,N1,u=1.45");
test.check(eq(N2[2], 0.334125),"P=3,N2,u=1.45");
test.check(eq(N2[3], 0.091125),"P=3,N3,u=1.45");


auto dN = Dune::IGA::Bspline<double>::basisFunctionDerivatives(1.45, knots, degree, 3);
//check values
test.check(eq(dN[0][0], 0.1109166666666667),"P=3,dN00,u=1.45");
test.check(eq(dN[0][1], 0.4638333333333333),"P=3,dN01,u=1.45");
test.check(eq(dN[0][2], 0.334125),"P=3,dN02,u=1.45");
test.check(eq(dN[0][3], 0.091125),"P=3,dN03,u=1.45");

//check first derivatives
test.check(eq(dN[1][0], -0.605),"P=3,dN10,u=1.45");
test.check(eq(dN[1][1], -0.88),"P=3,dN11,u=1.45");
test.check(eq(dN[1][2], 0.8775),"P=3,dN12,u=1.45");
test.check(eq(dN[1][3], 0.6075),"P=3,dN13,u=1.45");

//check second derivatives
test.check(eq(dN[2][0], 2.2),"P=3,dN20,u=1.45");
test.check(eq(dN[2][1], -2.8),"P=3,dN21,u=1.45");
test.check(eq(dN[2][2], -2.1),"P=3,dN22,u=1.45");
test.check(eq(dN[2][3], 2.7),"P=3,dN23,u=1.45");

//check third derivatives
test.check(eq(dN[3][0], -4.0),"P=3,dN30,u=1.45");
test.check(eq(dN[3][1], 16.0),"P=3,dN31,u=1.45");
test.check(eq(dN[3][2], -18.0),"P=3,dN32,u=1.45");
test.check(eq(dN[3][3], 6.0),"P=3,dN33,u=1.45");
//https://godbolt.org/z/Ta3fzW553
auto Nf = Dune::IGA::Bspline<double>(knots,degree);
std::vector<double> evalPoints={1, 0.1714677640603567, 0.001371742112482855,
                                  0.0740740740740741, 0.00274348422496571,
                                  0.4682213077274805, 0.1975308641975309, 0.05852766346593506,
                                  0.007315957933241894, 0};
for(int i = 0; i<evalPoints.size() ; ++i) {
  test.check(eq(Nf(i / (evalPoints.size() - 1.0) * 2.0)[0], evalPoints[i]));
}

std::array<double,2> xieta{0.2,0.25};
std::array<std::vector<double>,2> knots2 = {{{0,0,0,0.5,0.5,2,2,3,3,3},{0,0,0,2,2,2}}};
std::array<int,2> degree2{2,2};
const std::vector<std::vector<double >> weights2
    = {{{1,2,3,4,5,6,7 },{8,9,10,11,12,13,14},{15,16,17,18,19,20,21}}};
std::array< int, 2> dimsize
    = {static_cast< int>(weights2.size()), static_cast< int>(weights2[0].size())};
MultiDimensionNet<2,double> weightNet(dimsize,weights2);

auto N_Nurbs = Dune::IGA::Nurbs<double,2>::basisFunctions(xieta,knots2,degree2,weightNet);

test.check(N_Nurbs.size() == (degree2[0]+1)*(degree2[0]+1) );

test.check(eq(N_Nurbs[0], 0.04023722627737226),"Nurbs2d P=2,N0"); //check ansatzfunctions in domain
test.check(eq(N_Nurbs[1], 0.4291970802919708),"Nurbs2d P=2,N1");
test.check(eq(N_Nurbs[2], 0.2682481751824818),"Nurbs2d P=2,N2");
test.check(eq(N_Nurbs[3], 0.02299270072992701),"Nurbs2d P=2,N3");
test.check(eq(N_Nurbs[4], 0.137956204379562),"Nurbs2d P=2,N4");
test.check(eq(N_Nurbs[5],0.08175182481751825),"Nurbs2d P=2,N5");
test.check(eq(N_Nurbs[6], 0.002463503649635036),"Nurbs2d P=2,N6");
test.check(eq(N_Nurbs[7], 0.01094890510948905),"Nurbs2d P=2,N7");
test.check(eq(N_Nurbs[8], 0.006204379562043796),"Nurbs2d P=2,N8");
test.check(eq(std::accumulate(N_Nurbs.begin(),N_Nurbs.end(),0.0), 1.0),"partition of unity in domain"); //partition of unity in domain

xieta={0,0.1};
N_Nurbs = Dune::IGA::Nurbs<double,2>::basisFunctions(xieta,knots2,degree2,weightNet);
test.check(eq(N_Nurbs[0], 0.8204545454545455),"Nurbs P=2,N0"); //check ansatzfunctions on boundaries
test.check(eq(N_Nurbs[1], 0.0),"Nurbs P=2,N1");
test.check(eq(N_Nurbs[2], 0.0),"Nurbs P=2,N2");
test.check(eq(N_Nurbs[3], 0.1727272727272728),"Nurbs P=2,N3");
test.check(eq(N_Nurbs[4], 0.0),"Nurbs P=2,N4");
test.check(eq(N_Nurbs[5],0.0),"Nurbs P=2,N5");
test.check(eq(N_Nurbs[6], 0.00681818181818182),"Nurbs P=2,N6");
test.check(eq(N_Nurbs[7], 0.0),"Nurbs P=2,N7");
test.check(eq(N_Nurbs[8], 0.0),"Nurbs P=2,N8");

test.check(eq(std::accumulate(N_Nurbs.begin(),N_Nurbs.end(),0.0), 1.0),"partition of unity on boundary"); //partition of unity on boundary

//std::ranges::for_each(N_Nurbs,[](auto& Ni){std::cout<<Ni<<" ";});

}


int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

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
//  testNURBSGridSurfaceOfRevolution();
//  std::cout << "done with NURBS surface cylinder" << std::endl;

  testNurbsBasis();
  std::cout << "done with NURBS basis test "<< std::endl;

  testBsplineBasisFunctions();
  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
