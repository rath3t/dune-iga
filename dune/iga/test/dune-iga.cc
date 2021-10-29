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
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

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

  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
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

  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};

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


  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
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


  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
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

  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
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
  const int dim                    = 2;
  const int dimworld               = 3;
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

  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
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

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

  std::cout << "test NURBS grid surface" << std::endl;
  testNURBSGridSurface();
  std::cout << "done with NURBS grid surface" << std::endl;

  std::cout << "test NURBS grid Curve" << std::endl;
  testNURBSGridCurve();
  std::cout << "done with NURBS grid Curve" << std::endl;

  std::cout << "test NURBS grid Curve" << std::endl;
  testNURBSCurve();
  std::cout << "done with NURBS curve" << std::endl;

  testNURBSSurface();
  std::cout << "done with NURBS surface" << std::endl;

  testNurbsGridCylinder();
  std::cout << "done with NURBS surface cylinder" << std::endl;

  testNURBSGridSurfaceOfRevolution();
  std::cout << "done with NURBS surface cylinder" << std::endl;
  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
