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
  int subSampling = 20;

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};
//  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};

  const std::vector<FieldVector<double, dimworld> > controlPoints
      = {{1, 3, 4}, {2, 2, 2}, {3, 4, 5}, {5, 1, 7}, {4, 7, 2}, {8, 6, 2}, {2, 9, 9}, {1, 4, 3}, {1, 7, 1}};
//  const std::vector<FieldVector<double, dimworld> > controlPoints
//      = {{-1, 0, 0}, {1, 0, 0},};

          //    const std::vector<FieldVector<double, dimworld> > controlPoints
//        = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0}, {5, 0, 0}, {6, 0, 0}, {7, 0, 0}, {8, 0, 0}};
  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = Dune::IGA::MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

//  const std::vector<FieldVector<double, 1> > weight = {{1}, {1}};
  const std::vector<double > weight = {{1}, {3}, {1}, {2}, {1}, {1}, {7}, {1}, {5}};
  auto weightNet                                    = MultiDimensionNet<dim, double>(dimsize, weight);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, weightNet, order);
  grid.globalRefine(1);
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
  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPoints
      = {{{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}}};
  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
  //
  const std::vector<std::vector<double > > weight = {{{2}, {2}, {1}}, {{1}, {4}, {1}}, {{1}, {2}, {4}}};
  auto weightNet                                                  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, weightNet, order);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");
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
  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPoints
      = {{{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}}};
  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
  //
  const std::vector<std::vector<double > > weight = {{{2}, {2}, {1}}, {{1}, {4}, {1}}, {{1}, {2}, {4}}};

  auto weightNet  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, weightNet, order);

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

  const std::vector<FieldVector<double, dimworld> > controlPoints
      = {{1, 3, 4}, {2, 2, 2}, {3, 4, 5}, {5, 1, 7}, {4, 7, 2}, {8, 6, 2}, {2, 9, 9}, {1, 4, 3}, {1, 7, 1}};
  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

  const std::vector<double > weight = {{2}, {2}, {1}, {1}, {4}, {2}, {1}, {2}, {4}};

  auto weightNet = MultiDimensionNet<dim, double>(dimsize, weight);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, weightNet, order);

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

  const std::vector<FieldVector<double, dimworld> > controlPoints
      = {{1, 3, 4}, {2, 2, 2}, {3, 4, 5}, {5, 1, 7}, {4, 7, 2}, {8, 6, 2}, {2, 9, 9}, {1, 4, 3}, {1, 7, 1}};
  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

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
  const double cylLength                               = 10;
  const double cylRadius                               = 5;
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPoints
      = {{{0, 0, cylRadius}, {0, cylLength, cylRadius}},
         {{cylRadius, 0, cylRadius}, {cylRadius, cylLength, cylRadius}},
         {{cylRadius, 0, 0}, {cylRadius, cylLength, 0}}};

  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
  auto controlNet = MultiDimensionNetFVd<dim, dimworld>(dimsize, controlPoints);

  const std::vector<std::vector<double > > weight = {{{1}, {1}}, {{invsqr2}, {invsqr2}}, {{1}, {1}}};
  auto weightNet = MultiDimensionNet<dim, double>(dimsize, weight);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, weightNet, order);

  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, weightNet, order);

  testSuite.check(nurbsPatch.size(0) == 1);
  testSuite.check(nurbsPatch.size(1) == 4);
  testSuite.check(nurbsPatch.size(2) == 6);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("foo");

  grid.globalRefine(2);

  auto gridViewRefined  = grid.leafGridView();

  SubsamplingVTKWriter<decltype(gridViewRefined)> vtkWriter2(gridViewRefined, refinementIntervals1);
  vtkWriter.write("fooRefine");


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
  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
