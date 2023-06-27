// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include <iostream>

#include "dune/iga/gridcapabilities.hh"
#include "dune/iga/nurbsbasis.hh"
#include "dune/iga/nurbsgrid.hh"
#include "dune/iga/nurbspatch.hh"
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>

using namespace Dune;
using namespace Dune::IGA;

template <typename T, int worldDim, int Items>
struct Compare {
  constexpr bool operator()(const std::array<FieldVector<double, worldDim>, Items>& lhs,
                            const std::array<FieldVector<double, worldDim>, Items>& rhs) const {
    return std::ranges::lexicographical_compare(std::ranges::join_view(lhs), std::ranges::join_view(rhs));
  };
};

auto checkUniqueEdges(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  constexpr int gridDimension      = std::remove_cvref_t<decltype(gridView)>::dimension;
  std::set<std::array<FieldVector<double, gridDimensionworld>, 2>, Compare<double, gridDimensionworld, 2>>
      edgeVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    edgeVertexPairSet.clear();
    for (int edgeIndex = 0; edgeIndex < element.subEntities(gridDimension - 1); ++edgeIndex) {
      auto edge = element.template subEntity<gridDimension - 1>(edgeIndex);
      std::array<FieldVector<double, gridDimensionworld>, 2> pair;
      for (int c = 0; c < edge.geometry().corners(); ++c)
        pair[c] = edge.geometry().corner(c);
      bool inserted = edgeVertexPairSet.insert(pair).second;
      t.require(inserted) << "Duplicate edge detected in Element " << eleIndex << " Edges: " << pair[0] << ", "
                          << pair[1];
    }
    ++eleIndex;
  }
  return t;
}

auto checkUniqueSurfaces(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  std::set<std::array<FieldVector<double, gridDimensionworld>, 4>, Compare<double, gridDimensionworld, 4>>
      edgeVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    edgeVertexPairSet.clear();
    for (int edgeIndex = 0; edgeIndex < element.subEntities(2); ++edgeIndex) {
      auto edge = element.template subEntity<2>(edgeIndex);
      std::array<FieldVector<double, gridDimensionworld>, 4> tuple;
      for (int c = 0; c < edge.geometry().corners(); ++c)
        tuple[c] = edge.geometry().corner(c);

      t.require(edgeVertexPairSet.insert(tuple).second)
          << "Duplicate surface detected in Element " << eleIndex << " Surfaces: " << tuple[0] << ", " << tuple[1]
          << ", " << tuple[2] << ", " << tuple[3];
    }
    ++eleIndex;
  }
  return t;
}

auto testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const auto dim      = 1;
  const auto dimworld = 3;

  const std::array<int, dim> order                     = {3};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 3, 4, 4, 4, 5, 5, 5, 5}}};
  //  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints
      = {{.p = {1, 3, 4}, .w = 1}, {.p = {2, 2, 2}, .w = 3}, {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 2}, {.p = {4, 7, 2}, .w = 1}, {.p = {8, 6, 2}, .w = 1},
         {.p = {2, 9, 9}, .w = 7}, {.p = {1, 4, 3}, .w = 1}, {.p = {1, 7, 1}, .w = 5}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto additionalKnots    = std::vector<double>(2);
  additionalKnots[0]      = 0.5;
  additionalKnots[1]      = 3.5;
  patchData               = knotRefinement<dim>(patchData, additionalKnots, 0);
  patchData               = degreeElevate(patchData, 0, 1);

  TestSuite t;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  grid.globalRefine(3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  for (int eleIndex = 0;
       auto& ele :
       elements(gridView))  // This test also exists in grid check, but it is more convenient to debug it here
  {
    const int numCorners = ele.subEntities(dim);
    for (int c = 0; c < numCorners; ++c) {
      auto vertex = ele.template subEntity<dim>(c).geometry();
      auto elegeo = ele.geometry();
      t.check(Dune::FloatCmp::eq((elegeo.corner(c) - vertex.corner(0)).two_norm(), 0.0))
          << "Corner " << c << " " << elegeo.corner(c) << " Alt: " << vertex.corner(0);
    }
    ++eleIndex;
  }

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  //  vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-CurveNewFineResample-R");
  //  vtkWriter.write("NURBSGridTest-CurveNewFineResample_knotRefine");
  gridcheck(grid);
  return t;
}

void testNURBSGridSurface() {
  TestSuite t;
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2;
  const auto dimworld              = 3;
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

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");

  t.subTest(checkUniqueEdges(gridView));
}

auto test3DGrid() {
  constexpr std::size_t dim        = 3;
  constexpr std::size_t dimworld   = 3;
  const std::array<int, dim> order = {2, 2, 2};
  // quarter cylindrical hyperSurface
  const double lx = 2;
  const double ly = 1;
  const double lz = 1;

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  std::vector<std::vector<std::vector<ControlPoint>>> controlp;
  std::array<int, dim> dimSize = {3, 3, 3};
  for (int i = 0; i < dimSize[0]; ++i) {
    controlp.emplace_back();
    for (int j = 0; j < dimSize[1]; ++j) {
      controlp[i].emplace_back();
      for (int k = 0; k < dimSize[2]; ++k) {
        controlp[i][j].push_back(
            {.p = {i * i * lx / (dimSize[0] - 1) + 1, 2 * i * j * k * ly / (dimSize[1] - 1) + (k + 1) + (j + 1),
                   k * k * lz / (dimSize[2] - 1)},
             .w = 1});
      }
    }
  }

  nurbsPatchData.controlPoints = MultiDimensionNet(dimSize, controlp);
  nurbsPatchData.degree        = order;

  auto additionalKnots = std::vector<double>(2);
  additionalKnots[0]   = 0.1;
  additionalKnots[1]   = 0.3;
  //  additionalKnots[2] = 0.6;
  //  additionalKnots[1] = 3.5;
  //  nurbsPatchData = knotRefinement<dim>(nurbsPatchData, additionalKnots, 2);
  //  nurbsPatchData = degreeElevate(nurbsPatchData,0,1);
  //  nurbsPatchData = degreeElevate(nurbsPatchData, 1, 2);
  //  nurbsPatchData = degreeElevate(nurbsPatchData,2,1);
  IGA::NURBSGrid<3, 3> grid(nurbsPatchData);
  //  grid.globalRefine(1);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(0,1);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(1,2);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(2, 3);
  //  gridcheck(grid);

  auto gridView = grid.leafGridView();
  TestSuite t;
  t.subTest(checkUniqueEdges(gridView));
  t.subTest(checkUniqueSurfaces(gridView));

  const int subSampling = 10;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Solid2");

  Dune::GeometryChecker<decltype(grid)> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(grid, gridView, std::cout);

  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);

  return t;
}

void testFactory() {
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
  nurbsPatchData.controlPoints = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  nurbsPatchData.degree        = order;
};

auto testTorusGeometry() {
  const double R       = 2.0;
  const double r       = 1.0;
  auto circle          = makeCircularArc(r);
  auto nurbsPatchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 360.0);
  nurbsPatchData       = degreeElevate(nurbsPatchData, 0, 1);
  nurbsPatchData       = degreeElevate(nurbsPatchData, 1, 2);
  auto additionalKnots = std::vector<double>(1);
  additionalKnots[0]   = 0.1;
  nurbsPatchData       = knotRefinement<2>(nurbsPatchData, additionalKnots, 1);
  for (int refDirection = 0; refDirection < 2; ++refDirection)
    nurbsPatchData = degreeElevate(nurbsPatchData, refDirection, 1);

  IGA::NURBSGrid<2, 3> grid(nurbsPatchData);
  grid.globalRefine(1);
  grid.globalRefineInDirection(1, 1);
  grid.globalRefineInDirection(0, 2);
  grid.globalDegreeElevate(2);

  auto gridView = grid.leafGridView();

  const int subSampling = 2;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-SurfaceRevolutionFLAT");

  TestSuite test;
  Dune::GeometryChecker<IGA::NURBSGrid<2UL, 3UL>> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(grid, gridView, std::cout);

  double area = 0.0;
  for (auto&& ele : elements(gridView)) {
    area += ele.geometry().volume();
  }
  const double pi                        = std::numbers::pi_v<double>;
  const double referenceTorusSurfaceArea = 4.0 * pi * pi * r * R;
  test.check(area - referenceTorusSurfaceArea < 1e-4)
      << "The integrated area of the torus surface is wrong! Expected: " << referenceTorusSurfaceArea
      << " Computed: " << area;

  double gaussBonnet = 0.0;
  for (auto& ele : elements(gridView)) {
    const auto rule
        = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * (*std::ranges::max_element(grid.patchData().degree)));
    for (auto& gp : rule) {
      const auto Kinc = ele.geometry().impl().gaussianCurvature(gp.position());
      const auto Kmax = 1 / (r * (R + r));
      const auto Kmin = -1 / (r * (R - r));
      test.check(Kinc < Kmax && Kinc > Kmin)
          << "The Gaussian curvature should be within bounds " << Kmin << " < " << Kinc << " < " << Kmax << std::endl;
      gaussBonnet += Kinc * gp.weight() * ele.geometry().integrationElement(gp.position());
    }
  }

  test.check(std::abs(gaussBonnet) < 1e-5)
      << "Gauss-Bonnet theorem dictates a vanishing integrated Gaussian curvature for the torus! It is " << gaussBonnet;
  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);

  gridcheck(grid);
  test.subTest(checkUniqueEdges(gridView));
  return test;
}

auto testNURBSSurface() {
  // parameters
  int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2UL;
  const auto dimworld              = 3UL;
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

  return testSuite;
}

auto testNURBSCurve() {
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
      = {{.p = {1, 3, 4}, .w = 2}, {.p = {2, 2, 2}, .w = 2}, {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 1}, {.p = {4, 7, 2}, .w = 4}, {.p = {8, 6, 2}, .w = 2},
         {.p = {2, 9, 9}, .w = 1}, {.p = {1, 4, 3}, .w = 2}, {.p = {1, 7, 1}, .w = 4}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 5);
  testSuite.check(patch.size(1) == controlPoints.size());
  return testSuite;
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
  const double invsqr2             = 1.0 / std::sqrt(2.0);
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

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  grid.globalRefine(5);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, order);

  testSuite.check(nurbsPatch.size(0) == 1);
  testSuite.check(nurbsPatch.size(1) == 4);
  testSuite.check(nurbsPatch.size(2) == 4);

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  VTKWriter vtkWriter(gridView);

  vtkWriter.write("ZylRefine");

  testSuite.subTest(checkUniqueEdges(gridView));
}

// template <int dim>
auto testNurbsBasis() {
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
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  nurbsPatchData.controlPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};
  nurbsPatchData.degree = order;

  IGA::NURBSGrid<dim, dimworld> grid(nurbsPatchData);
  //  grid.globalRefine(1);
  grid.globalRefineInDirection(0, 2);
  grid.globalDegreeElevate(2);
  //  grid.globalRefineInDirection(1, 3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite test;

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");
  using GridView = decltype(gridView);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.impl().getPatchData());

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

  {
    using namespace Functions::BasisFactory;
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView, nurbs());
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView);
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check basis created via makeBasis
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, nurbs());
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check whether a B-Spline basis can be combined with other bases.
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, power<2>(nurbs()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    grid.globalDegreeElevate(1);
    auto gridViewNew = grid.leafGridView();
    // Check lower order basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridViewNew, gridViewNew.impl().lowerOrderPatchData());
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }
  return test;
}

auto testBsplineBasisFunctions() {
  std::vector<double> knots = {0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3};
  int degree                = 2;
  TestSuite test;

  auto N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.3, knots, degree);
  using Dune::FloatCmp::eq;
  test.check(eq(N[0], 0.16), "P=2,N0,u=0.3");
  test.check(eq(N[1], 0.48), "P=2,N1,u=0.3");
  test.check(eq(N[2], 0.36), "P=2,N2,u=0.3");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.5, knots, degree);  // try knot span boundary

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.5");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.5");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.5");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.0, knots, degree);  // try left end

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.0");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.0");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.0");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(3.0, knots, degree);  // try right end

  test.check(eq(N[0], 0.0), "P=2,N0,u=3");
  test.check(eq(N[1], 0.0), "P=2,N1,u=3");
  test.check(eq(N[2], 1.0), "P=2,N2,u=3");

  knots  = {0, 0, 0, 0.5, 1, 1, 1};
  degree = 2;

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.1, knots, degree);
  test.check(eq(N[0], 0.64), "P=2,N0,u=0.1");
  test.check(eq(N[1], 0.34), "P=2,N1,u=0.1");
  test.check(eq(N[2], 0.02), "P=2,N2,u=0.1");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.01, knots, degree);
  test.check(eq(N[0], 0.9604), "P=2,N0,u=0.01");
  test.check(eq(N[1], 0.0394), "P=2,N1,u=0.01");
  test.check(eq(N[2], 0.0002), "P=2,N2,u=0.01");

  knots  = {0, 0, 0, 0, 0.5, 1, 1, 2, 2, 2, 2};
  degree = 3;

  auto N2 = Dune::IGA::BsplineBasis1D<double>::basisFunctions(1.45, knots, degree);
  test.check(eq(N2[0], 0.1109166666666667), "P=3,N0,u=1.45");
  test.check(eq(N2[1], 0.4638333333333333), "P=3,N1,u=1.45");
  test.check(eq(N2[2], 0.334125), "P=3,N2,u=1.45");
  test.check(eq(N2[3], 0.091125), "P=3,N3,u=1.45");

  auto dN = Dune::IGA::BsplineBasis1D<double>::basisFunctionDerivatives(1.45, knots, degree, 3);
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
  auto Nf                           = Dune::IGA::BsplineBasis1D<double>(knots, degree);
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

  Dune::FieldVector<double, 2> xieta{0.2, 0.25};
  std::array<std::vector<double>, 2> knots2 = {{{0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3}, {0, 0, 0, 2, 2, 2}}};
  std::array<int, 2> degree2{2, 2};
  const std::vector<std::vector<double>> weights2
      = {{{1, 2, 3, 4, 5, 6, 7}, {8, 9, 10, 11, 12, 13, 14}, {15, 16, 17, 18, 19, 20, 21}}};
  std::array<int, 2> dimsize = {static_cast<int>(weights2.size()), static_cast<int>(weights2[0].size())};
  MultiDimensionNet<2UL, double> weightNet(dimsize, weights2);

  auto N_Nurbs = Dune::IGA::Nurbs<2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();

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
  N_Nurbs = Dune::IGA::Nurbs<2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();
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
  auto dN_Nurbs = Dune::IGA::Nurbs<2>::basisFunctionDerivatives(xieta, knots2, degree2, weightNet, 5);

  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 0}), 0.8204545454545455),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 1}), 0.1727272727272728), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 2}), 0.00681818181818182), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // first derivative in u dir
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 0}), -24.166115702479335),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 0}), 26.25454545454545), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 1}), -5.0876033057851231), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 1}), 3.1090909090909089), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 2}), -0.20082644628099172), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 2}), 0.090909090909090925), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in u dir
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 0}), 1236.8386175807659),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 0}), -1441.6132231404956), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 0}), 98.454545454545439), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 1}), 260.38707738542439), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 1}), -170.71735537190079), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 1}), 11.054545454545453), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 2}), 10.278437265214125), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 2}), -4.9917355371900829), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 2}), 0.30909090909090914), "Nurbs P=2,N8");

  // third derivative in u dir
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 0}), -94449.494433440283),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 0}), 110086.82794891056), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 0}), -7518.3471074380141), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 1}), -19884.104091250589), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 1}), 13036.598046581514), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 1}), -844.16528925619821), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 2}), -784.89884570726042), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 2}), 381.18707738542452), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 2}), -23.603305785123968), "Nurbs P=2,N8");

  // fourth derivative in u dir
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 0}), 9616675.7968593743),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 0}), -11208840.663889075), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 0}), 765504.432757325), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 1}), 2024563.3256546052), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 1}), -1327362.7101973905), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 1}), 85951.374906085635), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 2}), 79916.973381102871), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 2}), -38811.775151970498), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 2}), 2403.2456799398947), "Nurbs P=2,N8");

  // fifth derivative in u dir
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 0}), -1223940555.9639204),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 0}), 1426579720.8586094), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 0}), -97427836.896386802), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 1}), -257671695.99240425), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 1}), 168937072.20694059), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 1}), -10939265.897138171), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 2}), -10171251.15759491), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 2}), 4939680.4738871539), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 2}), -305867.63199235022), "Nurbs P=2,N8");

  // first derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 0}), -1.609504132231405),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 1}), 1.479338842975206), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 2}), 0.1301652892561984), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 0}), 3.380916604057099),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 1}), -4.507888805409466), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 2}), 1.126972201352367), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // third derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 0}), -9.220681647428448),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 1}), 12.29424219657127), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 2}), -3.073560549142818), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fourth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 0}), 33.52975144519434),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 1}), -44.70633526025914), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 2}), 11.17658381506479), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fifth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 0}), -152.4079611145197),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 1}), 203.2106148193597), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 2}), -50.80265370483993), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 1,1 mixed derivative
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 0}), 66.39293764087149),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 0}), -51.50413223140495), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 1}), -39.5762584522915), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 1}), 26.62809917355371), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 2}), -3.676183320811419), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 2}), 1.735537190082645), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 2,1 mixed derivative
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 0}), -4511.311932245063),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
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

  // 3,1 mixed derivative
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 0}), 430363.3606745686),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 0}), -408827.1566149851), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 0}), 21583.77160030052), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 1}), -118703.546081676), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 1}), 88813.60562803084), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 1}), -6462.50939143501), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 2}), -12947.75940540574), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 2}), 6609.384604876715), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 2}), -429.1510142749812), "Nurbs P=2,N8");
  return test;
}

#include "dune/iga/gridcapabilities.hh"
#include <dune/grid/test/checkindexset.hh>
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

template <typename C>
bool checkIfFinite(const C& v) {
  for (auto& vi : v) {
    if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vi)>>) {
      if (not std::isfinite(vi)) return false;
    } else {
      return checkIfFinite(vi);
    }
  }
  return true;
};

template <typename C>
bool checkIfZero(const C& v) {
  for (auto& vi : v) {
    if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vi)>>) {
      if (std::abs(vi) > 1e-13) return false;
    } else {
      return checkIfZero(vi);
    }
  }
  return true;
};

auto testPlate() {
  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 3;
  const std::array<int, gridDim> order = {2, 2};
  TestSuite t;

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0.5, 0, 0}, .w = 1}, {.p = {1, 0, 0}, .w = 1}},
         {{.p = {0, 0.5, 0}, .w = 1}, {.p = {0.5, 0.5, 0}, .w = 3}, {.p = {1, 0.5, 0}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {0.5, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<gridDim, dimworld>;

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto grid               = std::make_shared<Grid>(patchData);
  //  grid->globalRefine(1);
  //
  grid->globalDegreeElevate(4);
  //  grid->globalRefine(1);
  auto gridView = grid->leafGridView();
  Dune::GeometryChecker<decltype(grid)> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(*grid, gridView, std::cout);

  checkEntityLifetime(gridView, gridView.size(0));
  std::cout << "checkJacobians Start" << std::endl;
  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);
  std::cout << "checkJacobians End" << std::endl;
  std::cout << "checkIterators Start" << std::endl;
  checkIterators(gridView);
  std::cout << "checkIterators End" << std::endl;
  auto basis     = Dune::Functions::BasisFactory::makeBasis(gridView, Dune::Functions::BasisFactory::nurbs());
  auto localView = basis.localView();
  std::vector<Dune::FieldVector<double, 1>> N;
  std::vector<Dune::FieldVector<double, 1>> ddNi02;
  std::vector<Dune::FieldVector<double, 1>> ddNi20;
  std::vector<Dune::FieldVector<double, 1>> ddNi11;
  std::vector<Dune::FieldMatrix<double, 1, 2>> dN;
  for (auto& element : elements(gridView)) {
    localView.bind(element);
    const auto& fe     = localView.tree().finiteElement();
    auto& localBasis   = fe.localBasis();
    auto referenceEle  = referenceElement(element);
    auto geo           = element.geometry();
    auto localGeometry = referenceEle.template geometry<0>(0);
    for (int i = 0; i < localGeometry.corners(); ++i) {
      auto localV = localGeometry.corner(i);
      auto global = geo.global(localV);
      auto localT = geo.local(global);
      for (int j = 0; j < gridDim; ++j) {
        t.check(std::abs(localV[j] - localT[j]) < 1e-14)
            << "Dune::FloatCmp:eq(localV[i],localT[i]) failed Local " << localV << " Calculated: " << localT;
      }
      auto J  = geo.jacobianTransposed(localV);
      auto J2 = geo.impl().secondDerivativeOfPosition(localV);

      t.check(checkIfFinite(global)) << "checkIfFinite(global) failed";
      t.check(checkIfFinite(J)) << "checkIfFinite(J) failed";
      //      t.check(checkIfZero(J2))<<"checkIfZero(J2) "<<J2;  // J2 should be zero since the grid is linearly
      //      parametrized
      t.check(std::abs(J2[0][2]) < 1e-14)
          << "checkIfZero(J2[0][2])" << J2;  // The grid is plane, therefore, the curvature should be also inplane
      t.check(std::abs(J2[1][2]) < 1e-14)
          << "checkIfZero(J2[1][2])" << J2;  // The grid is plane, therefore, the curvature should be also inplane

      localBasis.evaluateFunction(localV, N);
      localBasis.evaluateJacobian(localV, dN);
      localBasis.partial({2, 0}, localV, ddNi20);
      localBasis.partial({1, 1}, localV, ddNi11);
      localBasis.partial({0, 2}, localV, ddNi02);

      t.check(checkIfFinite(N)) << "checkIfFinite(N) failed";
      t.check(checkIfFinite(dN)) << "checkIfFinite(dN) failed";
      t.check(checkIfFinite(ddNi20)) << "checkIfFinite(ddNi20) failed";
      t.check(checkIfFinite(ddNi11)) << "checkIfFinite(ddNi11) failed";
      t.check(checkIfFinite(ddNi02)) << "checkIfFinite(ddNi02) failed";
    }
  }
  t.subTest(checkUniqueEdges(gridView));
  return t;
}

auto testCurveHigherOrderDerivatives() {
  constexpr auto gridDim               = 1;
  constexpr auto dimworld              = 3;
  const std::array<int, gridDim> order = {2};
  TestSuite t;
  // parameters

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 0.5, 1, 1, 1}}};
  //  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints
      = {{.p = {0.086956521739130, -0.434782608695652, 0}, .w = 11.5},
         {.p = {0.200000000000000, 5.400000000000000, 0.200000000000000}, .w = 5},
         {.p = {1.857142857142857, 0.142857142857143, 0}, .w = 7},
         {.p = {3.714285714285714, 0.285714285714286, 2.000000000000000}, .w = 3.5}};

  std::array<int, gridDim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  std::vector<Dune::FieldVector<double, 3>> expectedControlPoints
      = {{0.086956521739130, -0.434782608695652, 0},
         {0.121212121212121, 1.333333333333333, 0.060606060606061},
         {0.320000000000000, 3.120000000000000, 0.120000000000000},
         {0.727272727272727, 3.727272727272727, 0.136363636363636},
         {1.538461538461539, 1.153846153846154, 0.038461538461538},
         {1.920000000000000, 0.506666666666667, 0.200000000000000},
         {2.476190476190476, 0.190476190476190, 0.666666666666667},
         {3.714285714285714, 0.285714285714286, 2.000000000000000}};

  std::vector<double> expectedWeights = {11.5, 8.25, 6.25, 5.5, 6.5, 6.25, 5.25, 3.5};

  auto patchDataP = std::make_shared<Dune::IGA::NURBSPatchData<gridDim, dimworld>>(patchData);

  NURBSPatchGeometry<gridDim, dimworld> geo(patchDataP);
  double vol  = geo.volume();
  int samples = 20;
  std::vector<Dune::FieldVector<double, 3>> evaluatedPoints;
  std::vector<Dune::FieldMatrix<double, 1, 3>> evaluatedJacobians;
  std::vector<Dune::FieldMatrix<double, 1, 3>> evaluatedHessians;

  for (int i = 0; i < samples + 1; ++i) {
    const Dune::FieldVector<double, 1> u = {i / static_cast<double>(samples)};
    std::cout << geo.global(u) << std::endl;
    evaluatedPoints.push_back(geo.global(u));
    const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
    evaluatedHessians.push_back(H);
    evaluatedJacobians.push_back(geo.jacobianTransposed(u));
  }

  for (int i = 1; i < 5; ++i) {
    auto patchDataN  = degreeElevate(patchData, 0, i);
    auto patchDataPN = std::make_shared<Dune::IGA::NURBSPatchData<gridDim, dimworld>>(patchDataN);
    NURBSPatchGeometry<gridDim, dimworld> geo2(patchDataPN);
    //    t.check(Dune::FloatCmp::eq(geo2.volume(),vol,1e-10))<<"vol "<<vol<<" is "<<geo2.volume();
    for (int j = 0; j < samples + 1; ++j) {
      const Dune::FieldVector<double, 1> u = {j / static_cast<double>(samples)};
      t.check(Dune::FloatCmp::eq(geo2.global(u), evaluatedPoints[j], 1e-10))
          << "evaluatedPoints[i] " << evaluatedPoints[j] << " is " << geo2.global(u);
      const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
      t.check(Dune::FloatCmp::eq(geo2.jacobianTransposed(u)[0], evaluatedJacobians[j][0], 1e-10))
          << "evaluatedJacobians[i] " << evaluatedJacobians[j] << " is " << geo2.jacobianTransposed(u);
      t.check(Dune::FloatCmp::eq(H[0], evaluatedHessians[j][0], 1e-10))
          << "evaluatedJacobians[i] " << evaluatedHessians[j] << " is " << H;
    }
  }

  patchData = degreeElevate(patchData, 0, 2);
  t.check(expectedControlPoints.size() == patchData.controlPoints.directGetAll().size())
      << "Size mismatch " << expectedControlPoints.size() << "is " << patchData.controlPoints.directGetAll().size();
  t.check(expectedWeights.size() == patchData.controlPoints.directGetAll().size())
      << "Size mismatch" << expectedWeights.size() << "is " << patchData.controlPoints.directGetAll().size();
  t.check(patchData.degree[0] == 4) << "Size mismatch";

  std::cout << "Knots:" << std::endl;
  for (int i = 0; auto kn : patchData.knotSpans[0]) {
    std::cout << kn << " ";
  }
  std::cout << std::endl;
  for (int i = 0; auto cp : patchData.controlPoints.directGetAll()) {
    std::cout << cp.p << " w: " << cp.w << std::endl;
    t.check(Dune::FloatCmp::eq(cp.p, expectedControlPoints[i], 1e-10))
        << "expectedControlPoints[i] " << expectedControlPoints[i] << " is " << cp.p;
    t.check(Dune::FloatCmp::eq(cp.w, expectedWeights[i], 1e-10))
        << "expectedWeights[i] " << expectedWeights[i] << " is " << cp.w;
    ++i;
  }

  return t;
}

auto testSurfaceHigherOrderDerivatives() {
  TestSuite t;

  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 3;
  const std::array<int, gridDim> order = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  //  const std::vector<std::vector<ControlPoint>> controlPoints
  //      = {{{.p = {0, 0,0}, .w = 1}, {.p = {0, 0.5,0}, .w = 1}},
  //         {{.p = {0.5, 0,0}, .w = 1}, {.p = {0.5, 0.5,0}, .w = 1}},
  //         {{.p = {1, 0,0}, .w = 1}, {.p = {1, 0.5,0}, .w = 1}}};

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0, 0.5, 0}, .w = 1}, {.p = {0, 1, 0}, .w = 1}},
         {{.p = {0.5, 0, 0}, .w = 1}, {.p = {0.5, 0.5, 0}, .w = 3}, {.p = {0.5, 1, 0}, .w = 1}},
         {{.p = {1, 0, 0}, .w = 1}, {.p = {1, 0.5, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 1}}};

  //  const std::vector<std::vector<ControlPoint>> controlPoints
  //      = {{{.p = {0, 0,0}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}},
  //         {{.p = {}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}},
  //         {{.p = {}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}}};
  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto additionalKnots    = std::vector<double>(1);
  additionalKnots[0]      = 0.5;
  patchData               = knotRefinement<2>(patchData, additionalKnots, 1);
  //  const double R       = 2.0;
  //  const double r       = 1.0;
  //  auto circle          = makeCircularArc(r);
  //  auto patchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 10.0);
  //  patchData       = degreeElevate(patchData, 0, 2);
  //  patchData       = degreeElevate(patchData, 1, 2);
  //  auto additionalKnots = std::vector<double>(1);
  ////  additionalKnots[0]   = 0.1;
  ////  patchData       = knotRefinement<2>(patchData, additionalKnots, 1);

  auto patchDataP = std::make_shared<Dune::IGA::NURBSPatchData<gridDim, dimworld>>(patchData);

  NURBSPatchGeometry<gridDim, dimworld> geo(patchDataP);
  double vol  = geo.volume();
  int samples = 2;
  std::vector<Dune::FieldVector<double, 3>> evaluatedPoints;
  std::vector<Dune::FieldMatrix<double, 2, 3>> evaluatedJacobians;
  std::vector<Dune::FieldMatrix<double, 3, 3>> evaluatedHessians;
  //  std::vector<double> evaluatedJacobians;
  for (int i = 0; i < samples + 1; ++i) {
    for (int j = 0; j < samples + 1; ++j) {
      const Dune::FieldVector<double, 2> u = {i / static_cast<double>(samples), j / static_cast<double>(samples)};
      std::cout << geo.global(u) << std::endl;
      evaluatedPoints.push_back(geo.global(u));
      const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
      evaluatedHessians.push_back(H);
      evaluatedJacobians.push_back(geo.jacobianTransposed(u));
      //    evaluatedJacobians.push_back(geo.impl()(u));
    }
  }

  for (int i = 1; i < 2; ++i) {
    std::cout << "degree elevate " << i << std::endl;
    //    auto patchDataN= degreeElevate(patchData,0,i);
    auto patchDataN = degreeElevate(patchData, 1, i);
    std::cout << "KnotsU:" << std::endl;
    for (auto kn : patchDataN.knotSpans[0])
      std::cout << kn << " ";
    std::cout << "KnotsV:" << std::endl;
    for (auto kn : patchDataN.knotSpans[1])
      std::cout << kn << " ";
    auto patchDataPN = std::make_shared<Dune::IGA::NURBSPatchData<gridDim, dimworld>>(patchDataN);
    for (auto cp : patchDataN.controlPoints.directGetAll())
      std::cout << cp.p << " w: " << cp.w << std::endl;
    NURBSPatchGeometry<gridDim, dimworld> geo2(patchDataPN);
    //    t.check(Dune::FloatCmp::eq(geo2.volume(),vol,1e-10))<<"vol "<<vol<<" is "<<geo2.volume();
    for (int j = 0, index = 0; j < samples + 1; ++j) {
      for (int k = 0; k < samples + 1; ++k) {
        const Dune::FieldVector<double, 2> u = {j / static_cast<double>(samples), k / static_cast<double>(samples)};
        t.check(Dune::FloatCmp::eq(geo2.global(u), evaluatedPoints[index], 1e-10))
            << "evaluatedPoints[i] " << evaluatedPoints[index] << " is " << geo2.global(u);
        const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
        t.check((geo2.jacobianTransposed(u)[0] - evaluatedJacobians[index][0]).two_norm() < 1e-8)
            << "evaluatedJacobians[i][0] \n"
            << evaluatedJacobians[index] << "\n is \n"
            << geo2.jacobianTransposed(u)
            << "\n norm:" << (geo2.jacobianTransposed(u)[0] - evaluatedJacobians[index][0]).two_norm();
        t.check((geo2.jacobianTransposed(u)[1] - evaluatedJacobians[index][1]).two_norm() < 1e-8)
            << "evaluatedJacobians[i][1] \n"
            << evaluatedJacobians[index] << "\n is \n"
            << geo2.jacobianTransposed(u)
            << "\n norm:" << (geo2.jacobianTransposed(u)[1] - evaluatedJacobians[index][1]).two_norm();
        t.check(Dune::FloatCmp::eq(H[0], evaluatedHessians[index][0], 1e-10))
            << "evaluatedJacobians[i][0] " << evaluatedHessians[index][0] << " is " << H;
        t.check(Dune::FloatCmp::eq(H[1], evaluatedHessians[index][1], 1e-10))
            << "evaluatedJacobians[i][1] " << evaluatedHessians[index] << " is " << H;
        t.check(Dune::FloatCmp::eq(H[2], evaluatedHessians[index][2], 1e-10))
            << "evaluatedJacobians[i][2] " << evaluatedHessians[index] << " is " << H;
        ++index;
      }
    }
  }

  return t;
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t;
  t.subTest(test3DGrid());
  t.subTest(testNURBSGridCurve());
  t.subTest(testPlate());
  testNurbsGridCylinder();
  t.subTest(testTorusGeometry());
  std::cout << " t.subTest(testNurbsBasis());" << std::endl;
  t.subTest(testNurbsBasis());
  t.subTest(testCurveHigherOrderDerivatives());
  t.subTest(testSurfaceHigherOrderDerivatives());

  gridCheck();
  t.subTest(testBsplineBasisFunctions());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
