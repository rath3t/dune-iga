// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include <iostream>

// #include "dune/iga/gridcapabilities.hh"
// #include "dune/iga/nurbsbasis.hh"
// #include "dune/iga/nurbsgrid.hh"
// #include "dune/iga/nurbspatch.hh"
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
// #include <dune/functions/functionspacebases/flatmultiindex.hh>
// #include <dune/functions/functionspacebases/powerbasis.hh>
// #include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/hierarchicpatch/hierachicpatchgrid.hh>
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
  Dune::TestSuite t;

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
  IGANEW::PatchGrid grid(patchData);
  grid.globalRefine(5);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;



  testSuite.subTest(checkUniqueEdges(gridView));
}

auto testHierarchicPatch() {
  TestSuite t;
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

  std::array<std::vector<double>, 2> uniqueKnotVector_;
  for (int i = 0; i < 2; ++i)  // create unique knotspan vectors
    std::ranges::unique_copy(nurbsPatchData.knotSpans[i], std::back_inserter(uniqueKnotVector_[i]),
                             [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

  Dune::TensorProductCoordinates<double, 2> coords(uniqueKnotVector_, std::array<int, 2>());
  std::cout << coords << std::endl;
  coords = coords.refine(std::bitset<2>(), std::bitset<2>(), 0, false);
  std::cout << coords << std::endl;

  // Dune::YaspGrid<2,Dune::TensorProductCoordinates<double,2>> g;
  // g.leafGridView().begin<0>()->geometry().impl().jacobian()
  // Dune::TensorProductCoordinates<double,2> coords();

  Dune::IGANEW::PatchGrid patch(nurbsPatchData);
  patch.globalRefine(3);
  auto gridView = patch.leafGridView();



  t.subTest(checkUniqueEdges(gridView));

  auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });
  for (auto&& elegeo : elements(gridView) | extractGeo)
    checkJacobians(elegeo);
  for (auto&& edgegeo : edges(gridView) | extractGeo)
    checkJacobians(edgegeo);

  for (auto&& vertGeo : vertices(gridView) | extractGeo)
    checkJacobians(vertGeo);


  checkEntityLifetime(gridView, gridView.size(0));
  gridcheck(patch);
  for (int i = 0; i < 3; ++i)
    checkIterators(patch.levelGridView(i));
    checkIterators(patch.leafGridView());
  //   Dune::TensorProductCoordinates<double,2> coords2(patch.uniqueKnotSpans,std::array<int,2>());
  //   std::cout<<coords2<<std::endl;
  // }
  checkEntityLifetime(gridView, gridView.size(0));



  return t;
}

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

  IGANEW::PatchGrid grid(nurbsPatchData);
  grid.globalRefine(1);
  // grid.globalRefineInDirection(1, 1);
  // grid.globalRefineInDirection(0, 2);
  // grid.globalDegreeElevate(2);

  auto gridView = grid.leafGridView();

  const int subSampling = 2;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-SurfaceRevolutionFLAT");

  TestSuite test;
  gridcheck(grid);

  // checkIntersectionIterator(grid,true);

  // for(int i=0; i<1; i++)
  // {
  //   std::cout << ">>> Refining grid and checking again..." << std::endl;
  //   grid.globalRefine( 1 );
  //   gridcheck(grid);
  //   checkIterators( grid.leafGridView() );
  //   checkIntersectionIterator(grid,true);
  //   checkTwists( grid.leafGridView(), NoMapTwist() );
  // }


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
  // for (auto& ele : elements(gridView)) {
  //   const auto rule
  //       = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * (*std::ranges::max_element(grid.patchData().degree)));
  //   for (auto& gp : rule) {
  //     const auto Kinc = ele.geometry().impl().gaussianCurvature(gp.position());
  //     const auto Kmax = 1 / (r * (R + r));
  //     const auto Kmin = -1 / (r * (R - r));
  //     test.check(Kinc < Kmax && Kinc > Kmin)
  //         << "The Gaussian curvature should be within bounds " << Kmin << " < " << Kinc << " < " << Kmax << std::endl;
  //     gaussBonnet += Kinc * gp.weight() * ele.geometry().integrationElement(gp.position());
  //   }
  // }
  //
  // test.check(std::abs(gaussBonnet) < 1e-5)
  //     << "Gauss-Bonnet theorem dictates a vanishing integrated Gaussian curvature for the torus! It is " << gaussBonnet;
  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);

  gridcheck(grid);
  test.subTest(checkUniqueEdges(gridView));
  return test;
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;
  t.subTest(testHierarchicPatch());
  // TestSuite t;
  // t.subTest(test3DGrid());
  // t.subTest(testNURBSGridCurve());
  // t.subTest(testPlate());
  // testNurbsGridCylinder();
  // t.subTest(testTorusGeometry());
  // std::cout << " t.subTest(testNurbsBasis());" << std::endl;
  // t.subTest(testNurbsBasis());
  // t.subTest(testCurveHigherOrderDerivatives());
  // t.subTest(testSurfaceHigherOrderDerivatives());
  //
  // gridCheck();
  // t.subTest(testBsplineBasisFunctions());
  //
  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
