// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include <iostream>

#include "dune/iga/gridcapabilities.hh"
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

  std::array<std::vector<double>, 2>uniqueKnotVector_;
  for (int i = 0; i < 2; ++i)  // create unique knotspan vectors
    std::ranges::unique_copy(nurbsPatchData.knotSpans[i], std::back_inserter(uniqueKnotVector_[i]),
                             [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

  Dune::TensorProductCoordinates<double,2> coords(uniqueKnotVector_,std::array<int,2>());
  std::cout<<coords<<std::endl;
  coords=coords.refine(std::bitset<2>(),std::bitset<2>(),0,false);
  std::cout<<coords<<std::endl;

  // Dune::YaspGrid<2,Dune::TensorProductCoordinates<double,2>> g;
  // g.leafGridView().begin<0>()->geometry().impl().jacobian()
  // Dune::TensorProductCoordinates<double,2> coords();

  Dune::IGANEW::PatchGrid patch(nurbsPatchData);
  patch.globalRefine(3);
  auto gridView= patch.leafGridView();
  Dune::GeometryChecker<decltype(patch)> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  for (int i = 0; i < 3; ++i)
  geometryChecker.checkGeometry(patch.levelGridView(i));

  t.subTest(checkUniqueEdges(gridView));

  gridcheck(patch);
  // for (int i = 0; i < 3; ++i) {
  //   Dune::TensorProductCoordinates<double,2> coords2(patch.uniqueKnotSpans,std::array<int,2>());
  //   std::cout<<coords2<<std::endl;
  // }

return t;
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
