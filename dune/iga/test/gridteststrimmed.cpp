// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif
#include "testhelper.hh"

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/tuplevector.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/concepts.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/subgrid/test/common.hh>

using namespace Dune;
using namespace Dune::IGA;

auto checkUniqueEdges(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  constexpr int gridDimension      = std::remove_cvref_t<decltype(gridView)>::dimension;
  if constexpr (gridDimension < 2)
    return t;
  else {
    std::set<std::array<FieldVector<double, gridDimensionworld>, 2>, Compare<double, gridDimensionworld, 2>>
        edgeVertexPairSet;
    for (int eleIndex = 0; auto&& element : elements(gridView)) {
      edgeVertexPairSet.clear();
      for (int edgeIndex = 0; edgeIndex < element.subEntities(gridDimension - 1); ++edgeIndex) {
        auto edge = element.template subEntity<gridDimension - 1>(edgeIndex);
        std::array<FieldVector<double, gridDimensionworld>, 2> pair;
        for (int c = 0; c < edge.geometry().corners(); ++c)
          pair[c] = edge.geometry().corner(c);
        const bool inserted = edgeVertexPairSet.insert(pair).second;
        t.require(inserted) << "Duplicate edge detected in Element " << eleIndex << " Edges: " << pair[0] << ", "
                            << pair[1];
      }
      ++eleIndex;
    }
    return t;
  }
}

template <typename GridView>
requires(GridView::dimension == 2)
auto checkUniqueVertices(const GridView& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = GridView::dimensionworld;
  constexpr int gridDimension      = GridView::dimension;

  std::set<std::array<FieldVector<double, gridDimensionworld>, 1>, Compare<double, gridDimensionworld, 1>>
      elementVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    elementVertexPairSet.clear();
    for (auto vertexIdx : Dune::range(element.subEntities(gridDimension))) {
      auto vertex = element.template subEntity<gridDimension>(vertexIdx);

      std::array<FieldVector<double, gridDimensionworld>, 1> pair;
      auto geo = vertex.geometry();
      for (auto c : Dune::range(geo.corners()))
        pair[c] = geo.corner(c);

      const bool inserted = elementVertexPairSet.insert(pair).second;
      t.require(inserted) << "Duplicate vertex detected in Element " << eleIndex << " Vertex: " << pair[0];
    }
    ++eleIndex;
  }
  return t;
}

auto checkUniqueSurfaces(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  constexpr int gridDimension      = std::remove_cvref_t<decltype(gridView)>::dimension;
  if constexpr (gridDimension < 3)
    return t;
  else {
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
}

template <typename G>
auto myGridCheck(G& grid) {
  TestSuite t;

  static_assert(G::dimension == 2);

  auto testGV = [&]<typename GV>(const GV& gv) {
    auto& indexSet = gv.indexSet();
    for (int eleIdx = 0; const auto& ele : elements(gv)) {
      const int numCorners  = ele.subEntities(2);
      const int numCorners2 = ele.geometry().corners();

      t.check(numCorners == numCorners2) << "Ele: " << eleIdx
                                         << " Corners from geometry not the same as subEntities(2)";

      // Check if conrers from corner and center from subentity are the same
      for (auto c : Dune::range(numCorners)) {
        auto corner  = ele.geometry().corner(c);
        auto corner2 = ele.template subEntity<2>(c).geometry().center();
        t.check(FloatCmp::eq(corner, corner2))
            << "Ele: " << eleIdx << "Position of corner(i) from the element is not the same as subentity(i)";
      }

      // Intersections
      int intersectionCount{0};
      for (const auto& intersection : intersections(gv, ele)) {
        auto geometry = intersection.geometry();
        ++intersectionCount;
      }
      const int numEdges = ele.subEntities(1);

      t.check(intersectionCount == numEdges) << "There should be as many edges as intersections";

      for (const auto i : Dune::range(numEdges)) {
        auto edge = ele.template subEntity<1>(i);
      }

      ++eleIdx;
      std::cout << std::endl;
    }
  };

  testGV(grid.levelGridView(grid.maxLevel()));
  testGV(grid.leafGridView());

  return t;
}

auto thoroughGridCheck(auto& grid) {
  TestSuite t("thoroughGridCheck");
  constexpr int gridDimension = std::remove_cvref_t<decltype(grid)>::dimension;

  auto gvTest = [&]<typename GV>(GV&& gv) {
    TestSuite tl;

    tl.subTest(checkUniqueEdges(gv));
    tl.subTest(checkUniqueSurfaces(gv));
    //
    if constexpr (GV::dimension == 2)
      tl.subTest(checkUniqueVertices(gv));

    auto extractGeo = std::views::transform([](const auto& ent) { return ent.geometry(); });
    for (auto&& elegeo : elements(gv) | extractGeo)
      checkJacobians(elegeo);

    for (auto&& vertGeo : vertices(gv) | extractGeo)
      checkJacobians(vertGeo);

    if constexpr (gridDimension > 1)
      for (auto&& edgegeo : edges(gv) | extractGeo)
        checkJacobians(edgegeo);

    if constexpr (gridDimension > 2)
      for (auto&& edgegeo : facets(gv) | extractGeo)
        checkJacobians(edgegeo);

    checkIterators(gv);
    checkEntityLifetime(gv);
    return tl;
  };

  for (int lvl = 0; lvl <= grid.maxLevel(); ++lvl) {
    t.subTest(gvTest(grid.levelGridView(lvl)));
  }
  t.subTest(gvTest(grid.leafGridView()));

  gridcheck(grid);
  t.subTest(myGridCheck(grid));

  try {
    // checkIntersectionIterator(grid);
    checkLeafIntersections(grid);
  } catch (const Dune::NotImplemented& e) {
    std::cout << e.what() << std::endl;
  }

  return t;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::Trimmer<typename GridFamily<2, 2, double>::Trimmer>
auto makeTestCases2d(TestSuite& t) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  const std::vector testCases{
      std::tuple<std::string, int, int, unsigned int>{"auxiliaryfiles/element_trim_xb.ibra", 0, 3, 100},
      {   "auxiliaryfiles/element_trim.ibra", 0, 3, 100},
      {    "auxiliaryfiles/trim_2edges.ibra", 0, 2, 100},
      {     "auxiliaryfiles/trim_multi.ibra", 0, 0, 100},
      {   "auxiliaryfiles/surface-hole.ibra", 1, 3, 150},
      // {"auxiliaryfiles/surface-hole-skew.ibra", 1, 3, 100},
      // {"auxiliaryfiles/surface-hole-square.ibra", 1, 3, 100}
  };

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();

  for (auto& [file_name, min, max, splitter] : testCases) {
    auto name = file_name.substr(file_name.find_last_of('/') + 1);

    gridFactory.insertJson(file_name, true, {min, min});
    gridFactory.insertTrimParameters(GridFactory::TrimParameterType{splitter});
    auto grid = gridFactory.createGrid();

    for (int i = min; i <= max; i++) {
      std::cout << "Testing now " << name << " (Refinement " << i << ")" << std::endl;
      if (i > min)
        grid->globalRefine(1);
      t.subTest(thoroughGridCheck(*grid));
    }
  }
}

template <template <int, int, typename> typename GridFamily>
auto testGrids() {
  TestSuite t("testTrimmedGrids", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  makeTestCases2d<GridFamily>(t);

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t("testTrimmedGrids", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  makeTestCases2d<DefaultTrim::PatchGridFamily>(t);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
