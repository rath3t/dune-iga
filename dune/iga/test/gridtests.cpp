// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#define CHECK_RESERVEDVECTOR
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
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/geometrykernel/makecirculararc.hh>
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/parameterspace/concepts.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/subgrid/test/common.hh>

using namespace Dune;
using namespace Dune::IGA;
using namespace Dune::IGA;

auto checkUniqueEdges(const auto& gridView) {
  Dune::TestSuite t;

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
        bool inserted = edgeVertexPairSet.insert(pair).second;
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
  Dune::TestSuite t;

  constexpr int gridDimensionworld = GridView::dimensionworld;
  constexpr int gridDimension      = GridView::dimension;

  std::set<std::array<FieldVector<double, gridDimensionworld>, 1>, Compare<double, gridDimensionworld, 1>>
      elementVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    elementVertexPairSet.clear();
    for (auto vertexIdx : Dune::range(element.subEntities(gridDimension))) {
      auto vertex = element.template subEntity<gridDimension>(vertexIdx);

      std::array<FieldVector<double, gridDimensionworld>, 1> pair;
      for (auto c : Dune::range(vertex.geometry().corners()))
        pair[c] = vertex.geometry().corner(c);

      bool inserted = elementVertexPairSet.insert(pair).second;
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

auto thoroughGridCheck(auto& grid) {
  TestSuite t;
  constexpr int gridDimension = std::remove_cvref_t<decltype(grid)>::dimension;
  auto gvTest                 = [&](auto&& gv) {
    TestSuite tl;

    using GV = std::remove_cvref_t<decltype(gv)>;

    tl.subTest(checkUniqueEdges(gv));
    tl.subTest(checkUniqueSurfaces(gv));

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

  try {
    checkIntersectionIterator(grid);
    checkLeafIntersections(grid);
  } catch (const Dune::NotImplemented& e) {
    std::cout << e.what() << std::endl;
  }

  return t;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<2, 3, double>::ParameterSpace>
auto testNurbsGridCylinder() {
  ////////////////////////////////////////////////////////////////
  // First test
  // A B-Spline surface of dimWorld 3
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
  // const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::array<std::vector<double>, dim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}
  };

  using ControlPoint                                         = NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {        {.p = {0, 0, rad}, .w = 1},         {.p = {0, l, rad}, .w = 1}},
      {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
      // {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
      {        {.p = {rad, 0, 0}, .w = 1},         {.p = {rad, l, 0}, .w = 1}}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
  auto controlNet              = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::PatchGrid<2, 3, GridFamily> grid(patchData);

  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  auto patchGeo = grid.patchGeometry(0);

  // Since this is only a quarter of a cylinder no glueing should be present
  auto [glued, edgeIndex] = patchGeo.isConnectedAtBoundary(0);
  testSuite.check(!glued) << "The edge " << 0 << " should be glued!";
  testSuite.check(edgeIndex == -1);
  auto [glued2, edgeIndex2] = patchGeo.isConnectedAtBoundary(0);
  testSuite.check(!glued2) << "The edge " << 2 << " shouldn't be glued!";
  testSuite.check(edgeIndex == -1);

  Dune::RefinementIntervals refinementIntervals1(3);
  SubsamplingVTKWriter<decltype(grid.leafGridView())> vtkWriter(grid.leafGridView(), refinementIntervals1);
  // vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-Zylinder");

  testSuite.subTest(thoroughGridCheck(grid));
  return testSuite;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<2, 3, double>::ParameterSpace>
auto testHierarchicPatch() {
  TestSuite t("testHierarchicPatch", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  const double R       = 2.0;
  const double r       = 1.0;
  auto circle          = makeCircularArc(r);
  auto nurbsPatchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 360.0);
  nurbsPatchData       = Splines::degreeElevate(nurbsPatchData, 0, 1);
  nurbsPatchData       = Splines::degreeElevate(nurbsPatchData, 1, 2);
  auto additionalKnots = std::vector<double>(1);
  additionalKnots[0]   = 0.1;
  nurbsPatchData       = Splines::knotRefinement<2>(nurbsPatchData, additionalKnots, 1);
  for (int refDirection = 0; refDirection < 2; ++refDirection)
    nurbsPatchData = Splines::degreeElevate(nurbsPatchData, refDirection, 1);

  Dune::IGA::PatchGrid<2, 3, GridFamily> patch(nurbsPatchData);

  auto glueIndicator = [](int i) {
    switch (i) {
      case 0:
        return 1;
      case 1:
        return 0;
      case 2:
        return 3;
      case 3:
        return 2;
      default:
        DUNE_THROW(RangeError, "Invalid index");
    }
  };
  auto patchGeo = patch.patchGeometry(0);
  for (auto i : Dune::range(4)) {
    auto [glued, edgeIndex] = patchGeo.isConnectedAtBoundary(i);
    t.check(glued) << "The edge " << i << " should be glued!";
    t.check(edgeIndex == glueIndicator(i))
        << "The edge " << i << " should be glued to edge " << glueIndicator(i) << " but is glued to " << edgeIndex;
  }

  Dune::RefinementIntervals refinementIntervals1(3);
  SubsamplingVTKWriter<decltype(patch.leafGridView())> vtkWriter(patch.leafGridView(), refinementIntervals1);
  // vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-Torus");

  t.subTest(thoroughGridCheck(patch));

  return t;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<2, 3, double>::ParameterSpace>
auto testTorusGeometry() {
  const double R       = 2.0;
  const double r       = 1.0;
  auto circle          = makeCircularArc(r);
  auto nurbsPatchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 360.0);
  nurbsPatchData       = Splines::degreeElevate(nurbsPatchData, 0, 1);
  nurbsPatchData       = Splines::degreeElevate(nurbsPatchData, 1, 2);
  auto additionalKnots = std::vector<double>(1);
  additionalKnots[0]   = 0.1;
  nurbsPatchData       = Splines::knotRefinement<2>(nurbsPatchData, additionalKnots, 1);
  for (int refDirection = 0; refDirection < 2; ++refDirection)
    nurbsPatchData = Splines::degreeElevate(nurbsPatchData, refDirection, 1);

  IGA::PatchGrid<2, 3, GridFamily> grid(nurbsPatchData);
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
  grid.globalRefine(1);

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
  test.check(abs(area - referenceTorusSurfaceArea) < 5e-1)
      << "The integrated area of the torus surface is wrong! Expected: " << referenceTorusSurfaceArea
      << " Computed: " << area;

  double gaussBonnet = 0.0;
  // for (auto& ele : elements(gridView)) {
  //   const auto rule
  //       = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 *
  //       (*std::ranges::max_element(grid.patchData().degree)));
  //   for (auto& gp : rule) {
  //     const auto Kinc = ele.geometry().impl().gaussianCurvature(gp.position());
  //     const auto Kmax = 1 / (r * (R + r));
  //     const auto Kmin = -1 / (r * (R - r));
  //     test.check(Kinc < Kmax && Kinc > Kmin)
  //         << "The Gaussian curvature should be within bounds " << Kmin << " < " << Kinc << " < " << Kmax <<
  //         std::endl;
  //     gaussBonnet += Kinc * gp.weight() * ele.geometry().integrationElement(gp.position());
  //   }
  // }
  //
  // test.check(std::abs(gaussBonnet) < 1e-5)
  //     << "Gauss-Bonnet theorem dictates a vanishing integrated Gaussian curvature for the torus! It is " <<
  //     gaussBonnet;
  test.subTest(thoroughGridCheck(grid));
  return test;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<1, 3, double>::ParameterSpace>
auto testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  // Second test
  // A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const auto dim      = 1;
  const auto dimworld = 3;

  const std::array<int, dim> order                     = {3};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 3, 4, 4, 4, 5, 5, 5, 5}}};
  // const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints = {
      {.p = {1, 3, 4}, .w = 1},
      {.p = {2, 2, 2}, .w = 3},
      {.p = {3, 4, 5}, .w = 1},
      {.p = {5, 1, 7}, .w = 2},
      {.p = {4, 7, 2}, .w = 1},
      {.p = {8, 6, 2}, .w = 1},
      {.p = {2, 9, 9}, .w = 7},
      {.p = {1, 4, 3}, .w = 1},
      {.p = {1, 7, 1}, .w = 5}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto additionalKnots    = std::vector<double>(2);
  additionalKnots[0]      = 0.5;
  additionalKnots[1]      = 3.5;
  patchData               = Dune::IGA::Splines::knotRefinement<dim>(patchData, additionalKnots, 0);
  patchData               = Dune::IGA::Splines::degreeElevate(patchData, 0, 1);

  TestSuite t;
  IGA::PatchGrid<dim, dimworld, GridFamily> grid(patchData);

  auto patchGeo = grid.patchGeometry(0);

  auto [glued, edgeIndex] = patchGeo.isConnectedAtBoundary(0);
  t.check(glued) << "The edge " << 0 << " should be glued!";
  t.check(edgeIndex == -1);

  grid.globalRefine(3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  for (int eleIndex = 0;
       const auto& ele :
       elements(gridView)) // This test also exists in grid check, but it is more convenient to debug it here
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
  // vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-CurveNewFineResample-R");
  // vtkWriter.write("NURBSGridTest-CurveNewFineResample_knotRefine");
  t.subTest(thoroughGridCheck(grid));
  return t;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<2, 3, double>::ParameterSpace>
auto testNURBSGridSurface() {
  TestSuite t;
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2;
  const auto dimworld              = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };
  // const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
  //     = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
  using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
      {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
      {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}
  };

  std::array dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::PatchGrid<dim, dimworld, GridFamily> grid(patchData);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  auto patchGeo = grid.patchGeometry(0);
  for (auto i : Dune::range(4)) {
    auto [glued, edgeIndex] = patchGeo.isConnectedAtBoundary(i);
    t.check(!glued) << "The edge " << i << " shouldn't be glued!";
    t.check(edgeIndex == -1);
  }

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");

  t.subTest(thoroughGridCheck(grid));
  return t;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<3, 3, double>::ParameterSpace>
auto test3DGrid() {
  constexpr std::size_t dim        = 3;
  constexpr std::size_t dimworld   = 3;
  const std::array<int, dim> order = {2, 2, 2};
  // quarter cylindrical hyperSurface
  const double lx = 2;
  const double ly = 1;
  const double lz = 1;

  using ControlPoint = NURBSPatchData<dim, dimworld>::ControlPointType;
  NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  std::vector<std::vector<std::vector<ControlPoint>>> controlp;
  std::array<int, dim> dimSize = {3, 3, 3};
  for (int i = 0; i < dimSize[0]; ++i) {
    controlp.emplace_back();
    for (int j = 0; j < dimSize[1]; ++j) {
      controlp[i].emplace_back();
      for (int k = 0; k < dimSize[2]; ++k) {
        controlp[i][j].push_back({
            .p = {i * i * lx / (dimSize[0] - 1) + 1, 2 * i * j * k * ly / (dimSize[1] - 1) + (k + 1) + (j + 1),
                  k * k * lz / (dimSize[2] - 1)},
            .w = 1
        });
      }
    }
  }

  nurbsPatchData.controlPoints = MultiDimensionalNet(dimSize, controlp);
  nurbsPatchData.degree        = order;

  auto additionalKnots = std::vector<double>(2);
  additionalKnots[0]   = 0.1;
  additionalKnots[1]   = 0.3;
  // additionalKnots[2] = 0.6;
  // additionalKnots[1] = 3.5;
  // nurbsPatchData = knotRefinement<dim>(nurbsPatchData, additionalKnots, 2);
  // nurbsPatchData = degreeElevate(nurbsPatchData,0,1);
  // nurbsPatchData = degreeElevate(nurbsPatchData, 1, 2);
  // nurbsPatchData = degreeElevate(nurbsPatchData,2,1);
  IGA::PatchGrid<3, 3, GridFamily> grid(nurbsPatchData);
  // grid.globalRefine(1);
  // gridcheck(grid);
  // grid.globalRefineInDirection(0,1);
  // gridcheck(grid);
  // grid.globalRefineInDirection(1,2);
  // gridcheck(grid);
  // grid.globalRefineInDirection(2, 3);
  // gridcheck(grid);

  auto gridView = grid.leafGridView();
  TestSuite t;
  auto patchGeo = grid.patchGeometry(0);

  for (auto i : Dune::range(6)) {
    auto [glued, edgeIndex] = patchGeo.isConnectedAtBoundary(i);
    t.check(!glued) << "The surface " << i << " shouldn't be glued!";
    t.check(edgeIndex == -1);
  }

  t.subTest(checkUniqueEdges(gridView));
  t.subTest(checkUniqueSurfaces(gridView));

  const int subSampling = 10;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Solid2");

  t.subTest(thoroughGridCheck(grid));

  return t;
}

auto testCurveHigherOrderDerivatives() {
  constexpr auto gridDim               = 1;
  constexpr auto dimworld              = 3;
  const std::array<int, gridDim> order = {2};
  TestSuite t;
  // parameters

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 0.5, 1, 1, 1}}};
  // const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints = {
      {               .p = {0.086956521739130, -0.434782608695652, 0}, .w = 11.5},
      {.p = {0.200000000000000, 5.400000000000000, 0.200000000000000},    .w = 5},
      {                .p = {1.857142857142857, 0.142857142857143, 0},    .w = 7},
      {.p = {3.714285714285714, 0.285714285714286, 2.000000000000000},  .w = 3.5}
  };

  std::array<int, gridDim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet                  = NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  std::vector<Dune::FieldVector<double, 3>> expectedControlPoints = {
      {0.086956521739130, -0.434782608695652,                 0},
      {0.121212121212121,  1.333333333333333, 0.060606060606061},
      {0.320000000000000,  3.120000000000000, 0.120000000000000},
      {0.727272727272727,  3.727272727272727, 0.136363636363636},
      {1.538461538461539,  1.153846153846154, 0.038461538461538},
      {1.920000000000000,  0.506666666666667, 0.200000000000000},
      {2.476190476190476,  0.190476190476190, 0.666666666666667},
      {3.714285714285714,  0.285714285714286, 2.000000000000000}
  };

  std::vector<double> expectedWeights = {11.5, 8.25, 6.25, 5.5, 6.5, 6.25, 5.25, 3.5};

  GeometryKernel::NURBSPatch<gridDim, dimworld> geo(patchData);
  double vol  = geo.volume();
  int samples = 20;
  std::vector<Dune::FieldVector<double, 3>> evaluatedPoints;
  std::vector<Dune::FieldMatrix<double, 1, 3>> evaluatedJacobians;
  std::vector<Dune::FieldMatrix<double, 1, 3>> evaluatedHessians;

  for (int i = 0; i < samples + 1; ++i) {
    const Dune::FieldVector<double, 1> u = {i / static_cast<double>(samples)};
    evaluatedPoints.push_back(geo.global(u));
    const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
    evaluatedHessians.push_back(H);
    evaluatedJacobians.push_back(geo.jacobianTransposed(u));
  }

  for (int i = 1; i < 5; ++i) {
    auto patchDataN = Dune::IGA::Splines::degreeElevate(patchData, 0, i);
    GeometryKernel::NURBSPatch<gridDim, dimworld> geo2(patchDataN);
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

  patchData = Dune::IGA::Splines::degreeElevate(patchData, 0, 2);
  t.check(expectedControlPoints.size() == patchData.controlPoints.directGetAll().size())
      << "Size mismatch " << expectedControlPoints.size() << "is " << patchData.controlPoints.directGetAll().size();
  t.check(expectedWeights.size() == patchData.controlPoints.directGetAll().size())
      << "Size mismatch" << expectedWeights.size() << "is " << patchData.controlPoints.directGetAll().size();
  t.check(patchData.degree[0] == 4) << "Size mismatch";

  for (int i = 0; auto cp : patchData.controlPoints.directGetAll()) {
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

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  // const std::vector<std::vector<ControlPoint>> controlPoints
  //     = {{{.p = {0, 0,0}, .w = 1}, {.p = {0, 0.5,0}, .w = 1}},
  //        {{.p = {0.5, 0,0}, .w = 1}, {.p = {0.5, 0.5,0}, .w = 1}},
  //        {{.p = {1, 0,0}, .w = 1}, {.p = {1, 0.5,0}, .w = 1}}};

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {  {.p = {0, 0, 0}, .w = 1},   {.p = {0, 0.5, 0}, .w = 1},   {.p = {0, 1, 0}, .w = 1}},
      {{.p = {0.5, 0, 0}, .w = 1}, {.p = {0.5, 0.5, 0}, .w = 3}, {.p = {0.5, 1, 0}, .w = 1}},
      {  {.p = {1, 0, 0}, .w = 1},   {.p = {1, 0.5, 0}, .w = 1},   {.p = {1, 1, 0}, .w = 1}}
  };

  // const std::vector<std::vector<ControlPoint>> controlPoints
  //     = {{{.p = {0, 0,0}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}},
  //        {{.p = {}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}},
  //        {{.p = {}, .w = 1}, {.p = {}, .w = 1}, {.p = {}, .w = 1}}};
  std::array dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto additionalKnots    = std::vector<double>(1);
  additionalKnots[0]      = 0.5;
  patchData               = Dune::IGA::Splines::knotRefinement<2>(patchData, additionalKnots, 1);
  // const double R       = 2.0;
  // const double r       = 1.0;
  // auto circle          = makeCircularArc(r);
  // auto patchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 10.0);
  // patchData       = degreeElevate(patchData, 0, 2);
  // patchData       = degreeElevate(patchData, 1, 2);
  // auto additionalKnots = std::vector<double>(1);
  ////  additionalKnots[0]   = 0.1;
  ////  patchData       = knotRefinement<2>(patchData, additionalKnots, 1);

  GeometryKernel::NURBSPatch<gridDim, dimworld> geo(patchData);
  double vol  = geo.volume();
  int samples = 2;
  std::vector<Dune::FieldVector<double, 3>> evaluatedPoints;
  std::vector<Dune::FieldMatrix<double, 2, 3>> evaluatedJacobians;
  std::vector<Dune::FieldMatrix<double, 3, 3>> evaluatedHessians;
  // std::vector<double> evaluatedJacobians;
  for (int i = 0; i < samples + 1; ++i) {
    for (int j = 0; j < samples + 1; ++j) {
      const Dune::FieldVector<double, 2> u = {i / static_cast<double>(samples), j / static_cast<double>(samples)};
      evaluatedPoints.push_back(geo.global(u));
      const auto [pos, J, H] = geo.zeroFirstAndSecondDerivativeOfPosition(u);
      evaluatedHessians.push_back(H);
      evaluatedJacobians.push_back(geo.jacobianTransposed(u));
      // evaluatedJacobians.push_back(geo.impl()(u));
    }
  }

  for (int i = 1; i < 2; ++i) {
    auto patchDataN = Dune::IGA::Splines::degreeElevate(patchData, 1, i);
    GeometryKernel::NURBSPatch<gridDim, dimworld> geo2(patchDataN);
    // t.check(Dune::FloatCmp::eq(geo2.volume(),vol,1e-10))<<"vol "<<vol<<" is "<<geo2.volume();
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

auto testNURBSCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  // Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint                            = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<ControlPoint> controlPoints = {
      {.p = {1, 3, 4}, .w = 2},
      {.p = {2, 2, 2}, .w = 2},
      {.p = {3, 4, 5}, .w = 1},
      {.p = {5, 1, 7}, .w = 1},
      {.p = {4, 7, 2}, .w = 4},
      {.p = {8, 6, 2}, .w = 2},
      {.p = {2, 9, 9}, .w = 1},
      {.p = {1, 4, 3}, .w = 2},
      {.p = {1, 7, 1}, .w = 4}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  GeometryKernel::NURBSPatch<dim, dimworld> patch({knotSpans, controlNet, order});

  TestSuite testSuite;
  testSuite.check(patch.numberOfControlPoints().size() == 1)
      << "Number of controlpoints net size should be " << 1 << " but is " << patch.numberOfControlPoints().size();
  testSuite.check(patch.numberOfControlPoints()[0] == 9)
      << "Number of controlpoints in dir " << 0 << " should be " << 9 << " but is " << patch.numberOfControlPoints()[0];
  testSuite.check(patch.numberOfSpans().size() == 1)
      << "Number of net size should be " << 1 << " but is " << patch.numberOfSpans().size();
  testSuite.check(patch.numberOfSpans()[0] == 5)
      << "Number of elements in dir " << 0 << " should be " << 5 << " but is " << patch.numberOfSpans()[0];

  // Make grid
  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::PatchGrid<dim, dimworld, IdentityParameterSpace::PatchGridFamily> grid(patchData);

  auto leafGridView  = grid.leafGridView();
  auto& leafIndexSet = leafGridView.indexSet();

  for (const auto& ele : elements(leafGridView))
    testSuite.checkNoThrow([&]() { auto index = leafIndexSet.index(ele); });

  auto levelGridView  = grid.levelGridView(grid.maxLevel());
  auto& levelIndexSet = levelGridView.indexSet();

  for (const auto& ele : elements(levelGridView))
    testSuite.checkNoThrow([&]() { auto index = levelIndexSet.index(ele); });

  // Test MCMGM
  MultipleCodimMultipleGeomTypeMapper mapper(leafGridView, mcmgElementLayout());

  testSuite.check(mapper.size() == leafGridView.size(0));
  for (const auto& ele : elements(leafGridView)) {
    decltype(mapper)::Index idx;
    testSuite.check(mapper.contains(ele, idx));
  }

  // Test IDSets
  auto& idSet = grid.globalIdSet();
  for (const auto& ele : elements(leafGridView)) {
    testSuite.checkNoThrow([&]() { auto id = idSet.id(ele); });
  }

  return testSuite;
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

  const std::array<std::vector<double>, dim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
      {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
      {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}
  };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
  //

  // auto weightNet  = MultiDimensionalNet<dim, double>(dimsize, weight);
  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::GeometryKernel::NURBSPatch<dim, dimworld> patch({knotSpans, controlNet, order});

  TestSuite testSuite;
  testSuite.check(patch.numberOfControlPoints().size() == 2)
      << "Number of controlpoints net size should be " << 2 << " but is " << patch.numberOfControlPoints().size();
  testSuite.check(patch.numberOfControlPoints()[0] == 3)
      << "Number of controlpoints in dir " << 0 << " should be " << 3 << " but is " << patch.numberOfControlPoints()[0];
  testSuite.check(patch.numberOfControlPoints()[1] == 3)
      << "Number of controlpoints in dir " << 1 << " should be " << 3 << " but is " << patch.numberOfControlPoints()[1];
  testSuite.check(patch.numberOfSpans().size() == 2)
      << "Number of net size should be " << 2 << " but is " << patch.numberOfSpans().size();
  testSuite.check(patch.numberOfSpans()[0] == 1)
      << "Number of elements in dir " << 0 << " should be " << 1 << " but is " << patch.numberOfSpans()[0];
  testSuite.check(patch.numberOfSpans()[1] == 1)
      << "Number of elements in dir " << 1 << " should be " << 1 << " but is " << patch.numberOfSpans()[1];

  return testSuite;
}

template <template <int, int, typename> typename GridFamily>
requires IGA::Concept::ParameterSpace<typename GridFamily<2, 2, double>::ParameterSpace>
auto testPlate() {
  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 2;
  const std::array<int, gridDim> order = {2, 2};
  TestSuite t;

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {  {.p = {0, 0}, .w = 1},   {.p = {0.5, 0}, .w = 1},   {.p = {1, 0}, .w = 1}},
      {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 1}, {.p = {1, 0.5}, .w = 1}},
      {  {.p = {0, 1}, .w = 1},   {.p = {0.5, 1}, .w = 1},   {.p = {1, 1}, .w = 1}}
  };

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::PatchGrid<gridDim, dimworld, GridFamily>;

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto grid               = std::make_shared<Grid>(patchData);
  grid->globalRefine(1);

  auto gridView = grid->leafGridView();
  VTKWriter writer(gridView);
  writer.write("platetest");

  t.subTest(thoroughGridCheck(*grid));

  return t;
}

template <template <int, int, typename> typename GridFamily>
auto testGrids() {
  TestSuite t("testGrids");

  std::cout << "testHierarchicPatch" << std::endl;
  t.subTest(testHierarchicPatch<GridFamily>());

  std::cout << "testNurbsGridCylinder" << std::endl;
  t.subTest(testNurbsGridCylinder<GridFamily>());

  std::cout << "testPlate" << std::endl;
  t.subTest(testPlate<GridFamily>());

  std::cout << "testTorusGeometry" << std::endl;
  t.subTest(testTorusGeometry<GridFamily>());

  if constexpr (GridFamily<2, 2, double>::ParameterSpace::isAlwaysTrivial)
    t.subTest(test3DGrid<GridFamily>());

  // Not sure, this fails
  // t.subTest(testNURBSGridCurve<GridFamily>());

  std::cout << "testPlate" << std::endl;
  t.subTest(testPlate<GridFamily>());

  std::cout << "testTorusGeometry" << std::endl;
  t.subTest(testTorusGeometry<GridFamily>());

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  std::cout << "==================================" << std::endl;
  std::cout << "===============TEST DefaultParameterSpace===" << std::endl;
  std::cout << "==================================" << std::endl;
  t.subTest(testGrids<DefaultParameterSpace::PatchGridFamily>());

  std::cout << "==================================" << std::endl;
  std::cout << "===============TEST IdentityParameterSpace===" << std::endl;
  std::cout << "==================================" << std::endl;
  t.subTest(testGrids<IdentityParameterSpace::PatchGridFamily>());

  std::cout << "testNURBSCurve" << std::endl;
  t.subTest(testNURBSCurve());
  std::cout << "testNURBSSurface" << std::endl;
  t.subTest(testNURBSSurface());

  t.subTest(testCurveHigherOrderDerivatives());
  t.subTest(testSurfaceHigherOrderDerivatives());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
