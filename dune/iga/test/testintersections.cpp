// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/tuplevector.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGANEW;
using namespace Dune;

std::vector<GeometryType> geometryTypes(bool trimmed, int refLevel) {
  assert(refLevel == 0);
  if (not trimmed)
    return {GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::cube(1)};
  else
    return {GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::none(1),
            GeometryTypes::cube(1)};
}

std::vector<Dune::FieldVector<double, 2>> unitOuterNormals(bool trimmed, int refLevel) {
  assert(refLevel == 0);
  if (not trimmed)
    return {
        {-1,  0},
        { 1,  0},
        { 0, -1},
        { 0,  1}
    };
  // /todo this has to be verified somehow
  return {
      {             0,              -1},
      {             1,               0},
      {             0,               1},
      {-0.62562913456, -0.780120622711},
      {            -1,               0}
  };
}

template <typename PatchGrid>
auto testIntersections(auto& grid, bool trimmed, int refLevel) {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  using LeafGridView  = typename PatchGrid::LeafGridView;
  using LevelGridView = typename PatchGrid::LevelGridView;

  using EntitySeed   = typename PatchGrid::template Codim<0>::EntitySeed;
  using NormalVector = Dune::FieldVector<double, 2>;

  // NumIntersection, UnitOuterNormal, centerUnitOuterNormal, outerNormal, geometryType, insideEntitySeed
  using ResultTuple = std::tuple<int, std::vector<NormalVector>, std::vector<NormalVector>, std::vector<NormalVector>,
                                 std::vector<GeometryType>, std::vector<EntitySeed>>;

  ResultTuple resLevel{};
  ResultTuple resLeaf{};

  auto gvs = Dune::TupleVector<LevelGridView, LeafGridView>(grid.levelGridView(grid.maxLevel()), grid.leafGridView());
  Hybrid::forEach(gvs, [&]<typename GV>(const GV& gridView) {
    ResultTuple resTuple{};

    int globalIntersectionCount{0};
    for (const auto& ele : elements(gridView)) {
      int eleIntersectionCount{};
      for (const auto& intersection : intersections(gridView, ele)) {
        ++eleIntersectionCount;
        ++globalIntersectionCount;

        std::get<1>(resTuple).push_back(intersection.unitOuterNormal({0.5}));

        std::get<2>(resTuple).push_back(intersection.centerUnitOuterNormal());
        std::get<3>(resTuple).push_back(intersection.outerNormal({0.5}));
        std::get<4>(resTuple).push_back(intersection.type());
        std::get<5>(resTuple).push_back(intersection.inside().seed());
      }
      t.check(eleIntersectionCount == ele.subEntities(1))
          << "There should be " << ele.subEntities(1) << " intersections, but there are " << eleIntersectionCount;
    }

    std::get<0>(resTuple) = globalIntersectionCount;
    // Save appropriate tuple
    if constexpr (std::is_same_v<GV, LeafGridView>)
      resLeaf = resTuple;
    else
      resLevel = resTuple;
  });

  // Cehck level == leaf
  t.check(std::get<0>(resLevel) == std::get<0>(resLeaf));
  //
  for (const auto i : Dune::range(std::get<0>(resLevel))) {
    // Check outerNormals
    auto lvlUnitOuterNormal  = std::get<1>(resLevel)[i];
    auto leafUnitOuterNormal = std::get<1>(resLeaf)[i];

    auto lvlCenterOuterNormal  = std::get<2>(resLevel)[i];
    auto leafCenterOuterNormal = std::get<2>(resLeaf)[i];

    auto lvlOuterNormal  = std::get<3>(resLevel)[i];
    auto leafOuterNormal = std::get<3>(resLeaf)[i];

    t.check(FloatCmp::eq(lvlUnitOuterNormal, leafUnitOuterNormal));
    t.check(FloatCmp::eq(lvlCenterOuterNormal, leafCenterOuterNormal));
    t.check(FloatCmp::eq(lvlOuterNormal, leafOuterNormal));

    if (refLevel == 0) {
      const auto expectedOuterNormals = unitOuterNormals(trimmed, refLevel);
      t.check(FloatCmp::eq(lvlUnitOuterNormal, expectedOuterNormals[i], 1e-8));
      t.check(FloatCmp::eq(lvlCenterOuterNormal, expectedOuterNormals[i], 1e-8));
    }

    // Check unit length
    t.check(FloatCmp::eq(lvlUnitOuterNormal.two_norm(), 1.0, 1e-8));

    // Check linear dependence of unitOuterNormal und outerNormal
    auto testJ      = Dune::FieldMatrix<double, 2>{lvlUnitOuterNormal, lvlOuterNormal};
    double detTestJ = testJ.determinant();
    t.check(FloatCmp::lt(detTestJ, 1e-8)) << "Determinant should be 0, but is " << detTestJ;

    // Check geometryTypes
    auto levelGeoType = std::get<4>(resLevel)[i];
    auto leafGeoType  = std::get<4>(resLeaf)[i];

    t.check(levelGeoType == leafGeoType);
    if (refLevel == 0) {
      auto expectedGeometryTypes = geometryTypes(trimmed, refLevel);
      t.check(levelGeoType == expectedGeometryTypes[i]);
    }

    // Seed
    auto levelEntityInside = grid.entity(std::get<5>(resLevel)[i]);
    auto leafEntityInside  = grid.entity(std::get<5>(resLeaf)[i]);

    t.check(levelEntityInside == leafEntityInside);
  }
  // Also for refLevel == 0
  if (refLevel == 0) {
    t.check(grid.entity(std::get<5>(resLevel)[0]) == grid.entity(std::get<5>(resLevel)[1]));
    t.check(grid.entity(std::get<5>(resLevel)[1]) == grid.entity(std::get<5>(resLevel)[2]));
    t.check(grid.entity(std::get<5>(resLevel)[2]) == grid.entity(std::get<5>(resLevel)[3]));
    t.check(grid.entity(std::get<5>(resLevel)[3]) == grid.entity(std::get<5>(resLevel)[0]));
  }

  return t;
}

auto testInsideOutside(auto& grid) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  Dune::TestSuite t("Inside Outisde Test", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  auto gridView  = grid.leafGridView();
  auto& indexSet = gridView.indexSet();

  for (const auto& element : elements(gridView)) {
    auto eleGeometry = element.geometry();
    auto elementIdx  = indexSet.index(element);

    for (int iIdx = 0; const auto& intersection : intersections(gridView, element)) {
      auto insideElement = intersection.inside();
      t.check(insideElement == element);

      if (intersection.boundary()) {
        ++iIdx;
        continue;
      }
      auto outsideElement   = intersection.outside();
      auto ousideElementIdx = indexSet.index(outsideElement);

      auto geometryInInside  = intersection.geometryInInside();
      auto geometryInOutside = intersection.geometryInOutside();

      auto centerInsideLocal  = geometryInInside.center();
      auto centerOutsideLocal = geometryInOutside.center();

      auto centerInside  = eleGeometry.global(centerInsideLocal);
      auto centerOutside = outsideElement.geometry().global(centerOutsideLocal);

      auto centerGlobal = intersection.geometry().center();

      t.check(FloatCmp::eq(centerInside, centerOutside, 1e-8));
      t.check(FloatCmp::eq(centerGlobal, centerInside, 1e-8));

      int indexInInside  = intersection.indexInInside();
      int indexInOutside = intersection.indexInOutside();

      auto it1 = std::ranges::find_if(
          gridView.ibegin(outsideElement), gridView.iend(outsideElement), [&](const auto& outsideI) {
            if (not outsideI.boundary() and indexSet.index(outsideI.outside()) == elementIdx)
              return outsideI.indexInOutside() == indexInInside;
            return false;
          });

      auto it2 = std::ranges::find_if(
          gridView.ibegin(outsideElement), gridView.iend(outsideElement), [&](const auto& outsideI) {
            if (not outsideI.boundary() and indexSet.index(outsideI.outside()) == elementIdx)
              return outsideI.indexInInside() == indexInOutside;
            return false;
          });

      t.require(it2 != gridView.iend(outsideElement));
      t.require(it1 != gridView.iend(outsideElement));

      t.check(it1->neighbor());
      t.check(it2->neighbor());

      t.check(it1->indexInInside() == it2->indexInInside());

      auto intersectionOnTheOutside = *it1;
      t.check(indexSet.index(intersectionOnTheOutside.inside()) == ousideElementIdx);
      t.check(indexSet.index(intersectionOnTheOutside.outside()) == elementIdx);

      auto insideCenterIntersectionOutsideLocal  = intersectionOnTheOutside.geometryInInside().center();
      auto outsideCenterIntersectionOutsideLocal = intersectionOnTheOutside.geometryInOutside().center();

      auto insideCenterIntersectionOutside  = outsideElement.geometry().global(insideCenterIntersectionOutsideLocal);
      auto outsideCenterIntersectionOutside = eleGeometry.global(outsideCenterIntersectionOutsideLocal);
      auto centerIntersectionOnTheOutside   = intersectionOnTheOutside.geometry().center();

      t.check(FloatCmp::eq(centerInside, centerIntersectionOnTheOutside, 1e-8));
      t.check(FloatCmp::eq(insideCenterIntersectionOutside, outsideCenterIntersectionOutside, 1e-8));

      ++iIdx;
    }
  }
  return t;
}

auto runIntersectionTests(Dune::TestSuite& t, const std::string& fileName, bool trimmed, int refLevel) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson(fileName, trimmed, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  t.subTest(testIntersections<PatchGrid>(*grid, trimmed, refLevel));
  t.subTest(testInsideOutside(*grid));
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  /******************
   * These are unit tests for the intersection implementation, we are mostly testing geometry, and not correct indices,
   * etc, as this is tested by the gridtests
   ******************/

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  runIntersectionTests(t, "auxiliaryfiles/element_trim.ibra", true, 0);
  runIntersectionTests(t, "auxiliaryfiles/element_trim.ibra", true, 1);
  runIntersectionTests(t, "auxiliaryfiles/element_trim.ibra", true, 2);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}