// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "testhelper.hh"
#include "testreferenceelement.hh"

#include <cfenv>
#include <thread>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/test/checkgeometry.hh>
#include <dune/iga/geometrykernel/nurbspatchtransform.hh>
#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGANEW;

using Grid = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;

struct PatchTrimDataResults
{
  double trimmingCurveTotalLength  = 0;
  double trimmingCurveCurvedLength = 0;
  int notAffineCounter             = 0;
  double straightLength            = 0;
};

template <Concept::Vector VectorType>
requires(VectorType::dimension == 2)
double cross(const VectorType& a, const VectorType& b) {
  return a[0] * b[1] - a[1] * b[0];
}

// template<typename Curve>
// auto windingNumberOfCurve(const Curve& curve, const Dune::FieldVector<double,2>&
// p=Dune::FieldVector<double,2>{0.5,0.5}) {
//
//   using HostGeo= Dune::MultiLinearGeometry<double,1,2>;
//
// auto degree = curve.degree();
//   const auto& rule = Dune::QuadratureRules<double, 1>::rule(
//       Dune::GeometryTypes::cube(1), degree[0]+20);
//   double windingNumber=0;
//   std::cout<<"Curve: "<<curve.corner(0)<<" "<<curve.corner(1)<<" Rule order: "<<rule.order()<<std::endl;
//   std::cout<<"Gps: ";
//   double volume=0;
//   for (auto& gp : rule) {
//     std::cout<<""<<gp.position()<<", "<<curve.global(gp.position())<<"; ";
//     auto fV = curve.global(gp.position())-p;
//     // std::cout<<"windingNumberOfCurve: "<<gp.position()<<" "<<fV<<" "<<std::endl;
//
//     auto j = curve.jacobianTransposed(gp.position());
//     auto intEle = curve.integrationElement(gp.position());
//     volume+=intEle*gp.weight();
//     fV=fV/fV.two_norm();
//     double pi =  std::numbers::pi_v<double>;
//     double tau = 2* pi;
//     double theta = atan2(fV[1],fV[0]);
//     windingNumber+= theta*intEle*gp.weight();
//   }
//   std::cout<<"windingNumber "<<windingNumber<<" Volume: "<<volume<<std::endl;
// return  windingNumber;
// }

// constexpr std::array unTrimmededgeVertexPairs{std::array<Dune::FieldVector<double, 2>,2>({{0, 0},{1,0}}),
// std::array<Dune::FieldVector<double, 2>,2>({{1, 0},{1,1}}),
//                                     std::array<Dune::FieldVector<double, 2>,2>({{1, 1},{0,1}}),
//                                     std::array<Dune::FieldVector<double, 2>,2>({{0, 1},{0,0}})};

struct GeomWrapper
{
  using CurveLocalViewType = Grid::Trimmer::TrimmingCurve;
  GeomWrapper(const CurveLocalViewType& curve)
      : variant(curve) {
  }
  GeomWrapper(const Dune::MultiLinearGeometry<double, 1, 2>& curve)
      : variant(curve) {
  }

  using Variant = std::variant<CurveLocalViewType, Dune::MultiLinearGeometry<double, 1, 2>>;

  auto domain() const {
    return std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<CurveLocalViewType, std::remove_cvref_t<decltype(var)>>)
            return var.domain();
          else
            return std::array<Utilities::Domain<double>, 1>{};
        },
        variant);
  }

  auto degree() const {
    return std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<CurveLocalViewType, std::remove_cvref_t<decltype(var)>>)
            return var.degree();
          else
            return std::array{1};
        },
        variant);
  }

  auto global(Dune::FieldVector<double, 1> u) const {
    return std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<CurveLocalViewType, std::remove_cvref_t<decltype(var)>>) {
            const std::array<Utilities::Domain<double>, 1> input{};
            // auto gpInSpan=Utilities::mapToRange(u,input,var.domain());
            // std::cout<<"TrimmingCurve: "<<gpInSpan<<" "<<var.uniqueKnotVector()[0].front()<<"
            // "<<var.uniqueKnotVector()[0].back()<<std::endl;

            // std::cout<<"Domain: "<<gpInSpan<<" Left: "<<var.domain()[0].left()<<" Right:
            // "<<var.domain()[0].right()<<std::endl; assert(var.domain()[0].checkInside(gpInSpan));
            return var.global(u);
          } else {
            // std::cout<<"Host: "<<u<<std::endl;
            return var.global(u);
          }
        },
        variant);
  }

  auto corner(int i) const {
    return std::visit([&](auto& var) { return var.corner(i); }, variant);
  }

  auto jacobianTransposed(Dune::FieldVector<double, 1> u) const {
    return std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<CurveLocalViewType, std::remove_cvref_t<decltype(var)>>) {
            return var.jacobianTransposed(u);
          } else {
            // std::cout<<"Host: "<<u<<std::endl;
            return var.jacobianTransposed(u);
          }
        },
        variant);
  }

  auto integrationElement(Dune::FieldVector<double, 1> u) const {
    return std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<CurveLocalViewType, std::remove_cvref_t<decltype(var)>>) {
            double volParamterSpace = 1.0;
            for (int i = 0; i < 1; ++i)
              volParamterSpace *= var.domain()[i].volume();
            return var.integrationElement(u) * volParamterSpace;
          } else {
            // std::cout<<"Host: "<<u<<std::endl;
            return var.integrationElement(u);
          }
        },
        variant);
  }
  std::variant<CurveLocalViewType, Dune::MultiLinearGeometry<double, 1, 2>> variant;
};

auto toDuneEdgeId(int idx) {
  switch (idx) {
    case 0:
      return 2;
    case 1:
      return 1;
    case 2:
      return 3;
    case 3:
      return 0;
    default:
      assert(false);
  }
}

template <typename Edge, typename HostGrid>
auto edgeGeometry(const Edge& edge, const HostGrid& grid, const typename Grid::Trimmer::ElementTrimData& eleTrimData) {
  if (edge.isTrimmed)
    return GeomWrapper(transform(edge.geometry.value()));

  std::array<Dune::FieldVector<double, 2>, 2> pair;

  auto hostEdgeGeometry = eleTrimData.hostEntity().subEntity<1>(toDuneEdgeId(edge.idx)).geometry();
  for (int i = 0; i < 2; ++i) {
    auto currentVertexCoord = edge.isTrimmed ? edge.geometry->corner(i) : hostEdgeGeometry.corner(i);
    pair[i]                 = currentVertexCoord;
  }
  // If edge.idx == 2 or 3 then switch
  if (edge.idx == 2 or edge.idx == 3)
    std::swap(pair[0], pair[1]);

  Dune::MultiLinearGeometry<double, 1, 2> edgeCube{
      Dune::GeometryTypes::line, std::vector{pair[0], pair[1]}
  };
  return GeomWrapper(edgeCube);
}

template <typename GridElement, typename GridView>
auto elementTrimDataObstacleCourse(const GridElement& ele, const Grid::Trimmer::ElementTrimData& eleTrimData,
                                   const GridView& gridView, const PatchTrimDataResults& resTrimPatch) {
  int eleIndex = gridView.indexSet().index(ele);
  Dune::TestSuite t("elementTrimDataObstacleCourse for element " + std::to_string(eleIndex) + std::string(" Flag: ") +
                        (eleTrimData.flag() == DefaultTrim::ElementTrimFlag::empty
                             ? "empty"
                             : (eleTrimData.flag() == DefaultTrim::ElementTrimFlag::full ? "full" : "trimmed")),
                    Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  if (eleTrimData.flag() != DefaultTrim::ElementTrimFlag::empty) {
    auto referenceElement = Grid::Trimmer::TrimmerTraits::ReferenceElementType(eleTrimData);
    t.subTest(checkReferenceElement(referenceElement, eleTrimData));
  }
  if (eleTrimData.flag() == DefaultTrim::ElementTrimFlag::empty) {
    t.check(eleTrimData.edges().empty()) << "Empty element should not have edges";
    t.check(eleTrimData.vertices().empty()) << "Empty element should not have vertices";
  }

  double lengthOfTrimmedCurves = 0;
  if (eleTrimData.flag() == DefaultTrim::ElementTrimFlag::trimmed) {
    int hostEdgesCounter = 0;
    // the edges should form a closed loop
    std::set<std::array<Dune::FieldVector<double, 2>, 2>, Compare<double, 2, 2>> edgeSet;
    std::array<Dune::FieldVector<double, 2>, 2> oldVertexCoordPair;
    double windingNumber = 0;

    t.check(eleTrimData.edges().size() == eleTrimData.vertices().size())
        << "The number of edges should be the number of vertices";
    for (int edgeIndex = 0, untrimmedEdgeIndex = 0; auto& edge : eleTrimData.edges()) {
      if (edge.isHost)
        ++hostEdgesCounter;
      auto edgeGeo = edgeGeometry(edge, gridView.grid(), eleTrimData);
      std::visit([&](auto& var) { t.check(Dune::checkGeometry(var)); }, edgeGeo.variant);
      if (edge.isTrimmed)
        t.check(edge.geometry.has_value()) << "Trimmed edge should have the geometry";
      else
        t.check(not edge.geometry.has_value() and edge.isHost) << "Host edge should not have the geometry";
      t.check(edge.geometry->corners() == 2) << "The edges should have 2 corners";
      std::array<Dune::FieldVector<double, 2>, 2> curveCorners;

      for (int i = 0; i < 2; ++i)
        curveCorners[i] = edgeGeo.corner(i);

      if (edgeIndex > 0) {
        std::stringstream ssE, ssS;
        ssE << "\n The end vertex (" << curveCorners[1] << ") of the new edge " << std::to_string(edgeIndex)
            << " is the same as end vertex (" << oldVertexCoordPair[1] << ") of the old edge "
            << std::to_string(edgeIndex - 1) << ".\n One of the two edges is in the wrong order.\n";
        ssS << "\n The start vertex (" << curveCorners[0] << ") of the new edge " << std::to_string(edgeIndex)
            << " is the same as start vertex (" << oldVertexCoordPair[0] << ") of the old edge "
            << std::to_string(edgeIndex - 1) << ".\n One of the two edges is in the wrong order.\n";

        t.check(Dune::FloatCmp::eq(curveCorners[0], oldVertexCoordPair[1]))
            << "The start vertex (" << curveCorners[0] << ") of edge " << edgeIndex
            << " should be the same as the end vertex (" << oldVertexCoordPair[1] << ") of edge " << edgeIndex - 1
            << "." << (Dune::FloatCmp::eq(curveCorners[1], oldVertexCoordPair[1]) ? ssE.str() : "")
            << (Dune::FloatCmp::eq(curveCorners[0], oldVertexCoordPair[0]) ? ssS.str() : "");
      }

      if (edgeIndex == eleTrimData.edges().size() - 1) {
        auto firstVertex = edgeGeometry(eleTrimData.edges()[0], gridView.grid(), eleTrimData).corner(0);
        t.check(Dune::FloatCmp::eq(curveCorners[1], firstVertex))
            << "The first vertex of the first edge should be the same as the end vertex of the last edge " << edgeIndex;
      }

      if (not edge.isHost) {
        const auto& rule =
            Dune::QuadratureRules<double, 1>::rule(Dune::GeometryTypes::cube(1), edgeGeo.degree()[0] + 20);
        for (auto& gp : rule)

          lengthOfTrimmedCurves += edgeGeo.integrationElement(gp.position()) * gp.weight();
      }

      oldVertexCoordPair[0] = curveCorners[0];
      oldVertexCoordPair[1] = curveCorners[1];

      // check that edges are unique
      const bool inserted = edgeSet.insert(curveCorners).second;
      t.require(inserted) << "Duplicate edge detected in Element " << eleIndex << " Edges: " << curveCorners[0] << ", "
                          << curveCorners[1];

      // compute winding number part

      // windingNumber+=windingNumberOfCurve( edgeGeo,ele.geometry().global({0.5,0.5}));

      ++edgeIndex;
      if (edge.isHost)
        ++untrimmedEdgeIndex;
    }

    std::set<Dune::FieldVector<double, 2>> vertexSet;
    int hostVertexCounter = 0;

    for (int vertexIndex = 0; auto& vertex : eleTrimData.vertices()) {
      if (vertex.isHost)
        ++hostVertexCounter;
      if (not vertex.isHost)
        t.check(vertex.geometry.has_value()) << "Trimmed edge should have the geometry";
      else
        t.check(not vertex.geometry.has_value()) << "Host edge should not have the geometry";
    }
    size_t totalEdges      = eleTrimData.edges().size();
    size_t totalVertices   = eleTrimData.vertices().size();
    size_t nonHostEdges    = totalEdges - hostEdgesCounter;
    size_t nonHostVertices = totalVertices - hostVertexCounter;

    // @todo This test is not true, i think, check trim_multi 1, 1 ele 3 (top right)
    // if (totalEdges > 2)
    //   t.check(nonHostEdges + 1 == nonHostVertices)
    //       << "Each non-host edge produces two non-host vertices."
    //       << "There are nonHostEdges " << nonHostEdges << " and nonHostVertices " << nonHostVertices;
    t.check(totalEdges > 1) << "There have to be at least 2 edges";

    // windingNumber /= 2*std::numbers::pi_v<double>;
    // t.check(Dune::FloatCmp::eq(windingNumber,1.0,1e-10) or Dune::FloatCmp::eq(windingNumber,1.0,1e-10))<<"The winding
    // number should be 1.0 or 0.0 but is "<<windingNumber;
  }

  return std::make_pair(t, lengthOfTrimmedCurves);
}

auto computeControlPointsFromMatrix(const Dune::DynamicMatrix<double>& cps, const std::vector<double>& weights) {
  int n_controlPoints = cps.N();
  using CP            = ControlPoint<Dune::FieldVector<double, 2>>;
  std::vector<CP> cpsLinear{};
  for (int i = 0; i < n_controlPoints; ++i) {
    Dune::FieldVector<double, 2> p;
    for (int k = 0; k < 2; ++k)
      p[k] = cps[i][k];

    double w = weights[i];

    cpsLinear.emplace_back(CP{.p = p, .w = w});
  }
  std::array dimSize{static_cast<int>(cpsLinear.size())};
  MultiDimensionalNet<1, CP> controlNet{dimSize, cpsLinear};
  return controlNet;
}

struct ExpectedValues
{
  double straightLength;
  double trimmingCurveCurvedLength;
  double trimmingCurveTotalLength;
  int firstLoopCurvesSize;
  int notAffineCounter;
};

template <typename ExecutionPolicy>
auto checkTrim(std::string filename, const ExpectedValues& expectedValues, ExecutionPolicy&& policy) {
  // Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  // Setup
  std::atomic<std::shared_ptr<Dune::TestSuite>> t =
      std::make_shared<Dune::TestSuite>(filename, Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  using GridFactory = Dune::GridFactory<Grid>;

  using Trimmer         = Grid::Trimmer;
  using ElementTrimData = Trimmer::ElementTrimData;

  auto range = Dune::range(4);
  for (auto refx : range) {
    // copy into vector sicne for_each does not work with iota_view and parallel execution
    std::vector<int> yR(range.begin(), range.end());

    std::for_each(policy, yR.begin(), yR.end(), [&](auto refy) {
      std::cout << "Thread: " << std::this_thread::get_id() << std::endl;
      auto gridFactory = GridFactory();
      auto brep        = readJson<2>(filename);
      std::cout << "Grid: " << filename << " Refinement level: (" << refx << ", " << refy << ")" << std::endl;
      gridFactory.insertJson(filename, true, {refx, refy});

      std::unique_ptr<Grid> igaGrid{};
      try {
        igaGrid = gridFactory.createGrid();
      } catch (Dune::GridError& e) {
        t.load(std::memory_order_relaxed)->check(false) << "Grid Creation failed ...\n";
        return;
      }

      auto patchTrimData           = igaGrid->trimmer().patchTrimData();
      const auto tensorCoordinates = GeometryKernel::NURBSPatch{gridFactory.patchData_}.uniqueKnotVector();
      Dune::YaspGrid gridH{tensorCoordinates};
      Dune::SubGrid<2, decltype(gridH)> grid(gridH);
      grid.createBegin();
      grid.insertLeaf();
      grid.createEnd();

      // const auto& patchTrimData = gridFactory.patchTrimData_.value();
      assert(patchTrimData.loops().size() == 1);
      const auto& firstLoop = patchTrimData.loops()[0];
      t.load(std::memory_order_relaxed)->check(firstLoop.curves().size() == expectedValues.firstLoopCurvesSize)
          << "firstLoop.curves().size() is " << firstLoop.curves().size() << " but should be "
          << expectedValues.firstLoopCurvesSize;
      PatchTrimDataResults resTrimPatch;

      std::vector<GeometryKernel::NURBSPatch<1, 2, double>> trimmingCurves;
      for (auto c : brep.curves2D) {
        auto knots = c.knots;
        knots[0].clear();
        knots[0].insert(knots[0].end(), c.knots[0].front());
        knots[0].insert(knots[0].end(), c.knots[0].begin(), c.knots[0].end());
        knots[0].insert(knots[0].end(), c.knots[0].back());
        auto cps = computeControlPointsFromMatrix(c.controlPoints, c.weights);
        NURBSPatchData<1, 2, double> curvePatchData(knots, cps, c.degree);
        GeometryKernel::NURBSPatch patch{curvePatchData};
        trimmingCurves.emplace_back(patch);
      }

      for (const auto& c : firstLoop.curves()) {
        const auto curveVolume = c.curveLength();
        resTrimPatch.trimmingCurveTotalLength += curveVolume;
        if (not c.affine()) {
          resTrimPatch.trimmingCurveCurvedLength += curveVolume;
          ++resTrimPatch.notAffineCounter;
        }
      }
      resTrimPatch.straightLength = resTrimPatch.trimmingCurveTotalLength - resTrimPatch.trimmingCurveCurvedLength;
      std::cout << std::setprecision(16);
      t.load(std::memory_order_relaxed)->check(resTrimPatch.notAffineCounter == expectedValues.notAffineCounter)
          << "notAffineCounter is " << resTrimPatch.notAffineCounter << " but should be "
          << expectedValues.notAffineCounter;
      // std::cout << resTrimPatch.trimmingCurveTotalLength << std::endl;
      // std::cout << resTrimPatch.trimmingCurveCurvedLength << std::endl;
      //  std::cout << resTrimPatch.straightLength << std::endl;

      double comparePrecision = 0.1;
      t.load(std::memory_order_relaxed)
              ->check(Dune::FloatCmp::eq(resTrimPatch.trimmingCurveTotalLength, expectedValues.trimmingCurveTotalLength,
                                         comparePrecision))
          << "trimmingCurveTotalLength is " << resTrimPatch.trimmingCurveTotalLength << " but should be "
          << expectedValues.trimmingCurveTotalLength;
      t.load(std::memory_order_relaxed)
              ->check(Dune::FloatCmp::eq(resTrimPatch.trimmingCurveCurvedLength,
                                         expectedValues.trimmingCurveCurvedLength, comparePrecision))
          << "trimmingCurveCurvedLength is " << resTrimPatch.trimmingCurveCurvedLength << " but should be "
          << expectedValues.trimmingCurveCurvedLength;
      t.load(std::memory_order_relaxed)
              ->check(Dune::FloatCmp::eq(resTrimPatch.straightLength, expectedValues.straightLength, comparePrecision))
          << "straightLength is " << resTrimPatch.straightLength << " but should be " << expectedValues.straightLength;

      auto gridView = grid.leafGridView();

      std::atomic<double> trimmedEdgeLengthsAccumulated{0};
      std::for_each(policy, gridView.template begin<0>(), gridView.template end<0>(), [&](const auto& ele) {
        ElementTrimData elementTrimData = DefaultTrim::TrimmerImpl<2, 2, double>::trimElement(ele, patchTrimData);
        auto [subTestEle, trimmedEdgeLength] =
            elementTrimDataObstacleCourse(ele, elementTrimData, gridView, resTrimPatch);
        trimmedEdgeLengthsAccumulated.fetch_add(trimmedEdgeLength, std::memory_order_relaxed);
        t.load(std::memory_order_relaxed)->subTest(subTestEle);
      });
      // for (const auto& ele : elements(gridView)) {
      //   ElementTrimData elementTrimData      = DefaultTrim::TrimmerImpl<2, 2, double>::trimElement(ele,
      //   patchTrimData); auto [subTestEle, trimmedEdgeLength] = elementTrimDataObstacleCourse(ele, elementTrimData,
      //   gridView, resTrimPatch); trimmedEdgeLengthsAccumulated += trimmedEdgeLength; t.subTest(subTestEle);
      // }

      t.load(std::memory_order_relaxed)
              ->check(Dune::FloatCmp::eq(trimmedEdgeLengthsAccumulated.load(std::memory_order_relaxed),
                                         resTrimPatch.trimmingCurveCurvedLength, comparePrecision))
          << "trimmedEdgeLengthsAccumulated is " << trimmedEdgeLengthsAccumulated << " but should be "
          << resTrimPatch.trimmingCurveCurvedLength;
    });
  }

  return *t.load(std::memory_order_relaxed);
}

#include <cfenv>

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  std::cout << std::thread::hardware_concurrency() << "\n";

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  /*
    t.subTest(checkTrim("auxiliaryfiles/element_trim.ibra",
                        ExpectedValues({.straightLength            = 2.63795970,
                                        .trimmingCurveCurvedLength = 0.97649438,
                                        .trimmingCurveTotalLength  = 3.61445407,
                                        .firstLoopCurvesSize       = 5,
                                        .notAffineCounter          = 1}),
                        std::execution::seq));
    t.subTest(checkTrim("auxiliaryfiles/element.ibra",
                        ExpectedValues({.straightLength            = 4,
                                        .trimmingCurveCurvedLength = 0,
                                        .trimmingCurveTotalLength  = 4,
                                        .firstLoopCurvesSize       = 4,
                                        .notAffineCounter          = 0}),
                        std::execution::seq));
                        */
  t.subTest(checkTrim("auxiliaryfiles/trim_2edges.ibra",
                      ExpectedValues({.straightLength            = 31.9472585,
                                      .trimmingCurveCurvedLength = 6.3246084,
                                      .trimmingCurveTotalLength  = 38.2718669,
                                      .firstLoopCurvesSize       = 6,
                                      .notAffineCounter          = 2}),
                      std::execution::seq)); /*

  t.subTest(checkTrim("auxiliaryfiles/trim_multi.ibra",
                      ExpectedValues({.straightLength            = 27.0662326,
                                      .trimmingCurveCurvedLength = 13.6037540,
                                      .trimmingCurveTotalLength  = 40.6699865,
                                      .firstLoopCurvesSize       = 5,
                                      .notAffineCounter          = 1}),
                      std::execution::seq));
  t.subTest(checkTrim("auxiliaryfiles/element_trim_xb.ibra",
                      ExpectedValues({.straightLength            = 2.79197791,
                                      .trimmingCurveCurvedLength = 1.57240470,
                                      .trimmingCurveTotalLength  = 4.36438261,
                                      .firstLoopCurvesSize       = 5,
                                      .notAffineCounter          = 1}),
                      std::execution::seq)); */

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
