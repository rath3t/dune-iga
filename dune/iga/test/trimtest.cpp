// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/grid/test/gridcheck.hh>

#include <dune/iga/geometrykernel/nurbspatchgeometry.hh>

#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/geometrykernel/makesurfaceofrevolution.hh>
#include <dune/iga/hierarchicpatch/gridcapabilities.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <clipper2/clipper.h>
#include <unordered_set>
#include <dune/iga/trimmer/defaulttrimmer/idset.hh>


using namespace Dune::IGANEW;


auto diagonalTrimmingCurve(double offset) {
  const std::array<std::vector<double>, 1> knotSpansCurve = {{
      {0, 0, 1, 1},
  }};
  using ControlPoint                                      = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointType;

  const std::vector<ControlPoint> controlPointsCurve = {{{.p = {-offset, offset}, .w = 1}, {.p = {1-offset, 1+offset}, .w = 1}}};
  const std::array orderCurve                        = {1};
  auto controlNetCurve = Dune::IGANEW::NURBSPatchData<1, 2>::ControlPointNetType(controlPointsCurve);
  Dune::IGANEW::NURBSPatchData<1, 2> patchDataCurve;
  patchDataCurve.knotSpans     = knotSpansCurve;
  patchDataCurve.degree        = orderCurve;
  patchDataCurve.controlPoints = controlNetCurve;
  return Dune::IGANEW::GeometryKernel::NURBSPatch(patchDataCurve);
}

using UntrimmedParameterSpaceGrid
            = YaspGrid<2,
                       TensorProductCoordinates<double, 2>>;
  using YaspGridIdType = typename UntrimmedParameterSpaceGrid::Traits::GlobalIdSet::IdType;
  struct IndexType {
    YaspGridIdType id;
    auto touint()const {return id.touint();}
  };

using IndexVariant = std::variant<YaspGridIdType, IndexType>;


bool operator<(const IndexVariant& lhs, const IndexVariant& rhs) {
    return std::visit([](auto&& arg1, auto&& arg2) {
      if constexpr(std::is_same_v<std::remove_cvref_t<decltype(arg1)>,YaspGridIdType> and std::is_same_v<std::remove_cvref_t<decltype(arg2)>,YaspGridIdType>)
      {
      return arg1 < arg2;
      }
        else if constexpr(std::is_same_v<std::remove_cvref_t<decltype(arg1)>,IndexType> and std::is_same_v<std::remove_cvref_t<decltype(arg2)>,IndexType>)
        {
            return arg1.id < arg2.id; }
            else if constexpr(std::is_same_v<std::remove_cvref_t<decltype(arg1)>,IndexType> and std::is_same_v<std::remove_cvref_t<decltype(arg2)>,YaspGridIdType>)
        {

            return true;
            }
      else
    {
        return false;
        }
    }, lhs, rhs);
}


auto testFactoryWithPlateWithTriangularTrim2D() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int gridDim   = 2;
  constexpr auto dimworld = 2;
  using Grid              = Dune::IGANEW::PatchGrid<gridDim, dimworld, DefaultTrim::Trimmer>;
  const std::array order  = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double Lx = 2;
  const double Ly = 3;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {Lx / 2, 0}, .w = 1}, {.p = {Lx, 0}, .w = 1}},
         {{.p = {0, Ly / 2}, .w = 1}, {.p = {Lx / 2, Ly / 2}, .w = 1}, {.p = {Lx, Ly / 2}, .w = 1}},
         {{.p = {0, Ly}, .w = 1}, {.p = {Lx / 2, Ly}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGANEW::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  const auto trimCurve = diagonalTrimmingCurve(0.5);
  using PatchTrimData= typename Grid::TrimmerType::PatchTrimData;
  PatchTrimData patchTrimData;
  patchTrimData.insertTrimCurve(trimCurve);

  Dune::IGANEW::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  Grid grid(patchData,patchTrimData);
  grid.globalRefine(1);
  auto& parameterGrid = grid.parameterSpaceGrid();
  auto gv =parameterGrid.leafGridView();

  auto& globalIdSet = parameterGrid.globalIdSet();
  auto& globalIdSetNURBS = grid.globalIdSet();

  

//   using namespace Clipper2Lib;
//
//
//   PathsD edgesPath(1);
//
//   struct ElementIndices {
//     std::array<int,2> indices;
//   };
//   std::map<int,ElementIndices> interSectionIndicesOfElement;
//
//   std::unordered_set<int> usedIndicesByUnTrimmedParameterGrid;
//
//
//
//     std::map<IndexType, std::string> myMap;
//     std::map<IndexVariant, IndexType> myIndexMapping;
//
// PathsD clipCurve(1);
//  for(auto v: Utilities::linspace(trimCurve.domain()[0],10))
// { std::cout<<v<<std::endl;
// auto point = trimCurve.global(v);
// clipCurve[0].emplace_back(point[0], point[1]);
// }
// std::cout<<"clipCurve[0]"<<clipCurve[0].size()<<std::endl;
// for(int i =0;i<clipCurve[0].size();++i)
// {
// std::cout<<clipCurve[0][i]<<std::endl;
// }
//   PathsD edges(1);
//   edges[0].resize(5);
//   auto trimCubeReferenceElement=[&](const auto& element, const auto& clip) {
//
//     auto geo  = element.geometry();
//     auto pos0 = geo.corner(0);
//     auto pos1 = geo.corner(1);
//     auto pos2 = geo.corner(3);  // see dune book page 127 Figure 5.12
//     auto pos3 = geo.corner(2);
//     std::cout<<"C0"<<pos0<<std::endl;
//     std::cout<<"C1"<<pos1<<std::endl;
//     std::cout<<"C2"<<pos2<<std::endl;
//     std::cout<<"C3"<<pos3<<std::endl;
//     edges.front()[0]={pos0[0], pos0[1]};
//     edges.front()[1]={pos1[0], pos1[1]};
//     edges.front()[2]={pos2[0], pos2[1]}; //swap order of 2 and 3
//     edges.front()[3]={pos3[0], pos3[1]};
//       edges.front()[4]={pos0[0], pos0[1]};
//     Clipper2Lib::ClipperD clipper(5);
//     Clipper2Lib::PathsD clippedEdges;
//     clipper.AddSubject(edges);
//     clipper.AddClip(clip);
//     clipper.Execute(ClipType::Intersection, FillRule::NonZero, clippedEdges);
//     return clippedEdges;
//   };
//
//     // Function to get or generate a stable index for the own ids and host grid ids
//     auto getStableIndex = [&myIndexMapping](IndexVariant thirdPartyIndex) -> IndexType {
//         auto it = myIndexMapping.find(thirdPartyIndex);
//         if (it == myIndexMapping.end()) {
//             // If not found, generate a new index
//             IndexType newIndex{myIndexMapping.size()};
//             myIndexMapping[thirdPartyIndex] = newIndex;
//             return myIndexMapping[thirdPartyIndex];
//         }
//             // If found, return the existing index
//             return it->second;
//     };
//
//
//
//   std::map<Grid::Traits::GlobalIdSet::IdType,
//                                           std::vector<IndexType>> globalEdgesIdOfElementsMap;
//   using ParameterSpaceEdge = typename Grid::ParameterSpaceGrid::template Codim<1>::Entity;
//   using Edge = typename Grid::ParameterSpaceGrid::template Codim<1>::Entity;
//   //std::map< IndexType,Edge> globalEdgeMap;
//   std::map< IndexType,std::string> globalEdgeMap;
//
//   IndexType myTrimmedEdgeIndex;
//   myTrimmedEdgeIndex.id=0;
//
//   for (const auto& ele: elements(gv)) {
//   auto trimedEdges= trimCubeReferenceElement(ele,clipCurve);
//   std::cout<<"T"<<trimedEdges.size()<<std::endl;
//   for(auto v: trimedEdges)
//   std::cout<<v<<std::endl;
//   std::cout<<"Element"<<std::endl;
//     auto elementId = globalIdSet.id(ele);
//     auto& elementEdgeIndices= globalEdgesIdOfElementsMap[elementId];
//
//
//     for(const auto& intersection: intersections(gv,ele)){
//     auto localEdgeIndex= intersection.indexInInside();
//       std::cout<<"Edge "<<localEdgeIndex<<std::endl;
//
//       auto edge= intersection.inside().subEntity<1>(localEdgeIndex);
//       if (true /* untrimmed */) {
//         auto edgeId = globalIdSet.subId(ele,localEdgeIndex,1);
//         IndexType id =getStableIndex(edgeId);
//         elementEdgeIndices.push_back(id);
//        // globalEdgeMap.insert({id,Edge(edge)});
//         globalEdgeMap.insert({id,"Untrimmed edge " + std::to_string(id.id.touint())});
//       }
//       else /* trimmed */
//         {
//       if(true/* edge is trimmed but has a rest on the original intersection*/)
//       //if edge is trimmed but contains a part of the original intersection,
//       // then this part of the edge that belongs to the original intersection will be inserted using the id of the host
//         // the rest of the edge will be created using a new index
//       {
//         auto edgeId = globalIdSet.subId(ele,localEdgeIndex,1);
//         IndexType id =getStableIndex(edgeId);
//         elementEdgeIndices.push_back(id);
//         globalEdgeMap.insert({id,"Part of trimmed edge is part of untrimmed edge" + std::to_string(id.id.touint())});
//         IndexType id2 = getStableIndex(myTrimmedEdgeIndex);
//         globalEdgeMap.insert({id2,"Part of Trimmed edge with ID, that is not completely new but this part is"+std::to_string(myTrimmedEdgeIndex.id.touint())});
//         elementEdgeIndices.push_back(id);
//
//         ++myTrimmedEdgeIndex.id;
//
//       }
//           else if (false /* the edge is completly gone?*/) {
//             //
//           }
//             //
//               IndexType id =getStableIndex(myTrimmedEdgeIndex);
//         elementEdgeIndices.push_back(id);
//               //  globalEdgeMap.insert({id,Edge(edge)});
//
//                               globalEdgeMap.insert({getStableIndex(myTrimmedEdgeIndex),"Trimmed edge with ID "+std::to_string(myTrimmedEdgeIndex.id.touint())});
//  ++myTrimmedEdgeIndex.id;
//       }
//     //  auto edgeGeo = intersection.geometry();
//     //  auto firstCorner = edgeGeo.corner(0);
//      // auto secondCorner = edgeGeo.corner(1);
//      // edges[0][indexMapper(i)](edgeGeo.corner(0), pos0[1]);
//
//     }
//
//   }
//
//
//              for (const auto& [key, value] : globalEdgeMap)
//                     std::cout << "Index value: " << key.id<< " Edge: " << value << std::endl;

  return t;
}


#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  t.subTest(testFactoryWithPlateWithTriangularTrim2D());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
