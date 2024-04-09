// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "dune/iga/trimmer/defaulttrimmer/trimmer.hh"
#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
namespace Dune::IGANEW::DefaultTrim {
// template <int dim, int dimworld, typename ScalarType>
// void TrimmerImpl<dim, dimworld, ScalarType>::createLevel(GridImp& grid, int lvl) {
//   using IdType                 = typename GridImp::Traits::GlobalIdSet::IdType;
//   using EdgeHostType           = typename UntrimmedParameterSpaceGrid::template Codim<1>::Entity;
//   using EdgeGridType           = typename GridImp::template Codim<1>::Entity;
//   using EleGridType            = typename GridImp::template Codim<0>::Entity;
//   using EdgeParameterSpaceType = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
//   using EleParameterSpaceType  = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;
//
//   entityContainer_.push_back();
//
//   auto& entityContainer           = entityContainer_;
//   auto& globalIdSet               = grid.globalIdSet_;
//   auto gv                         = parameterSpaceGrid_->levelGridView(lvl);
//   auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
//
//   for (const auto& ele : elements(gv)) {
//     // auto trimedEdges= trimCubeReferenceElement(ele,clipCurve);
//     // std::cout<<"T"<<trimedEdges.size()<<std::endl;
//     // for(auto v: trimedEdges)
//     //   std::cout<<v<<std::endl;
//     std::cout << "Element" << std::endl;
//     // auto elementId = globalIdSet->getStableIndex();
//     IdType elementId = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = globalIdSetParameterSpace.id(ele)};
//
//     auto& elementEdgeIndices   = entityContainer.globalEdgesIdOfElementsMap_[elementId];
//     auto& elementVertexIndices = entityContainer.globalVerticesIdOfElementsMap[elementId];
//     std::set<IdType> elementVertexSet;
//     elementVertexSet.clear();
//     for (const auto& intersection : intersections(gv, ele)) {
//       auto localEdgeIndex = intersection.indexInInside();
//       std::cout << "Edge " << localEdgeIndex << std::endl;
//
//       auto edge    = intersection.inside().template subEntity<1>(localEdgeIndex);
//       auto vertex0 = intersection.inside().template subEntity<2>(0);
//       auto vertex1 = intersection.inside().template subEntity<2>(1);
//       auto refEle  = Dune::ReferenceElements<ctype, mydimension>::cube();
//
//       if (true /* untrimmed */) {
//         auto eleParameterSpace = EleParameterSpaceType(&grid, ele, elementId);
//         // auto realEntity = typename EleGridType::Implementation(&grid,eleParameterSpace);
//         entityContainer.globalElementMap.insert({elementId, eleParameterSpace});
//         auto edgeId = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
//         // IndexType id =globalIdSet->g(edgeId);
//         IdType id = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = edgeId};
//         elementEdgeIndices.push_back(id);
//         auto localVertexid0 = refEle.subEntity(localEdgeIndex, 1, 0, 2);
//         auto localVertexid1 = refEle.subEntity(localEdgeIndex, 1, 1, 2);
//         auto vertexId0      = globalIdSetParameterSpace.subId(ele, localVertexid0, 2);
//         auto vertexId1      = globalIdSetParameterSpace.subId(ele, localVertexid1, 2);
//
//         IdType idVert0 = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = vertexId0};
//         IdType idVert1 = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = vertexId1};
//
//         elementVertexSet.insert(idVert0);
//         elementVertexSet.insert(idVert1);
//         // globalEdgeMap.insert({id,Edge(edge)});
//         auto paraEnt = EdgeParameterSpaceType(&grid, edge, id);
//         // auto realEdge = typename EdgeGridType::Implementation(&grid,paraEnt);
//         entityContainer.globalEdgeMap.insert({id, paraEnt});
//       } else /* trimmed */
//       {
//         if (true /* edge is trimmed but has a rest on the original intersection*/)
//         // if edge is trimmed but contains a part of the original intersection,
//         //  then this part of the edge that belongs to the original intersection will be inserted using the id of
//         the
//         //  host the rest of the edge will be created using a new index
//         {
//           // auto edgeId = globalIdSet.subId(ele,localEdgeIndex,1);
//           // IndexType id =getStableIndex(edgeId);
//           // elementEdgeIndices.push_back(id);
//           // globalEdgeMap.insert({id,"Part of trimmed edge is part of untrimmed edge" +
//           // std::to_string(id.id.touint())}); IndexType id2 = getStableIndex(myTrimmedEdgeIndex);
//           // globalEdgeMap.insert({id2,"Part of Trimmed edge with ID, that is not completely new but this part
//           // is"+std::to_string(myTrimmedEdgeIndex.id.touint())}); elementEdgeIndices.push_back(id);
//           //
//           // ++myTrimmedEdgeIndex.id;
//
//         } else if (false /* the edge is completly gone?*/) {
//           //
//         }
//         //
//         // IndexType id =getStableIndex(myTrimmedEdgeIndex);
//         // elementEdgeIndices.push_back(id);
//         // //  globalEdgeMap.insert({id,Edge(edge)});
//         //
//         // globalEdgeMap.insert({getStableIndex(myTrimmedEdgeIndex),"Trimmed edge with ID
//         // "+std::to_string(myTrimmedEdgeIndex.id.touint())});
//         // ++myTrimmedEdgeIndex.id;
//       }
//       //  auto edgeGeo = intersection.geometry();
//       //  auto firstCorner = edgeGeo.corner(0);
//       // auto secondCorner = edgeGeo.corner(1);
//       // edges[0][indexMapper(i)](edgeGeo.corner(0), pos0[1]);
//     }
//     elementVertexIndices.resize(4);
//     std::ranges::copy(elementVertexSet, elementVertexIndices.begin());
//   }
//
//   for (const auto& [key, value] : entityContainer.globalEdgeMap)
//     if (key.host_or_trimmed == IdType::HostOrTrimmed::host)
//       std::cout << "Host: Index value: " << key.id << std::endl;
//     else
//       std::cout << "Trimmed: Index value: " << key.id << std::endl;
//
//   for (const auto& [key, value_vec] : entityContainer.globalEdgesIdOfElementsMap_) {
//     std::cout << "Host: Index value: " << key.id << std::endl;
//     for (auto value : value_vec)
//       std::cout << " Edge:  id: " << value.id << std::endl;
//   }
//
//   for (const auto& [key, value_vec] : entityContainer.globalVerticesIdOfElementsMap) {
//     std::cout << "Host: Index value: " << key.id << std::endl;
//     for (auto value : value_vec)
//       std::cout << " Vertex: " << value.id << std::endl;
//   }
// }

/**
 * \brief Create the paramter grid levels
 * \param grid
 */
template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::refineParameterSpaceGrid(int refCount, bool initFlag) {
  using IdType = typename GridFamily::TrimmerTraits::GlobalIdSetId;

  using EdgeHostType           = typename UntrimmedParameterSpaceGrid::template Codim<1>::Entity;
  using EdgeGridType           = typename GridImp::template Codim<1>::Entity;
  using EleGridType            = typename GridImp::template Codim<0>::Entity;
  using EdgeParameterSpaceType = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
  using EleParameterSpaceType  = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;
  const int oldLevel           = untrimmedParameterSpaceGrid_->maxLevel();
  untrimmedParameterSpaceGrid_->globalRefine(refCount);
  auto gvu = untrimmedParameterSpaceGrid_->leafGridView();
  parameterSpaceGrid_->createBegin();
  parameterSpaceGrid_->insertLeaf();
  parameterSpaceGrid_->createEnd();
  assert((initFlag and oldLevel == 0 and refCount == 0) or
         !initFlag && "If we initialize the grid, we untrimmedParameterSpaceGrid_ should only have one level");
  // if we init we start at 0 otherwise at 1
  // if the grid is refined once later this would yield int i = 1; i < 2;
  for (int i = !initFlag; i < refCount + 1; ++i) {
    const int newLevel = oldLevel + i;
    entityContainer_.entityImps_.emplace_back();
    auto& entityContainer  = entityContainer_;
    auto& elementContainer = std::get<0>(entityContainer_.entityImps_.back());
    auto& edgeContainer    = std::get<1>(entityContainer_.entityImps_.back());
    auto& vertexContainer  = std::get<2>(entityContainer_.entityImps_.back());

    auto gv                         = parameterSpaceGrid_->levelGridView(newLevel);
    auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
    int unTrimmedElementIndex       = 0;
    int trimmedElementIndex         = 0;
    int edgeIndex                   = 0;
    int vertexIndex                 = 0;

    auto figure = matplot::figure(true);

    for (const auto& ele : elements(gv)) {
      ElementTrimFlag eleTrimFlag{ElementTrimFlag::full};
      ElementTrimData eleTrimData(eleTrimFlag, ele);
      if (trimData_.has_value()) {
        eleTrimData = trimElement(ele, trimData_.value());
        eleTrimFlag = eleTrimData.flag();

        // Testing Purposes
        auto resName = "ele_" + std::to_string(unTrimmedElementIndex + trimmedElementIndex);
        // eleTrimData.drawResult(resName, true);
        // eleTrimData.drawResult(resName, false);
        eleTrimData.drawResult(resName, false, false);
      }

      auto hostId = globalIdSetParameterSpace.id(ele);

      IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
      if (ele.hasFather()) {
        auto fatherHostId = globalIdSetParameterSpace.id(ele.father());
        // if we know that we (the element) are untrimmed we know that our father is also untrimmed
        if (eleTrimFlag == ElementTrimFlag::full) {
          IdType fatherId = {.entityIdType = IdType::EntityIdType::host, .id = fatherHostId};
          // insert new element info
          EntityInfo<0> elementInfo{.indexInLvlStorage   = trimmedElementIndex + unTrimmedElementIndex,
                                    .unTrimmedIndexInLvl = unTrimmedElementIndex,
                                    .lvl                 = newLevel,
                                    .id                  = elementId,
                                    .fatherId            = std::make_optional<IdType>(fatherId)};
          ++unTrimmedElementIndex;
          entityContainer.idToElementInfoMap.insert({elementId, elementInfo});
          elementContainer.emplace_back(grid_, ele, elementInfo);

          // since we have a father we have to add us as his son, this can be faster, we can store in decendantIds,
          // the indexInLvlStorage and lvl, which would provide faster access
          entityContainer.idToElementInfoMap.at(fatherId).decendantIds.push_back(elementId);
          entityContainer.template entity<0>(fatherId, newLevel - 1).entityInfo_.decendantIds.push_back(elementId);
        } else if (eleTrimFlag == ElementTrimFlag::trimmed) {
          // trimmed
          ++trimmedElementIndex;
        } else {
          // empty
        }
      } else /* no father */ {
        if (eleTrimFlag == ElementTrimFlag::full) {
          EntityInfo<0> elementInfo{.indexInLvlStorage   = trimmedElementIndex + unTrimmedElementIndex,
                                    .unTrimmedIndexInLvl = unTrimmedElementIndex,
                                    .lvl                 = newLevel,
                                    .id                  = elementId};
          elementContainer.emplace_back(grid_, ele, elementInfo);
          entityContainer.idToElementInfoMap.insert({elementId, elementInfo});
          ++unTrimmedElementIndex;
        } else if (eleTrimFlag == ElementTrimFlag::trimmed) {
          // trimmed
          ++trimmedElementIndex;
        } else {
          // empty
        }
      }
      auto& elementEdgeIndices   = entityContainer.globalEdgesIdOfElementsMap_[elementId];
      auto& elementVertexIndices = entityContainer.globalVerticesIdOfElementsMap[elementId];

      for (int localEdgeIndex = 0; localEdgeIndex < ele.subEntities(1); ++localEdgeIndex) {
        // setup all edge indices for given element
        if (true /* untrimmed */) {
          auto hostEdgeId = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
          IdType edgeId   = {.entityIdType = IdType::EntityIdType::host, .id = hostEdgeId};
          elementEdgeIndices.emplace_back(edgeId);

          // store vertex ids of edges, for subIndex method of indeSet
          auto& edgeVertexIndices = entityContainer.globalVertexIdOfEdgesMap_[edgeId];
          if (edgeVertexIndices.size() < 2) // if this edge already has two vertices we don't visit again
          {
            const auto& cube = Dune::ReferenceElements<ctype, mydimension>::cube();
            for (auto vertexLocalIndexWRTElement : cube.subEntities(localEdgeIndex, 1, 2)) {
              auto hostVertexId = globalIdSetParameterSpace.subId(ele, vertexLocalIndexWRTElement, 2);
              IdType vertexId   = {.entityIdType = IdType::EntityIdType::host, .id = hostVertexId};
              edgeVertexIndices.emplace_back(vertexId);
            }
          }
        }
      }

      // store vertices ids of element
      for (int localVertexId = 0; localVertexId < ele.subEntities(2); ++localVertexId) {
        // setup all vertex indices for given element
        if (true /* untrimmed */) {
          auto hostVertexId = globalIdSetParameterSpace.subId(ele, localVertexId, 2);
          IdType vertexId   = {.entityIdType = IdType::EntityIdType::host, .id = hostVertexId};
          elementVertexIndices.emplace_back(vertexId);
        }
      }
    }
    matplot::save("grid", "gif");

    // save numbers of untrimmed and trimmed elements per level
    entityContainer.numberOfTrimmedElements.push_back(trimmedElementIndex);
    entityContainer.numberOfUnTrimmedElements.push_back(unTrimmedElementIndex);
    for (const auto& edge : edges(gv)) {
      auto edgeHostId = globalIdSetParameterSpace.id(edge);
      IdType edgeId   = {.entityIdType = IdType::EntityIdType::host, .id = edgeHostId};
      EntityInfo<1> edgeInfo{.indexInLvlStorage = edgeIndex++, .lvl = newLevel, .stemFromTrim = false, .id = edgeId};

      entityContainer.idToEdgeInfoMap.insert({edgeId, edgeInfo});

      edgeContainer.emplace_back(grid_, edge, edgeInfo);
    }
    entityContainer.idToVertexInfoMap.emplace_back();
    for (const auto& vertex : vertices(gv)) {
      auto vertexHostId = globalIdSetParameterSpace.id(vertex);
      IdType vertexId   = {.entityIdType = IdType::EntityIdType::host, .id = vertexHostId};
      EntityInfo<2> vertexInfo{
          .indexInLvlStorage = vertexIndex++, .lvl = newLevel, .stemFromTrim = false, .id = vertexId};

      entityContainer.idToVertexInfoMap.back().insert({vertexId, vertexInfo});
      vertexContainer.emplace_back(grid_, vertex, vertexInfo);
    }
  }
}
} // namespace Dune::IGANEW::DefaultTrim
