// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <ranges>

#include <dune/iga/geometrykernel/geohelper.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/cliputils.hh>

namespace Dune::IGANEW::DefaultTrim {

namespace Util {
  auto sameCorner(const auto& corner1, const auto& corner2, double precision = 1e-8) -> bool {
    return std::ranges::all_of(Dune::range(2),
                               [&](auto i) -> bool { return Dune::FloatCmp::eq(corner1[i], corner2[i], precision); });
  };

  bool isSameEdgeGeometry(const auto& geo1, const auto& geo2) {
    return sameCorner(geo1.corner(0), geo2.corner(0)) and sameCorner(geo1.corner(1), geo2.corner(1)) ||
           sameCorner(geo1.corner(0), geo2.corner(1)) and sameCorner(geo1.corner(1), geo2.corner(0)) ||
           sameCorner(geo1.corner(1), geo2.corner(0)) and sameCorner(geo1.corner(0), geo2.corner(1));
  }
} // namespace Util

template <int dim, int dimworld, typename ScalarType>
auto TrimmerImpl<dim, dimworld, ScalarType>::makeElementID(const HostEntity<0>& ele) -> GlobalIdType {
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  auto hostId                     = globalIdSetParameterSpace.id(ele);
  return {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostId};
}

template <int dim, int dimworld, typename ScalarType>
// auto TrimmerImpl<dim, dimworld, ScalarType>::idForTrimmedVertex(const typename TrimmerTraits::template
// Codim<2>::TrimmedParameterSpaceGeometry::PatchGeometry& vertex)
auto TrimmerImpl<dim, dimworld, ScalarType>::idForTrimmedVertex(const FieldVector<double, 2>& vertex) -> GlobalIdType {
  using PersistentIndexType = typename TrimmerTraits::PersistentIndexType;

  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  auto& vertexIdContainer         = entityContainer_.trimmedVertexIds_.back();
  bool alreadyThere               = false;
  PersistentIndexType foundIndex{};

  for (auto& [vertexId, pos] : vertexIdContainer) {
    if (Util::sameCorner(vertex, pos)) {
      alreadyThere = true;
      foundIndex   = vertexId.id;
    }
  }
  if (not alreadyThere)
    foundIndex = globalIdSet_->newFreeIndex();

  GlobalIdType vertexId{.entityIdType = GlobalIdType::EntityIdType::newId, .id = foundIndex};
  if (not alreadyThere)
    vertexIdContainer.emplace(vertexId, vertex);
  return vertexId;
}

template <int dim, int dimworld, typename ScalarType>
auto TrimmerImpl<dim, dimworld, ScalarType>::idForTrimmedHostEdge(
    typename TrimmerTraits::PersistentIndexType hostEdgeId,
    const typename ElementTrimData::EdgeInfo& trimmedEdge) -> GlobalIdType {
  using HostEdge            = typename TrimmerTraits::template Codim<1>::HostParameterSpaceGridEntity;
  using PersistentIndexType = typename TrimmerTraits::PersistentIndexType;

  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  // We have to check weather the edgeId is already there
  bool alreadyThere = false;
  PersistentIndexType foundIndex{};
  for (auto& [edgeId, edgeInfo] : entityContainer_.idToEdgeInfoMap) {
    if (edgeId.hostId.has_value() and edgeInfo.id.hostId.value() == hostEdgeId) {
      // We have the same hostId, this doens't neccesarily mean that its the same edge
      auto edgeGeo = edgeInfo.trimmedEntityGeometries.front().geometry;
      if (Util::isSameEdgeGeometry(edgeGeo, trimmedEdge.geometry.value())) {
        alreadyThere = true;
        foundIndex   = edgeId.id;
        break;
      }
    }
  }
  if (not alreadyThere)
    foundIndex = globalIdSet_->newFreeIndex();
  GlobalIdType edgeId = {
      .entityIdType = GlobalIdType::EntityIdType::newId, .id = foundIndex, .hostId = std::make_optional(hostEdgeId)};
  return edgeId;
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::createAndSaveElementInfo(
    const std::tuple<unsigned int, unsigned int, int>& indices, const HostEntity<0>& ele, bool trimmed) {
  const unsigned int unTrimmedElementIndex = std::get<0>(indices);
  const unsigned int trimmedElementIndex   = std::get<1>(indices);
  const int newLevel                       = std::get<2>(indices);

  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  auto& indexSet                  = parameterSpaceGrid_->levelGridView(newLevel).indexSet();

  std::optional<GlobalIdType> fatherId{};
  if (ele.hasFather())
    fatherId = {.entityIdType = GlobalIdType::EntityIdType::host, .id = globalIdSetParameterSpace.id(ele.father())};

  GlobalIdType elementId = makeElementID(ele);

  EntityInfo<0> elementInfo = {
      .indexInLvlStorage   = trimmedElementIndex + unTrimmedElementIndex,
      .unTrimmedIndexInLvl = unTrimmedElementIndex,
      .trimmedIndexInLvl   = trimmedElementIndex,
      .hostIndexInLvl      = indexSet.index(ele),
      .lvl                 = newLevel,
      .stemFromTrim        = trimmed,
      .id                  = elementId,
      .hostSeed            = ele.seed(),
      .fatherId            = fatherId,
  };

  entityContainer_.idToElementInfoMap.insert({elementId, elementInfo});

  // If we have a father we have to add us as his son, this can be faster, we can store in decendantIds,
  // the indexInLvlStorage and lvl, which would provide faster access
  if (fatherId.has_value()) {
    entityContainer_.idToElementInfoMap.at(fatherId.value()).decendantIds.push_back(elementId);
    entityContainer_.template entity<0>(fatherId.value(), newLevel - 1).entityInfo_.decendantIds.push_back(elementId);
  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::collectElementEdges(int level, const HostEntity<0>& ele,
                                                                 const ElementTrimData& eleTrimData) {
  auto& indexSet         = parameterSpaceGrid_->levelGridView(level).indexSet();
  GlobalIdType elementId = makeElementID(ele);

  auto& elementEdgeIndices        = entityContainer_.globalEdgesIdOfElementsMap_[elementId];
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  // Define some lambdas (can be moved to seperate functions later)
  auto addHostEdge = [&](const unsigned int localEdgeIndex) {
    auto hostEdgeId     = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
    GlobalIdType edgeId = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostEdgeId};
    elementEdgeIndices.emplace_back(edgeId);

    if (entityContainer_.idToEdgeInfoMap.contains(edgeId))
      return;

    auto edge = ele.template subEntity<1>(localEdgeIndex);
    EntityInfo<1> edgeInfo{.indexInLvlStorage = indexSet.index(edge),
                           .lvl               = level,
                           .stemFromTrim      = false,
                           .id                = edgeId,
                           .hostSeed          = edge.seed()};
    entityContainer_.idToEdgeInfoMap.insert({edgeId, edgeInfo});

    // store vertex ids of edges, for subIndex method of indexSet
    auto& edgeVertexIndices = entityContainer_.globalVertexIdOfEdgesMap_[edgeId];

    const auto& cube = ReferenceElements<ctype, mydimension>::cube();
    for (auto vertexLocalIndexWRTElement : cube.subEntities(localEdgeIndex, 1, 2)) {
      auto hostVertexId     = globalIdSetParameterSpace.subId(ele, vertexLocalIndexWRTElement, 2);
      GlobalIdType vertexId = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostVertexId};
      edgeVertexIndices.emplace_back(vertexId);
    }
  };

  auto addTrimmedHostEdge = [&](int localEdgeIndex, const typename ElementTrimData::EdgeInfo& edgeOfTrimmedElement) {
    // To get neighborhoud information we save the original hostId
    auto hostEdgeId     = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
    GlobalIdType edgeId = idForTrimmedHostEdge(hostEdgeId, edgeOfTrimmedElement);
    elementEdgeIndices.emplace_back(edgeId);

    if (entityContainer_.idToEdgeInfoMap.contains(edgeId)) {
      entityContainer_.idToEdgeInfoMap.at(edgeId).trimmedEntityGeometries.emplace_back(
          indexSet.index(ele), edgeOfTrimmedElement.geometry.value());
      return;
    }

    auto edge = ele.template subEntity<1>(localEdgeIndex);
    EntityInfo<1> edgeInfo{.indexInLvlStorage = entityContainer_.edgeCount.back()++,
                           .lvl               = level,
                           .stemFromTrim      = true,
                           .id                = edgeId,
                           .hostSeed          = edge.seed(),
                           .trimInfo          = edgeOfTrimmedElement};
    edgeInfo.trimmedEntityGeometries.emplace_back(indexSet.index(ele), edgeOfTrimmedElement.geometry.value());
    entityContainer_.idToEdgeInfoMap.insert({edgeId, edgeInfo});

    // store vertex ids of edges, for subIndex method of indexSet
    auto& edgeVertexIndices = entityContainer_.globalVertexIdOfEdgesMap_[edgeId];

    const auto& cube    = ReferenceElements<ctype, mydimension>::cube();
    auto subEntityRange = cube.subEntities(localEdgeIndex, 1, 2);
    std::vector<std::size_t> subEntityVector{subEntityRange.begin(), subEntityRange.end()};

    // We have one vertex that belongs to the hostgrid, and one new one
    assert(edgeOfTrimmedElement.direction == ElementTrimData::TrimmedHostEdgeDirection::HostNew or
           edgeOfTrimmedElement.direction == ElementTrimData::TrimmedHostEdgeDirection::NewHost);

    std::size_t vertexLocalIndexWRTElement =
        edgeOfTrimmedElement.direction == ElementTrimData::TrimmedHostEdgeDirection::HostNew ? subEntityVector[0]
                                                                                             : subEntityVector[1];
    auto hostVertexId = globalIdSetParameterSpace.subId(ele, vertexLocalIndexWRTElement, 2);
    edgeVertexIndices.emplace_back(GlobalIdType{.entityIdType = GlobalIdType::EntityIdType::host, .id = hostVertexId});

    // Now the new one
    auto edgeGeo   = edgeOfTrimmedElement.geometry.value();
    auto newVertex = edgeOfTrimmedElement.direction == ElementTrimData::TrimmedHostEdgeDirection::HostNew
                         ? edgeGeo.corner(1)
                         : edgeGeo.corner(0);
    edgeVertexIndices.emplace_back(idForTrimmedVertex(newVertex));
  };

  auto addTrimmedEdge = [&](const typename ElementTrimData::EdgeInfo& edgeOfTrimmedElement) {
    GlobalIdType edgeId = {.entityIdType = GlobalIdType::EntityIdType::newId, .id = globalIdSet_->newFreeIndex()};
    elementEdgeIndices.emplace_back(edgeId);

    EntityInfo<1> edgeInfo{.indexInLvlStorage = entityContainer_.edgeCount.back()++,
                           .lvl               = level,
                           .stemFromTrim      = true,
                           .id                = edgeId,
                           .trimInfo          = edgeOfTrimmedElement};
    entityContainer_.idToEdgeInfoMap.insert({edgeId, edgeInfo});

    auto& edgeVertexIndices = entityContainer_.globalVertexIdOfEdgesMap_[edgeId];
    for (auto c : Dune::range(2)) {
      auto vertexId = idForTrimmedVertex(edgeOfTrimmedElement.geometry.value().corner(c));
      edgeVertexIndices.emplace_back(vertexId);
    }
  };

  if (eleTrimData.flag() == ElementTrimFlag::full) {
    for (auto localEdgeIndex : Dune::range(ele.subEntities(1))) {
      addHostEdge(localEdgeIndex);
    }
  } else /* trimmed */ {
    for (const auto& edgeOfTrimmedElement : eleTrimData.edges()) {
      if (edgeOfTrimmedElement.isHost and not edgeOfTrimmedElement.isTrimmed) {
        addHostEdge(Util::edgeIndexMapping[edgeOfTrimmedElement.idx]);
      } else if (edgeOfTrimmedElement.isHost and edgeOfTrimmedElement.isTrimmed) {
        // This is a host Edge which is partially trimmed
        auto localEdgeIndex = Util::edgeIndexMapping[edgeOfTrimmedElement.idx];
        addTrimmedHostEdge(localEdgeIndex, edgeOfTrimmedElement);
      } else /* new Edge*/ {
        addTrimmedEdge(edgeOfTrimmedElement);
      }
    }
  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::collectElementVertices(int level, const HostEntity<0>& ele,
                                                                    const ElementTrimData& eleTrimData) {
  using VertexGeo = typename TrimmerTraits::template Codim<2>::TrimmedParameterSpaceGeometry::PatchGeometry;

  GlobalIdType elementId = makeElementID(ele);
  auto& indexSet         = parameterSpaceGrid_->levelGridView(level).indexSet();

  auto& elementVertexIndices      = entityContainer_.globalVerticesIdOfElementsMap[elementId];
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  auto addHostVertex = [&](unsigned int localVertexIndex) {
    auto hostVertexId     = globalIdSetParameterSpace.subId(ele, localVertexIndex, 2);
    GlobalIdType vertexId = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostVertexId};
    elementVertexIndices.emplace_back(vertexId);

    auto vertex = ele.template subEntity<2>(localVertexIndex);
    EntityInfo<2> vertexInfo{.indexInLvlStorage = indexSet.index(vertex),
                             .lvl               = level,
                             .stemFromTrim      = false,
                             .id                = vertexId,
                             .hostSeed          = vertex.seed()};
    entityContainer_.idToVertexInfoMap.back().insert({vertexId, vertexInfo});
  };

  auto addNewVertex = [&](const typename ElementTrimData::VertexInfo& vertex) {
    GlobalIdType vertexId = idForTrimmedVertex(vertex.geometry.value());
    elementVertexIndices.emplace_back(vertexId);
    EntityInfo<2> vertexInfo{.indexInLvlStorage = entityContainer_.vertexCount.back()++,
                             .lvl               = level,
                             .stemFromTrim      = true,
                             .id                = vertexId,
                             .trimInfo          = vertex};
    entityContainer_.idToVertexInfoMap.back().insert({vertexId, vertexInfo});
  };

  if (eleTrimData.flag() == ElementTrimFlag::full) {
    // Save eleVertices in local order
    for (auto localVertexIndex : Dune::range(ele.subEntities(2))) {
      addHostVertex(localVertexIndex);
    }

  } else /* trimmed */ {
    for (auto& vertexInfo : eleTrimData.vertices()) {
      if (vertexInfo.isHost) {
        addHostVertex(Util::vertexIndexMapping[vertexInfo.idx]);
      } else /* new */ {
        addNewVertex(vertexInfo);
      }
    }
  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::createElements(int level, const std::vector<ElementTrimData>& trimDatas) {
  using ElementEntity = typename TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

  auto& elementMap       = entityContainer_.idToElementInfoMap;
  auto& elementContainer = std::get<0>(entityContainer_.entityImps_.back());

  elementContainer.resize(entityContainer_.numberOfTrimmedElements.back() +
                          entityContainer_.numberOfUnTrimmedElements.back());
  for (auto& [eleId, eleInfo] :
       elementMap | std::ranges::views::filter([level](const auto& pair) { return pair.second.lvl == level; })) {
    auto ele = parameterSpaceGrid_->entity(eleInfo.hostSeed);
    if (not eleInfo.stemFromTrim) {
      elementContainer[eleInfo.indexInLvlStorage] = ElementEntity{grid_, ele, eleInfo};
    } else {
      elementContainer[eleInfo.indexInLvlStorage] =
          ElementEntity{grid_, ele, eleInfo, trimDatas[eleInfo.hostIndexInLvl]};
    }
  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::createSubEntities(int level) {
  using EdgeEntity   = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
  using VertexEntity = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;

  // auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  // auto& indexSet                  = parameterSpaceGrid_->levelGridView(level).indexSet();

  auto& vertexMap = entityContainer_.idToVertexInfoMap.back();
  auto& edgeMap   = entityContainer_.idToEdgeInfoMap;

  auto& vertexContainer = std::get<2>(entityContainer_.entityImps_.back());
  auto& edgeContainer   = std::get<1>(entityContainer_.entityImps_.back());

  // we are resizing the containers, so we can add the entities in the indexSet order, the index is obtained from
  // the infos (also this is a bit more efficient as pushing back all the time)

  vertexContainer.resize(entityContainer_.vertexCount.back());
  for (auto& [vertexId, vertexInfo] : vertexMap) {
    if (not vertexInfo.stemFromTrim) {
      auto vertex                                   = parameterSpaceGrid_->entity(vertexInfo.hostSeed);
      vertexContainer[vertexInfo.indexInLvlStorage] = VertexEntity{grid_, vertex, vertexInfo};
    } else {
      vertexContainer[vertexInfo.indexInLvlStorage] = VertexEntity(grid_, vertexInfo);
    }
  }

  edgeContainer.resize(entityContainer_.edgeCount.back());
  for (auto& [edgeId, edgeInfo] :
       edgeMap | std::ranges::views::filter([level](const auto& pair) { return pair.second.lvl == level; })) {
    if (not edgeInfo.stemFromTrim) {
      auto edge                                 = parameterSpaceGrid_->entity(edgeInfo.hostSeed);
      edgeContainer[edgeInfo.indexInLvlStorage] = EdgeEntity{grid_, edge, edgeInfo};
    } else {
      edgeContainer[edgeInfo.indexInLvlStorage] = EdgeEntity{grid_, edgeInfo};
    }
  }
}
} // namespace Dune::IGANEW::DefaultTrim