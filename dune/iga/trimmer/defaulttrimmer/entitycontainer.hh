// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <map>
#include <set>

#include <dune/common/reservedvector.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/indextransformations.hh>

namespace Dune::IGA::DefaultTrim {
template <typename GridImp>
struct VectorEntityContainer
{
  using Trimmer = typename GridImp::Trimmer;

  static constexpr int gridDim        = GridImp::dimension;
  static constexpr int dimensionworld = GridImp::dimensionworld;
  using ctype                         = typename GridImp::ctype;
  template <int cd>
  using EntityInfo = typename Trimmer::TrimmerTraits::template Codim<cd>::EntityInfo;

  using GeoTypes = typename GridImp::GridFamily::GeometryTypes;

  using IdType = typename Trimmer::TrimmerTraits::GlobalIdSetId;
  struct StringAndIndex
  {
    std::string msg;
    IdType id;
  };
  template <int codim>
  using Entity = typename Trimmer::TrimmerTraits::template Codim<codim>::ParameterSpaceGridEntity;
  // template <int codim>
  // using EntityMap    = std::map<IdType, ElementEntity<codim>>;
  // using EdgeEntity   = typename GridImp::template Codim<1>::Entity;
  // using VertexEntity = typename GridImp::template Codim<2>::Entity;

  template <int codim>
  auto begin(int lvl) const {
    return std::get<codim>(entityImps_[lvl]).begin();
  }

  template <int codim>
  auto end(int lvl) const {
    return std::get<codim>(entityImps_[lvl]).end();
  }

  template <int codim>
  auto begin(int lvl) {
    return std::get<codim>(entityImps_[lvl]).begin();
  }

  template <int codim>
  auto end(int lvl) {
    return std::get<codim>(entityImps_[lvl]).end();
  }

  const IdType& subId(const IdType& elementId, int localSubIndex, int codim) const {
    if (codim == 0)
      return elementId;
    if (codim == 1)
      return globalEdgesIdOfElementsMap_.at(elementId)[localSubIndex];
    if (codim == 2)
      return globalVerticesIdOfElementsMap.at(elementId)[localSubIndex];
    assert(codim >= 0 and codim <= 2);
    __builtin_unreachable();
  }

  size_t subIds(const IdType& elementId, int codim) const {
    if (codim == 0)
      return 0;
    if (codim == 1)
      return globalEdgesIdOfElementsMap_.at(elementId).size();
    if (codim == 2)
      return globalVerticesIdOfElementsMap.at(elementId).size();
    assert(codim >= 0 and codim <= 2);
    __builtin_unreachable();
  }

  int outsideIntersectionIndex(const IdType& insideElementId, const IdType& outsideElementId, int indexInInside) const {
    const IdType& insideSubId = subId(insideElementId, indexInInside, 1);

    for (const auto i : Dune::range(subIds(outsideElementId, 1))) {
      const IdType& outsideSubId = subId(outsideElementId, i, 1);
      if (outsideSubId == insideSubId) {
        return i;
      }
    }
    DUNE_THROW(GridError, "outsideIntersectionIndex not successfull");
  }

  bool isElementTrimmed(const IdType& elementId) const {
    return infoFromId<0>(elementId).trimmed;
  }

  template <int codim>
  requires(codim >= 0 and codim <= 2)
  const auto& entity(int lvl, int indexInLvlStorage) const {
    return std::get<codim>(entityImps_[lvl])[indexInLvlStorage];
  }

  template <int codim>
  requires(codim >= 0 and codim <= 2)
  auto& entity(int lvl, int indexInLvlStorage) {
    return std::get<codim>(entityImps_[lvl])[indexInLvlStorage];
  }

  template <int codim>
  auto& entity(const EntityInfo<codim>& info) {
    return entity<codim>(info.lvl, info.indexInLvlStorage);
  }

  template <int codim>
  const auto& entity(const EntityInfo<codim>& info) const {
    return entity<codim>(info.lvl, info.indexInLvlStorage);
  }

  template <int codim>
  requires(codim >= 0 and codim <= 2)
  auto& entity(const IdType& id, int lvl) {
    return entity<codim>(infoFromId<codim>(id, lvl));
  }

  template <int codim>
  requires(codim >= 0 and codim <= 2)
  const auto& entity(const IdType& id, int lvl) const {
    return entity<codim>(infoFromId<codim>(id, lvl));
  }

  // todo specialize for codim 0 and 2 to not need lvl
  template <int codim>
  requires(codim >= 0 and codim <= 2)
  const auto& infoFromId(const IdType& id, int lvl = std::numeric_limits<int>::max()) const {
    if constexpr (codim == 0)
      return idToElementInfoMap.at(id);
    else if constexpr (codim == 1)
      return idToEdgeInfoMap.at(id);
    else
      return idToVertexInfoMap[lvl].at(id);
  }

  template <int cc>
  auto subIndexFromId(const IdType& id, int i, int codim, int lvl) const {
    assert(codim > 0);
    if constexpr (cc == 0) {
      if (codim == 1)
        return idToEdgeInfoMap.at(globalEdgesIdOfElementsMap_.at(id)[i]).indexInLvlStorage;
      else if (codim == 2)
        return idToVertexInfoMap[lvl].at(globalVerticesIdOfElementsMap.at(id)[i]).indexInLvlStorage;
      __builtin_unreachable();
    } else if constexpr (cc == 1) {
      if (codim == 1)
        return idToEdgeInfoMap.at(id).indexInLvlStorage;
      else if (codim == 2) // get index of vertex from edge id
        return idToVertexInfoMap[lvl].at(globalVertexIdOfEdgesMap_.at(id)[i]).indexInLvlStorage;
      __builtin_unreachable();
    } else {
      assert(codim == 2 and i == 0);
      return idToVertexInfoMap[lvl].at(id).indexInLvlStorage;
    }
  }

  // for each level we have codim+1 vectors of entities
  template <std::size_t... codim>
  static auto makeEntityImpsImpl_(std::index_sequence<codim...>, std::size_t numLevels) {
    return std::vector<std::tuple<std::vector<Entity<codim>>...>>(numLevels);
  }

  // Create the lists of vertices, edges, elements for each level
  static auto makeEntityImps_(std::size_t numLevels = 1) {
    return makeEntityImpsImpl_(std::make_index_sequence<gridDim + 1>{}, numLevels);
  }

  GeoTypes types(int codim, int level) const {
    if (codim == 0) {
      if (numberOfTrimmedElements[level] == 0)
        return GeoTypes{GeometryTypes::cube(gridDim)};

      if (numberOfUnTrimmedElements[level] == 0)
        return GeoTypes{GeometryTypes::none(gridDim)};

      return GeoTypes{GeometryTypes::cube(gridDim), GeometryTypes::none(gridDim)};
    }
    return GeoTypes{GeometryTypes::cube(gridDim - codim)};
  }

  template <typename EntityType>
  bool contains(const EntityType& e, int lvl) const {
    auto& entities = std::get<EntityType::codimension>(entityImps_[lvl]);
    return std::ranges::find(entities, e.impl().getLocalEntity()) != entities.end();
  }

  std::size_t size(int codim, int lvl) const {
    if (codim == 0)
      return std::get<0>(entityImps_[lvl]).size();
    else if (codim == 1)
      return std::get<1>(entityImps_[lvl]).size();
    else if (codim == 2)
      return std::get<2>(entityImps_[lvl]).size();
    assert(codim >= 0 and codim <= 2);
    __builtin_unreachable();
  }

  std::size_t size(GeometryType type, int lvl) const {
    if (type == GeometryTypes::quadrilateral)
      return numberOfUnTrimmedElements[lvl];
    if (type == GeometryTypes::none(gridDim))
      return numberOfTrimmedElements[lvl];
    if (type == GeometryTypes::line)
      return std::get<1>(entityImps_[lvl]).size();
    if (type == GeometryTypes::vertex)
      return std::get<2>(entityImps_[lvl]).size();
    return 0;
  }

  std::size_t sizeOfInfos(int codim, int lvl) const {
    if (codim == 2)
      return idToVertexInfoMap[lvl].size();

    auto count = [=](auto&& range) -> std::size_t {
      return std::ranges::count_if(range, [=](const auto& entityInfoPair) { return entityInfoPair.second.lvl == lvl; });
    };
    if (codim == 0)
      return count(idToElementInfoMap);
    if (codim == 1)
      return count(idToEdgeInfoMap);
    __builtin_unreachable();
  }

  // The vector type for tuple of lists of vertices, edges, elements for each level
  using EntityImps = std::decay_t<decltype(makeEntityImps_())>;

  // The tuple type of lists of vertices, edges, elements
  using EntityTuple = typename EntityImps::value_type;

  template <int codim>
  using EntityConstInteratorImpl = typename std::tuple_element_t<codim, EntityTuple>::const_iterator;

  template <int codim>
  using EntityInteratorImpl = typename std::tuple_element_t<codim, EntityTuple>::iterator;

  // The lists of vertices, edges, elements for each level
  EntityImps entityImps_;

  std::map<IdType, Dune::ReservedVector<IdType, 8>> globalEdgesIdOfElementsMap_;
  std::map<IdType, Dune::ReservedVector<IdType, 2>> globalVertexIdOfEdgesMap_;
  std::map<IdType, Dune::ReservedVector<IdType, 8>> globalVerticesIdOfElementsMap;

  std::vector<std::map<IdType, EntityInfo<2>>>
      idToVertexInfoMap; // vertices are repeated on different levels and share the same id, to differentiate them we
                         // warp it inisdew a vector which runs over the levels
  std::map<IdType, EntityInfo<1>> idToEdgeInfoMap;
  std::map<IdType, EntityInfo<0>> idToElementInfoMap;
  // store information to know what geometry types we have to return.
  std::vector<int> numberOfTrimmedElements{};
  std::vector<int> numberOfUnTrimmedElements{};

  // Inserted by original Element Idx
  std::vector<std::vector<ElementTrimFlag>> trimFlags_;

  // We store trimmed vertexIds for each level
  std::vector<std::map<IdType, FieldVector<double, 2>>> trimmedVertexIds_;
  // Store current amound of vertices and edges (untrimmed configuration + n) per lvl
  std::vector<unsigned int> edgeCount;
  std::vector<unsigned int> vertexCount;
};
} // namespace Dune::IGA::DefaultTrim
