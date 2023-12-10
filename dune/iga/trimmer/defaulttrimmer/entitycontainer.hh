// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <map>

#include <dune/common/reservedvector.hh>
namespace Dune::IGANEW::DefaultTrim {
  template <typename GridImp>
  struct VectorEntityContainer {
    using Trimmer = typename GridImp::Trimmer;

    static constexpr int gridDim        = GridImp::dimension;
    static constexpr int dimensionworld = GridImp::dimensionworld;
    using ctype                         = typename GridImp::ctype;
    template <int cd>
    using EntityInfo = typename Trimmer::TrimmerTraits::template Codim<cd>::EntityInfo;

    using IdType = typename Trimmer::TrimmerTraits::GlobalIdSetId;
    struct StringAndIndex {
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
      else if (codim == 1)
        return globalEdgesIdOfElementsMap_.at(elementId)[localSubIndex];
      else if (codim == 2)
        return globalVerticesIdOfElementsMap.at(elementId)[localSubIndex];
      assert(codim >= 0 and codim <= 2);
      __builtin_unreachable();
    }

    template <int codim>
    requires(codim >= 0 and codim <= 2) const auto& entity(int lvl, int indexInLvlStorage) const {
      return std::get<codim>(entityImps_[lvl])[indexInLvlStorage];
    }

    template <int codim>
    requires(codim >= 0 and codim <= 2) auto& entity(int lvl, int indexInLvlStorage) {
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
    requires(codim >= 0 and codim <= 2) auto& entity(const IdType& id, int lvl) { return entity<codim>(infoFromId<codim>(id,lvl)); }

    template <int codim>
    requires(codim >= 0 and codim <= 2) const auto& entity(const IdType& id, int lvl) const {
      return entity<codim>(infoFromId<codim>(id,lvl));
    }

    template <int codim>
    requires(codim >= 0 and codim <= 2) const auto& infoFromId(const IdType& id, int lvl) const {
      if constexpr (codim == 0)
        return idToElementInfoMap.at(id);
      else if constexpr (codim == 1)
        return idToEdgeInfoMap.at(id);
      else
        return idToVertexInfoMap[lvl].at(id);
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

    // The vector type for tuple of lists of vertices, edges, elements for each level
    using EntityImps = std::decay_t<decltype(makeEntityImps_())>;

    // The tuple type of lists of vertices, edges, elements
    using EntityTuple = typename EntityImps::value_type;

    template <int codim>
    using EntityConstInteratorImpl = typename std::tuple_element_t<codim, EntityTuple>::const_iterator;

    template <int codim>
    using EntityInteratorImpl = typename std::tuple_element_t<codim, EntityTuple>::iterator;

    //! The lists of vertices, edges, elements for each level
    EntityImps entityImps_;

    std::map<IdType, Dune::ReservedVector<IdType, 8>> globalEdgesIdOfElementsMap_;
    std::map<IdType, Dune::ReservedVector<IdType, 8>> globalVerticesIdOfElementsMap;

    std::vector<std::map<IdType, EntityInfo<2>>> idToVertexInfoMap; // vertices are repeated on different levels and share the same id, to differentiate them we warp it inisdew a vector which runs over the levels
    std::map<IdType, EntityInfo<1>> idToEdgeInfoMap;
    std::map<IdType, EntityInfo<0>> idToElementInfoMap;
  };
}  // namespace Dune::IGANEW::DefaultTrim
