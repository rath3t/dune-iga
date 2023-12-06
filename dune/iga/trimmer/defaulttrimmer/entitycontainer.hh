// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

namespace Dune::IGANEW::DefaultTrim {
  template <typename GridImp>
  struct VectorEntityContainer {
    using Trimmer = typename GridImp::Trimmer;

    static constexpr int gridDim = GridImp::dimension;
    static constexpr int dimensionworld = GridImp::dimensionworld;
    using ctype = typename GridImp::ctype;

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

    // EntityMap<0> globalElementMap;
    // EntityMap<1> globalEdgeMap;
    // EntityMap<2> globalVertexMap;


    //for each level we have codim+1 vectors of entities
    template<std::size_t... codim>
    static auto makeEntityImpsImpl_(std::index_sequence<codim...>, std::size_t numLevels)
    { return std::vector<std::tuple<std::vector<Entity<codim>>...>>(numLevels); }

    // Create the lists of vertices, edges, elements for each level
    static auto makeEntityImps_(std::size_t numLevels = 1)
    { return makeEntityImpsImpl_(std::make_index_sequence<gridDim+1>{}, numLevels); }

    // The vector type for tuple of lists of vertices, edges, elements for each level
    using EntityImps = std::decay_t<decltype(makeEntityImps_())>;

    // The tuple type of lists of vertices, edges, elements
    using EntityTuple = typename EntityImps::value_type;

    template <int codim>
    using EntityConstInteratorImpl = typename std::tuple_element<codim,EntityTuple>::const_iterator;

    template <int codim>
    using EntityInteratorImpl = typename std::tuple_element<codim,EntityTuple>::iterator;

    //! The lists of vertices, edges, elements for each level
    EntityImps entityImps_;

    std::map<IdType, std::vector<IdType>> globalEdgesIdOfElementsMap_;
    std::map<IdType, std::vector<IdType>> globalVerticesIdOfElementsMap;
  };
}  // namespace Dune::IGANEW::DefaultTrim
