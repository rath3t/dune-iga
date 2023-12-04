// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

namespace Dune::IGANEW::DefaultTrim
{
template<typename GridImp>
struct EntityContainer
{
 using TrimmerType= typename GridImp::TrimmerType;


  using IdType = typename TrimmerType::GlobalIdSetIdType;
  struct StringAndIndex {
    std::string msg;
    IdType id;
  };
  template<int codim>
  using ElementEntity= typename TrimmerType::template TrimmedParameterSpaceGridEntity<codim,GridImp>;
  template<int codim>
  using EntityMap =   std::map< IdType,ElementEntity<codim>> ;
  using EdgeEntity= typename GridImp::template Codim<1>::Entity;
  using VertexEntity= typename GridImp::template Codim<2>::Entity;

  template<int codim>
  using EntityConstInteratorImpl = typename EntityMap<codim>::const_iterator;

  template<int codim>
using EntityInteratorImpl = typename EntityMap<codim>::const_iterator;

  template<int codim>
  auto begin() const {
     if constexpr (codim==0)
      return globalElementMap.begin();
     else      if constexpr (codim==1)
       return globalEdgeMap.begin();
     else      if constexpr (codim==2)
       return globalVertexMap.begin();
  }

  template<int codim>
auto end() const {
    if constexpr (codim==0)
      return globalElementMap.end();
    else      if constexpr (codim==1)
      return globalEdgeMap.end();
    else      if constexpr (codim==2)
      return globalVertexMap.end();
  }

  template<int codim>
auto begin()  {
    if constexpr (codim==0)
      return globalElementMap.begin();
    else      if constexpr (codim==1)
      return globalEdgeMap.begin();
    else      if constexpr (codim==2)
      return globalVertexMap.begin();
  }

  template<int codim>
auto end()  {
    if constexpr (codim==0)
      return globalElementMap.end();
    else      if constexpr (codim==1)
      return globalEdgeMap.end();
    else      if constexpr (codim==2)
      return globalVertexMap.end();
  }

  auto subId(const IdType& idEle, int localId, int codim) {
    assert(codim>0 and codim<3);
    if (codim==1)
      return globalEdgesIdOfElementsMap_.at(idEle)[localId];
    else     if (codim==2)
      return globalVerticesIdOfElementsMap.at(idEle)[localId];
  }

  EntityMap<0> globalElementMap;
  EntityMap<1> globalEdgeMap;
  EntityMap<2> globalVertexMap;
  std::map<IdType, std::vector<IdType>> globalEdgesIdOfElementsMap_;
  std::map<IdType, std::vector<IdType>> globalVerticesIdOfElementsMap;
};
}
