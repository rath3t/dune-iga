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
  std::map< IdType,StringAndIndex> globalEdgeMap;
  std::map<IdType, std::vector<IdType>> globalEdgesIdOfElementsMap_;
  std::map<IdType, std::vector<IdType>> globalVerticesIdOfElementsMap;
};
}
