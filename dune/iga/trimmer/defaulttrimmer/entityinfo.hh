// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/reservedvector.hh>

namespace Dune::IGA::DefaultTrim {

template <typename HostIdType>
struct IdType;

/**
 *
 * @tparam Traits The trimmer traits
 * @tparam codim
 */
template <typename Traits, int codim>
struct EntityInfoImpl
{
};

template <typename Traits>
struct EntityInfoImpl<Traits, 2>
{
  static constexpr int codimension = 2;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;
  using TrimInfo       = typename Traits::ElementTrimData::VertexInfo;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id;
  EntitySeedType hostSeed{};

  std::optional<TrimInfo> trimInfo{};

  bool isTrimmed() const {
    return trimmed;
  }

  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }
};

template <typename Traits>
struct EntityInfoImpl<Traits, 1>
{
  static constexpr int codimension = 1;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;
  using TrimmedEntityGeometry =
      typename Traits::template Codim<codimension>::TrimmedParameterSpaceGeometry::PatchGeometry;

  using TrimInfo = typename Traits::ElementTrimData::EdgeInfo;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id;
  EntitySeedType hostSeed{};

  struct GeometryMap
  {
    unsigned int indexOfInsideElementinLvl;
    TrimmedEntityGeometry geometry;
  };

  /**
   *
   * @return an oriented geometry (if trimmed Host)
   */
  TrimmedEntityGeometry intersectionGeometry() const {
    assert(isTrimmed());
    if (isTrimmedHost())
      return edgeGeometries.front().geometry;

    // Check if intersection is vertical or horizontally orientated
    auto jac = edgeGeometries.front().geometry.jacobian({0.5});

    return edgeGeometries.front().geometry;
  }

  std::vector<GeometryMap> edgeGeometries{};
  std::optional<TrimInfo> trimInfo{};

  bool isTrimmed() const {
    return trimmed;
  }
  bool isTrimmedHost() const {
    return hostSeed.isValid() and isTrimmed();
  }
  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }
};

template <typename Traits>
struct EntityInfoImpl<Traits, 0>
{
  static constexpr int codimension = 0;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  unsigned int unTrimmedIndexInLvl{std::numeric_limits<unsigned int>::max()};
  unsigned int trimmedIndexInLvl{std::numeric_limits<unsigned int>::max()};
  unsigned int hostIndexInLvl{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id{};
  EntitySeedType hostSeed{};

  auto isTrimmed() const {
    return trimmed;
  }
  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }

  std::optional<IdType<HostIdType>> fatherId;
  ReservedVector<IdType<HostIdType>, 4> decendantIds;
};

} // namespace Dune::IGA::DefaultTrim