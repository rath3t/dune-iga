// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune::IGA {

template <class GridImp>
struct PatchGridLevelGridView;

template <class GridImp>
struct PatchGridLeafGridView;

template <class GridImp>
struct PatchGridLevelGridViewTraits : DefaultLevelGridViewTraits<const GridImp>
{
  typedef PatchGridLevelGridView<GridImp> GridViewImp;
};

template <class GridImp>
struct PatchGridLevelGridView : DefaultLevelGridView<const GridImp>
{
  typedef PatchGridLevelGridView ThisType;
  using ParameterSpaceType = typename GridImp::ParameterSpace;

  PatchGridLevelGridView(const typename DefaultLevelGridView<const GridImp>::Grid& grid, int level)
      : DefaultLevelGridView<const GridImp>(grid, level) {}

  const auto& patchData() const {
    return this->grid().patchGeometries[this->level_].patchData();
  }

  const auto& unTrimmedPatch() const {
    return this->grid().patchGeometries[this->level_];
  }

  auto untrimmedElementNumbers() const {
    return this->grid().untrimmedElementNumbers(this->level_);
  }

  const auto& tensorProductCoordinates() const {
    return this->grid().tensorProductCoordinates(this->level_);
  }

  int level() const {
    return this->level_;
  }
};

template <class GridImp>
struct PatchGridLeafGridViewTraits : public DefaultLeafGridViewTraits<const GridImp>
{
  typedef PatchGridLeafGridView<GridImp> GridViewImp;
};

template <class GridImp>
struct PatchGridLeafGridView : public DefaultLeafGridView<const GridImp>
{
  typedef PatchGridLeafGridView ThisType;

  PatchGridLeafGridView(const typename DefaultLeafGridView<const GridImp>::Grid& grid)
      : DefaultLeafGridView<const GridImp>(grid) {}

  using ParameterSpaceType = typename GridImp::ParameterSpace;
  const auto& patchData() const {
    return this->grid().patchGeometries_[this->grid().maxLevel()].patchData();
  }

  const auto& unTrimmedPatch() const {
    return this->grid().patchGeometries[this->grid().maxLevel()];
  }

  const auto& tensorProductCoordinates() const {
    return this->grid().tensorProductCoordinates(this->grid().maxLevel());
  }

  std::array<int, GridImp::dimension> untrimmedElementNumbers() const {
    return this->grid().untrimmedElementNumbers(this->grid().maxLevel());
  }

  int level() const {
    return this->grid().maxLevel();
  }
};

} // namespace Dune::IGA
