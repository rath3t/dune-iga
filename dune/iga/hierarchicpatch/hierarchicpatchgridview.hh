// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/defaultgridview.hh>

namespace Dune::IGANEW {

  template <class GridImp>
  struct PatchGridLevelGridView;

  template <class GridImp>
  struct PatchGridLeafGridView;

  template <class GridImp>
  struct PatchGridLevelGridViewTraits : DefaultLevelGridViewTraits<const GridImp> {
    typedef PatchGridLevelGridView<GridImp> GridViewImp;
  };

  template <class GridImp>
  struct PatchGridLevelGridView : public DefaultLevelGridView<const GridImp> {
    typedef PatchGridLevelGridView ThisType;
    static constexpr Trimming trim = GridImp::trim;

    PatchGridLevelGridView(const typename DefaultLevelGridView<const GridImp>::Grid& grid, int level)
        : DefaultLevelGridView<const GridImp>(grid, level) {}

    const auto& patchData() const { return this->grid().patchGeometries[this->level_].patchData(); }

    const auto& unTrimmedPatch() const {
      // TODO Trim
      return this->grid().patchGeometries[this->level_];
    }

    auto untrimmedElementNumbers() const {
      // TODO Trim this should be the quantity from the untrimmed grid
      if constexpr (trim == Trimming::Disabled)
        return this->grid().getHostGrid().levelSize(this->level_);
      else
        DUNE_THROW(Dune::NotImplemented, "This needs to be implemented");
    }
  };

  template <class GridImp>
  struct PatchGridLeafGridViewTraits : public DefaultLeafGridViewTraits<const GridImp> {
    typedef PatchGridLeafGridView<GridImp> GridViewImp;
  };

  template <class GridImp>
  struct PatchGridLeafGridView : public DefaultLeafGridView<const GridImp> {
    typedef PatchGridLeafGridView ThisType;

    PatchGridLeafGridView(const typename DefaultLeafGridView<const GridImp>::Grid& grid)
        : DefaultLeafGridView<const GridImp>(grid) {}

    static constexpr Trimming trim = GridImp::trim;
    const auto& patchData() const { return this->grid().patchGeometries[this->grid().maxLevel()].patchData(); }

    const auto& unTrimmedPatch() const {
      // TODO Trim
      return this->grid().patchGeometries[this->grid().maxLevel()];
    }

    auto untrimmedElementNumbers() const { return this->grid().getHostGrid().levelSize(this->grid().maxLevel()); }
  };

}  // namespace Dune::IGANEW
