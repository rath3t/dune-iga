// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/**
 * \file
 * @brief The PatchGridEntitySeed class
 */

namespace Dune::IGANEW::DefaultTrim {

/**
 * @brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
 * @ingroup PatchGrid
 *
 */
template <int codim, class GridImp>
class PatchGridEntitySeed
{
protected:
  using Trimmer = typename GridImp::Trimmer;
  friend Trimmer;
  // Entity type of the hostgrid
  using EntityImp = typename Trimmer::TrimmerTraits::template Codim<codim>::EntityImp;
  using Entity    = typename GridImp::template Codim<codim>::Entity;

  // EntitySeed type of the hostgrid
  using ParameterSpaceGridEntitySeed = typename Trimmer::template Codim<codim>::ParameterSpaceGridEntitySeed;
  using EntityInfo                   = typename Trimmer::TrimmerTraits::template Codim<codim>::EntityInfo;

public:
  constexpr static int codimension = codim;

  /**
   * @brief Construct an empty (i.e. isValid() == false) seed.
   */
  PatchGridEntitySeed() = default;

  /**
   * @brief Create EntitySeed from hostgrid Entity
   *
   * We call hostEntity.seed() directly in the constructor
   * of PatchGridEntitySeed to allow for return value optimization.
   */
  explicit PatchGridEntitySeed(const EntityImp& ent)
      : lvl_(ent.getHostEntity().entityInfo_.lvl),
        indexInLvlStorage_{ent.getHostEntity().entityInfo_.indexInLvlStorage} {
  }

  /**
   * @brief Get stored ParameterSpaceGridEntitySeed
   */
  // const ParameterSpaceGridEntitySeed& hostEntitySeed() const { return hostEntitySeed_; }

  /**
   * @brief Check whether it is safe to create an Entity from this Seed
   */
  bool isValid() const {
    return indexInLvlStorage_ != std::numeric_limits<unsigned int>::max();
  }

private:
  auto data() const {
    return std::make_pair(lvl_, indexInLvlStorage_);
  }

  int lvl_{};
  unsigned int indexInLvlStorage_{std::numeric_limits<unsigned int>::max()};
};

} // namespace Dune::IGANEW::DefaultTrim

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
