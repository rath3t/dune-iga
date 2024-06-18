// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

/**
 * \file
 * @brief The PatchGridEntitySeed class
 */

namespace Dune::IGA::DefaultParameterSpace {

/**
 * @brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
 * @ingroup PatchGrid
 *
 */
template <int codim, class GridImp>
class PatchGridEntitySeed
{
protected:
  using ParameterSpace = typename GridImp::ParameterSpace;
  friend ParameterSpace;
  // Entity type of the hostgrid
  using EntityImp = typename ParameterSpace::ParameterSpaceTraits::template Codim<codim>::EntityImp;
  using Entity    = typename GridImp::template Codim<codim>::Entity;

  // EntitySeed type of the hostgrid
  using ParameterSpaceGridEntitySeed = typename ParameterSpace::template Codim<codim>::ParameterSpaceGridEntitySeed;
  using EntityInfo                   = typename ParameterSpace::ParameterSpaceTraits::template Codim<codim>::EntityInfo;

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
      : lvl_(ent.getLocalEntity().entityInfo_.lvl),
        indexInLvlStorage_{ent.getLocalEntity().entityInfo_.indexInLvlStorage} {
    assert(isValid());
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

} // namespace Dune::IGA::DefaultParameterSpace

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
