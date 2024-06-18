// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

/**
 * \file
 * @brief The PatchGridEntitySeed class
 */

namespace Dune::IGA::IdentityParameterSpace {

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
  // Entity type of the hostgrid
  using ParameterSpaceGridEntity = typename ParameterSpace::template Codim<codim>::ParameterSpaceGridEntity;

  // EntitySeed type of the hostgrid
  using ParameterSpaceGridEntitySeed = typename ParameterSpace::template Codim<codim>::ParameterSpaceGridEntitySeed;

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
  explicit PatchGridEntitySeed(const ParameterSpaceGridEntity& hostEntity)
      : hostEntitySeed_(hostEntity.seed()) {}

  /**
   * @brief Get stored ParameterSpaceGridEntitySeed
   */
  const ParameterSpaceGridEntitySeed& hostEntitySeed() const {
    return hostEntitySeed_;
  }

  /**
   * @brief Check whether it is safe to create an Entity from this Seed
   */
  bool isValid() const {
    return hostEntitySeed_.isValid();
  }

private:
  ParameterSpaceGridEntitySeed hostEntitySeed_;
};

} // namespace Dune::IGA::IdentityParameterSpace

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
