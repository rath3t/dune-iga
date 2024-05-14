// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/**
 * \file
 * @brief The PatchGridEntitySeed class
 */

namespace Dune::IGA::IdentityTrim {

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
  // Entity type of the hostgrid
  using ParameterSpaceGridEntity = typename Trimmer::template Codim<codim>::ParameterSpaceGridEntity;

  // EntitySeed type of the hostgrid
  using ParameterSpaceGridEntitySeed = typename Trimmer::template Codim<codim>::ParameterSpaceGridEntitySeed;

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

} // namespace Dune::IGA::IdentityTrim

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
