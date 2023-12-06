// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
  class PatchGridEntitySeed {
   protected:
    using Trimmer = typename GridImp::Trimmer;
    friend Trimmer;
    // Entity type of the hostgrid
    using Entity= typename GridImp::template Codim<codim>::Entity ;

    // EntitySeed type of the hostgrid
        using ParameterSpaceGridEntitySeed =typename Trimmer::template Codim<codim>::ParameterSpaceGridEntitySeed ;

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
    explicit PatchGridEntitySeed(const Entity& ent) : entity_(&ent) {}

    /**
     * @brief Get stored ParameterSpaceGridEntitySeed
     */
    // const ParameterSpaceGridEntitySeed& hostEntitySeed() const { return hostEntitySeed_; }

    /**
     * @brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const { return entity_!=nullptr; }

   private:
    /** \brief Access to the underlying FoamGrid data structure */
    const Entity* target() const
    {
      return entity_;
    }

    const Entity* entity_{nullptr};
  };

}  // namespace Dune::IGANEW

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
