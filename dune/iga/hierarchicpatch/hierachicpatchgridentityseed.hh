// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITY_GRID_ENTITY_SEED_HH
#define DUNE_IDENTITY_GRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The PatchGridEntitySeed class
 */


namespace Dune::IGANEW {


  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup PatchGrid
   *
   */
  template<int codim, class GridImp>
  class PatchGridEntitySeed
  {
  protected:

    // Entity type of the hostgrid
    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::Entity HostEntity;

    // EntitySeed type of the hostgrid
    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::EntitySeed HostEntitySeed;

  public:

    constexpr static int codimension = codim;

    /**
     * \brief Construct an empty (i.e. isValid() == false) seed.
     */
    PatchGridEntitySeed()
    {}

    /**
     * \brief Create EntitySeed from hostgrid Entity
     *
     * We call hostEntity.seed() directly in the constructor
     * of PatchGridEntitySeed to allow for return value optimization.
     */
    PatchGridEntitySeed(const HostEntity& hostEntity) :
      hostEntitySeed_(hostEntity.seed())
    {}

    /**
     * \brief Get stored HostEntitySeed
     */
    const HostEntitySeed& hostEntitySeed() const
    {
      return hostEntitySeed_;
    }

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const
    {
      return hostEntitySeed_.isValid();
    }
  private:

    HostEntitySeed hostEntitySeed_;
  };

} // namespace Dune

#endif  // #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
