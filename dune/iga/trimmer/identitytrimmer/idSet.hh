// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
namespace Dune::IGA::IdentityTrim {
template <class GridImp>
class PatchGridGlobalIdSet : public IdSet<GridImp, PatchGridGlobalIdSet<GridImp>,
                                          typename std::remove_const_t<GridImp>::Traits::GlobalIdSet::IdType>
{
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

public:
  // constructor stores reference to a grid
  explicit PatchGridGlobalIdSet(const GridImp& g)
      : grid_(&g) {}

  // define the type used for persistent indices
  typedef typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType IdType;

  // get id of an entity
  /*
     We use the remove_const to extract the Type from the mutable class,
     because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return grid_->parameterSpaceGrid().globalIdSet().id(e.impl().getLocalEntity());
  }

  // get id of subEntity
  /*

   */
  IdType subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i,
               int codim) const {
    return grid_->parameterSpaceGrid().globalIdSet().subId(e.impl().getLocalEntity(), i, codim);
  }

  void update() {}

  const GridImp* grid_;
};

template <class GridImp>
class PatchGridLocalIdSet
    : public IdSet<GridImp, PatchGridLocalIdSet<GridImp>,
                   typename std::remove_const<GridImp>::type::ParameterSpaceGrid::Traits::LocalIdSet::IdType>
{
private:
  typedef typename std::remove_const<GridImp>::type::ParameterSpaceGrid ParameterSpaceGrid;

public:
  // define the type used for persistent local ids
  typedef typename ParameterSpaceGrid::Traits::LocalIdSet::IdType IdType;

  // constructor stores reference to a grid
  PatchGridLocalIdSet(const GridImp& g)
      : grid_(&g) {}

  // get id of an entity
  /*
      We use the remove_const to extract the Type from the mutable class,
      because the const class is not instantiated yet.
   */
  template <int cd>
  IdType id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const {
    // Return id of the host entity
    return grid_->parameterSpaceGrid().localIdSet().id(e.impl().getLocalEntity());
  }

  // get id of subEntity
  /*
   * We use the remove_const to extract the Type from the mutable class,
   * because the const class is not instantiated yet.
   */
  IdType subId(const typename std::remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const {
    // Return sub id of the host entity
    /* @todo Trim, the sub indices are wrong!!! */

    return grid_->parameterSpaceGrid().localIdSet().subId(e.impl().getLocalEntity(), i, codim);
  }

  /** @todo Should be private */
  void update() {}

  const GridImp* grid_;
};
} // namespace Dune::IGA::IdentityTrim
