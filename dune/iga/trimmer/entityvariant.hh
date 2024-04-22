// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/** \file
 * @brief The GridEntityVariant class
 */
#include <variant>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"

namespace Dune::IGANEW::Trim {

//**********************************************************************
//
// --GridEntityVariant
// --Entity
//
/** @brief The implementation of entities in a PatchGrid
 *   @ingroup PatchGrid
 *
 *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
 *  An entity of codimension c in dimension d is a d-c dimensional object.
 *
 */
// template <int codim, typename TrimmerType_, class... Implementations>
// class GridEntityVariant  {
//
//  public:
//   using TrimmerType = TrimmerType_;
//
//   auto visit(auto&& lambda) const { return std::visit(lambda, impl_); }
//
//   template <class Implementation>
//   GridEntityVariant(const Implementation& impl) : impl_(impl) {}
//
//   GridEntityVariant()                                  = default;
//   GridEntityVariant(const GridEntityVariant& other) = default;
//   template <class Implementation>
//   requires(!std::is_same_v<Implementation, GridEntityVariant>) GridEntityVariant& operator=(
//       const Implementation& impl) {
//     impl_ = impl;
//     return *this;
//   };
//   GridEntityVariant(GridEntityVariant&& other) noexcept = default;
//   GridEntityVariant& operator=(const GridEntityVariant& other) = default;
//   GridEntityVariant& operator=(GridEntityVariant&& other) noexcept = default;
//
//   bool equals(const GridEntityVariant& other) const {
//
//     return hostEntity_ == other.hostEntity_;
// }
//
//   //! returns true if father entity exists
//   bool hasFather() const {
//     return visit([](const auto& impl) { return impl.hasFather(); });
// }
//
//   //! Create EntitySeed
//   auto seed() const {
//     return visit([](const auto& impl) { return impl.seed(); });
// }
//
//   //! level of this element
//   int level() const {
//     return visit([](const auto& impl) { return impl.level(); });
// }
//
//   /** @brief The partition type for parallel computing
//    */
//   PartitionType partitionType() const {
//     return visit([](const auto& impl) { return impl.partitionType(); });
// }
//
//   /** @brief Return the number of subEntities of codimension codim.
//    */
//   unsigned int subEntities(unsigned int cc) const {
//     return visit([](const auto& impl) { return impl.subEntities(cc); });
// }
//
//     using ParameterSpaceGeometry = typename TrimmerType::template LocalParameterSpaceGeometry<codim>;
//
//
//   //! geometry of this entity
//   auto geometry() const {
//     return visit([](const auto& impl) { return impl.geometry(); });
//   }
//
//   ParameterSpaceGridEntity hostEntity_;
//
//  private:
//   const GridImp* patchGrid_;
// };

//***********************
//
//  --GridEntityVariant
//
//***********************
/** @brief Specialization for codim-0-entities.
 * @ingroup PatchGrid
 *
 * This class embodies the topological parts of elements of the grid.
 * It has an extended interface compared to the general entity class.
 * For example, Entities of codimension 0  allow to visit all neighbors.
 */
template <int codim_, typename TrimmerType_, class HostImplementation, class TrimmedImplementation>
class ParameterSpaceGridEntityVariant
{
public:
  using TrimmerType                 = TrimmerType_;
  using LocalParameterSpaceGeometry = typename TrimmerType::template LocalParameterSpaceGeometry<codim_>;

  auto visit(auto&& lambda) const {
    return std::visit(lambda, impl_);
  }

  template <class Implementation>
  ParameterSpaceGridEntityVariant(const Implementation& impl)
      : impl_(impl) {
  }

  ParameterSpaceGridEntityVariant()                                             = default;
  ParameterSpaceGridEntityVariant(const ParameterSpaceGridEntityVariant& other) = default;
  template <class Implementation>
  requires(!std::is_same_v<Implementation, ParameterSpaceGridEntityVariant>)
  ParameterSpaceGridEntityVariant& operator=(const Implementation& impl) {
    impl_ = impl;
    return *this;
  };
  ParameterSpaceGridEntityVariant(ParameterSpaceGridEntityVariant&& other) noexcept            = default;
  ParameterSpaceGridEntityVariant& operator=(const ParameterSpaceGridEntityVariant& other)     = default;
  ParameterSpaceGridEntityVariant& operator=(ParameterSpaceGridEntityVariant&& other) noexcept = default;

  [[nodiscard]] bool equals(const ParameterSpaceGridEntityVariant& other) const {
    return *this == other;
  }

  //! returns true if father entity exists
  [[nodiscard]] bool hasFather() const {
    return visit([](const auto& impl) { return impl.hasFather(); });
  }

  //! Create EntitySeed
  [[nodiscard]] auto seed() const {
    return visit([](const auto& impl) { return impl.seed(); });
  }

  //! Level of this element
  [[nodiscard]] int level() const {
    return visit([](const auto& impl) { return impl.level(); });
  }

  /** @brief The partition type for parallel computing */
  [[nodiscard]] PartitionType partitionType() const {
    return visit([](const auto& impl) { return impl.partitionType(); });
  }

  //! Geometry of this entity
  [[nodiscard]] decltype(auto) geometry() const {
    return visit([](const auto& impl) { return LocalParameterSpaceGeometry(impl.geometry()); });
  }

  /** @brief Return the number of subEntities of codimension codim.
   */
  [[nodiscard]] unsigned int subEntities(unsigned int codim) const {
    return visit([&](const auto& impl) { return impl.subEntities(codim); });
  }

  /** @brief Provide access to sub entity i of given codimension. Entities
   *  are numbered 0 ... subEntities(cc)-1
   */
  template <int cc>
  requires(codim_ == 0)
  [[nodiscard]] decltype(auto) subEntity(int i) const {
    return visit([&](const auto& impl) { return impl.template subEntity<cc>(i); });
  }

  //! First level intersection
  template <typename = void>
  requires(codim_ == 0)
  [[nodiscard]] decltype(auto) ilevelbegin() const {
    return visit([](const auto& impl) { return impl.ilevelbegin(); });
  }

  //! Reference to one past the last neighbor
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ilevelend() const {
    return visit([](const auto& impl) { return impl.ilevelend(); });
  }

  //! First leaf intersection
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ileafbegin() const {
    return visit([](const auto& impl) { return impl.ileafbegin(); });
  }

  //! Reference to one past the last leaf intersection
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) ileafend() const {
    return visit([](const auto& impl) { return impl.ileafend(); });
  }

  //! returns true if Entity has NO children
  template <typename = void>
  requires(codim_ == 0)
  bool isLeaf() const {
    return visit([](const auto& impl) { return impl.isLeaf(); });
  }

  //! Inter-level access to father element on coarser grid.
  //! Assumes that meshes are nested.
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) father() const {
    return visit([](const auto& impl) { return impl.father(); });
  }

  /** @brief Location of this element relative to the reference element element of the father.
   * This is sufficient to interpolate all dofs in conforming case.
   * Nonconforming may require access to neighbors of father and
   * computations with local coordinates.
   * On the fly case is somewhat inefficient since dofs  are visited several times.
   * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
   * implementation of numerical algorithms is only done for simple discretizations.
   * Assumes that meshes are nested.
   */
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) geometryInFather() const {
    return visit([](const auto& impl) { return impl.geometryInFather(); });
  }

  /** @brief Inter-level access to son elements on higher levels<=maxlevel.
   * This is provided for sparsely stored nested unstructured meshes.
   * Returns iterator to first son.
   */
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) hbegin(int maxLevel) const {
    return visit([&](const auto& impl) { return impl.hbegin(maxLevel); });
  }

  //! Returns iterator to one past the last son
  template <typename = void>
  requires(codim_ == 0)
  decltype(auto) hend(int maxLevel) const {
    return visit([&](const auto& impl) { return impl.hend(maxLevel); });
  }

  //! @todo Please doc me !
  template <typename = void>
  requires(codim_ == 0)
  bool wasRefined() const {
    return visit([](const auto& impl) { return impl.wasRefined(); });
  }

  //! @todo Please doc me !
  template <typename = void>
  requires(codim_ == 0)

  bool mightBeCoarsened() const {
    return visit([](const auto& impl) { return impl.mightBeCoarsened(); });
  }

  // /////////////////////////////////////////
  //   Internal stuff
  // /////////////////////////////////////////

  // private:
  const auto& untrimmedHostEntity() const {
    if (std::holds_alternative<HostImplementation>(impl_))
      return std::get<HostImplementation>(impl_);
    else
      return std::get<TrimmedImplementation>(impl_).untrimmedHostEntity();
  }

  std::variant<HostImplementation, TrimmedImplementation> impl_;
};

template <int cd, int dim, int dimworld, typename ScalarType, template <int, int, typename> typename TrimmerType,
          template <int, int, class> class GridEntityVariant>
auto referenceElement(
    const GridEntityVariant<cd, dim, const PatchGrid<dim, dimworld, TrimmerType, ScalarType>>& entity) {
  return TrimmerType<dim, dimworld, ScalarType>::referenceElement(entity);
}

template <int cd, int dim, int dimworld, typename ScalarType, template <int, int, typename> typename TrimmerType,
          template <int, int, class> class GridEntityVariant>
auto referenceElement(
    const Entity<cd, dim, const PatchGrid<dim, dimworld, TrimmerType, ScalarType>, GridEntityVariant>& entity) {
  return referenceElement(entity.impl());
}

} // namespace Dune::IGANEW::Trim
