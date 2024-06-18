// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <variant>

namespace Dune::IGA::Trim {

template <typename ParameterSpaceType_, class... Implementations>
class LocalGeometryVariant
{
public:
  using Variant                       = std::variant<Implementations...>;
  using ctype                         = std::common_type_t<typename Implementations::ctype...>;
  using FirstElement                  = std::tuple_element_t<0, std::tuple<Implementations...>>;
  static constexpr int mydimension    = FirstElement::mydimension;
  static constexpr int coorddimension = FirstElement::coorddimension;
  using ParameterSpaceType            = ParameterSpaceType_;
  // using TrimDataType               = DefaultElementTrimData<mydimension, ctype>;

  using LocalCoordinate  = std::common_type_t<typename Implementations::LocalCoordinate...>;
  using GlobalCoordinate = std::common_type_t<typename Implementations::GlobalCoordinate...>;
  // using JacobianTransposed            = std::common_type_t<typename Implementations::JacobianTransposed...>;
  // using Hessian                       = std::common_type_t<typename Implementations::JacobianTransposed>;
  using Jacobian                  = std::common_type_t<typename Implementations::Jacobian...>;
  using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
  using JacobianInverseTransposed = std::common_type_t<typename Implementations::JacobianInverseTransposed...>;
  using JacobianInverse           = std::common_type_t<typename Implementations::JacobianInverse...>;
  using Volume                    = ctype;

  auto visit(auto&& lambda) const {
    return std::visit(lambda, impl_);
  }

  template <class Implementation>
  LocalGeometryVariant(const Implementation& impl)
      : impl_(impl) {}

  LocalGeometryVariant()                                  = default;
  LocalGeometryVariant(const LocalGeometryVariant& other) = default;
  template <class Implementation>
  requires(!std::is_same_v<Implementation, LocalGeometryVariant>)
  LocalGeometryVariant& operator=(const Implementation& impl) {
    impl_ = impl;
    return *this;
  };
  LocalGeometryVariant(LocalGeometryVariant&& other) noexcept            = default;
  LocalGeometryVariant& operator=(const LocalGeometryVariant& other)     = default;
  LocalGeometryVariant& operator=(LocalGeometryVariant&& other) noexcept = default;

  /** @brief Return the element type identifier
   */
  [[nodiscard]] GeometryType type() const {
    return visit([](const auto& impl) { return impl.type(); });
  }

  // return whether we have an affine mapping
  [[nodiscard]] bool affine() const {
    return visit([](const auto& impl) { return impl.affine(); });
  }

  // return the number of corners of this element. Corners are numbered 0...n-1
  [[nodiscard]] int corners() const {
    return visit([](const auto& impl) { return impl.corners(); });
  }

  // access to coordinates of corners. Index is the number of the corner
  auto corner(int i) const {
    return visit([&](const auto& impl) { return impl.corner(i); });
  }

  /** @brief Maps a local coordinate within reference element to
   * global coordinate in element  */
  auto global(const auto& local) const {
    return visit([&](const auto& impl) { return impl.global(local); });
  }

  auto center() const {
    return visit([](const auto& impl) { return impl.center(); });
  }

  /** @brief Return the transposed of the Jacobian
   */
  JacobianTransposed jacobianTransposed(const auto& local) const {
    return visit([&](const auto& impl) -> JacobianTransposed { return impl.jacobianTransposed(local); });
  }

  /** @brief Maps a global coordinate within the element to a
   * local coordinate in its reference element */
  auto local(const auto& global) const {
    return visit([&](const auto& impl) { return impl.local(global); });
  }

  // Returns true if the point is in the current element
  bool checkInside(const auto& local) const {
    return visit([&](const auto& impl) { return impl.checkInside(local); });
  }

  [[nodiscard]] auto integrationElement(const auto& local) const {
    return visit([&](const auto& impl) { return impl.integrationElement(local); });
  }

  // The Jacobian matrix of the mapping from the reference element to this element
  [[nodiscard]] auto jacobianInverseTransposed(const auto& local) const {
    return visit([&](const auto& impl) { return impl.jacobianInverseTransposed(local); });
  }

  auto getQuadratureRule(const std::optional<int>& p_order = std::nullopt,
                         const QuadratureType::Enum qt     = QuadratureType::GaussLegendre) const {
    // The parameterSpace has linear geometry, so if no order is specified, we assume it to be myDimension
    return visit([&](const auto& impl) {
      if constexpr (requires { impl.getQuadratureRule(p_order, qt); })
        return impl.getQuadratureRule(p_order, qt);
      else
        return QuadratureRules<double, mydimension>::rule(type(), p_order.value_or(mydimension), qt);
    });
  }

  bool isTrimmed() const {
    return visit([&](const auto& impl) {
      if constexpr (requires { impl.isTrimmed(); })
        return impl.isTrimmed();
      return false;
    });
  }

private:
  Variant impl_;
};
} // namespace Dune::IGA::Trim
