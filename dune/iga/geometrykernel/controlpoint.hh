// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/iga/hierarchicpatch/concepts.hh>

namespace Dune::IGANEW {

/**
 * @brief Represents a NURBS control point.
 *
 * This struct stores the coordinates and weight of a NURBS control point. The coordinates are
 * represented by the position vector `p` of type `VectorType`, and the weight is represented by
 * the scalar `w`.
 *
 * @tparam VT The type of coordinates of the control point position.
 * @note The vector type must satisfy the Concept::Vector concept.
 */
template <Concept::Vector VT>
struct ControlPoint
{
  using VectorType = VT;
  VectorType p;
  typename VectorType::value_type w{1.0};

  /**
   * @brief Resets the coordinates and the weight of the control point to zero.
   */
  void setZero() {
    p = typename VectorType::value_type{0};
    w = typename VectorType::value_type{0};
  }

  /**
   * @brief Adds the coordinates and the weights of two control points.
   *
   * @param cp The control point to be added.
   * @return A reference to the modified control point.
   */
  ControlPoint& operator+=(const ControlPoint& cp) {
    p += cp.p;
    w += cp.w;
    return *this;
  }

  /**
   * @brief Multiplies and assigns the coordinates and the weight by a scalar.
   *
   * @param v The scalar value for multiplication.
   * @return A reference to the modified control point.
   */
  ControlPoint& operator*=(const typename VectorType::value_type& v) {
    p *= v;
    w *= v;
    return *this;
  }
};

/**
 * @brief Multiplies the coordinates and the weight of a control point by a scalar from the right.
 *
 * @tparam VectorType The vector type satisfying the Concept::Vector concept.
 * @param cp The control point.
 * @param scal The scalar value for multiplication.
 * @return A new control point with the scaled coordinates and weight.
 */
template <Concept::Vector VectorType>
ControlPoint<VectorType> operator*(const ControlPoint<VectorType>& cp, const typename VectorType::value_type& scal) {
  return {.p = cp.p * scal, .w = cp.w * scal};
}

/**
 * @brief Multiplies the coordinates and the weight of a control point by a scalar from the left.
 *
 * @tparam VectorType The vector type satisfying the Concept::Vector concept.
 * @param scal The scalar value for multiplication.
 * @param cp The control point.
 * @return A new control point with the scaled coordinates and weight.
 */
template <Concept::Vector VectorType>
ControlPoint<VectorType> operator*(const typename VectorType::value_type& scal, const ControlPoint<VectorType>& cp) {
  return cp * scal;
}

/**
 * @brief Adds the coordinates and the weight of two control points.
 *
 * @tparam VectorType The vector type satisfying the Concept::Vector concept.
 * @param cpL The left control point.
 * @param cpR The right control point.
 * @return A new control point with the summed coordinates and weight.
 */
template <Concept::Vector VectorType>
ControlPoint<VectorType> operator+(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
  return {.p = cpL.p + cpR.p, .w = cpL.w + cpR.w};
}

/**
 * @brief Negates the coordinates and the weight of a control point.
 *
 * @tparam VectorType The vector type satisfying the Concept::Vector concept.
 * @param cpL The control point to be negated.
 * @return A new control point with negated coordinates and weight.
 */
template <Concept::Vector VectorType>
ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL) {
  return {.p = -cpL.p, .w = -cpL.w};
}

/**
 * @brief Subtracts the coordinates and the weight of two control points.
 *
 * @tparam VectorType The vector type satisfying the Concept::Vector concept.
 * @param cpL The left control point.
 * @param cpR The right control point.
 * @return A new control point with the subtracted coordinates and weight.
 */
template <Concept::Vector VectorType>
ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
  return {.p = cpL.p - cpR.p, .w = cpL.w - cpR.w};
}

} // namespace Dune::IGANEW
