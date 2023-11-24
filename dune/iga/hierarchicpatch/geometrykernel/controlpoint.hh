// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/iga/hierarchicpatch/concepts.hh>
namespace Dune::IGANEW {

  /** \brief The class which stores a Nurbs ControlPoint
   * \tparam VT The type of coordinates of the controlpoint position
   * The constrolpoint has two public member variables p and w where p is the position and w is the weight of point
   */
  template <Concept::Vector VT>
  struct ControlPoint {
    using VectorType = VT;
    VectorType p;
    typename VectorType::value_type w{1.0};

    /** \brief Reset the coordinates and the weight of the controlpoint to zero */
    void setZero() {
      p = typename VectorType::value_type{0};
      w = typename VectorType::value_type{0};
    }

    /** \brief Adds the coordinates and the weights of two control points */
    ControlPoint& operator+=(const ControlPoint& cp) {
      p += cp.p;
      w += cp.w;
      return *this;
    }

    /** \brief Multiplies and assign the coordinates and the weight by a Scalar */
    ControlPoint& operator*=(const typename VectorType::value_type& v) {
      p *= v;
      w *= v;
      return *this;
    }
  };

  /** \brief Multiplies the coordinates and the weight by a Scalar from the right */
  template <Concept::Vector VectorType>
  ControlPoint<VectorType> operator*(const ControlPoint<VectorType>& cp, const typename VectorType::value_type& scal) {
    return {.p = cp.p * scal, .w = cp.w * scal};
  }

  /** \brief Multiplies the coordinates and the weight by a Scalar from the left */
  template <Concept::Vector VectorType>
  ControlPoint<VectorType> operator*(const typename VectorType::value_type& scal, const ControlPoint<VectorType>& cp) {
    return cp * scal;
  }

  /** \brief Adds the coordinates and the weight  of two controlpoints */
  template <Concept::Vector VectorType>
  ControlPoint<VectorType> operator+(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
    return {.p = cpL.p + cpR.p, .w = cpL.w + cpR.w};
  }

  /** \brief Negates the coordinates and the weight  of a control points */
  template <Concept::Vector VectorType>
  ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL) {
    return {.p = -cpL.p, .w = -cpL.w};
  }

  /** \brief Substracts the coordinates and the weight  of two control points */
  template <Concept::Vector VectorType>
  ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
    return {.p = cpL.p - cpR.p, .w = cpL.w - cpR.w};
  }

}  // namespace Dune::IGANEW
