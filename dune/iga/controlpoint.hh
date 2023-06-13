// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "dune/iga/utils/concepts.hh"
namespace Dune::IGA {

  /** \brief The class which stored a Nurbs ControlPoint
   *
   * @tparam VT The type of coordinates of the controlpoint position
   *
   * The constrolpoint has two public member variables p and w where p is the position and w is the weight of point
   */
  template <Vector VT>
  struct ControlPoint {
    using VectorType = VT;
    VectorType p;
    typename VectorType::value_type w{1.0};

    /** \brief Reset the coordinates and the weight of the controlpoint to zero */
    void setZero() {
      p = typename VectorType::value_type{0};
      w = typename VectorType::value_type{0};
    }

    ControlPoint& operator+=(const ControlPoint& cp) {
      p += cp.p;
      w += cp.w;
      return *this;
    }

    ControlPoint& operator*=(const typename VectorType::value_type& v) {
      p *= v;
      w *= v;
      return *this;
    }
  };

  template <Vector VectorType>
  ControlPoint<VectorType> operator*(const ControlPoint<VectorType>& cp, const typename VectorType::value_type& scal) {
    return {.p = cp.p * scal, .w = cp.w * scal};
  }

  template <Vector VectorType>
  ControlPoint<VectorType> operator*(const typename VectorType::value_type& scal, const ControlPoint<VectorType>& cp) {
    return cp * scal;
  }

  template <Vector VectorType>
  ControlPoint<VectorType> operator+(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
    return {.p = cpL.p + cpR.p, .w = cpL.w + cpR.w};
  }

  template <Vector VectorType>
  ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL) {
    return {.p = -cpL.p, .w = -cpL.w};
  }

  template <Vector VectorType>
  ControlPoint<VectorType> operator-(const ControlPoint<VectorType>& cpL, const ControlPoint<VectorType>& cpR) {
    return {.p = cpL.p - cpR.p, .w = cpL.w - cpR.w};
  }

}  // namespace Dune::IGA
