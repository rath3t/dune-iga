//
// Created by lex on 27.10.21.
//

#pragma once
#include <dune/iga/concepts.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
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
      Dune::IGA::setZero(p);
      w = typename VectorType::value_type{0};
    }

    ControlPoint& operator+=(const ControlPoint& cp) {
      p += cp.p;
      w += cp.w;
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

}  // namespace Dune::IGA