//
// Created by lex on 27.10.21.
//

#pragma once
#include <dune/iga/concepts.hh>
namespace Dune::IGA {
  template <Vector VT>
  struct ControlPoint {
    using VectorType = VT;
    VectorType p;
    typename VectorType::value_type w{1.0};

    ControlPoint& operator= (const typename VectorType::value_type& a)
    {
      p=a;
      w=a;
      return *this;
    }

    ControlPoint& operator+= (const ControlPoint& cp)
    {
      p+=cp.p;
      w+=cp.w;
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