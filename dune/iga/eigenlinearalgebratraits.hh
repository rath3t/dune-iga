//
// Created by lex on 27.10.21.
//

#pragma once


#include <Eigen/Core>

namespace Dune::IGA {

  /** \brief Traits which satisfy the concept LinearAlgebra with dune types */
  template <typename ScalarType>
  struct EigenLinearAlgebraTraits {
    using value_type = ScalarType;
    template <int rows, int cols>
    using FixedMatrixType = Eigen::Matrix<ScalarType, rows, cols>;
    template <int rows>
    using FixedVectorType   =  Eigen::Vector<ScalarType, rows>;
    using DynamicMatrixType = Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>;
    template <int cols>
    using RowFixedMatrix = Eigen::Matrix<ScalarType,cols,Eigen::Dynamic>;
    template <int rows>
    using ColumnFixedMatrix =  Eigen::Matrix<ScalarType,Eigen::Dynamic,rows>;
    using DynamicVectorType = Eigen::Vector<ScalarType,Eigen::Dynamic>;
  };

  template <int rows, typename ScalarType>
  inline auto two_norm(const Eigen::Vector<ScalarType, rows>& v) {
    return v.norm();
  }

  template <int rows, typename ScalarType>
  inline auto dot(const Eigen::Vector<ScalarType, rows>& v1, const Eigen::Vector<ScalarType, rows>& v2) {
    return v1.dot(v2);
  }

  template <int rows, typename ScalarType>
  inline void setZero(Eigen::Vector<ScalarType, rows>& v) {
    v.setZero();
  }
}  // namespace Dune::IGA
