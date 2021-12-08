//
// Created by lex on 27.10.21.
//

#pragma once

#include <span>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune::IGA {

  /** \brief Traits which satisfy the concept LinearAlgebra with dune types */
  template <typename ScalarType>
  struct DuneLinearAlgebraTraits {
    using value_type = ScalarType;
    template <int rows, int cols>
    using FixedMatrixType = FieldMatrix<ScalarType, rows, cols>;
    template <int rows>
    using FixedVectorType   = FieldVector<ScalarType, rows>;
    using DynamicMatrixType = DynamicMatrix<ScalarType>;
    template <int cols>
    using RowFixedMatrix = DynamicMatrix<ScalarType>;
    template <int rows>
    using ColumnFixedMatrix = DynamicMatrix<ScalarType>;
    using DynamicVectorType = DynamicVector<ScalarType>;
  };

  template <int rows, typename ScalarType>
  inline auto two_norm(const Dune::FieldVector<ScalarType, rows>& v) {
    return v.two_norm();
  }

  template <typename ScalarType>
  inline auto two_norm(const Dune::DynamicVector<ScalarType>& v) {
    return v.two_norm();
  }

  template <int rows, typename ScalarType>
  inline void setZero(Dune::FieldVector<ScalarType, rows>& v) {
    v = ScalarType(0);
  }
}  // namespace Dune::IGA
