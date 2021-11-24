//
// Created by lex on 27.10.21.
//

#pragma once

#include <span>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune::IGA {
  template <typename ScalarType, std::integral auto dim, std::integral auto dimworld>
  struct LinearAlgebraTraits {
    using GlobalCoordinateType      = FieldVector<ScalarType, dimworld>;
    using LocalCoordinateType       = FieldVector<ScalarType, dim>;
    using JacobianTransposedType    = FieldMatrix<ScalarType, dim, dimworld>;
    using JacobianInverseTransposed = FieldMatrix<ScalarType, dimworld, dim>;
    using value_type                = typename FieldMatrix<ScalarType, dimworld, dim>::value_type;
    template <int rows, int cols>
    using FixedMatrixType = FieldMatrix<ScalarType, rows, cols>;
    template <int rows>
    using FixedVectorType = FieldVector<ScalarType, rows>;
  };
  }  // namespace Dune::IGA