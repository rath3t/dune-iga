//
// Created by lex on 27.10.21.
//

#pragma once

namespace Dune::IGA {
  template <typename ScalarType, int dim, int dimworld>
  struct LinearAlgebraTraits {
    using GlobalCoordinateType      = FieldVector<ScalarType, dimworld>;
    using LocalCoordinateType       = FieldVector<ScalarType, dim>;
    using JacobianTransposedType    = FieldMatrix<ScalarType, dim, dimworld>;
    using JacobianInverseTransposed = FieldMatrix<ScalarType, dimworld, dim>;
  };
}