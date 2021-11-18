//
// Created by lex on 27.10.21.
//

#pragma once

#include <span>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune::IGA {
  template <typename ScalarType, std::integral auto  dim, std::integral auto  dimworld>
  struct LinearAlgebraTraits {
    using GlobalCoordinateType      = FieldVector<ScalarType, dimworld>;
    using LocalCoordinateType       = FieldVector<ScalarType, dim>;
    using JacobianTransposedType    = FieldMatrix<ScalarType, dim, dimworld>;
    using JacobianInverseTransposed = FieldMatrix<ScalarType, dimworld, dim>;
    using value_type = typename FieldMatrix<ScalarType, dimworld, dim>::value_type;
    template<int rows,int cols>
    using FixedMatrixType = FieldMatrix<ScalarType, rows, cols>;
    template<int rows>
    using FixedVectorType = FieldVector<ScalarType, rows>;

  };

  template<std::floating_point ScalarType,template<std::floating_point> typename DynamicType,template<std::floating_point, int> typename FixedType, int degreeT= -1>
  struct fixedOrDynamic
  {
    static constexpr long unsigned int degree = (degreeT==-1 ) ? 0 : degreeT;
    using type = typename std::conditional<degreeT==-1,DynamicType<ScalarType>,FixedType<ScalarType,degree>>::type;
  };

  template<std::floating_point ScalarType,int degreeT= -1>
  using fixedOrDynamicSpan = typename fixedOrDynamic<ScalarType,std::span,std::span,degreeT>::type;


}