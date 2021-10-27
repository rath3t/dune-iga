//
// Created by lex on 27.10.21.
//

#pragma once

#include <concepts>

namespace Dune::IGA {


  template< typename VectorType>
  concept Vector =  requires(VectorType v, double a, std::size_t index) {
    typename VectorType::value_type;
    { a*v } -> std::same_as<VectorType >;
    { v*a } -> std::same_as<VectorType >;
    { v+=v } -> std::same_as<VectorType&>;
    { v-=v } -> std::same_as<VectorType&>;
    { v*=a } -> std::same_as<VectorType&>;
    { v/=a } -> std::same_as<VectorType&>;
  };

  template< typename MatrixType>
  concept Matrix =  requires(MatrixType A, double a, int index) {
    typename MatrixType::value_type;
    { a*A } -> std::same_as<MatrixType >;
    { A*a } -> std::same_as<MatrixType >;
    { A+=A } -> std::same_as<MatrixType&>;
    { A-=A } -> std::same_as<MatrixType&>;
    { A*=a } -> std::same_as<MatrixType&>;
    { A/=a } -> std::same_as<MatrixType&>;
  };

  template< typename LinearAlgebraTraits>
  concept NurbsGridLinearAlgebra = Matrix<typename LinearAlgebraTraits::JacobianTransposedType> &&
      Matrix<typename LinearAlgebraTraits::JacobianInverseTransposed> &&
      Vector<typename LinearAlgebraTraits::GlobalCoordinateType> &&
      Vector<typename LinearAlgebraTraits::LocalCoordinateType> && requires()
  {
    typename LinearAlgebraTraits::JacobianTransposedType;
    typename LinearAlgebraTraits::JacobianInverseTransposed;
    typename LinearAlgebraTraits::GlobalCoordinateType;
    typename LinearAlgebraTraits::LocalCoordinateType;
  };


}
