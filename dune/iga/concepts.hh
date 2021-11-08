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

//  template< template<std::floating_point,int,int> typename MatrixType, typename ScalarType, int rows, int cols>
//  concept FixedMatrix =  Matrix<MatrixType<ScalarType,rows,cols>>;
//
//  template < class > struct checkFixedMatrixtemplate : std::false_type {};
//
//  // Specialize for template classes
//  template <template<std::floating_point,int,int> typename MatrixType, typename ScalarType, int rows, int cols>
//  struct checkFixedMatrixtemplate< MatrixType<ScalarType,rows,cols> > : std::true_type {};


  template< typename LinearAlgebraTraits, int a=1>
  concept NurbsGridLinearAlgebra = Matrix<typename LinearAlgebraTraits::JacobianTransposedType> &&
      Matrix<typename LinearAlgebraTraits::JacobianInverseTransposed> &&
      Vector<typename LinearAlgebraTraits::GlobalCoordinateType> &&
      Vector<typename LinearAlgebraTraits::LocalCoordinateType>
     && requires()
  {
    typename LinearAlgebraTraits::JacobianTransposedType;
    typename LinearAlgebraTraits::JacobianInverseTransposed;
    typename LinearAlgebraTraits::GlobalCoordinateType;
    typename LinearAlgebraTraits::LocalCoordinateType;
    typename LinearAlgebraTraits::value_type;

    typename LinearAlgebraTraits::template FixedMatrixType<a,a>;


  };


}
