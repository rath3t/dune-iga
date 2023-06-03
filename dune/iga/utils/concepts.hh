// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <concepts>

namespace Dune::IGA {

  template <typename VectorType>
  concept Vector = requires(VectorType v, double a, int index) {
    typename VectorType::value_type;
    { v[index] } -> std::same_as<typename VectorType::value_type&>;
    { v += v } -> std::same_as<VectorType&>;
    { v -= v } -> std::same_as<VectorType&>;
    { v *= a } -> std::same_as<VectorType&>;
    { v /= a } -> std::same_as<VectorType&>;
  };

  template <typename MatrixType>
  concept Matrix = requires(MatrixType A, double a) {
    typename MatrixType::value_type;
    { A += A } -> std::same_as<MatrixType&>;
    { A -= A } -> std::same_as<MatrixType&>;
    { A *= a } -> std::same_as<MatrixType&>;
    { A /= a } -> std::same_as<MatrixType&>;
  };

  template <typename ControlPointType>
  concept ControlPointConcept = Vector<typename ControlPointType::VectorType> && requires(ControlPointType cp) {
    typename ControlPointType::VectorType;
    typename ControlPointType::VectorType::value_type;
    { cp.p } -> std::same_as<typename ControlPointType::VectorType&>;
    { cp.w } -> std::same_as<typename ControlPointType::VectorType::value_type&>;
  };  // namespace Dune::IGA

  template <typename L, typename R>
  concept MultiplyAble = requires(L x, R y) {
    x* y;
  };

  template <typename L, typename R>
  concept AddAble = requires(L x, R y) {
    x + y;
  };

  template <typename L, typename R>
  concept SubstractAble = requires(L x, R y) {
    x - y;
  };

  template <typename L, typename R>
  concept MultiplyAssignAble = requires(L x, R y) {
    x *= y;
  };

  template <typename L, typename R>
  concept DivideAssignAble = requires(L x, R y) {
    x /= y;
  };

  template <typename L, typename R>
  concept DivideAble = requires(L x, R y) {
    x / y;
  };

  template <typename V>
  concept StdVectorLikeContainer = requires(V v, int a) {
    typename V::value_type;
    { v.resize(a) } -> std::same_as<void>;
    { v.back() } -> std::same_as<typename V::value_type&>;
    { v.front() } -> std::same_as<typename V::value_type&>;
  };

}  // namespace Dune::IGA
