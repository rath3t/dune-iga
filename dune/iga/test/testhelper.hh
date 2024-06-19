// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <array>
#include <filesystem>
#include <iostream>
#include <ranges>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

template <typename T, int worldDim, int Items>
struct Compare
{
  constexpr bool operator()(const std::array<Dune::FieldVector<double, worldDim>, Items>& lhs,
                            const std::array<Dune::FieldVector<double, worldDim>, Items>& rhs) const {
    return std::ranges::lexicographical_compare(std::ranges::join_view(lhs), std::ranges::join_view(rhs));
  };
};

inline void createOutputFolder(const std::string& folderName) {
  if (!std::filesystem::exists(folderName)) {
    try {
      std::filesystem::create_directory(folderName);
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "Error creating folder: " << e.what() << std::endl;
    }
  }
}

template <typename TestSuiteType, typename ScalarType>
requires std::is_integral_v<ScalarType>
void checkScalars(TestSuiteType& t, const ScalarType val, const ScalarType expectedVal,
                  const std::string& messageIfFailed = "") {
  if constexpr (std::is_integral_v<ScalarType>)
    t.check(val == expectedVal) << std::setprecision(16) << "Incorrect Scalar integer:\t" << expectedVal << " Actual:\t"
                                << val << messageIfFailed;
}

template <typename TestSuiteType, typename ScalarType>
requires(not std::is_integral_v<ScalarType>)
void checkScalars(TestSuiteType& t, const ScalarType val, const ScalarType expectedVal,
                  const std::string& messageIfFailed = "",
                  double tol                         = Dune::FloatCmp::DefaultEpsilon<ScalarType>::value()) {
  t.check(Dune::FloatCmp::eq(val, expectedVal, tol))
      << std::setprecision(16) << "Incorrect Scalar floating point:\t" << expectedVal << " Actual:\t" << val
      << ". The used tolerance was " << tol << messageIfFailed;
}