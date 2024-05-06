// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <array>
#include <filesystem>
#include <iostream>
#include <ranges>

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