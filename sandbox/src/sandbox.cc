// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "config.h"

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);

  std::cout << "Hello this is dune-iga" << std::endl;

  using Grid = Dune::YaspGrid<2>;

  auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({0, 0}, {1, 1}, {1, 1});

  for (const auto& v : Dune::vertices(grid->leafGridView())) {
    auto vGeo = v.geometry();

    std::cout << "Center: " << vGeo.center() << std::endl;
    std::cout << "Corner: " << vGeo.corner(0) << std::endl;
    std::cout << "Local: " << vGeo.local({}) << std::endl;
    std::cout << "Global: " << vGeo.global({}) << std::endl;
    std::cout << std::endl;
  }

  return 0;
}
