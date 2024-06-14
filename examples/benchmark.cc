// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/timer.hh>
#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>

using namespace Dune::IGA;

template <typename T>
auto mean(const std::vector<T>& v) -> T {
  T sum = std::accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}

auto runBenchmark(const auto& lambda, int n = 10) -> double {
  std::vector<double> times(n);
  for (const auto i : Dune::range(n)) {
    Dune::Timer timer{};
    lambda();
    times[i] = timer.stop();
  }
  return mean(times);
}

int main(int argc, char** argv) {
  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);

  using PatchGridI   = PatchGrid<2, 2>;
  using GridFactoryI = Dune::GridFactory<PatchGridI>;

  using PatchGridD   = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;
  using GridFactoryD = Dune::GridFactory<PatchGridD>;

  std::string input_file = "../dune/iga/test/auxiliaryfiles/quarter_plate.ibra";

  // Grid Creation

  auto gridFactoryI = GridFactoryI();
  gridFactoryI.insertJson(input_file, false, {8, 8});
  const double timeCreationI = runBenchmark([&]() { gridFactoryI.createGrid(); });

  auto gridFactoryD = GridFactoryD();
  gridFactoryD.insertJson(input_file, false, {8, 8});
  const double timeCreationD = runBenchmark([&]() { gridFactoryD.createGrid(); });

  auto gridFactoryDT = GridFactoryD();
  gridFactoryDT.insertJson(input_file, true, {8, 8});
  const double timeCreationDT = runBenchmark([&]() { gridFactoryDT.createGrid(); }, 3);

  std::cout << "Grid Creation: \n"
               "Identity (Untrimmed): \t\t "
            << timeCreationI
            << "\n"
               "Default (Untrimmed): \t\t "
            << timeCreationD
            << "\n"
               "Default (Trimmed): \t\t\t "
            << timeCreationDT << "\n";

  // Integration Rule generator

  gridFactoryDT.insertJson(input_file, true, {6, 6});
  auto grid6 = gridFactoryDT.createGrid();

  gridFactoryDT.insertJson(input_file, true, {7, 7});
  auto grid7 = gridFactoryDT.createGrid();

  gridFactoryDT.insertJson(input_file, true, {8, 8});
  auto grid8 = gridFactoryDT.createGrid();

  auto makeQR = []<typename G>(G* grid) {
    auto gv = grid->leafGridView();
    for (auto&& ele : elements(gv))
      ele.impl().getQuadratureRule(2 * 2);
  };

  DefaultTrim::Preferences::getInstance().boundaryDivisions(4);
  DefaultTrim::Preferences::getInstance().targetAccuracy(1);

  // auto timeQR6 = runBenchmark([&]() { makeQR(grid6.get()); }, 5);
  // auto timeQR7 = runBenchmark([&]() { makeQR(grid7.get()); }, 5);
  // auto timeQR8 = runBenchmark([&]() { makeQR(grid8.get()); }, 4);

  // std::cout << "Integration Rule Creation (BoundaryDivisions 4): \n"
  //              "6x6: \t\t "
  //           << timeQR6
  //           << "\n"
  //              "7x7:"
  //              " \t\t "
  //           << timeQR7
  //           << "\n"
  //              "8x8: \t\t "
  //           << timeQR8 << "\n";
  //
  // DefaultTrim::Preferences::getInstance().boundaryDivisions(1);
  // DefaultTrim::Preferences::getInstance().targetAccuracy(1e-5);
  //
  // timeQR6 = runBenchmark([&]() { makeQR(grid6.get()); }, 5);
  // timeQR7 = runBenchmark([&]() { makeQR(grid7.get()); }, 5);
  // timeQR8 = runBenchmark([&]() { makeQR(grid8.get()); }, 4);
  //
  // std::cout << "Integration Rule Creation (Target Accuracy 1e-5): \n"
  //              "6x6: \t\t "
  //           << timeQR6
  //           << "\n"
  //              "7x7:"
  //              " \t\t "
  //           << timeQR7
  //           << "\n"
  //              "8x8: \t\t "
  //           << timeQR8 << "\n";

  // Basis creation
  using namespace Dune::Functions::BasisFactory;

  gridFactoryI.insertJson(input_file, false, {8, 8});
  auto grid8I           = gridFactoryI.createGrid();
  auto gridView8I       = grid8I->leafGridView();
  const auto timeBasisI = runBenchmark([&]() { makeBasis(gridView8I, nurbs()); }, 15);

  gridFactoryD.insertJson(input_file, false, {8, 8});
  auto grid8D           = gridFactoryD.createGrid();
  auto gridView8D       = grid8D->leafGridView();
  const auto timeBasisD = runBenchmark([&]() { makeBasis(gridView8D, nurbs()); }, 15);

  auto gridView8DT = grid8->leafGridView();
  auto timeBasisDT = runBenchmark([&]() { makeBasis(gridView8DT, nurbs()); }, 15);

  std::cout << "Basis Creation: \n"
               "Identity (Untrimmed): \t\t "
            << timeBasisI << "\tSize: " << makeBasis(gridView8I, nurbs()).size() << "\n"
            << "Default (Untrimmed): \t\t " << timeBasisD << "\tSize: " << makeBasis(gridView8D, nurbs()).size()
            << "\n"
               "Default (Trimmed): \t\t\t "
            << timeBasisDT << "\tSize: " << makeBasis(gridView8DT, nurbs()).size() << "\n";
}