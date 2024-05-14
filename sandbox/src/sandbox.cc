// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);

  std::cout << "Hello this is dune-iga" << std::endl;

  return 0;
}
