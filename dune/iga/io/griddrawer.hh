// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/parameterspace/default/parameterspace.hh>

namespace Dune::IGA {

template <typename PatchGrid>
void drawGrid(PatchGrid* grid, std::string&& file_name) {
  const typename PatchGrid::ParameterSpace& parameterspace = grid->parameterSpace();
  auto& nonConstParameterSpace = const_cast<typename PatchGrid::ParameterSpace&>(parameterspace);

  auto eleTrimDatas = nonConstParameterSpace.trimElements();

  auto figure = matplot::figure(true);
  figure->size(1000, 1000);

  for (auto& eleTrimData : eleTrimDatas)
    eleTrimData.drawResult("out/ele", false, false);

  matplot::save(file_name, "gif");
}

} // namespace Dune::IGA