// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

namespace Dune::IGANEW {

template <typename PatchGrid>
void drawGrid(PatchGrid* grid, std::string&& file_name) {
  const typename PatchGrid::Trimmer& trimmer = grid->trimmer();
  auto eleTrimDatas                          = trimmer.trimElements();

  auto figure = matplot::figure(true);
  figure->size(1000, 1000);

  for (auto& eleTrimData : eleTrimDatas)
    eleTrimData.drawResult("resName", false, false);

  matplot::save(file_name, "gif");
}

} // namespace Dune::IGANEW