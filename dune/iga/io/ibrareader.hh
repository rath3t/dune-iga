// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "ibrageometry.hh"
#include "ibrajsonreader.hh"

#include <fstream>
#include <nlohmann/json.hpp>

namespace Dune::IGANEW {

template <int dim, int dimworld, typename PatchGrid>
requires(dim == 2) && (dimworld == 2 || dimworld == 3) && (dim <= dimworld)
class IbraReader
{
  using GridFamily    = typename PatchGrid::GridFamily;
  using PatchData     = NURBSPatchData<dim, dimworld, typename GridFamily::ctype>;
  using PatchTrimData = typename GridFamily::TrimmerTraits::PatchTrimData;
  using TrimmingCurve = typename GridFamily::TrimmerTraits::TrimmingCurve;

  using ControlPointType    = typename PatchData::ControlPointType;
  using ControlPointNetType = typename PatchData::ControlPointNetType;

public:
  static auto read(const std::string& fileName, const bool trim = true,
                   std::array<int, 2> preKnotRefine = {0, 0}) -> std::tuple<PatchData, std::optional<PatchTrimData>> {
    Ibra::Brep brep = readJson<dimworld>(fileName);
    // Each surface in a brep is one Patch, as for now only one brep is supported in NURBSGrid
    assert(brep.surfaces.size() == 1);

    std::cout << "IbraReader successfully read in 1 patch with " << brep.trims.size() << " edges." << std::endl;

    // Reader has done its job, now NURBSGrid can be constructed

    Ibra::Surface<dimworld> _surface                     = brep.surfaces[0];
    const std::array<std::vector<double>, dim> knotSpans = _surface.compileKnotVectors();
    const std::vector<std::vector<ControlPointType>> controlPoints =
        _surface.template transformControlPoints<ControlPointType>();
    std::array<int, dim> dimsize = _surface.n_controlPoints;

    auto controlNet = ControlPointNetType(dimsize, controlPoints);
    PatchData _patchData{knotSpans, controlNet, _surface.degree};

    // Optional preKnot refinement
    for (const auto i : std::views::iota(0, dim)) {
      if (preKnotRefine[i] > 0) {
        auto newKnots = Splines::generateRefinedKnots(knotSpans, i, preKnotRefine[i]);
        _patchData    = Splines::knotRefinement(_patchData, newKnots, i);
      }
    }

    if (trim) {
      PatchTrimData trimData{};
      constructTrimmingCurves(brep, trimData);
      return std::make_tuple(_patchData, trimData);
    }
    return std::make_tuple(_patchData, std::nullopt);
  }

private:
  static void constructTrimmingCurves(const Ibra::Brep<dimworld>& brep, PatchTrimData& trimData) {
    const std::vector<Ibra::BrepLoop> loops = brep.loops;
    for (int i = 0; const Ibra::BrepLoop& loop : loops) {
      trimData.addLoop();
      for (const Ibra::BrepTrim& trim : loop.trims)
        trimData.insertTrimCurve(trim.asCurve<TrimmingCurve>(), i);
      i++;
    }
  }
};

} // namespace Dune::IGANEW
