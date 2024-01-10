// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <fstream>
#include <nlohmann/json.hpp>

#include "ibrageometry.hh"

namespace Dune::IGANEW {

  template <int dim, int dimworld, typename PatchGrid>
  requires(dim == 2) && (dimworld == 2 || dimworld == 3) && (dim <= dimworld) class IbraReader {
    using GridFamily    = typename PatchGrid::GridFamily;
    using PatchData     = NURBSPatchData<dim, dimworld, typename GridFamily::ctype>;
    using PatchTrimData = typename GridFamily::TrimmerTraits::PatchTrimData;
    using TrimmingCurve = typename GridFamily::TrimmerTraits::TrimmingCurve;

    using ControlPointType    = typename PatchData::ControlPointType;
    using ControlPointNetType = typename PatchData::ControlPointNetType;

  public:
    static auto read(const std::string& fileName, const bool trim = true, std::array<int, 2> preKnotRefine = {0, 0}) {
      std::ifstream ibraInputFile;
      ibraInputFile.open(fileName);
      return read(ibraInputFile, trim, preKnotRefine);
    }

    template <typename InputStringType>
    requires(
        not std::convertible_to<
            std::string,
            InputStringType> and not std::convertible_to<InputStringType, const char*>) static auto read(InputStringType&
                                                                                                             ibraInputFile,
                                                                                                         const bool trim
                                                                                                         = true,
                                                                                                         std::array<int,
                                                                                                                    2>
                                                                                                             preKnotRefine
                                                                                                         = {0, 0}) {
      using json = nlohmann::json;

      std::vector<Ibra::Surface<dimworld>> surfaces;
      std::vector<Ibra::Curve2D> curves2D;

      std::vector<Ibra::BrepRepresentation> brepRepresentations;
      std::vector<Ibra::BrepLoopRepresentation> brepLoopRepresentations;
      std::vector<Ibra::BrepTrimRepresentation> brepTrimRepresentations;

      json ibraJson;
      try {
        ibraJson = json::parse(ibraInputFile);

        for (auto& j : ibraJson) {
          auto geo = j.get<Ibra::IbraBase>();

          switch (geo.type) {
            case Ibra::Type::NurbsSurfaceGeometry3D:
              surfaces.push_back(j.get<Ibra::Surface<dimworld>>());
              break;
            case Ibra::Type::NurbsCurveGeometry2D:
              curves2D.push_back(j.get<Ibra::Curve2D>());
              break;
            case Ibra::Type::BrepType:
              brepRepresentations.push_back(j.get<Ibra::BrepRepresentation>());
              break;
            case Ibra::Type::BrepLoopType:
              brepLoopRepresentations.push_back(j.get<Ibra::BrepLoopRepresentation>());
              break;
            case Ibra::Type::BrepTrimType:
              brepTrimRepresentations.push_back(j.get<Ibra::BrepTrimRepresentation>());
              break;
            default:
              // Do nothing -- keep the compiler happy
              break;
          }
        }
      } catch (json::parse_error& ex) {
        DUNE_THROW(Dune::IOError, "Error in reading input stream: "
                                      << ", parse error at byte " << ex.byte << " What: " << ex.what());
      }

      // Make Connections
      auto brepRepresentation = brepRepresentations[0];

      Ibra::Brep brep{brepRepresentation, curves2D, surfaces, brepLoopRepresentations, brepTrimRepresentations};
      // Each surface in a brep is one Patch, as for now only one brep is supported in NURBSGrid
      assert(brep.surfaces.size() == 1);

      std::cout << "IbraReader successfully read in 1 patch with " << brep.trims.size() << " edges." << std::endl;

      // Reader has done its job, now NURBSGrid can be constructed

      Ibra::Surface<dimworld> _surface                     = brep.surfaces[0];
      const std::array<std::vector<double>, dim> knotSpans = _surface.compileKnotVectors();
      const std::vector<std::vector<ControlPointType>> controlPoints
          = _surface.template transformControlPoints<ControlPointType>();
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

      PatchTrimData trimData{};
      if (trim) constructTrimmingCurves(brep, trimData);

      return std::make_tuple(_patchData, trimData);
    }

  private:
    static void constructTrimmingCurves(const Ibra::Brep<dimworld>& brep, PatchTrimData& trimData) {
      const std::vector<Ibra::BrepLoop> loops = brep.loops;
      assert(!loops.empty() && "Only one boundary loop is currently supported");
      for (int i = 0; const Ibra::BrepLoop& loop : loops) {
        trimData.addLoop();
        for (const Ibra::BrepTrim& trim : loop.trims)
          trimData.insertTrimCurve(trim.asCurve<TrimmingCurve>(), i);
        i++;
      }
    }
  };

}  // namespace Dune::IGANEW
