// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once


#include <fstream>
#include <dune/iga/trimmer/identitytrimmer/trimmer.hh>
#include <nlohmann/json.hpp>

#include "ibrageometry.hh"


namespace Dune::IGANEW {

  template <int dim, int dimworld, typename PatchGrid>
  class IbraReader {

    using GridFamily = typename PatchGrid::GridFamily;
    using PatchData = NURBSPatchData<dim, dimworld, typename GridFamily::ctype>;
    using PatchTrimData = typename GridFamily::TrimmerTraits::PatchTrimData;

    using ControlPoint        = typename PatchData::ControlPointType;
    using ControlPointNetType = typename PatchData::ControlPointNetType;
  public:
    static auto read(const std::string& fileName) {
      std::ifstream ibraInputFile;
      ibraInputFile.open(fileName);
      return read(ibraInputFile);
    }

    template <typename InputStringType>
    requires(not std::convertible_to<
                 std::string, InputStringType> and not std::convertible_to<InputStringType, const char*>) static
    auto read(InputStringType& ibraInputFile, const bool trim = true,
                              std::array<int, 2> elevateDegree = {0, 0}, std::array<int, 2> preKnotRefine = {0, 0},
                              std::array<int, 2> postKnotRefine = {0, 0}) {
      using json = nlohmann::json;

      // Result
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

      Ibra::Surface<dimworld> _surface                           = brep.surfaces[0];
      const std::array<std::vector<double>, dim> knotSpans   = _surface.compileKnotVectors();
      const std::vector<std::vector<ControlPoint>> controlPoints = _surface.transformControlPoints();
      std::array<int, dim> dimsize                           = _surface.n_controlPoints;

      auto controlNet = ControlPointNetType(dimsize, controlPoints);
      PatchData _patchData{knotSpans, controlNet, _surface.degree};


      // Make the trimData and pass them into the grid
      // So the grid can figure out what to do with it
      // auto trimData = constructGlobalBoundaries(brep);

      return std::make_tuple(_patchData, PatchTrimData());
    }

   private:
    // static auto constructGlobalBoundaries(const Ibra::Brep<worldDim>& brep) -> std::shared_ptr<TrimData> {
    //   auto data = std::make_shared<TrimData>();
    //   std::vector<Boundary> boundaries;
    //   std::vector<Ibra::BrepLoop> loops = brep.loops;
    //   for (Ibra::BrepLoop& loop : loops) {
    //     boundaries.clear();
    //     for (Ibra::BrepTrim& trim : loop.trims) {
    //       boundaries.emplace_back(trim);
    //     }
    //     data->addLoop(boundaries);
    //   }
    //
    //   return data;
    // }
  };

}  // namespace Dune::IGA
