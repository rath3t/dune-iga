//
// Created by Henri on 11.03.2023.
//
#pragma once

#ifndef DUNE_IGA_IBRAREADER_HH
#  define DUNE_IGA_IBRAREADER_HH

#  include <clipper2/clipper.h>
#  include <nlohmann/json.hpp>

#  include <dune/iga/ibraGeometry.hh>
#  include <dune/iga/nurbsgrid.hh>
#  include <dune/iga/nurbstrimboundary.hh>

namespace Dune::IGA {

  // At the moment only gridDim und worldDim == 2 or 3 supported
  template <int gridDim, int worldDim>
    requires(gridDim == 2) && (worldDim == 2 || worldDim == 3) && (gridDim <= worldDim)
  class IbraReader {
   public:
    using Grid                = Dune::IGA::NURBSGrid<gridDim, worldDim>;
    using ControlPoint        = Dune::IGA::NURBSPatchData<gridDim, worldDim>::ControlPointType;
    using PatchData           = Dune::IGA::NURBSPatchData<gridDim, worldDim>;
    using ControlPointNetType = Dune::IGA::MultiDimensionNet<gridDim, ControlPoint>;

    static std::shared_ptr<Grid> read(const std::string& fileName) {
      // Get Standard Parameter for trimming
      Dune::IGA::Utilities::setStandardParameters();

      using json = nlohmann::json;

      // Result
      std::vector<Ibra::Surface<worldDim>> surfaces;
      std::vector<Ibra::Curve2D> curves2D;

      std::vector<Ibra::BrepRepresentation> brepRepresentations;
      std::vector<Ibra::BrepLoopRepresentation> brepLoopRepresentations;
      std::vector<Ibra::BrepTrimRepresentation> brepTrimRepresentations;

      std::ifstream ibraInputFile;
      ibraInputFile.open(fileName);

      json ibraJson;
      try {
        ibraJson = json::parse(ibraInputFile);

        for (auto& j : ibraJson) {
          auto geo = j.get<Ibra::IbraBase>();

          switch (geo.type) {
            case Ibra::Type::NurbsSurfaceGeometry3D:
              surfaces.push_back(j.get<Ibra::Surface<worldDim>>());
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
        DUNE_THROW(Dune::IOError, "Error in file: " << fileName << ", parse error at byte " << ex.byte);
      }

      // Make Connections
      auto brepRepresentation = brepRepresentations[0];

      Ibra::Brep brep{brepRepresentation, curves2D, surfaces, brepLoopRepresentations, brepTrimRepresentations};
      // Each surface in a brep is one Patch, as for now only one brep is supported in NURBSGrid
      assert(brep.surfaces.size() == 1);

      std::cout << "IbraReader successfully read in 1 patch with " << brep.trims.size() << " edges." << std::endl;

      // Reader has done its job, now NURBSGrid can be constructed

      Ibra::Surface<worldDim> _surface = brep.surfaces[0];

      const std::array<std::vector<double>, gridDim> knotSpans   = _surface.compileKnotVectors();
      const std::vector<std::vector<ControlPoint>> controlPoints = _surface.transformControlPoints();

      std::array<int, gridDim> dimsize = _surface.n_controlPoints;

      auto controlNet = ControlPointNetType(dimsize, controlPoints);
      PatchData _patchData {knotSpans, controlNet, _surface.degree};

      // Make the trimData, and pass them into the grid
      // So the grid can figure out what to do with it
      auto trimData = constructGlobalBoundaries(brep);

      auto grid = std::make_unique<Grid>(_patchData, trimData);
      return std::move(grid);
    }

    // Deleted Default Constructor, use static method read
    IbraReader() = delete;

   private:
    static auto constructGlobalBoundaries(Ibra::Brep<worldDim>& brep) -> std::shared_ptr<TrimData> {
      TrimData data{};

      std::vector<Boundary> boundaries;
      auto loops = brep.loops;
      for (Ibra::BrepLoop& loop : loops) {
        boundaries.clear();
        for (auto& trim : loop.trims) {
          boundaries.emplace_back(trim);
        }
        data.addLoop(boundaries);
      }

      return std::make_shared<TrimData>(data);
    }
  };

}  // namespace Dune::IGA

#endif  // DUNE_IGA_IBRAREADER_HH
