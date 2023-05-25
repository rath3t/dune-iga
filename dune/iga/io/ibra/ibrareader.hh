//
// Created by Henri on 11.03.2023.
//
#pragma once

#include "ibrageometry.hh"

#include <clipper2/clipper.h>
#include <nlohmann/json.hpp>

#include "dune/iga/nurbsalgorithms.hh"
#include "dune/iga/nurbsgrid.hh"
#include "dune/iga/trim/nurbstrimboundary.hh"

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

    static std::shared_ptr<Grid> read(const std::string& fileName, const bool trim = true,
                                      std::array<int, 2> elevateDegree = {0, 0}, std::array<int, 2> preKnotRefine = {0, 0}) {
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
        DUNE_THROW(Dune::IOError,
                   "Error in file: " << fileName << ", parse error at byte " << ex.byte << " What: " << ex.what());
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
      PatchData _patchData{knotSpans, controlNet, _surface.degree};

      // Optional Pre Refinement of knots
//      for (int i = 0; i < gridDim; ++i) {
//        if (preKnotRefine[i] > 0) {
//          auto additionalKnots = generateRefinedKnots(knotSpans, i, preKnotRefine[i]);
//          _patchData = knotRefinement<gridDim>(_patchData, additionalKnots, i);
//        }
//      }

      // Optional Degree Elevate
      for (int i = 0; i < gridDim; ++i)
        if (elevateDegree[i] > 0) _patchData = degreeElevate(_patchData, i, elevateDegree[i]);

      // Make the trimData, and pass them into the grid
      // So the grid can figure out what to do with it
      auto trimData = constructGlobalBoundaries(brep);

      auto grid = (trim) ? std::make_unique<Grid>(_patchData, trimData) : std::make_unique<Grid>(_patchData);
      return std::move(grid);
    }

    // Deleted Default Constructor, use static method read
    IbraReader() = delete;

   private:

    static auto constructGlobalBoundaries(const Ibra::Brep<worldDim>& brep) -> std::shared_ptr<TrimData> {
      auto data =std::make_shared<TrimData>();
      std::vector<Boundary> boundaries;
      std::vector<Ibra::BrepLoop> loops = brep.loops;
      for (Ibra::BrepLoop& loop : loops) {
        boundaries.clear();
        for (Ibra::BrepTrim& trim : loop.trims) {
          boundaries.emplace_back(trim);
        }
        data->addLoop(boundaries);
      }

      return data;
    }
  };

}  // namespace Dune::IGA
