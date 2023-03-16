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

  // At the moment only gridDim und worldDim == 2 supported
  template <int gridDim, int worldDim>
    requires(gridDim == 2 && worldDim == 2)
  class IbraReader {
   public:
    using Grid                = Dune::IGA::NURBSGrid<gridDim, worldDim>;
    using ControlPoint        = Dune::IGA::NURBSPatchData<gridDim, worldDim>::ControlPointType;
    using PatchData           = Dune::IGA::NURBSPatchData<gridDim, worldDim>;
    using ControlPointNetType = Dune::IGA::MultiDimensionNet<gridDim, ControlPoint>;

    static std::shared_ptr<Grid> read(const std::string& fileName) {
      static_assert(gridDim == 2);
      assert(gridDim <= worldDim);

      using json = nlohmann::json;

      // Result
      std::vector<Ibra::Surface> surfaces;
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
          auto geo = j.get<Ibra::Geometry>();

          // Surface
          if (geo.type == Ibra::Type::NurbsSurfaceGeometry3D) surfaces.push_back(j.get<Ibra::Surface>());
          // Curve
          else if (geo.type == Ibra::Type::NurbsCurveGeometry2D)
            curves2D.push_back(j.get<Ibra::Curve2D>());
          // BrepType Representation
          else if (geo.type == Ibra::Type::BrepType)
            brepRepresentations.push_back(j.get<Ibra::BrepRepresentation>());
          else if (geo.type == Ibra::Type::BrepLoopType)
            brepLoopRepresentations.push_back(j.get<Ibra::BrepLoopRepresentation>());
          else if (geo.type == Ibra::Type::BrepTrimType)
            brepTrimRepresentations.push_back(j.get<Ibra::BrepTrimRepresentation>());
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

      auto _surface = brep.surfaces[0];

      const std::array<std::vector<double>, gridDim> knotSpans   = _surface.compileKnotVectors();
      const std::vector<std::vector<ControlPoint>> controlPoints = _surface.compileControlPoints<ControlPoint>();

      std::array<int, gridDim> dimsize = _surface.n_controlPoints;

      PatchData _patchData;
      auto controlNet = ControlPointNetType(dimsize, controlPoints);

      _patchData.knotSpans     = knotSpans;
      _patchData.degree        = _surface.degree;
      _patchData.controlPoints = controlNet;

      // Make the boundaries, and pass them into the grid
      // So the grid can figure out what to do with it
      std::vector<Boundary> globalBoundaries = constructGlobalBoundaries(brep);

      // Create Grid as unique_ptr and move to its destination
      auto grid = std::make_unique<Grid>(_patchData, globalBoundaries);
      return std::move(grid);
    }

    // Deleted Default Constructor, use static method read
    IbraReader() = delete;

   private:
    static std::vector<Boundary> constructGlobalBoundaries(Ibra::Brep& brep) {
      auto trims = brep.trims;

      // Get Boundarys
      std::vector<Boundary> boundaries;
      for (auto& trim : trims) {
        boundaries.emplace_back(trim);
      }

      return boundaries;
    }
  };

}  // namespace Dune::IGA

#endif  // DUNE_IGA_IBRAREADER_HH
