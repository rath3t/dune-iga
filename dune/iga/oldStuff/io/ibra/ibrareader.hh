// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "ibrageometry.hh"

#include <clipper2/clipper.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <utility>

#include "dune/iga/nurbsalgorithms.hh"
#include "dune/iga/trim/nurbstrimboundary.hh"
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

namespace Dune::IGA {

template <int dim, int dimworld, typename ScalarType>
class NURBSGrid;

// At the moment only gridDim und worldDim == 2 or 3 supported
template <int gridDim, int worldDim, typename ScalarType = double>
requires(gridDim == 2) && (worldDim == 2 || worldDim == 3) && (gridDim <= worldDim)
class IbraReader
{
public:
  using Grid                = Dune::IGA::NURBSGrid<gridDim, worldDim, ScalarType>;
  using ControlPoint        = Dune::IGA::NURBSPatchData<gridDim, worldDim>::ControlPointType;
  using PatchData           = Dune::IGA::NURBSPatchData<gridDim, worldDim>;
  using ControlPointNetType = Dune::IGA::MultiDimensionalNet<gridDim, ControlPoint>;

  static std::unique_ptr<Grid> read(const std::string& fileName, const bool trim = true,
                                    std::array<int, 2> elevateDegree  = {0, 0},
                                    std::array<int, 2> preKnotRefine  = {0, 0},
                                    std::array<int, 2> postKnotRefine = {0, 0}) {
    std::ifstream ibraInputFile;
    ibraInputFile.open(fileName);
    return read(ibraInputFile, trim, elevateDegree, preKnotRefine, postKnotRefine);
  }

  template <typename InputStringType>
  requires(not std::convertible_to<std::string, InputStringType> and
           not std::convertible_to<InputStringType, const char*>)
  static std::unique_ptr<Grid> read(InputStringType& ibraInputFile, const bool trim = true,
                                    std::array<int, 2> elevateDegree  = {0, 0},
                                    std::array<int, 2> preKnotRefine  = {0, 0},
                                    std::array<int, 2> postKnotRefine = {0, 0}) {
    using json = nlohmann::json;

    // Result
    std::vector<Ibra::Surface<worldDim>> surfaces;
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
                 "Error in reading input stream: " << ", parse error at byte " << ex.byte << " What: " << ex.what());
    }

    // Make Connections
    auto brepRepresentation = brepRepresentations[0];

    Ibra::Brep brep{brepRepresentation, curves2D, surfaces, brepLoopRepresentations, brepTrimRepresentations};
    // Each surface in a brep is one Patch, as for now only one brep is supported in NURBSGrid
    assert(brep.surfaces.size() == 1);

    std::cout << "IbraReader successfully read in 1 patch with " << brep.trims.size() << " edges." << std::endl;

    // Reader has done its job, now NURBSGrid can be constructed

    Ibra::Surface<worldDim> _surface                           = brep.surfaces[0];
    const std::array<std::vector<double>, gridDim> knotSpans   = _surface.compileKnotVectors();
    const std::vector<std::vector<ControlPoint>> controlPoints = _surface.transformControlPoints();
    std::array<int, gridDim> dimsize                           = _surface.n_controlPoints;

    auto controlNet = ControlPointNetType(dimsize, controlPoints);
    PatchData _patchData{knotSpans, controlNet, _surface.degree};
    // Optional Pre Refinement of knots
    for (int i = 0; i < gridDim; ++i) {
      if (preKnotRefine[i] > 0) {
        auto additionalKnots = generateRefinedKnots(knotSpans, i, preKnotRefine[i]);
        _patchData           = knotRefinement<gridDim>(_patchData, additionalKnots, i);
      }
    }

    // Optional Degree Elevate
    for (int i = 0; i < gridDim; ++i)
      if (elevateDegree[i] > 0)
        _patchData = degreeElevate(_patchData, i, elevateDegree[i]);

    // Optional Post Refinement of knots
    for (int i = 0; i < gridDim; ++i) {
      if (postKnotRefine[i] > 0) {
        auto additionalKnots = generateRefinedKnots(knotSpans, i, postKnotRefine[i]);
        _patchData           = knotRefinement<gridDim>(_patchData, additionalKnots, i);
      }
    }

    // Make the trimData and pass them into the grid
    // So the grid can figure out what to do with it
    auto trimData = constructGlobalBoundaries(brep);

    if (trim)
      return std::make_unique<Grid>(_patchData, trimData);
    else
      return std::make_unique<Grid>(_patchData);
  }

private:
  static auto constructGlobalBoundaries(const Ibra::Brep<worldDim>& brep) -> std::shared_ptr<TrimData> {
    auto data = std::make_shared<TrimData>();
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

} // namespace Dune::IGA

namespace Dune {

template <int gridDim, int worldDim, typename ScalarType>
struct DGFGridInfo<Dune::IGA::NURBSGrid<gridDim, worldDim, ScalarType>>
{
  static int refineStepsForHalf() {
    return 1;
  }
  static double refineWeight() {
    return std::pow(0.5, gridDim);
  }
};

template <typename Grid_>
struct JSONGridFactory
{
  static constexpr std::integral auto gridDim  = Grid_::dim;
  static constexpr std::integral auto worldDim = Grid_::dimworld;
  using ScalarType                             = typename Grid_::ScalarType;

  using Grid = Grid_;
  JSONGridFactory(std::string p_fileName) {
    grid_ = IGA::IbraReader<gridDim, worldDim, ScalarType>::read(p_fileName);
  }
  JSONGridFactory(std::istream& p_istreamGrid) {
    grid_ = IGA::IbraReader<gridDim, worldDim, ScalarType>::read(p_istreamGrid);
  }

  std::unique_ptr<Grid> grid_;

  [[nodiscard]] Grid* grid() const {
    return grid_.release();
  }
};

template <int gridDim, int worldDim, typename ScalarType>
struct DGFGridFactory<Dune::IGA::NURBSGrid<gridDim, worldDim, ScalarType>>
{
  using Grid = Dune::IGA::NURBSGrid<gridDim, worldDim, ScalarType>;
  DGFGridFactory(std::string p_fileName) {
  }
  DGFGridFactory(std::istream& p_istreamGrid) {
  }

  Grid* grid() const {
    DUNE_THROW(Dune::NotImplemented, "DGFGridFactory not implemtented use JSONGridFactory");
    return nullptr;
  }
};
} // namespace Dune
