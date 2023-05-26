// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>



#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

#include "dune/iga/io/ibra/ibrareader.hh"
#include "dune/iga/io/igaDataCollector.h"
#include <dune/common/parametertreeparser.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/vtk/vtkwriter.hh>
#include "timer.h"


int main(int argc, char **argv) {
  Ikarus::init(argc, argv);

  constexpr int gridDim  = 2;
  constexpr int worldDim = 2;
  double lambdaLoad      = 1;

  /// Read Parameter
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree &gridParameters        = parameterSet.sub("GridParameters");
  const Dune::ParameterTree &materialParameters    = parameterSet.sub("MaterialParameters");
  const Dune::ParameterTree &postProcessParameters = parameterSet.sub("PostProcessParameters");

  const auto gridFileName       = gridParameters.get<std::string>("filename");
  const bool trimGrid           = gridParameters.get<bool>("trim");
  const auto u_degreeElevate    = gridParameters.get<int>("u_degreeElevate");
  const auto v_degreeElevate    = gridParameters.get<int>("v_degreeElevate");
  const auto globalRefine       = gridParameters.get<int>("globalRefine");
  const auto refineInUDirection = gridParameters.get<int>("u_refine");
  const auto refineInVDirection = gridParameters.get<int>("v_refine");

  const auto E  = materialParameters.get<double>("E");
  const auto nu = materialParameters.get<double>("nu");

  const int subsample = postProcessParameters.get<int>("subsample");

  // Log the Paramaters
  spdlog::info(
      "Filename: {} \n The following parameters were used: \nMaterial: E {}, nu {} \nRefinements: global {}, u {}, v {}", gridFileName, E, nu, globalRefine, refineInUDirection, refineInVDirection);

  /// Instantiate a timer
  Timer timer;
  auto outputFileName = Timer<>::makeUniqueName(argv[0]);

  /// Create Grid

  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate}, {globalRefine, globalRefine});

  grid->globalMultiRefine(globalRefine, refineInUDirection, refineInVDirection);

  GridView gridView = grid->leafGridView();

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  // Calculate Area
  Dune::QuadratureRule<double, 2> rule;
  double area = 0.0;
  for (auto& ele : elements(gridView)) {
    ele.impl().fillQuadratureRule(rule, 2);
    auto geo = ele.geometry();
    for (auto& ip : rule)
      area += geo.integrationElement(ip.position()) * ip.weight();
  }
  auto targetArea = 4 * 4  - 0.25 * (std::numbers::pi * std::pow(1, 2));
  spdlog::info("Area {}, (target {})", area, targetArea);

  vtkWriter.write(outputFileName);

  return 0;
}
