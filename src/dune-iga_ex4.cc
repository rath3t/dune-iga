// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include "kirchhoffPlate.hh"
#include "stressEvaluator.h"
#include "timer.h"

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

#include "dune/iga/ibra/ibraReader.hh"
#include "dune/iga/io/igaDataCollector.h"
#include "dune/iga/utils/igaHelpers.h"
#include <dune/common/parametertreeparser.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/vtk/vtkwriter.hh>

int main(int argc, char **argv) {
  Ikarus::init(argc, argv);

  constexpr int gridDim  = 2;
  constexpr int worldDim = 2;

  /// Read Parameter
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree &gridParameters        = parameterSet.sub("GridParameters");
  const Dune::ParameterTree &materialParameters    = parameterSet.sub("MaterialParameters");
  const Dune::ParameterTree &postProcessParameters = parameterSet.sub("PostProcessParameters");

  const auto gridFileName    = gridParameters.get<std::string>("filename");
  const bool trimGrid        = gridParameters.get<bool>("trim");
  const auto u_degreeElevate = gridParameters.get<int>("u_degreeElevate");
  const auto v_degreeElevate = gridParameters.get<int>("v_degreeElevate");
  const auto globalRefine    = gridParameters.get<int>("globalRefine");

  const auto E  = materialParameters.get<double>("E");
  const auto nu = materialParameters.get<double>("nu");
  const auto thk = materialParameters.get<double>("thk");

  double lambdaLoad      = 1 * std::pow(thk, 3);

  const int subsample = postProcessParameters.get<int>("subsample");

  // Log the Paramaters
  spdlog::info(
      "Filename: {} \n The following parameters were used: \nMaterial: E {}, nu {} \nRefinements: global {}, u {}, v {}", gridFileName, E, nu, globalRefine, refineInUDirection, refineInVDirection);

  /// Instantiate a timer
  Timer timer;
  timer.startTimer("all");
  auto outputFileName = Timer<>::makeUniqueName(argv[0]);

  /// Create Grid
  timer.startTimer("grid");

  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate});
  grid->globalRefine(globalRefine);
  GridView gridView = grid->leafGridView();

  const auto& patchData = grid->getPatch().getPatchData();
  spdlog::info("Degree: u {}, v {}", patchData->degree[0], patchData->degree[1]);

  spdlog::info("Loading and trimming the grid took {} milliseconds ", timer.stopTimer("grid").count());
  timer.startTimer("basis");

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, gridView.impl().getPreBasis());

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(basis_, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
      dirichletFlags[localView.index(localIndex)] = true;
    });
  });
  spdlog::info("Creating a basis and fixing dofs took {} milliseconds, fixed {} dofs",
               timer.stopTimer("basis").count(), dirichletValues.fixedDOFsize());

  /// Declare a vector "fes" of kirchhoff plate 2D elements
  using LinearElasticType = Ikarus::KirchhoffPlate<decltype(basis)>;
  std::vector<LinearElasticType> fes;

  /// Add the linear elastic 2D planar solid elements to the vector "fes"
  timer.startTimer("createElement");
  for (auto &element : elements(gridView)) {
    auto localView = basis.flat().localView();
    fes.emplace_back(basis, element, E, nu, thk);
  }
  spdlog::info("Creating the {} Finite Elements took {} milliseconds ", gridView.size(0), timer.stopTimer("createElement").count());

  /// Create a sparse assembler
  auto sparseAssembler = Ikarus::SparseFlatAssembler(fes, dirichletValues);

  /// Define "elastoStatics" affordances and create functions for stiffness matrix and residual calculations
  auto req = Ikarus::FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto KFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getMatrix(req);
  };

  auto residualFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getVector(req);
  };

  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.flat().size());

  /// Create a non-linear operator
  timer.startTimer("assemble");

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::linearAlgebraFunctions(residualFunction, KFunction),
                                            Ikarus::parameter(D_Glob, lambdaLoad));

  spdlog::info("The assembly took {} milliseconds with {} dofs", timer.stopTimer("assemble").count(),
               basis.flat().size());

  const auto &K    = nonLinOp.derivative();
  const auto &Fext = nonLinOp.value();

  /// solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  timer.startTimer("solve");

  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);

  spdlog::info("The solver took {} milliseconds ", timer.stopTimer("solve").count());
  spdlog::info("Total time spent until solution: {} milliseconds  ", timer.stopTimer("all").count());

  /// Postprocess
  auto dispGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 1>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 1>>(basis.flat(), Fext);

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::scalar, 1));


  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::kirchhoff_moments>>(
          gridView, &fes, D_Glob, lambdaLoad)));


  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += f;
  std::cout << "Total Force: " << totalForce << std::endl;

  vtkWriter.write(gridFileName);

  return 0;
}
