// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>


#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/finiteElements/feRequirements.hh>

#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/ibraReader.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/basis.hh>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/iga/igaDataCollector.h>
#include <dune/vtk/vtkwriter.hh>

#include "igaHelpers.h"
#include "stressEvaluator.h"
#include "timer.h"
#include "ReissnerMindlinPlate.hh"

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

  const auto gridFileName    = gridParameters.get<std::string>("filename");
  const bool trimGrid        = gridParameters.get<bool>("trim");
  const auto u_degreeElevate = gridParameters.get<int>("u_degreeElevate");
  const auto v_degreeElevate = gridParameters.get<int>("v_degreeElevate");
  const auto globalRefine    = gridParameters.get<int>("globalRefine");

  const auto E  = materialParameters.get<double>("E");
  const auto nu = materialParameters.get<double>("nu");
  const auto thk = materialParameters.get<double>("thk");

  const int subsample = postProcessParameters.get<int>("subsample");

  /// Instantiate a timer
  Timer timer;
  timer.startTimer("all");

  /// Create Grid
  timer.startTimer("grid");

  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate});
  grid->globalRefine(globalRefine);
  GridView gridView = grid->leafGridView();

  spdlog::info("Loading and trimming the grid took {} milliseconds ", timer.stopTimer("grid").count());
  timer.startTimer("basis");

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<3>(gridView.impl().getPreBasis(), FlatInterleaved()));

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(basis_, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
      dirichletFlags[localView.index(localIndex)] = true;
    });
  });
  spdlog::info("Creating a basis and fixing dofs took {} milliseconds, fixed {} dofs ",
               timer.stopTimer("basis").count(), dirichletValues.fixedDOFsize());

  /// Declare a vector "fes" of kirchhoff plate 2D elements
  using LinearElasticType = Ikarus::ReissnerMindlinPlate<decltype(basis)>;
  std::vector<LinearElasticType> fes;

  auto volumeLoad = [](auto &globalCoord, auto &lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    fext[0] = 1;
    return fext;
  };

  /// Add the linear elastic 2D planar solid elements to the vector "fes"
  for (auto &element : elements(gridView)) {
    auto localView = basis.flat().localView();
    fes.emplace_back(basis, element, E, nu, thk, volumeLoad);
  }
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
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(basis.flat(), Fext);

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 3));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 3));


  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::RM_moments>>(
          gridView, &fes, D_Glob, lambdaLoad)));

  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::RM_forces>>(
          gridView, &fes, D_Glob, lambdaLoad)));

  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += f;
  std::cout << "Total Force: " << totalForce << std::endl;


  vtkWriter.write(gridFileName);

  return 0;
}
