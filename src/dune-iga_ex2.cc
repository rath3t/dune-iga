// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

#include "dune/functions/gridfunctions/analyticgridviewfunction.hh"
#include "dune/geometry/refinement/base.cc"
#include "dune/iga/test/plotFunctionality.h"
#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/iga/ibraReader.hh>
#include <dune/iga/nurbsgrid.hh>

#include "linearElasticTrimmed.h"
#include "stressEvaluator.h"
#include "helpers.h"

#include <ikarus/utils/observer/controlVTKWriter.hh>

#include <dune/grid/io/file/vtk/function.hh>
#include <dune/iga/igaDataCollector.h>
#include <dune/vtk/vtkwriter.hh>

#include <dune/iga/test/plotFunctionality.h>

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

  /// Create Grid
  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate});

  grid->globalRefine(globalRefine);
  grid->globalRefineInDirection(0, refineInUDirection);
  grid->globalRefineInDirection(1, refineInVDirection);

  GridView gridView = grid->leafGridView();
  Plot::plotParametricGridAndPhysicalGrid(grid);

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(gridView.impl().getPreBasis(), FlatInterleaved()));

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(basis_, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
      if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
    });
  });
  std::cout << dirichletValues.fixedDOFsize() << " Dofs fixed" << std::endl;

  /// Declare a vector "fes" of linear elastic 2D planar solid elements
  using LinearElasticType = Ikarus::LinearElasticTrimmed<decltype(basis)>;
  std::vector<LinearElasticType> fes;

  /// function for volume load- here: returns zero
  auto volumeLoad = [](auto &globalCoord, auto &lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    return fext;
  };

  /// neumann boundary load in horizontal direction
  auto neumannBoundaryLoad = [&](auto &globalCoord, auto &lamb) {
    Eigen::Vector2d F = Eigen::Vector2d::Zero();
    F[1]              = -lamb;
    return F;
  };

  const auto &indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto neumannPredicate = [](auto &vertex) -> bool { return Dune::FloatCmp::eq<double>(vertex[0], 15); };
  for (auto &&vertex : vertices(gridView)) {
    auto coords                             = vertex.geometry().corner(0);
    neumannVertices[indexSet.index(vertex)] = neumannPredicate(coords);
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView);

  // create property for insertion by vertex vector
  BoundaryPatchEnclosingVerticesPropertyTrimmed<GridView, 1> prop(gridView, neumannVertices);
  neumannBoundary.insertFacesByProperty(prop);

  /// Add the linear elastic 2D planar solid elements to the vector "fes"
  for (auto &element : elements(gridView)) {
    auto localView = basis.flat().localView();
    fes.emplace_back(basis, element, E, nu, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
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
  auto startAssembly    = std::chrono::high_resolution_clock::now();
  auto nonLinOp         = Ikarus::NonLinearOperator(Ikarus::linearAlgebraFunctions(residualFunction, KFunction),
                                                    Ikarus::parameter(D_Glob, lambdaLoad));
  auto stopAssembly     = std::chrono::high_resolution_clock::now();
  auto durationAssembly = duration_cast<std::chrono::milliseconds>(stopAssembly - startAssembly);
  spdlog::info("The assembly took {:>6d} milliseconds with {:>7d} dofs", durationAssembly.count(), basis.flat().size());

  const auto &K    = nonLinOp.derivative();
  const auto &Fext = nonLinOp.value();

  /// solve the linear system
  auto linSolver   = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  auto startSolver = std::chrono::high_resolution_clock::now();

  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);
  auto stopSolver     = std::chrono::high_resolution_clock::now();
  auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
  spdlog::info("The solver took {} milliseconds ", durationSolver.count());

  /// Stresses
  Ikarus::ResultRequirements resultRequirements = Ikarus::ResultRequirements()
                                                      .insertGlobalSolution(Ikarus::FESolutions::displacement, D_Glob)
                                                      .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                                                      .addResultRequest(Ikarus::ResultType::cauchyStress);

  /// Postprocess
  auto dispGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), Fext);

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));

  vtkWriter.addPointData(dispGlobalFunc,
                         Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));

  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::normalStress>>(D_Glob, lambdaLoad, fes)));
  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::shearStress>>(D_Glob, lambdaLoad, fes)));
  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::vonMieses>>(D_Glob, lambdaLoad, fes)));

  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += f;

  std::cout << "Total Force: " << totalForce << std::endl;

  vtkWriter.write(gridFileName);

  return 0;
}
