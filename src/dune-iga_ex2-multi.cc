// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include "anaSolution.h"
#include "linearElasticTrimmed.h"
// #include "stressEvaluator.h"
#include "timer.h"

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/io/resultFunction.hh>

#include "dune/iga/io/ibra/ibrareader.hh"
#include "dune/iga/io/igaDataCollector.h"
#include "dune/iga/utils/igaHelpers.h"
#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/vtk/vtkwriter.hh>

constexpr int gridDim  = 2;
constexpr int worldDim = 2;
double lambdaLoad      = 1;

void analyse(auto argv, int u_degreeElevate, int v_degreeElevate, int globalRefine, bool writeVTK) {
  /// Read Parameter
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree &gridParameters        = parameterSet.sub("GridParameters");
  const Dune::ParameterTree &materialParameters    = parameterSet.sub("MaterialParameters");
  const Dune::ParameterTree &postProcessParameters = parameterSet.sub("PostProcessParameters");

  const auto gridFileName       = gridParameters.get<std::string>("filename");
  const bool trimGrid           = gridParameters.get<bool>("trim");
  //const auto u_degreeElevate    = gridParameters.get<int>("u_degreeElevate");
  //const auto v_degreeElevate    = gridParameters.get<int>("v_degreeElevate");
  //const auto globalRefine       = gridParameters.get<int>("globalRefine");
  const auto refineInUDirection = gridParameters.get<int>("u_refine");
  const auto refineInVDirection = gridParameters.get<int>("v_refine");

  const auto E  = materialParameters.get<double>("E");
  const auto nu = materialParameters.get<double>("nu");

  const int subsample = postProcessParameters.get<int>("subsample");

  // Log the Paramaters
  spdlog::info(
      "Filename: {} \nThe following parameters were used: \nMaterial: E {}, nu {} \nRefinements: global {}, u {}, v {}", gridFileName, E, nu, globalRefine, refineInUDirection, refineInVDirection);

  /// Instantiate a timer
  Timer timer;
  timer.startTimer("all");
  auto outputFileName = Timer<>::makeUniqueName(argv[0]);

  /// Create Grid
  timer.startTimer("grid");

  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate}, {globalRefine, globalRefine});

  const auto& patchData = grid->getPatch().getPatchData();
  spdlog::info("Degree: u {}, v {}", patchData->degree[0], patchData->degree[1]);

  grid->globalMultiRefine(globalRefine, refineInUDirection, refineInVDirection);

  GridView gridView = grid->leafGridView();

  spdlog::info("Loading and trimming the grid took {} milliseconds ", timer.stopTimer("grid").count());
  timer.startTimer("basis");

  // Instantiate Analytical Solution
  AnalyticalSolution analyticalSolution{E, nu, lambdaLoad, 1.0, Dune::FieldVector<double, 2>{0, 0}};

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(gridView.impl().getPreBasis(), FlatInterleaved()));
  auto basisSize = basis.flat().size();
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(
        Dune::Functions::subspaceBasis(basis_, 1), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
          if (std::fabs(intersection.geometry().center()[1]) < 1e-8)
            dirichletFlags[localView.index(localIndex)] = true;
        });
    Dune::Functions::forEachUntrimmedBoundaryDOF(
        Dune::Functions::subspaceBasis(basis_, 0), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
          if (std::fabs(intersection.geometry().center()[0]) < 1e-8)
            dirichletFlags[localView.index(localIndex)] = true;
        });
  });
  spdlog::info("Creating a basis and fixing dofs took {} milliseconds, fixed {} dofs",
               timer.stopTimer("basis").count(), dirichletValues.fixedDOFsize());

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
    auto stresses = analyticalSolution.stressSolution(globalCoord);

    // left side
    if (globalCoord(0, 0) > 4 - 1e-8)
      F << stresses[0], stresses[2];
    // Top side
    else if (globalCoord(1, 0) > 4 - 1e-8)
      F << stresses[2], stresses[1];

    return F;
  };

  const auto &indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto neumannPredicate = [](auto &vertex) -> bool {
    return std::fabs(vertex[0]) > 4 - 1e-8 or std::fabs(vertex[1]) > 4 - 1e-8;
  };
  for (auto &&vertex : vertices(gridView)) {
    auto coords                             = vertex.geometry().corner(0);
    neumannVertices[indexSet.index(vertex)] = neumannPredicate(coords);
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView);

  // create property for insertion by vertex vector
  BoundaryPatchEnclosingVerticesPropertyTrimmed<GridView, 1> prop(gridView, neumannVertices);
  neumannBoundary.insertFacesByProperty(prop);

  /// Add the linear elastic 2D planar solid elements to the vector "fes"
  timer.startTimer("createElement");
  for (auto &element : elements(gridView))
    fes.emplace_back(basis, element, E, nu, volumeLoad, &neumannBoundary, neumannBoundaryLoad);

  spdlog::info("Creating the {} Finite Elements took {} milliseconds ", gridView.size(0), timer.stopTimer("createElement").count());

  /// Create a sparse assembler
  auto sparseAssembler = Ikarus::SparseFlatAssembler(fes, dirichletValues);

  /// Define "elastoStatics" affordances and create functions for stiffness matrix and residual calculations
  auto req = Ikarus::FErequirements().insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad).addAffordance(Ikarus::AffordanceCollections::elastoStatics);

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

  const auto& K    = nonLinOp.derivative();
  const auto& Fext = nonLinOp.value();

  /// solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  timer.startTimer("solve");

  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);

  spdlog::info("The solver took {} milliseconds ", timer.stopTimer("solve").count());
  spdlog::info("Total time spent until solution: {} milliseconds  ", timer.stopTimer("all").count());

  /// Postprocess
  auto dispGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), Fext);

  // Estimate error
  auto resReq = Ikarus::ResultRequirements()
                    .insertGlobalSolution(Ikarus::FESolutions::displacement, D_Glob)
                    .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad)
                    .addResultRequest(Ikarus::ResultType::linearStress);
  auto stressFunction = Dune::Vtk::Function<GridView>(std::make_shared<Ikarus::ResultFunction<LinearElasticType>>(&fes, resReq));

  auto l2_error = analyticalSolution.estimateError(gridView, stressFunction, dispGlobalFunc, fes.front().integrationOrder());

  spdlog::info("Stress l2 error: {}", l2_error[0]);
  spdlog::info("Displacement l2 error: {}", l2_error[1]);

  // Amount of quadrature Points
  int nQuadraturePoints = 0;
  std::ranges::for_each(fes, [&](const auto& fe) {nQuadraturePoints += fe.numOfQuadraturePoints();});
  spdlog::info("Num Quadrature Points: {}", nQuadraturePoints);
  spdlog::info("RefinementTarget: {}", Grid::Traits::TrimmedElementRepresentationType::targetTolerance);

  if (not writeVTK)
    return;
  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));


  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(stressFunction));

  auto stressGridFunction =  Dune::Functions::makeAnalyticGridViewFunction(analyticalSolution.stressLambda(), gridView);
  auto displacementGridFunction =  Dune::Functions::makeAnalyticGridViewFunction(analyticalSolution.displacementLambda(), gridView);

  vtkWriter.addPointData(stressGridFunction, Dune::VTK::FieldInfo("stress solution", Dune::VTK::FieldInfo::Type::vector, 3));
  vtkWriter.addPointData(displacementGridFunction, Dune::VTK::FieldInfo("displacement solution", Dune::VTK::FieldInfo::Type::vector, 2));

  vtkWriter.write(outputFileName);
}



int main(int argc, char **argv) {
  Ikarus::init(argc, argv);
  for (auto i : std::views::iota(0, 5))
    for (auto j : std::views::iota(2, 9))
      analyse(argv, i, i, j, false);



  return 0;
}
