//
// Created by Henri on 30.04.2023.
//
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
#include <ikarus/utils/observer/controlVTKWriter.hh>

#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/ibraReader.hh>
#include <dune/iga/igaDataCollector.h>
#include <dune/iga/nurbsgrid.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "igaHelpers.h"
#include "linearElasticTrimmed.h"
#include "stressEvaluator.h"
#include "timer.h"
#include "kirchhoffPlate.hh"

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

  /// Instantiate a timer
  Timer timer;
  timer.startTimer("all");

  /// Create Grid
  timer.startTimer("grid");

  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate});
  // Refine if neccesary
  grid->globalMultiRefine(globalRefine, refineInUDirection, refineInVDirection);
  GridView gridView = grid->leafGridView();

  spdlog::info("Loading and trimming the grid took {} milliseconds ", timer.stopTimer("grid").count());
  timer.startTimer("basis");

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(gridView.impl().getPreBasis(), FlatInterleaved()));

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(
        Dune::Functions::subspaceBasis(basis_, 0), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
          if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
        });

    Dune::Functions::forEachUntrimmedBoundaryDOF(
        Dune::Functions::subspaceBasis(basis_, 1), [&](auto &&localIndex, auto &&localView, auto &&intersection) {
          const auto& localFE = localView.tree().finiteElement();
          const auto& ele = localView.element();

          auto coeff = localFE.localCoefficients();
          auto localKey = coeff.localKey(localIndex - localFE.size());

          if (localKey.codim() == 2) {
            auto vertex = ele.template subEntity<2>(localKey.subEntity());
            auto vgeo = vertex.geometry().center();

            if (vgeo.two_norm() < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
          }
        });
  });

  spdlog::info("Creating a basis and fixing dofs took {} milliseconds, fixed {} dofs ",
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
    F[0]              = lamb;
    return F;
  };

  const auto &indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto neumannPredicate = [](auto &vertex) -> bool { return Dune::FloatCmp::eq(vertex[0], 20.0); };
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
    fes.emplace_back(basis, element, E, nu, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
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
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), Fext);

  // Analytical solution
  auto mapToRange = []<std::floating_point T>(T value, T inputMin, T inputMax, T outputMin, T outputMax) -> T {
    return (value - inputMin) * (outputMax - outputMin) / (inputMax - inputMin) + outputMin;
  };

  auto toPolar = [&](const auto& pos) -> std::pair<double, double> {
    auto x = pos[0] - 10;
    auto y = pos[1] - 1;

    auto r = std::hypot(x, y);
    auto theta = std::atan2(y, x);

    auto pi = std::numbers::pi;
    if (theta <= pi and theta > pi/2)
      theta = mapToRange(theta, pi, pi/2, 0.0, pi/2);
    else if (theta < 0 and theta > -pi/2)
      theta *= -1;
    else if (theta <= -pi/2 and theta >= -pi)
      theta = mapToRange(theta, -pi, -pi/2, 0.0, pi/2);

    return {r, theta};
  };

  double Tx = 1;
  double R = 0.3;


  auto analyticalSolution = [&](const auto& pos) -> Dune::FieldVector<double, 3>{
    auto [r, theta] = toPolar(pos);

    auto th  = Tx / 2;
    auto rr2 = std::pow(R, 2) / std::pow(r, 2);
    auto rr4 = std::pow(R, 4) / std::pow(r, 4);
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    Dune::FieldVector<double, 3> sigma_thetar;
    sigma_thetar[0] = th * (1 - rr2) + th * (1 + 3 * rr4 - 4 * rr2) * cos2t;
    sigma_thetar[1] = th * (1 + rr2) - th * (1 + 3 * rr4) * cos2t;
    sigma_thetar[2] = - th * (1 - 3 * rr4 + 2 * rr2) * sin2t;

    return sigma_thetar;
  };

  auto coordinateTransform = [&](const Dune::FieldVector<double, 3>& stress, double theta) -> Dune::FieldVector<double, 3> {
    auto cos2t = std::cos(2 * theta);
    auto sin2t = std::sin(2 * theta);

    const auto s_x = stress[0];
    const auto s_y = stress[1];
    const auto s_xy = stress[2];

    auto hpl = 0.5 * (s_x + s_y);
    auto hmi = 0.5 * (s_x - s_y);

    Dune::FieldVector<double, 3> sigma_transformed;
    sigma_transformed[0] = hpl + hmi * cos2t + s_xy * sin2t;
    sigma_transformed[1] = hpl - hmi * cos2t - s_xy * sin2t;
    sigma_transformed[2] = - hmi  * sin2t + s_xy * cos2t;

    return sigma_transformed;
  };

  auto analyticalSolutionBackTransformed = [&](const auto& pos) -> Dune::FieldVector<double, 3> {
    auto sigma_thetar = analyticalSolution(pos);

    auto [r, theta] = toPolar(pos);
    theta *= -1;
    return coordinateTransform(sigma_thetar, theta);
  };

  using namespace Dune::Functions;


  auto allSigmas = [](const auto& sigma) {
    Dune::FieldVector<double, 3> res;
    res[0] = sigma(0, 0);
    res[1] = sigma(1, 0);
    res[2] = sigma(2, 0);
    return res;
  };

  auto localStress = localFunction(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::user_function>>(
          gridView, &fes, D_Glob, allSigmas, "all_sigmas", lambdaLoad)));

  auto analyticalSolutionFunction =  Dune::Functions::makeAnalyticGridViewFunction(analyticalSolutionBackTransformed, gridView);
  auto localStressAnalytic = localFunction(analyticalSolutionFunction);

  double l2_error = 0.0;
  for (auto& ele : elements(gridView)) {
    localStress.bind(ele);
    localStressAnalytic.bind(ele);

    auto center = ele.geometry().center();
    if (center[0] < 10 and center[1] < 1)
      continue;

    const auto rule = ele.impl().getQuadratureRule();

    for (auto& gp : rule) {
      const auto stress_ana = localStressAnalytic(gp.position());
      const auto stress_fe = localStress(gp.position());

      l2_error += Dune::power(stress_ana[0] - stress_fe[0], 2) * ele.geometry().integrationElement(gp.position()) * gp.weight();
    }
  }
  spdlog::info("l2 error: {}", l2_error);


  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));

  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::normalStress>>(
          gridView, &fes, D_Glob, lambdaLoad)));

//  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::shearStress>>(
//          gridView, &fes, D_Glob, lambdaLoad)));
//
//  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::vonMises>>(
//          gridView, &fes, D_Glob, lambdaLoad)));
//
//  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::principalStress>>(
//          gridView, &fes, D_Glob, lambdaLoad)));

//  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::normalStress>>(
//          gridView, &fes, D_Glob, lambdaLoad)));
//  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::shearStress>>(
//          gridView, &fes, D_Glob, lambdaLoad)));
//  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::vonMises>>(
//          gridView, &fes, D_Glob, lambdaLoad)));
//  vtkWriter.addCellData(Dune::Vtk::Function<GridView>(
//      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::principalStress>>(
//          gridView, &fes, D_Glob, lambdaLoad)));




  vtkWriter.addPointData(analyticalSolutionFunction, Dune::VTK::FieldInfo("stress solution", Dune::VTK::FieldInfo::Type::vector, 3));


  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += f;
  std::cout << "Total Force: " << totalForce << std::endl;

  vtkWriter.write(gridFileName);

  return 0;
}
