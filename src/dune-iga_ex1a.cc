// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>


#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/ibraReader.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/basis.hh>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>

#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/iga/igaDataCollector.h>
#include <dune/vtk/vtkwriter.hh>

#include "linearElasticTrimmed.h"
#include "igaHelpers.h"
#include "stressEvaluator.h"
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

  const auto gridFileName    = gridParameters.get<std::string>("filename");
  const bool trimGrid        = gridParameters.get<bool>("trim");
  const auto u_degreeElevate = gridParameters.get<int>("u_degreeElevate");
  const auto v_degreeElevate = gridParameters.get<int>("v_degreeElevate");
  const auto globalRefine    = gridParameters.get<int>("globalRefine");

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
  grid->globalRefine(globalRefine);
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
          const auto& localFE = localView.tree().finiteElement();
          const auto& ele = localView.element();

          auto coeff = localFE.localCoefficients();
          auto localKey = coeff.localKey(localIndex);

          if (localKey.codim() == 2) {
            auto vertex = ele.template subEntity<2>(localKey.subEntity());
            auto vgeo = vertex.geometry().center();

            if (Dune::FloatCmp::eq(vgeo, {5, 0}) or Dune::FloatCmp::eq(vgeo, {5, 10})) dirichletFlags[localView.index(localIndex)] = true;
          }
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

            if (Dune::FloatCmp::eq(vgeo, {0, 0})) dirichletFlags[localView.index(localIndex)] = true;
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
    if (globalCoord(0, 0) < 5)
      F[0]              = -lamb;
    else
      F[0]              = lamb;
    return F;
  };

  const auto &indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto neumannPredicate = [](auto &vertex) -> bool {
    return vertex[0] > 10 - 1e-8 or std::fabs(vertex[0]) < 1e-8;
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
  for (auto &element : elements(gridView)) {
    auto localView = basis.flat().localView();
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
    auto x = pos[0] - 5;
    auto y = pos[1] - 5;

    auto r = std::hypot(x, y);
    auto theta = std::atan2(y, x);

//    auto pi = std::numbers::pi;
//    if (theta <= pi and theta > pi/2)
//      theta = mapToRange(theta, pi, pi/2, 0.0, pi/2);
//    else if (theta < 0 and theta > -pi/2)
//      theta *= -1;
//    else if (theta <= -pi/2 and theta >= -pi)
//      theta = mapToRange(theta, -pi, -pi/2, 0.0, pi/2);

    return {r, theta};
  };


  double Tx = lambdaLoad;
  double R = 1;

  auto analyticalSolutionStress = [&](const auto& pos) -> Dune::FieldVector<double, 3>{
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

  auto analyticalSolutionStressBackTransformed = [&](const auto& pos) -> Dune::FieldVector<double, 3> {
    auto sigma_thetar = analyticalSolutionStress(pos);

    auto [r, theta] = toPolar(pos);
    theta *= -1;
    return coordinateTransform(sigma_thetar, theta);
  };

//  auto kappa = (3 - nu) / (1 + nu);
//  auto G = E / (1 - 2 * nu);

  auto kappa = (3 - nu) / (1 + nu);
  auto G = E / (2 * (1 + nu));

  auto factor = Tx * R / (8 * G);

  auto analyticalSolutionDisplacements = [&](const auto& pos) -> Dune::FieldVector<double, 2> {
    auto [r, theta] = toPolar(pos);

    auto costh = std::cos(theta);
    auto sinth = std::sin(theta);
    auto cos3th = std::cos(3 * theta);
    auto sin3th = std::sin(3 * theta);

    auto ra = r / R;
    auto ar = R / r;
    auto ar3 = std::pow(R, 3) / std::pow(r, 3);

    Dune::FieldVector<double, 2> res;
    res[0] = factor * (
                 ra * (kappa + 1) * costh +
                 2 * ar * ((1 + kappa) * costh + cos3th) -
                 2 * ar3 * cos3th
                 );
    res[1] = factor * (
                 ra * (kappa - 3) * sinth +
                 2 * ar * ((1 - kappa) * sinth + sin3th) -
                 2 * ar3 * sin3th
                 );
    return res;
  };

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView, subsample);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));

  auto analyticalSolutionFunctionStress =  Dune::Functions::makeAnalyticGridViewFunction(analyticalSolutionStressBackTransformed, gridView);
  auto analyticalSolutionFunctionDisplacements =  Dune::Functions::makeAnalyticGridViewFunction(analyticalSolutionDisplacements, gridView);

  vtkWriter.addPointData(analyticalSolutionFunctionDisplacements, Dune::VTK::FieldInfo("displacement solution", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(analyticalSolutionFunctionStress, Dune::VTK::FieldInfo("stress solution", Dune::VTK::FieldInfo::Type::vector, 3));

  // Dirichlet Flag func
  std::vector<int> flags(dirichletValues.size());
  for (size_t i : std::views::iota(0u, dirichletValues.size()))
    flags[i] = dirichletValues.isConstrained(i);
  auto dirichletFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), flags);

  vtkWriter.addPointData(dirichletFunc,
                         Dune::VTK::FieldInfo("dirichlet BC", Dune::VTK::FieldInfo::Type::vector, 2));

  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::normalStress>>(
          gridView, &fes, D_Glob, lambdaLoad)));

  vtkWriter.addPointData(Dune::Vtk::Function<GridView>(
      std::make_shared<StressEvaluator2D<GridView, LinearElasticType, StressEvaluatorComponents::vonMises>>(
          gridView, &fes, D_Glob, lambdaLoad)));


  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += std::fabs(f);
  std::cout << "Total Force: " << totalForce << std::endl;

  vtkWriter.write(gridFileName + "-ex1");

  return 0;
}
