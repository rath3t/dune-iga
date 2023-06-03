// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "kirchhoffplate.hh"

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

#include "dune/iga/io/ibra/ibrareader.hh"
#include "dune/iga/io/igadatacollector.hh"
#include "dune/iga/nurbsgrid.hh"
#include "dune/iga/utils/igahelpers.hh"
#include <dune/common/parametertreeparser.hh>
#include <dune/vtk/vtkwriter.hh>

int main(int argc, char **argv) {
  std::cout << "ALUGRID_VERBOSITY_LEVEL0: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
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
  const auto u_refine        = gridParameters.get<int>("u_refine");
  const auto v_refine        = gridParameters.get<int>("v_refine");

  const auto E   = materialParameters.get<double>("E");
  const auto nu  = materialParameters.get<double>("nu");
  const auto thk = materialParameters.get<double>("thk");

  double lambdaLoad = 1 * std::pow(thk, 3);

  const int subsample = postProcessParameters.get<int>("subsample");
  std::cout << "ALUGRID_VERBOSITY_LEVEL1: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
  // Log the Parameters
  spdlog::info(
      "Filename: {} \n The following parameters were used: \nMaterial: E {}, nu {} \nRefinements: global {}, u {}, v "
      "{}",
      gridFileName, E, nu, globalRefine, u_refine, v_refine);
  std::cout << "ALUGRID_VERBOSITY_LEVEL2: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
  using Grid     = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;
  std::cout << "ALUGRID_VERBOSITY_LEVEL2: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read(
      "auxiliaryFiles/" + gridFileName, trimGrid, {u_degreeElevate, v_degreeElevate});

  std::cout << "ALUGRID_VERBOSITY_LEVEL4: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
  grid->globalRefine(globalRefine);
  GridView gridView = grid->leafGridView();
  std::cout << "ALUGRID_VERBOSITY_LEVEL5: " << getenv("ALUGRID_VERBOSITY_LEVEL") << std::endl;
  const auto &patchData = grid->getPatch().getPatchData();
  spdlog::info("Degree: u {}, v {}", patchData->degree[0], patchData->degree[1]);

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, gridView.impl().getPreBasis());

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(basis_, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
      dirichletFlags[localView.index(localIndex)] = true;
    });
  });

  /// Declare a vector "fes" of kirchhoff plate 2D elements
  using LinearElasticType = Ikarus::KirchhoffPlate<decltype(basis)>;
  std::vector<LinearElasticType> fes;

  /// Add the linear elastic 2D planar solid elements to the vector "fes"
  for (auto &element : elements(gridView)) {
    auto localView = basis.flat().localView();
    fes.emplace_back(basis, element, E, nu, thk);
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
  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::linearAlgebraFunctions(residualFunction, KFunction),
                                            Ikarus::parameter(D_Glob, lambdaLoad));

  const auto &K    = nonLinOp.derivative();
  const auto &Fext = nonLinOp.value();

  /// solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);

  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);

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

  double totalForce = 0.0;
  for (auto &f : Fext)
    totalForce += f;
  std::cout << "Total Force: " << totalForce << std::endl;

  vtkWriter.write(gridFileName);

  return 0;
}
