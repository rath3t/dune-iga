// SPDX-FileCopyrightText: 2022 Alexander Müller mueller@ibb.uni-stuttgart.de
//
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <config.h>

#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>

#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/ibraReader.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/basis.hh>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>

// #include <Eigen/Eigen>

#include <ikarus/utils/observer/controlVTKWriter.hh>

#include <dune/iga/igaDataCollector.h>
#include <dune/vtk/vtkwriter.hh>

#include "linearElasticTrimmed.h"

template <class GridView, int ncomp = 1>
class BoundaryPatchEnclosingVerticesPropertyTrimmed
{
  typedef typename GridView::IndexSet IndexSet;
  static const int dim = GridView::dimension;
 public:
  typedef typename GridView::Intersection Intersection;

  /** \brief Create property from a marker vector
   */
  BoundaryPatchEnclosingVerticesPropertyTrimmed(const GridView& gridView, const Dune::BitSetVector<ncomp>& vertices) :
                                                                                                                indexSet_(gridView.indexSet()),
                                                                                                                vertices_(vertices)
  {}

  /** \brief Check if intersection is enclosed by vertices in the vector
   */
  bool operator() (const Intersection& i) const
  {
    const auto inside = i.inside();
    int localFaceIndex = i.indexInInside();

    if (inside.impl().isTrimmed())
      return false;

    auto refElement = Dune::ReferenceElements<double, dim>::general(inside.type());

    // Get global node ids
    int n = refElement.size(localFaceIndex, 1, dim);

    // Using ReferenceElement::subEntity is OK here, because we loop
    // over all subEntities (i.e. sub-vertices) and just return false
    // if _any_ of them is not marked.
    for (int i=0; i<n; i++) {
      int localVertexIndex = refElement.subEntity(localFaceIndex, 1, i, dim);
      if (not(vertices_[indexSet_.subIndex(inside, localVertexIndex, dim)].any()))
        return false;
    }
    return true;
  }

 private:
  const IndexSet& indexSet_;
  const Dune::BitSetVector<ncomp>& vertices_;
};

namespace Dune::Functions {
  template <class Basis, class F>
    requires(std::same_as<typename Basis::GridView, Dune::IGA::NURBSGrid<2, 2>::LeafGridView>)
  void forEachUntrimmedBoundaryDOF(const Basis &basis, F &&f) {
    auto localView       = basis.localView();
    auto seDOFs          = subEntityDOFs(basis);
    const auto &gridView = basis.gridView();
    for (auto &&element : elements(gridView))
      if (element.hasBoundaryIntersections() and !element.impl().isTrimmed()) {
        localView.bind(element);
        for (const auto &intersection : intersections(gridView, element))
          if (intersection.boundary())
            for (auto localIndex : seDOFs.bind(localView, intersection))
              f(localIndex, localView, intersection);
      }
  }
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);

  constexpr int gridDim = 2;
  constexpr int worldDim = 2;
  double lambdaLoad = 1;

  /// Read Parameter
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree &gridParameters = parameterSet.sub("GridParameters");
  const Dune::ParameterTree &materialParameters = parameterSet.sub("MaterialParameters");

  const auto gridFileName = gridParameters.get<std::string>("filename");
  const bool trimGrid = gridParameters.get<bool>("trim");
  const auto u_degreeElevate = gridParameters.get<int>("u_degreeElevate");
  const auto v_degreeElevate = gridParameters.get<int>("v_degreeElevate");
  const auto globalRefine = gridParameters.get<int>("globalRefine");

  const auto E = materialParameters.get<double>("E");
  const auto nu = materialParameters.get<double>("nu");

  /// Create Grid

  using Grid = Dune::IGA::NURBSGrid<gridDim, worldDim>;
  using GridView = Dune::IGA::NURBSGrid<gridDim, worldDim>::LeafGridView;

  std::shared_ptr<Grid> grid = Dune::IGA::IbraReader<gridDim, worldDim>::read("auxiliaryFiles/"+gridFileName,
                                                                              trimGrid, {u_degreeElevate, v_degreeElevate});
  grid->globalRefine(globalRefine);
  GridView gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(gridView.impl().getPreBasis(), FlatInterleaved()));

  // Clamp the left boundary
  Ikarus::DirichletValues dirichletValues(basis.flat());

  dirichletValues.fixDOFs([](auto &basis_, auto &dirichletFlags) {
    Dune::Functions::forEachUntrimmedBoundaryDOF(basis_,
                                        [&](auto &&localIndex, auto &&localView, auto &&intersection) {
                                          if (std::abs(intersection.geometry().center()[0]) < 1e-8)
                                            dirichletFlags[localView.index(localIndex)] = true;
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
    F[0]              = lamb;
    return F;
  };

  const auto &indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto neumannPredicate = [](auto &vertex) -> bool {
    return std::isgreaterequal(vertex[0], 10 - 1e-8);
  };
  for (auto &&vertex: vertices(gridView)) {
    auto coords = vertex.geometry().corner(0);
    neumannVertices[indexSet.index(vertex)] = neumannPredicate(coords);
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView);

  // create property for insertion by vertex vector
  BoundaryPatchEnclosingVerticesPropertyTrimmed<GridView, 1> prop(gridView, neumannVertices);
  neumannBoundary.insertFacesByProperty(prop);


  /// Add the linear elastic 2D planar solid elements decorated with EAS to the vector "fes"
  for (auto &element: elements(gridView)) {
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
  auto startAssembly = std::chrono::high_resolution_clock::now();
  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::linearAlgebraFunctions(residualFunction, KFunction),
                                                    Ikarus::parameter(D_Glob, lambdaLoad));
  auto stopAssembly = std::chrono::high_resolution_clock::now();
  auto durationAssembly = duration_cast<std::chrono::milliseconds>(stopAssembly - startAssembly);
  spdlog::info("The assembly took {:>6d} milliseconds with {:>7d} dofs",
               durationAssembly.count(), basis.flat().size());

  const auto &K = nonLinOp.derivative();
  const auto &Fext = nonLinOp.value();

  /// solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  auto startSolver = std::chrono::high_resolution_clock::now();

  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);
  auto stopSolver = std::chrono::high_resolution_clock::now();
  auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
  spdlog::info("The solver took {} milliseconds ",
               durationSolver.count());

//  /// Stresses
//  Ikarus::ResultRequirements resultRequirements = Ikarus::ResultRequirements().
//                                                  insertGlobalSolution(Ikarus::FESolutions::displacement, D_Glob).
//                                                  insertParameter(Ikarus::FEParameter::loadfactor, lambdaLoad).
//                                                  addResultRequest(Ikarus::ResultType::cauchyStress);
//
//
//  auto von_mieses2d = [](const auto& sigma){
//    const auto s_x = sigma(0, 0);
//    const auto s_y = sigma(1, 0);
//    const auto s_xy = sigma(2, 0);
//
//    return std::sqrt(std::pow(s_x, 2) + std::pow(s_y, 2) - s_x * s_y + 3 * std::pow(s_xy, 2));
//  };
//
//  Ikarus::ResultTypeMap<double> res;
//
//  auto  elementMapper = Dune::MultipleCodimMultipleGeomTypeMapper(gridView, Dune::mcmgElementLayout());
//  const int numEle = gridView.size(0);
//  std::vector<double> sig_x(numEle);
//  std::vector<double> sig_y(numEle);
//  std::vector<double> sig_xy(numEle);
//  std::vector<double> sig_vm(numEle);
//
//  for (const auto &ele: elements(gridView)) {
//    auto eleID = ele.impl().getIndex();
//    fes[eleID].calculateAt(resultRequirements, {0.5, 0.5}, res);
//    auto sigma = res.getResult(Ikarus::ResultType::cauchyStress);
//
//    sig_x[eleID] = sigma(0, 0);
//    sig_y[eleID] = sigma(1, 0);
//    sig_xy[eleID] = sigma(2, 0);
//    sig_vm[eleID] = von_mieses2d(sigma);
//  }

  /// Postprocess
  auto dispGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob);

  auto forceGlobalFunc
      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), Fext);

  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector(gridView);
  Dune::VtkUnstructuredGridWriter vtkWriter(dataCollector, Dune::Vtk::FormatTypes::ASCII);

  vtkWriter.addPointData(dispGlobalFunc,
                          Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  vtkWriter.addPointData(forceGlobalFunc,
                         Dune::VTK::FieldInfo("external force", Dune::VTK::FieldInfo::Type::vector, 2));


  double totalForce = 0.0;
  for (auto& f : Fext)
    totalForce += f;
  std::cout << "Total Force: " << totalForce << std::endl;

//  vtkWriter.addCellData(sig_x, "sig_x");
//  vtkWriter.addCellData(sig_y, "sig_y");
//  vtkWriter.addCellData(sig_xy, "sig_xy");
//  vtkWriter.addCellData(sig_vm, "sig_vm");

  vtkWriter.write(gridFileName);


  return 0;
}
