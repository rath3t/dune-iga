// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#define CHECK_RESERVEDVECTOR
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>

using namespace Dune;

template <typename LocalView>
auto checkJacobianAndPartial(const LocalView& localView,
                             Dune::FieldVector<double, LocalView::Element::Geometry::mydimension> pos, int i) {
  TestSuite t("Check the point with number" + std::to_string(i));
  auto& fe               = localView.tree().finiteElement();
  const auto& localBasis = fe.localBasis();

  static constexpr int gridDim = LocalView::Element::Geometry::mydimension;

  std::vector<Dune::FieldMatrix<double, 1, gridDim>> referenceGradients;
  localBasis.evaluateJacobian(pos, referenceGradients);

  std::array<unsigned int, gridDim> dir{};
  constexpr double tol = 1e-14;
  for (int d = 0; d < gridDim; ++d) {
    std::vector<Dune::FieldVector<double, 1>> dNdxi;
    dir[d] = 1;
    localBasis.partial(dir, pos, dNdxi);
    for (int j = 0; j < referenceGradients.size(); ++j)
      t.check(Dune::FloatCmp::eq(referenceGradients[j][0][d], dNdxi[j][0], tol))
          << "The failed direction is " << d << " with the values " << referenceGradients[j][0][d] << " " << dNdxi[j][0]
          << "\n The pos is " << pos;
    dir[d] = 0;
  }

  return t;
}

template <typename Basis>
auto checkJacobianAndPartialConsistency(const Basis& basis) {
  TestSuite t;
  constexpr int numTestPos = 13;
  constexpr double posTol  = 1e-8;

  static_assert(Basis::LocalView::Element::Geometry::mydimension == 2);

  // 4 random Gauss points, center position, 0th, 1st, 2nd and 3rd geometry corners and with tolerances
  std::array<Dune::FieldVector<double, 2>, numTestPos> gpPos = {
      {{0.2113248654051871, 0.7886751345948129},
       {0.7886751345948129, 0.2113248654051871},
       {0.2113248654051871, 0.2113248654051871},
       {0.7886751345948129, 0.7886751345948129},
       {0.5, 0.5},
       {0.0, 0.0},
       {1.0, 0.0},
       {0.0, 1.0},
       {1.0, 1.0},
       {0.0 + posTol, 0.0 + posTol},
       {1.0 - posTol, 0.0 + posTol},
       {0.0 + posTol, 1.0 - posTol},
       {1.0 - posTol, 1.0 - posTol}}
  };

  auto localView = basis.localView();

  for (const auto& e : elements(basis.gridView())) {
    localView.bind(e);
    for (int i = 0; i < numTestPos; ++i)
      t.subTest(checkJacobianAndPartial(localView, gpPos[i], i));
  }

  return t;
}

template <template <int, int, typename> typename GridFamily>
auto testNurbsBasis() {
  ////////////////////////////////////////////////////////////////
  // First test
  // A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  // const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {
      {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}
  };

  nurbsPatchData.controlPoints = {
      {        {.p = {0, 0, rad}, .w = 1},         {.p = {0, l, rad}, .w = 1}},
      {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
      // {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
      {        {.p = {rad, 0, 0}, .w = 1},         {.p = {rad, l, 0}, .w = 1}}
  };
  nurbsPatchData.degree = order;

  IGA::PatchGrid<dim, dimworld, GridFamily> grid(nurbsPatchData);
  // grid.globalRefine(1);
  grid.globalRefineInDirection({2, 0});
  grid.degreeElevateOnAllLevels({2, 2});
  // grid.globalRefineInDirection(1, 3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::TestSuite test(TestSuite::ThrowPolicy::AlwaysThrow);

  // Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");
  using GridView = decltype(gridView);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.impl().patchData());

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

  {
    using namespace Functions::BasisFactory;
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView, nurbs());
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    Dune::Functions::forEachBoundaryDOF(basis2, [](auto&& localIndex) {});
    test.subTest(checkJacobianAndPartialConsistency(basis2));
  }

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView);
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkJacobianAndPartialConsistency(basis2));
  }

  {
    // Check basis created via makeBasis
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, nurbs());
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkJacobianAndPartialConsistency(basis2));
  }

  {
    // Check whether a B-Spline basis can be combined with other bases.
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, power<2>(nurbs()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkJacobianAndPartialConsistency(Functions::subspaceBasis(basis2, 0)));
    test.subTest(checkJacobianAndPartialConsistency(Functions::subspaceBasis(basis2, 1)));
  }

  {
    grid.degreeElevateOnAllLevels({1, 1});
    auto gridViewNew = grid.leafGridView();
    // Check lower order basis created via its constructor
    using namespace Functions::BasisFactory;
    Functions::NurbsBasis<GridView> basis2(gridViewNew, nurbs(degreeElevate(1, 1)));
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
    test.subTest(checkJacobianAndPartialConsistency(basis2));
  }

  {
    grid.degreeElevateOnAllLevels({0, 1});
    auto gridViewNew = grid.leafGridView();
    // Check lower order basis created via its constructor
    using namespace Functions::BasisFactory;
    Functions::NurbsBasis<GridView> basis2(gridViewNew, nurbs(degreeElevate(1, 0)));
    test.subTest(checkBasis(basis2, EnableContinuityCheck()));
  }

  return test;
}

template <template <int, int, typename> typename GridFamily>
auto testPrePostDegreeRefinement() {
  TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);
  using namespace Functions::BasisFactory;

  using PatchGrid   = IGA::PatchGrid<2, 2, GridFamily>;
  using GridView    = typename PatchGrid::LeafGridView;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {0, 0}, {0, 0}, {0, 0});
  const auto gridNoRefine = gridFactory.createGrid();

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {1, 1}, {0, 0}, {0, 0});
  const auto gridPreRefine = gridFactory.createGrid();

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {1, 1}, {1, 1}, {0, 0});
  const auto gridPreRefineAndDegree = gridFactory.createGrid();

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {1, 1}, {1, 1}, {1, 1});
  const auto gridPrePostRefineAndDegree = gridFactory.createGrid();

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {0, 0}, {1, 1}, {1, 1});
  const auto gridPostRefineAndDegree = gridFactory.createGrid();

  // Pre Knot Refine gridNoRefine, no it has to have the same basis size as gridPreReinfe
  gridNoRefine->globalRefine(1);

  Functions::NurbsBasis<GridView> basisNoRefine(gridNoRefine->leafGridView(), nurbs());
  Functions::NurbsBasis<GridView> basisPreRefine(gridPreRefine->leafGridView(), nurbs());

  t.check(basisNoRefine.size() == basisPreRefine.size());

  // Now degree elevate gridNoRefine, this should now have the same effect as gridPreRefineAndDegree
  gridNoRefine->degreeElevateOnAllLevels({1, 1});

  Functions::NurbsBasis<GridView> basisNoRefine2(gridNoRefine->leafGridView(), nurbs());
  Functions::NurbsBasis<GridView> basisPreRefineAndDegree(gridPreRefineAndDegree->leafGridView(), nurbs());

  t.check(basisNoRefine2.size() == basisPreRefineAndDegree.size());

  // Now the same with gridPreRefine but without grid elevation but basis elevation
  Functions::NurbsBasis<GridView> basisPreRefine2(gridPreRefine->leafGridView(), nurbs(degreeElevate(1, 1)));

  t.check(basisNoRefine2.size() == basisPreRefine2.size());
  t.check(basisPreRefineAndDegree.size() == basisPreRefine2.size());

  // Now we post refine gridNoRefine, same as gridPrePostRefineAndDegree
  gridNoRefine->globalRefine(1);

  Functions::NurbsBasis<GridView> basisNoRefine3(gridNoRefine->leafGridView(), nurbs());
  Functions::NurbsBasis<GridView> basisPrePostRefineAndDegree(gridPrePostRefineAndDegree->leafGridView(), nurbs());

  t.check(basisNoRefine3.size() == basisPrePostRefineAndDegree.size());

  // Now make k refinement without elevating the grid geometry (this is key)
  auto kRefine = [](const auto& patchData, std::array<int, 2> refinement, std::array<int, 2> degreeElevate) {
    auto newPatchData = patchData;
    for (const auto i : Dune::range(2))
      if (degreeElevate[i] > 0)
        newPatchData = IGA::Splines::degreeElevate(newPatchData, i, degreeElevate[i]);
    for (const auto i : Dune::range(2)) {
      if (refinement[i] > 0) {
        auto newKnots = IGA::Splines::generateRefinedKnots(newPatchData.knotSpans, i, refinement[i]);
        newPatchData  = IGA::Splines::knotRefinement(newPatchData, newKnots, i);
      }
    }
    return newPatchData;
  };

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {0, 0}, {0, 0}, {0, 0});
  const auto freshGrid = gridFactory.createGrid();

  t.check(freshGrid->patchGeometryAtBack().degree()[0] == 1);

  // Now refine both grid and basis

  auto newPatchData = kRefine(freshGrid->patchGeometryAtBack().patchData(), {1, 1}, {1, 1});
  freshGrid->globalRefine(1);
  Functions::NurbsBasis<GridView> freshBasisElevated(freshGrid->leafGridView(), nurbs(newPatchData));
  Functions::NurbsBasis<GridView> basisPostRefineAndDegree(gridPostRefineAndDegree->leafGridView(), nurbs());

  t.check(freshBasisElevated.size() == basisPostRefineAndDegree.size());

  t.check(freshGrid->patchGeometryAtBack().degree()[0] == 1);
  t.check(gridPostRefineAndDegree->patchGeometryAtBack().degree()[0] == 2);

  // test the created basis

  t.subTest(checkBasis(basisNoRefine3, EnableContinuityCheck()));
  t.subTest(checkBasis(basisPostRefineAndDegree, EnableContinuityCheck()));
  t.subTest(checkBasis(freshBasisElevated, EnableContinuityCheck()));
  t.subTest(checkBasis(basisPrePostRefineAndDegree, EnableContinuityCheck()));

  // Now do one more with higher refinements

  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {0, 0}, {0, 0}, {0, 0});
  const auto grid = gridFactory.createGrid();

  auto newPatchData1 = kRefine(grid->patchGeometryAtBack().patchData(), {2, 2}, {1, 1});
  freshGrid->globalRefine(2);
  Functions::NurbsBasis<GridView> basis1(grid->leafGridView(), nurbs());

  t.subTest(checkBasis(basis1, EnableContinuityCheck()));

  return t;
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

  TestSuite t;
  std::cout << "==================================" << std::endl;
  std::cout << "===============TEST IdentityParameterSpace==" << std::endl;
  std::cout << "==================================" << std::endl;

  // t.subTest(testNurbsBasis<IGA::IdentityParameterSpace::PatchGridFamily>());

  std::cout << std::endl;
  std::cout << "==================================" << std::endl;
  std::cout << "===============TEST DefaultParameterSpace===" << std::endl;
  std::cout << "==================================" << std::endl;

  // t.subTest(testNurbsBasis<IGA::DefaultParameterSpace::PatchGridFamily>());

  t.subTest(testPrePostDegreeRefinement<IGA::IdentityParameterSpace::PatchGridFamily>());
  t.subTest(testPrePostDegreeRefinement<IGA::DefaultParameterSpace::PatchGridFamily>());

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
