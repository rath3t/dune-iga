// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "resultcollection.hh"
#include "testhelper.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/parameterspace/identity/parameterspace.hh>
#include <dune/iga/patchgrid.hh>

using namespace Dune::IGA;
using namespace Dune;

template <typename LocalView>
void getShapeFunctionsAndDerivativesWithLocal(const LocalView& localView, const auto& pos,
                                              Eigen::VectorXd& shapeFunctionValues, Eigen::Matrix2Xd& dNdXi,
                                              Eigen::Matrix2Xd& dNdXiPartial, Eigen::Matrix2Xd& dNdX,
                                              Eigen::Matrix3Xd& ddNddXi, Eigen::Matrix3Xd& ddNddX) {
  auto& ele              = localView.element();
  auto& fe               = localView.tree().finiteElement();
  const auto geo         = localView.element().geometry();
  const auto& localBasis = fe.localBasis();
  const auto JInvT       = geo.jacobianInverseTransposed(pos);
  const auto ddXddXi     = toEigen(geo.impl().secondDerivativeOfPosition(pos)).eval();

  std::vector<Dune::FieldMatrix<double, 1, 2>> referenceGradients;
  localBasis.evaluateJacobian(pos, referenceGradients);
  std::vector<Dune::FieldVector<double, 2>> gradients(referenceGradients.size());

  /// copying to an Eigen::Vector
  dNdXi.setZero(Eigen::NoChange, fe.size());
  dNdXiPartial.setZero(Eigen::NoChange, fe.size());
  for (size_t i = 0; i < referenceGradients.size(); ++i) {
    dNdXi(0, i) = referenceGradients[i][0][0];
    dNdXi(1, i) = referenceGradients[i][0][1];
  }

  for (size_t i = 0; i < gradients.size(); i++)
    JInvT.mv(referenceGradients[i][0], gradients[i]);

  std::vector<Dune::FieldVector<double, 1>> dNdxi;
  std::vector<Dune::FieldVector<double, 1>> dNdeta;
  localBasis.partial({1, 0}, pos, dNdxi);
  localBasis.partial({0, 1}, pos, dNdeta);

  /// copying to an Eigen::Vector
  for (auto i = 0U; i < fe.size(); ++i) {
    dNdXiPartial(0, i) = dNdxi[i];
    dNdXiPartial(1, i) = dNdeta[i];
  }

  std::vector<Dune::FieldVector<double, 1>> shapeFunctionValuesVec;
  localBasis.evaluateFunction(pos, shapeFunctionValuesVec);

  /// copying to an Eigen::Vector
  shapeFunctionValues.setZero(fe.size());
  for (int i = 0; i < fe.size(); ++i)
    shapeFunctionValues[i] = shapeFunctionValuesVec[i][0];

  dNdX.setZero(Eigen::NoChange, fe.size());
  ddNddX.setZero(Eigen::NoChange, fe.size());

  /// copying to an Eigen::Vector
  for (size_t i = 0; i < fe.size(); i++) {
    dNdX(0, i) = gradients[i][0];
    dNdX(1, i) = gradients[i][1];
  }

  using Dune::power;
  const auto Jacobian = toEigen(geo.jacobianTransposed(pos)).eval();
  Eigen::Matrix3d JJ;
  JJ << power(Jacobian(0, 0), 2), power(Jacobian(0, 1), 2), 2 * Jacobian(0, 0) * Jacobian(0, 1),
      power(Jacobian(1, 0), 2), power(Jacobian(1, 1), 2), 2 * Jacobian(1, 0) * Jacobian(1, 1),
      Jacobian(0, 0) * Jacobian(1, 0), Jacobian(0, 1) * Jacobian(1, 1),
      Jacobian(0, 0) * Jacobian(1, 1) + Jacobian(0, 1) * Jacobian(1, 0);
  auto JJInv = JJ.inverse();
  ddNddXi.setZero(Eigen::NoChange, fe.size());
  std::vector<Dune::FieldVector<double, 1>> dNdXiXi;
  std::vector<Dune::FieldVector<double, 1>> dNdEtaEta;
  std::vector<Dune::FieldVector<double, 1>> dNdXiEta;

  localBasis.partial({2, 0}, pos, dNdXiXi);
  localBasis.partial({0, 2}, pos, dNdEtaEta);
  localBasis.partial({1, 1}, pos, dNdXiEta);

  /// copying to an Eigen::Vector
  for (auto i = 0U; i < fe.size(); ++i) {
    ddNddXi(0, i) = dNdXiXi[i];
    ddNddXi(1, i) = dNdEtaEta[i];
    ddNddXi(2, i) = dNdXiEta[i];
  }

  ddNddX = JJInv * (ddNddXi - ddXddXi * dNdX);
}

template <typename LocalView>
void evaluateJacobianAndPartial(const LocalView& localView, const auto& pos, Eigen::Matrix2Xd& dNdXiJacobian,
                                Eigen::Matrix2Xd& dNdXiPartial) {
  auto& fe               = localView.tree().finiteElement();
  const auto& localBasis = fe.localBasis();

  std::vector<Dune::FieldMatrix<double, 1, 2>> referenceGradients;
  localBasis.evaluateJacobian(pos, referenceGradients);

  std::vector<Dune::FieldVector<double, 1>> dNdxi;
  std::vector<Dune::FieldVector<double, 1>> dNdeta;
  localBasis.partial({1, 0}, pos, dNdxi);
  localBasis.partial({0, 1}, pos, dNdeta);

  /// copying to an Eigen::Vector
  dNdXiJacobian.setZero(Eigen::NoChange, fe.size());
  for (size_t i = 0; i < referenceGradients.size(); ++i) {
    dNdXiJacobian(0, i) = referenceGradients[i][0][0];
    dNdXiJacobian(1, i) = referenceGradients[i][0][1];
  }

  /// copying to an Eigen::Vector
  dNdXiPartial.setZero(Eigen::NoChange, fe.size());
  for (auto i = 0U; i < fe.size(); ++i) {
    dNdXiPartial(0, i) = dNdxi[i];
    dNdXiPartial(1, i) = dNdeta[i];
  }
}

template <typename TestSuiteType>
void checkValueAndDerivative(TestSuiteType& t, const auto& gpPos, const auto& val, const auto& expectVal,
                             int eleCounter, const std::string& msg) {
  constexpr double tol = 1e-14;
  for (int j = 0; j < val.size(); ++j) {
    checkScalars(t, val[j], expectVal[j],
                 " Mismatch in " + msg + " at pos (" + std::to_string(gpPos[0]) + ", " + std::to_string(gpPos[1]) +
                     ") for element " + std::to_string(eleCounter) + " at j = " + std::to_string(j),
                 tol);
  }
  std::cout << std::endl;
}

template <int refinementLevel>
static auto ShapeFunctionsAndDerivativesIGA() {
  static_assert(refinementLevel == 0 or refinementLevel == 1,
                "ShapeFunctionsAndDerivativesIGA only implemented for refinement levels 0 and 1");
  TestSuite t("ShapeFunctionsAndDerivativesIGA");
  constexpr int gridDim  = 2;
  constexpr int dimWorld = 2;

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 1, 1}, {0, 0, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0}, .w = 1}, {.p = {0, 1}, .w = 1}},
      {{.p = {1, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}
  };

  std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::PatchGrid<gridDim, dimWorld>;

  Dune::IGA::NURBSPatchData<gridDim, dimWorld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = Dune::IGA::Splines::degreeElevate(patchData, i, 1);
  auto grid = std::make_shared<Grid>(patchData);
  grid->globalRefine(refinementLevel);
  auto gridView = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;

  auto CPs = grid->patchGeometryAtBack().patchData().controlPoints.directGetAll();
  auto ks  = grid->patchGeometryAtBack().patchData().knotSpans;
  std::cout << "Control Points" << std::endl;
  for (const auto& cp : CPs)
    std::cout << cp.p << "\t \t";
  std::cout << std::endl;
  std::cout << "Knot Spans" << std::endl;
  for (const auto& knotS : ks) {
    for (const auto& knot : knotS)
      std::cout << knot << "  ";
    std::cout << std::endl;
  }

  auto basis     = makeBasis(gridView, nurbs());
  auto localView = basis.localView();

  constexpr int numTestPos = 6;

  // random Gauss point, center position, 0th, 1st, 2nd and 3rd geometry corners
  std::array<Dune::FieldVector<double, 2>, numTestPos> gpPos = {
      {{0.2113248654051871, 0.7886751345948129}, {0.5, 0.5}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}}
  };

  if constexpr (refinementLevel == 0)
    t.check(gridView.size(0) == 1, "gridView.size(0) is not equal to one.");
  else
    t.check(gridView.size(0) == 4, "gridView.size(0) is not equal to four.");

  int eleCounter = 0;
  for (const auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto [expectN, expectdNdxi, expectdNdeta, expectddNdxixi, expectddNdetaeta, expectddNdxieta] =
        shapeFunctionsAndDerivativesIGAResults<refinementLevel, numTestPos>(eleCounter);
    Eigen::VectorXd shapeFunctionValues;
    Eigen::Matrix2Xd dNdX, dNdXi, dNdXiPartial;
    Eigen::Matrix3Xd ddNddX, ddNddXi;
    for (int i = 0; i < numTestPos; ++i) {
      getShapeFunctionsAndDerivativesWithLocal(localView, gpPos[i], shapeFunctionValues, dNdXi, dNdXiPartial, dNdX,
                                               ddNddXi, ddNddX);
      t.check(shapeFunctionValues.size() == 9) << "shapeFunctionValues.size() = " << shapeFunctionValues.size();
      checkValueAndDerivative(t, gpPos[i], shapeFunctionValues, expectN[i], eleCounter, "shapeFunctions");
      checkValueAndDerivative(t, gpPos[i], dNdXi.row(0), expectdNdxi[i], eleCounter, "dNdxi");
      checkValueAndDerivative(t, gpPos[i], dNdXi.row(1), expectdNdeta[i], eleCounter, "dNdeta");
      checkValueAndDerivative(t, gpPos[i], dNdXi.row(0), dNdXiPartial.row(0), eleCounter, "partial dNdxi");
      checkValueAndDerivative(t, gpPos[i], dNdXi.row(1), dNdXiPartial.row(1), eleCounter, "partial dNdeta");
      checkValueAndDerivative(t, gpPos[i], ddNddXi.row(0), expectddNdxixi[i], eleCounter, "ddNdxixi");
      checkValueAndDerivative(t, gpPos[i], ddNddXi.row(1), expectddNdetaeta[i], eleCounter, "ddNdetaeta");
      checkValueAndDerivative(t, gpPos[i], ddNddXi.row(2), expectddNdxieta[i], eleCounter, "ddNdxieta");
    }
    eleCounter++;
  }
  return t;
}

template <int refinementLevel>
static auto EvaluateJacobianAndPartialTest(int degreeElevated) {
  TestSuite t("EvaluateJacobianAndPartialTest");
  constexpr int gridDim  = 2;
  constexpr int dimWorld = 2;

  const std::array<std::vector<double>, gridDim> knotSpans = {
      {{0, 0, 1, 1}, {0, 0, 1, 1}}
  };

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints = {
      {{.p = {0, 0}, .w = 1}, {.p = {0, 1}, .w = 1}},
      {{.p = {1, 0}, .w = 1}, {.p = {1, 1}, .w = 1}}
  };

  std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::PatchGrid<gridDim, dimWorld>;

  Dune::IGA::NURBSPatchData<gridDim, dimWorld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = Dune::IGA::Splines::degreeElevate(patchData, i, degreeElevated);
  auto grid = std::make_shared<Grid>(patchData);

  for (int ref = 0; ref < refinementLevel; ++ref) {
    auto gridView = grid->leafGridView();
    using namespace Dune::Functions::BasisFactory;

    auto basis     = makeBasis(gridView, nurbs());
    auto localView = basis.localView();

    constexpr int numTestPos = 13;
    constexpr double posTol  = 1e-8;

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

    int eleCounter = 0;
    for (const auto& ele : elements(gridView)) {
      localView.bind(ele);
      Eigen::Matrix2Xd dNdXiJacobian, dNdXiPartial;
      for (int i = 0; i < numTestPos; ++i) {
        evaluateJacobianAndPartial(localView, gpPos[i], dNdXiJacobian, dNdXiPartial);
        checkValueAndDerivative(t, gpPos[i], dNdXiJacobian.row(0), dNdXiPartial.row(0), eleCounter, "partial dNdxi");
        checkValueAndDerivative(t, gpPos[i], dNdXiJacobian.row(1), dNdXiPartial.row(1), eleCounter, "partial dNdeta");
      }
      eleCounter++;
    }
    grid->globalRefine(1);
  }
  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  constexpr int maxDegreeElevate   = 4;
  constexpr int maxRefinementLevel = 5;
  for (int p = 0; p < maxDegreeElevate; ++p)
    t.subTest(EvaluateJacobianAndPartialTest<maxRefinementLevel>(p));

  t.report();
  std::cout << "All derivatives work as expected\n";
  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
