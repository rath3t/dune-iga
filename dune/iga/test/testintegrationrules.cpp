// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "testhelper.hh"

#include <cfenv>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/integrationrules/simplexintegrationrulegenerator.hh>
#include <dune/iga/io/vtk/igadatacollector.hh>
#include <dune/iga/parameterspace/default/parameterspace.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;
using namespace Dune::IGA;

template <bool useEle, bool useBoundaryDivisions = false>
requires(!(!useEle and useBoundaryDivisions))
auto testAreaIntegration(Dune::TestSuite& t, const std::string& file_name, int refLevel, double referenceArea,
                         unsigned int splitter = 100) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultParameterSpace::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  using IntegrationRuleGenerator = SimplexIntegrationRuleGenerator<const PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{splitter});
  gridFactory.insertJson(file_name, true, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  auto gv         = grid->leafGridView();

  if constexpr (useEle and not useBoundaryDivisions) {
    DefaultParameterSpace::Preferences::getInstance().boundaryDivisions(20);
    DefaultParameterSpace::Preferences::getInstance().targetAccuracy(1);
  }
  if constexpr (useEle and useBoundaryDivisions) {
    DefaultParameterSpace::Preferences::getInstance().boundaryDivisions(5);
    DefaultParameterSpace::Preferences::getInstance().targetAccuracy(1e-6);
  }

  auto parameters = IntegrationRuleGenerator::Parameters{.boundaryDivisions = 20};

  double area{0};
  for (const auto& ele : elements(gv)) {
    auto qR = [&]() {
      if constexpr (useEle)
        return ele.impl().getQuadratureRule(2 * gridDim);
      else
        return IntegrationRuleGenerator::createIntegrationRule(ele.impl().getLocalEntity(), 2 * gridDim, parameters);
    }();
    auto geo = ele.geometry();

    for (auto [gp, w] : qR) {
      area += geo.integrationElement(gp) * w;
    }
  }
  t.check(FloatCmp::eq(area, referenceArea, power(10.0, -(refLevel + refLevel < 3 ? 1 : 0))))
      << "Area is " << area << ", but should be " << referenceArea << " (" << file_name << ", " << refLevel << ")";
}

auto testBoundaryDivisionsPreference() {
  TestSuite t("", TestSuite::ThrowPolicy::AlwaysThrow);

  DefaultParameterSpace::Preferences::getInstance().targetAccuracy(2);

  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultParameterSpace::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {2, 2});

  const auto grid = gridFactory.createGrid();
  auto gridView   = grid->leafGridView();

  auto getNQP = [&]() {
    return std::accumulate(gridView.begin<0>(), gridView.end<0>(), 0ul,
                           [](size_t lhs, const auto& ele) { return lhs += ele.impl().getQuadratureRule().size(); });
  };

  auto lastnQP = 0ul;
  DefaultParameterSpace::Preferences::getInstance().targetAccuracy(2); // So we go into boundaryDivisions
  for (const auto i : Dune::range(5)) {
    DefaultParameterSpace::Preferences::getInstance().boundaryDivisions(i);
    const auto newQP = getNQP();
    t.check(newQP > lastnQP)
        << "There have to be more Quadrature Points for a higher number of prescribed boundary divisions";
    lastnQP = newQP;
  }

  return t;
}

auto testTargetAccuracyPreference() {
  TestSuite t("", TestSuite::ThrowPolicy::AlwaysThrow);

  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultParameterSpace::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {2, 2});

  const auto grid = gridFactory.createGrid();
  auto gridView   = grid->leafGridView();

  auto getNQP = [&]() {
    return std::accumulate(gridView.begin<0>(), gridView.end<0>(), 0ul,
                           [](size_t lhs, const auto& ele) { return lhs += ele.impl().getQuadratureRule().size(); });
  };

  std::array testAccuracies{1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

  auto lastnQP = 0ul;
  for (const auto acc : testAccuracies) {
    DefaultParameterSpace::Preferences::getInstance().targetAccuracy(acc);
    const auto newQP = getNQP();
    t.check(newQP >= lastnQP)
        << "There have to be more Quadrature Points for a higher or same number of prescribed target accuracy.";
    if (lastnQP == newQP)
      std::cout << "Same amount of QP for " << acc << " as for the accuracy before that.\n";
    lastnQP = newQP;
  }

  return t;
}

template <typename GridImp>
struct AlternativetIntegrationRuleGenerator
{
  bool wasCalled{};
  using Generator = SimplexIntegrationRuleGenerator<GridImp>;
  auto integrationRule() {
    return [&](const auto& element, int order, QuadratureType::Enum /* qt */) {
      wasCalled = true;
      std::cout << "This is the alternative IntegrationRuleGenerator" << std::endl;
      return Generator::createIntegrationRule(element, order);
    };
  }
};

auto testIntegrationRulePolicy() {
  TestSuite t("", TestSuite::ThrowPolicy::AlwaysThrow);

  using PatchGrid   = PatchGrid<2, 2, DefaultParameterSpace::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {2, 2});

  const auto grid = gridFactory.createGrid();
  auto gridView   = grid->leafGridView();

  auto ele =
      *std::ranges::find_if(gridView.begin<0>(), gridView.end<0>(), [](const auto& e) { return e.impl().isTrimmed(); });

  auto qr1 = ele.impl().getQuadratureRule();

  auto alternativeRuleGenerator = AlternativetIntegrationRuleGenerator<PatchGrid>{};
  grid->integrationRule(alternativeRuleGenerator.integrationRule());

  auto qr2 = ele.impl().getQuadratureRule();
  t.check(alternativeRuleGenerator.wasCalled);

  // Pass lambda directly
  auto emptyIntegrationRule = [](const auto& /* ele */, int /* order */, QuadratureType::Enum /* qt */) {
    return QuadratureRule<double, 2>{};
  };

  grid->integrationRule(emptyIntegrationRule);
  auto qr3 = ele.impl().getQuadratureRule();

  t.check(qr3.empty());

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim.ibra", 0, 0.73688393);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim.ibra", 1, 0.73688393);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim.ibra", 2, 0.73688393);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim.ibra", 3, 0.73688393);

  testAreaIntegration<true, true>(t, "auxiliaryfiles/element_trim.ibra", 2, 0.73688393);
  testAreaIntegration<true, true>(t, "auxiliaryfiles/element_trim.ibra", 3, 0.73688393);

  testAreaIntegration<false>(t, "auxiliaryfiles/element_trim.ibra", 2, 0.73688393);
  testAreaIntegration<false>(t, "auxiliaryfiles/element_trim.ibra", 3, 0.73688393);

  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim_xb.ibra", 0, 0.464634714);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim_xb.ibra", 1, 0.464634714);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim_xb.ibra", 2, 0.464634714);
  testAreaIntegration<true>(t, "auxiliaryfiles/element_trim_xb.ibra", 3, 0.464634714);

  testAreaIntegration<true, true>(t, "auxiliaryfiles/element_trim_xb.ibra", 2, 0.464634714);
  testAreaIntegration<true, true>(t, "auxiliaryfiles/element_trim_xb.ibra", 3, 0.464634714);

  testAreaIntegration<false>(t, "auxiliaryfiles/element_trim_xb.ibra", 2, 0.464634714);
  testAreaIntegration<false>(t, "auxiliaryfiles/element_trim_xb.ibra", 3, 0.464634714);

  testAreaIntegration<true>(t, "auxiliaryfiles/trim_multi.ibra", 0, 65.9490597);
  testAreaIntegration<true>(t, "auxiliaryfiles/trim_multi.ibra", 1, 65.9490597);
  testAreaIntegration<true>(t, "auxiliaryfiles/trim_multi.ibra", 2, 65.9490597);
  testAreaIntegration<true>(t, "auxiliaryfiles/trim_multi.ibra", 3, 65.9490597);

  testAreaIntegration<true>(t, "auxiliaryfiles/surface-hole-skew.ibra", 1, 49.069565);
  testAreaIntegration<true>(t, "auxiliaryfiles/surface-hole-skew.ibra", 2, 49.069565);
  testAreaIntegration<true>(t, "auxiliaryfiles/surface-hole-skew.ibra", 3, 49.069565);

  createOutputFolder("out");

  t.subTest(testBoundaryDivisionsPreference());
  t.subTest(testTargetAccuracyPreference());

  t.subTest(testIntegrationRulePolicy());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}