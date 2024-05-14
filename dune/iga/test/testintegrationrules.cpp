// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/integrationrules/simplexintegrationrulegenerator.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune;
using namespace Dune::IGA;

template <bool useEle, bool useBoundaryDivisions = false>
requires(!(!useEle and useBoundaryDivisions))
auto testAreaIntegration(Dune::TestSuite& t, const std::string& file_name, int refLevel, double referenceArea) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  using IntegrationRuleGenerator = DefaultTrim::SimplexIntegrationRuleGenerator<const PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson(file_name, true, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  auto gv         = grid->leafGridView();

  if constexpr (useEle and not useBoundaryDivisions) {
    Preferences::getInstance().boundaryDivisions(20);
    Preferences::getInstance().targetAccuracy(1);
  }
  if constexpr (useEle and useBoundaryDivisions) {
    Preferences::getInstance().boundaryDivisions(5);
    Preferences::getInstance().targetAccuracy(1e-4);
  }

  auto parameters = IntegrationRuleGenerator::Parameters{.boundaryDivisions = 20};

  double area{0};
  for (const auto& ele : elements(gv)) {
    auto qR = [&]() {
      if constexpr (useEle)
        return ele.impl().getQuadratureRule(2 * gridDim);
      else
        return IntegrationRuleGenerator::createIntegrationRule(ele, 2 * gridDim, parameters);
    }();
    auto geo = ele.geometry();

    for (auto [gp, w] : qR) {
      area += geo.integrationElement(gp) * w;
    }
  }
  t.check(FloatCmp::eq(area, referenceArea, power(10.0, -(refLevel + refLevel < 3 ? 1 : 0))))
      << "Area is " << area << ", but should be " << referenceArea << " (" << file_name << ", " << refLevel << ")";
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

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

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}