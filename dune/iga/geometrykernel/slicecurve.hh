// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/splines/bsplinealgorithms.hh>
#include <dune/iga/splines/nurbsalgorithms.hh>

namespace Dune::IGANEW {

  template <typename GeoCurve, typename ScalarType = typename GeoCurve::ctype>
  auto sliceCurve(const GeoCurve& geoCurve, std::array<ScalarType, 2> t) -> GeoCurve {
    assert(GeoCurve::mydimension == 1);

    auto findMultiplicity = [](auto t_, auto& knotVec) {
      return std::ranges::count_if(knotVec, [=](double val) {
        return FloatCmp::eq(t_, val, 1e-8);
      });
    };

    auto& patchData = geoCurve.patchData();
    auto domain = geoCurve.domain();
    auto degree = patchData.degree[0];

    if (FloatCmp::eq(domain[0][0], t[0]) and FloatCmp::eq(domain[0][1], t[1]))
      return geoCurve;

    auto newPatchData = typename GeoCurve::PatchData();
    size_t dropCtrptsFront = 0;
    size_t dropCtrptsBack = 0;

    bool skippedFirst = false;
    for (int i = 0; i < 2; ++i) {
      if (FloatCmp::eq(domain[0][i], t[i])) {
        skippedFirst = true;
        continue;
      }

      auto span = Splines::findSpan(degree, t[i], patchData.knotSpans.front());
      auto ks = span - degree + 1;

      auto r = degree + 1 - findMultiplicity(t[i], patchData.knotSpans.front());
      std::vector<ScalarType> newKnots(r);
      std::ranges::fill(newKnots, t[i]);

      newPatchData = Splines::knotRefinement(i == 0 or (i==1 and skippedFirst) ? patchData : newPatchData, newKnots, 0);

      auto& newKnotSpans = newPatchData.knotSpans.front();

      auto span2 = Splines::findSpan(degree, t[i], newKnotSpans) + 1;
      if (i == 1) {
        auto firstNView = std::ranges::take_view(newKnotSpans, span2);
        newPatchData.knotSpans.front() = std::vector<ScalarType>(firstNView.begin(), firstNView.end());
        dropCtrptsBack = patchData.knotSpans.front().size() - (ks + r - 1);
      } else {
        auto lastNView = std::ranges::subrange(std::ranges::end(newKnotSpans) - std::min<long>(span2, newKnotSpans.size()),
        std::ranges::end(newKnotSpans));
        newPatchData.knotSpans.front() = std::vector<ScalarType>(lastNView.begin(), lastNView.end());
        dropCtrptsFront = ks + r - 1;
      }
    }

    // Get rid of excess control points
    auto& controlPoints = newPatchData.controlPoints.directGetAll();

    auto trimmedView = std::ranges::subrange(
       std::ranges::begin(controlPoints) + std::min(dropCtrptsFront, controlPoints.size()),
       std::ranges::end(controlPoints) - std::min(dropCtrptsBack, controlPoints.size())
   );


    newPatchData.controlPoints =
      typename GeoCurve::ControlPointNetType(std::vector<typename GeoCurve::ControlPointType> (trimmedView.begin(), trimmedView.end()));

    return GeoCurve(newPatchData);
  }

}