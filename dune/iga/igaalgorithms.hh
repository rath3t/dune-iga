//
// Created by lex on 21.10.21.
//

#pragma once

#include <dune/iga/NURBSpatch.hh>
namespace Dune::IGA {

  auto findUpperSpanIndex(const int p, const double u, const std::vector<double>& U) {
    auto it = std::upper_bound(U.begin() + p, U.end() - p, u);
    return std::distance(U.begin(), it) - 1;
  }


  template <int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto curveKnotRefinement(const NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData,
                           const std::array<std::vector<double>, 1>& newKnots) {
    using NurbsPatchData = NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using namespace std::ranges;
    auto oldKnotVec = oldData.getKnots()[0];
    oldKnotVec.erase(oldKnotVec.begin());
    oldKnotVec.erase(oldKnotVec.end() - 1);

    std::array<std::vector<double>, 1> newKnotsArray;
    auto& newKnotVec = newKnotsArray[0];
    newKnotVec.reserve(oldKnotVec.size() + newKnots[0].size());
    merge(oldKnotVec, newKnots[0], std::back_inserter(newKnotVec));

    const auto& oldCPs     = line(oldData.getControlPoints());
    const auto& oldweights = oldCPs | views::transform([](auto& cp) -> auto& { return cp.w; });

    const auto newKSize = newKnots[0].size();
    const int p         = oldData.getOrder()[0];

    const int a = findUpperSpanIndex(p, newKnots[0].front(), oldKnotVec);
    const int b = findUpperSpanIndex(p, newKnots[0].back(), oldKnotVec);

    const unsigned int numberOfnewControlPoints = oldData.getControlPoints().directSize() + newKSize;

    typename NurbsPatchData::ControlPointNetType newCPv({numberOfnewControlPoints});

    auto newCPs = line(newCPv);

    using ControlPoint = typename NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };
    std::ranges::transform(oldCPs | views::take(a), newCPs.begin(), scaleCPWithW);

    std::ranges::transform(oldCPs.begin() + b + p - 1, oldCPs.end(), reverse_view(newCPs).begin(), scaleCPWithW);

    auto newCp = newCPs.begin() + b + newKSize;
    auto oldCp = oldCPs.begin() + b;

    auto currentNewKnot = newKnotVec.end();
    auto currentOldKnot = oldKnotVec.end();

    for (auto const& currentAdditionalKnot : reverse_view(newKnots[0])) {
      while (currentAdditionalKnot <= *(currentOldKnot - 1)
             && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1) {
        (*newCp).p = (*oldCp).p * (*oldCp).w;
        (*newCp).w = (*oldCp).w;
        --currentNewKnot;
        --currentOldKnot;
        --newCp;
        --oldCp;
      }

      *newCp = *(newCp + 1);
      ++newCp;
      currentOldKnot -= p;
      for ([[maybe_unused]] int _ : iota_view{0, p}) {
        ++newCp;
        auto& curCPValue       = *newCp;
        auto& leftOfcurCPValue = *(newCp - 1);

        auto alpha = *currentNewKnot - currentAdditionalKnot;
        using std::abs;
        if (abs(alpha) < 1e-7)
          leftOfcurCPValue = curCPValue;
        else {
          alpha            = alpha / (*currentNewKnot - *(currentOldKnot));
          leftOfcurCPValue = alpha * leftOfcurCPValue + curCPValue * (1.0 - alpha);
        }
        ++currentOldKnot;
        ++currentNewKnot;
      }

      currentNewKnot -= p + 1;
      newCp -= p + 2;
    }

    transform(newCPs, newCPs.begin(), [](auto& cp) -> ControlPoint { return {.p = cp.p / cp.w}; });

    newKnotsArray[0].insert(newKnotsArray[0].begin(), newKnotsArray[0].front());
    newKnotsArray[0].push_back(newKnotsArray[0].back());

    return NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.getOrder());
  }
}  // namespace Dune::IGA