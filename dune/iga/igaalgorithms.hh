//
// Created by lex on 21.10.21.
//

#pragma once

#include <dune/iga/NURBSpatch.hh>
namespace Dune::IGA {

  auto findLowerSpanIndex(const int p, const double u, const std::vector<double>& U) {
    auto it = std::lower_bound(U.begin() + p, U.end() - p, u);
    return std::distance(U.begin(), it) - 1;
  }

  auto findUpperSpanIndex(const int p, const double u, const std::vector<double>& U) {
    auto it = std::upper_bound(U.begin() + p, U.end() - p, u);
    return std::distance(U.begin(), it) - 1;
  }

  // see https://link.springer.com/book/10.1007/978-3-642-59223-2 A5.1
  template <int dimworld>
  auto curveKnotRefinement(const NURBSPatchData<1, dimworld>& oldData,
                           const std::array<std::vector<double>, 1>& additionalKnots) {
    auto oldKnots = oldData.getKnots()[0];
    oldKnots.erase(oldKnots.begin());
    oldKnots.erase(oldKnots.end() - 1);

    std::array<std::vector<double>, 1> newKnotsArray;
    auto& newKnots = newKnotsArray[0];
    newKnots.reserve(oldKnots.size() + additionalKnots[0].size());
    std::set_un(oldKnots.begin(), oldKnots.end(), additionalKnots[0].begin(), additionalKnots[0].end(),
               std::back_inserter(newKnots));

    const auto& oldCPs     = line(oldData.getControlPoints());
    const auto& oldweights = line(oldData.getWeights());

    const int newKSize            = additionalKnots[0].size();
    const int p                   = oldData.getOrder()[0];

    const int a = findUpperSpanIndex(p, additionalKnots[0].front(), oldKnots);
    const int b = findUpperSpanIndex(p, additionalKnots[0].back(), oldKnots);

    const int numberOfnewControlPoints = oldData.getControlPoints().directSize() + newKSize;

    MultiDimensionNetFVd<1, dimworld> newCPv({numberOfnewControlPoints});
    MultiDimensionNet<1, double> newWeightsv({numberOfnewControlPoints});

    auto newCPs     = line(newCPv);
    auto newWeights = line(newWeightsv);

    std::ranges::transform((oldCPs | std::views::take(a)), oldweights, newCPs.begin(), std::multiplies());
    std::ranges::copy_n(oldweights.begin(), a, newWeights.begin());

    std::transform(oldCPs.begin() + b + p - 1, oldCPs.end(), oldweights.begin() + b + p - 1,
                   newCPs.begin() + newKSize + b + p - 1, std::multiplies());
    std::copy(oldweights.begin() + b + p - 1, oldweights.end(), newWeights.begin() + newKSize + b + p - 1);

    auto newCp     = newCPs.begin() + b + newKSize ;
    auto newWeight = newWeights.begin() + b  + newKSize;
    auto oldCp     = oldCPs.begin() + b;
    auto oldWeight = oldweights.begin() + b;

    auto currentNewKnot = newKnots.end();
    auto currentOldKnot = oldKnots.end();

    for (auto const& currentAdditionalKnot : std::ranges::reverse_view(additionalKnots[0])) {
      while (currentAdditionalKnot <= *(currentOldKnot-1) && std::distance(oldKnots.begin(),currentOldKnot) > a + 1) {
        *newCp     = *(oldCp) * (*oldWeight);
        *newWeight = *oldWeight;
        --currentNewKnot;
        --currentOldKnot;
        --newCp;
        --newWeight;
        --oldCp;
        --oldWeight;
      }

      *newCp     = *(newCp + 1);
      *newWeight = *(newWeight + 1);

      ++newCp;
      ++newWeight;
      currentOldKnot -= p ;
      for (int jj : std::ranges::iota_view{0, p}) {
        ++newCp;
        ++newWeight;
        auto& curCPValue           = *newCp;
        auto& leftOfcurCPValue     = *(newCp - 1);
        auto& curWeightValue       = *newWeight;
        auto& leftOfcurWeightValue = *(newWeight - 1);

        auto alpha = *currentNewKnot - currentAdditionalKnot;
        if (std::abs(alpha) < 1e-7) {
          leftOfcurCPValue     = curCPValue;
          leftOfcurWeightValue = curWeightValue;
        } else {
          alpha                = alpha / (*currentNewKnot - *(currentOldKnot));
          leftOfcurCPValue     = alpha * leftOfcurCPValue + curCPValue * (1.0 - alpha);
          leftOfcurWeightValue = alpha * leftOfcurWeightValue + curWeightValue * (1 - alpha);
        }
        ++currentOldKnot;
        ++currentNewKnot;
      }

      currentNewKnot-=p+1;
      newCp -= p + 2;
      newWeight -= p + 2;
    }

    std::ranges::transform(newCPs, newWeights, newCPs.begin(), std::divides());

    newKnotsArray[0].insert(newKnotsArray[0].begin(), newKnotsArray[0].front());
    newKnotsArray[0].push_back(newKnotsArray[0].back());

    return NURBSPatchData<1, dimworld>(newKnotsArray, newCPv, newWeightsv, oldData.getOrder());
  }
}  // namespace Dune::IGA