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
    return it-std::begin(U)-1;
  }

  // see https://link.springer.com/book/10.1007/978-3-642-59223-2 A5.1
  template <int dimworld>
  auto curveKnotRefinement(const NURBSPatchData<1, dimworld>& oldData,
                           const std::array<std::vector<double>, 1>& additionalKnots) {
    auto oldKnots = oldData.getKnots()[0];
    oldKnots.erase(oldKnots.begin());
    oldKnots.erase(oldKnots.end()-1);

    std::array<std::vector<double>, 1> newKnotsArray;
    auto& newKnots = newKnotsArray[0];
    newKnots.reserve(oldKnots.size() + additionalKnots[0].size());
    std::merge(oldKnots.begin(), oldKnots.end(), additionalKnots[0].begin(), additionalKnots[0].end(),
               std::back_inserter(newKnots));

    const auto& oldCPs      = line(oldData.getControlPoints());
    const auto& oldweights = line(oldData.getWeights());

    const int newKSize            = additionalKnots[0].size();
    const int p                   = oldData.getOrder()[0];
    const int oldControlPointSize = oldData.getControlPoints().directSize();
//    const int nb_knots            = oldKnots.size();

    const int a = findUpperSpanIndex(p, additionalKnots[0].front(), oldKnots);
    const int b = findUpperSpanIndex(p, additionalKnots[0].back(), oldKnots);

    const int numberOfnewControlPoints = oldControlPointSize + newKSize;

    MultiDimensionNetFVd<1, dimworld> newCPv({numberOfnewControlPoints});
    MultiDimensionNet<1, double> newWeightsv({numberOfnewControlPoints});

    auto newCPs = line(newCPv);
    auto newWeights = line(newWeightsv);

    std::ranges::transform((oldCPs | std::views::take(a)),oldweights,newCPs.begin(),std::multiplies());
    std::ranges::copy_n(oldweights.begin(),a,newWeights.begin());

    std::transform(oldCPs.begin()+b+p-1,oldCPs.end() ,oldweights.begin()+b+p-1,newCPs.begin()+newKSize+b+p-1,std::multiplies());
    std::copy(oldweights.begin()+b+p-1,oldweights.end(),newWeights.begin()+newKSize+b+p-1);

//    for (int i = b+p-1 ; i < oldControlPointSize; i++) {
//      newCP(newKSize + i)      = oldCP(i) * oldweights(i);
//      newWeights(newKSize + i) = oldweights(i);
//    }
    for (auto&& cpOnLine : newCPs)
      std::cout<<cpOnLine<<std::endl;
    for (auto&& cpOnLine : newWeights)
      std::cout<<cpOnLine<<std::endl;

    const int r = newKSize - 1;

    int i = b+2+p-1;
    int k = b+2+p+r;
    for (auto const& currentKnot : std::ranges::reverse_view(additionalKnots[0])) {
      while (currentKnot <= oldKnots[i - 1] && i > a + 1) {
        newCPs[k - p - 1]         = oldCPs[i - p - 1]* oldweights[i - p - 1];
        newWeights[k - p - 1] = oldweights[i - p - 1];
        k -= 1;
        i -= 1;
      }

      newCPs[k - p - 1]         = newCPs[k - p];
      newWeights[k - p - 1] = newWeights[k - p];

      for (std::size_t l = 1; l < p+1; ++l) {
        const std::size_t index = k - p + l;
        auto alpha              = newKnots[k + l - 1] - currentKnot;
        if (std::abs(alpha) < 1e-7) {
          newCPs[index - 1]         = newCPs[index] ;
          newWeights[index - 1] = newWeights[index];
        } else {
          alpha = alpha / (newKnots[k + l - 1] - oldKnots[i + l - p - 1]);
          newCPs[index - 1]         = alpha * newCPs[index - 1] +  newCPs[index] * (1.0 - alpha) ;
          newWeights[index - 1] = alpha* newWeights[index - 1]  + newWeights[index] * (1 - alpha);
        }
      }

      k -= 1;
    }

//    for (int ii = 0; ii < newCP.size()[0]; ++ii)
//      newCPs[ii] /= newWeights[ii];
    std::ranges::transform(newCPs ,newWeights,newCPs.begin(),std::divides());

    newKnotsArray[0].insert(newKnotsArray[0].begin(),newKnotsArray[0].front());
    newKnotsArray[0].push_back(newKnotsArray[0].back());
    std::cerr << "=====newKnotsArray[0]" << std::endl;
    for (int kk = 0; kk < newKnotsArray[0].size(); ++kk)
      std::cerr << newKnotsArray[0][kk] << " ";
    std::cerr << "=======================" << std::endl;
    for (auto&& cpOnLine : newCPs)
      std::cout<<cpOnLine<<std::endl;
    return NURBSPatchData<1, dimworld>(newKnotsArray, newCPv, newWeightsv, oldData.getOrder());
  }
}  // namespace Dune::IGA