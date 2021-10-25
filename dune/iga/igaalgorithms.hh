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
    return std::distance(U.begin(), it);
  }

  // see https://link.springer.com/book/10.1007/978-3-642-59223-2 A5.1
  template <int dimworld>
  auto curveKnotRefinement(const NURBSPatchData<1, dimworld>& oldData,
                           const std::array<std::vector<double>, 1>& additionalKnots) {
    const auto& oldKnots = oldData.getKnots()[0];

    std::array<std::vector<double>, 1> newKnotsArray;
    auto& newKnots = newKnotsArray[0];
    newKnots.resize(oldKnots.size() + additionalKnots[0].size());
//    std::merge(oldKnots.begin(), oldKnots.end(), additionalKnots[0].begin(), additionalKnots[0].end(),
//               std::back_inserter(newKnots));

    const auto& oldCP      = oldData.getControlPoints();
    const auto& oldweights = oldData.getWeights();

    const int addituonalKnotSize  = additionalKnots[0].size();
    const int p                   = oldData.getOrder()[0];
    const int oldControlPointSize = oldData.getControlPoints().directSize();
    const int nb_knots            = oldKnots.size();

    const int a = findLowerSpanIndex(p, additionalKnots[0].front(), oldKnots);
    const int b = findUpperSpanIndex(p, additionalKnots[0].back(), oldKnots)-1;

    const int numberOfnewControlPoints = oldControlPointSize + addituonalKnotSize;
//    const int numberOfNewKnots         = nb_knots + addituonalKnotSize + 2;

    MultiDimensionNetFVd<1, dimworld> newCP({numberOfnewControlPoints});
    MultiDimensionNetFVd<1, 1> newWeights({numberOfnewControlPoints});

    for (int i = 0; i <= a - p; i++) {
      newCP(i)         = oldCP(i) * oldweights(i)[0];
      newWeights(i)[0] = oldweights(i)[0];
    }

    for (int i = b ; i < oldControlPointSize; i++) {
      newCP(addituonalKnotSize + i)      = oldCP(i) * oldweights(i)[0];
      newWeights(addituonalKnotSize + i) = oldweights(i)[0];
    }
    std::cerr << "============OLDWeight before algo===========" << std::endl;
    for (int kk = 0; kk < oldweights.size()[0]; ++kk)
      std::cerr << oldweights(kk)[0]  << std::endl;
    std::cerr << "============NEWWeight before algo===========" << std::endl;
    for (int kk = 0; kk < newWeights.size()[0]; ++kk)
      std::cerr << newWeights(kk)[0]  << std::endl;
    std::cerr << "=============OLDKNOTS==========" << std::endl;
    for (auto& addKnot : oldKnots)
      std::cerr << addKnot << " ";
    std::cerr << "=======================" << std::endl;
    std::cerr << "============CPs before algo===========" << std::endl;
    for (int kk = 0; kk < newCP.size()[0]; ++kk)
      std::cerr << newCP(kk)[0] << " " << newCP(kk)[1] << " " << newCP(kk)[2] << std::endl;
    std::cerr << "============oldCPs before algo===========" << std::endl;
    for (int kk = 0; kk < oldCP.size()[0]; ++kk)
      std::cerr << oldCP(kk)[0] << " " << oldCP(kk)[1] << " " << oldCP(kk)[2] << std::endl;
    std::cerr << "=======================" << std::endl;
    std::cerr << "numberOfnewControlPoints:" << numberOfnewControlPoints << std::endl;
    //    const int n = oldControlPointSize - 1;
    //    const int m = n + p + 1;
    const int r = addituonalKnotSize - 1;

    int i = b+2+p-1;
    int k = b+2+p+r-1;
    std::cout<<"newKnots.size"<<newKnots.size()<<std::endl;
    std::cout<<"k"<<k<<std::endl;
    std::cout<<"newKnots.size"<<oldKnots.size()<<std::endl;
    std::cout<<"i"<<i<<std::endl;
    //    const std::size_t r = nb_knots_u_to_insert - 1;
    // Create new knot span vector
    for (int ii = 0; ii < a + 1; ++ii) {
      newKnots[ii] = oldKnots[ii];
    }
    for (int ii = b + p - 1; ii < oldKnots.size(); ++ii) {
      newKnots[ii + addituonalKnotSize] = oldKnots[ii];
    }



    for (auto const& currentKnot : std::ranges::reverse_view(additionalKnots[0])) {
      std::cerr << "i: " << i << std::endl;
      //      std::cerr<<"j: "<<j<<std::endl;
      std::cerr << "r: " << r << std::endl;
      std::cerr << "a: " << a << std::endl;
      std::cerr << "b: " << b << std::endl;
      std::cerr << "k: " << k << std::endl;
      while (currentKnot <= oldKnots[i - 1] && i > a + 1) {
        newCP(k - p - 1)         = oldCP(i - p - 1)* oldweights(i - p - 1)[0];
        newWeights(k - p - 1)[0] = oldweights(i - p - 1)[0];
        newKnots[k - 1] = oldKnots[i - 1];
        k -= 1;
        i -= 1;
      }

      newCP(k - p - 1)         = newCP(k - p);
      newWeights(k - p - 1)[0] = newWeights(k - p)[0];
      std::cerr << "=====FinalCPs" << std::endl;
      for (int kk = 0; kk < newCP.size()[0]; ++kk)
        std::cerr << newCP(kk)[0] << " ";

      for (std::size_t l = 1; l < p+1; ++l) {
        const std::size_t index = k - p + l;
        auto alpha              = newKnots[k + l - 1] - currentKnot;
        std::cerr << "alpha: " << alpha << std::endl;
        std::cerr << "index: " << index << std::endl;
        std::cerr << "currentKnot: " << currentKnot << std::endl;
        std::cerr << "newKnots[k + l - 1]: " << newKnots[k + l - 1] << std::endl;
        if (std::abs(alpha) < 1e-7) {
          newCP(index - 1)         = newCP(index) ;
          newWeights(index - 1)[0] = newWeights(index)[0];
        } else {
          alpha = alpha / (newKnots[k + l - 1] - oldKnots[i + l - p - 1]);
          std::cerr << "alpha2: " << alpha << std::endl;
          std::cerr << "oldKnots[i + l - p - 1]: " << oldKnots[i + l - p - 1] << std::endl;
          std::cerr << "newWeights(index - 1)[0]: " << newWeights(index - 1)[0] << std::endl;
          std::cerr << "newWeights(index)[0]: " << newWeights(index)[0] << std::endl;
          std::cerr << "newCP(index)[0]: " << newCP(index)[0] << std::endl;
          std::cerr << "newCP(index-1)[0]: " << newCP(index-1)[0] << std::endl;
          newCP(index - 1)         = alpha * newCP(index - 1) + (1.0 - alpha) * newCP(index);
          newWeights(index - 1)[0] = newWeights(index - 1)[0] * alpha + newWeights(index)[0] * (1 - alpha);

        }
      }
//      break;
      newKnots[k - 1] = currentKnot;
      k -= 1;
    }
    std::cerr << "=====FinalWeights" << std::endl;
    for (int kk = 0; kk < newWeights.size()[0]; ++kk)
      std::cerr << newWeights(kk)[0] << " ";

    std::cerr << "=====FinalCPsBeforeDiv" << std::endl;
    for (int kk = 0; kk < newCP.size()[0]; ++kk)
      std::cerr << newCP(kk)[0] << " ";
    std::cerr << "=======================" << std::endl;

    for (int ii = 0; ii < newCP.size()[0]; ++ii)
      newCP(ii) /= newWeights(ii)[0];

    std::cerr << "=====FinalCPs" << std::endl;
    for (int kk = 0; kk < newCP.size()[0]; ++kk)
      std::cerr << newCP(kk)[0] << " " << newCP(kk)[1] << " " << newCP(kk)[2] << std::endl;
    std::cerr << "=======================" << std::endl;

    std::cerr << "=====newKnotsArray[0]" << std::endl;
    for (int kk = 0; kk < newKnotsArray[0].size(); ++kk)
      std::cerr << newKnotsArray[0][kk] << " ";
    std::cerr << "=======================" << std::endl;
    return NURBSPatchData<1, dimworld>(newKnotsArray, newCP, newWeights, oldData.getOrder());
  }
}  // namespace Dune::IGA