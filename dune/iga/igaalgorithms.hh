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
    //    oldKnotVec.erase(oldKnotVec.begin());
    //    oldKnotVec.erase(oldKnotVec.end() - 1);

    std::array<std::vector<double>, 1> newKnotsArray;
    auto& newKnotVec = newKnotsArray[0];
    newKnotVec.reserve(oldKnotVec.size() + newKnots[0].size());
    merge(oldKnotVec, newKnots[0], std::back_inserter(newKnotVec));

    auto oldCPs = line(oldData.getControlPoints());

    const auto newKSize = newKnots[0].size();

    const unsigned int numberOfnewControlPoints = oldData.getControlPoints().directSize() + newKSize;
    typename NurbsPatchData::ControlPointNetType newCPv({numberOfnewControlPoints});

    using ControlPoint = typename NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

    auto newCPs = line(newCPv);

    const int p = oldData.getOrder()[0];
    const int a = findUpperSpanIndex(p, newKnots[0].front(), oldKnotVec);
    const int b = findUpperSpanIndex(p, newKnots[0].back(), oldKnotVec);
    std::ranges::transform(oldCPs | views::take(a - p + 2), newCPs.begin(), scaleCPWithW);

    std::ranges::transform(oldCPs.begin() + b + 1, oldCPs.end(), newCPs.begin() + newKSize + b + 1, scaleCPWithW);

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
        const auto& curCPValue = *newCp;
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

    std::ranges::for_each(newCPs, [](auto& cp) {
      std::cout << "cp.p: " << cp.p << " "
                << "cp.w: " << cp.w << std::endl;
    });
    std::ranges::for_each(newKnotsArray[0], [](auto& kn) { std::cout << kn << " "; });
    //    newKnotsArray[0].insert(newKnotsArray[0].begin(), newKnotsArray[0].front());
    //    newKnotsArray[0].push_back(newKnotsArray[0].back());

    return NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.getOrder());
  }

  template <int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto surfaceKnotRefinement(const NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData,
                             const std::array<std::vector<double>, 2>& newKnots) {
    using NurbsPatchData = NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using namespace std::ranges;
    std::array<std::vector<double>, 2> newKnotsArray;
    typename NurbsPatchData::ControlPointNetType newCPv(oldData.getControlPoints().size());

    auto oldCPv = oldData.getControlPoints();
    for (int refDirection = 0; refDirection < 2; ++refDirection) {
      auto oldKnotVec  = oldData.getKnots()[refDirection];
      auto& newKnotVec = newKnotsArray[refDirection];

      newKnotVec.reserve(oldKnotVec.size() + newKnots[refDirection].size());
      merge(oldKnotVec, newKnots[refDirection], std::back_inserter(newKnotVec));

      const auto newKSize = newKnots[refDirection].size();
      std::array<unsigned int, 2> numberOfCPperDir{};

      numberOfCPperDir[refDirection] = oldCPv.size()[refDirection] + newKSize;

      for (int i = 0; i < 2; ++i) {
        if (i == refDirection) continue;
        numberOfCPperDir[i] = oldCPv.size()[i];
      }
      newCPv.enlarge(numberOfCPperDir);

      using ControlPoint = typename NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;

      auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

      const int degree = oldData.getOrder()[refDirection];
      const int a      = findUpperSpanIndex(degree, newKnots[refDirection].front(), oldKnotVec);
      const int b      = findUpperSpanIndex(degree, newKnots[refDirection].back(), oldKnotVec);

      const int otherLine = (refDirection == 0) ? 1 : 0;
      for (unsigned int i = 0; i < a - degree + 1; ++i)
        std::ranges::transform(line(oldCPv, otherLine, at(i)), line(newCPv, otherLine, at(i)).begin(), scaleCPWithW);

      for (unsigned int i = b; i < oldCPv.size()[refDirection]; ++i)
        std::ranges::transform(line(oldCPv, otherLine, at(i)), line(newCPv, otherLine, at(i + newKSize)).begin(),
                               scaleCPWithW);

      int k           = b + newKSize;
      int i           = b;
      auto newVLineAt = [&newCPv, &otherLine](const int index) { return line(newCPv, otherLine, at(index)); };
      auto oldVLineAt = [&oldCPv, &otherLine](const int index) { return line(oldCPv, otherLine, at(index)); };
      auto newCpLine  = newVLineAt(k);
      auto oldCpLine  = oldVLineAt(i);

      auto currentNewKnot = newKnotVec.end();
      auto currentOldKnot = oldKnotVec.end();

      for (auto const& currentAdditionalKnot : reverse_view(newKnots[refDirection])) {
        while (currentAdditionalKnot <= *(currentOldKnot - 1)
               && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1) {
          std::ranges::transform(oldCpLine, newCpLine.begin(), scaleCPWithW);
          --currentNewKnot;
          --currentOldKnot;
          newCpLine = newVLineAt(--k);
          oldCpLine = oldVLineAt(--i);
        }

        std::ranges::copy(newVLineAt(k + 1), newCpLine.begin());
        newCpLine = newVLineAt(++k);
        currentOldKnot -= degree;
        for ([[maybe_unused]] int _ : iota_view{0, degree}) {
          newCpLine = newVLineAt(++k);

          auto alpha = *currentNewKnot - currentAdditionalKnot;
          using std::abs;
          if (abs(alpha) < 1e-7)
            std::ranges::copy(newCpLine, newVLineAt(k - 1).begin());
          else {
            alpha = alpha / (*currentNewKnot - *(currentOldKnot));
            std::ranges::transform(newCpLine, newVLineAt(k - 1), newVLineAt(k - 1).begin(),
                                   [&alpha](auto& cp, auto& cpL) { return alpha * cpL + cp * (1.0 - alpha); });
          }
          ++currentOldKnot;
          ++currentNewKnot;
        }

        currentNewKnot -= degree + 1;
        k -= degree + 2;
        newCpLine = newVLineAt(k);
      }

      std::ranges::transform(newCPv.directGetAll(), newCPv.directGetAll().begin(),
                             [](auto& cp) -> ControlPoint { return {.p = cp.p / cp.w, .w = cp.w}; });

      oldCPv = newCPv;
    }
    //    std::ranges::for_each(newCPv.directGetAll(), [](auto& cp) {
    //      std::cout << "cp.p: " << cp.p << " "
    //                << "cp.w: " << cp.w << std::endl;
    //    });
    return NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.getOrder());
  }

}  // namespace Dune::IGA