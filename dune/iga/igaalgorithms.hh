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

    std::ranges::transform(oldCPs.begin() + b + 1, oldCPs.end(), newCPs.begin() +newKSize+ b + 1, scaleCPWithW);

    auto newCp = newCPs.begin() + b + newKSize;
    auto oldCp = oldCPs.begin() + b;

    auto currentNewKnot = newKnotVec.end();
    auto currentOldKnot = oldKnotVec.end();

    for (auto const& currentAdditionalKnot : reverse_view(newKnots[0])) {
      while (currentAdditionalKnot <= *(currentOldKnot - 1)
             && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1)
      {
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

std::ranges::for_each(newCPs,[](auto& cp){std::cout<<"cp.p: "<<cp.p<<" "<<"cp.w: "<<cp.w<<std::endl;});
std::ranges::for_each( newKnotsArray[0],[](auto& kn){std::cout<<kn<<" ";});
//    newKnotsArray[0].insert(newKnotsArray[0].begin(), newKnotsArray[0].front());
//    newKnotsArray[0].push_back(newKnotsArray[0].back());

    return NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.getOrder());
  }

  template <int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto surfaceKnotRefinement(const NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData,
                             const std::array<std::vector<double>, 2>& newKnots) {
    using NurbsPatchData = NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using namespace std::ranges;
    auto oldKnotVecU = oldData.getKnots()[0];
//    oldKnotVecU.erase(oldKnotVecU.begin());
//    oldKnotVecU.erase(oldKnotVecU.end() - 1);

    auto oldKnotVecV = oldData.getKnots()[1];
//    oldKnotVecV.erase(oldKnotVecV.begin());
//    oldKnotVecV.erase(oldKnotVecV.end() - 1);

    std::array<std::vector<double>, 2> newKnotsArray;
    auto& newKnotVecU = newKnotsArray[0];
    auto& newKnotVecV = newKnotsArray[1];
    newKnotVecV       = oldKnotVecV;
    newKnotVecU.reserve(oldKnotVecU.size() + newKnots[0].size());
    merge(oldKnotVecU, newKnots[0], std::back_inserter(newKnotVecU));
    const auto& oldCPv = oldData.getControlPoints();
    //    const auto& oldCPs     = line(oldData.getControlPoints());

    const auto newKSize = newKnots[0].size();

    const unsigned int numberOfnewControlPointsU = oldData.getControlPoints().size()[0] + newKSize;
    const unsigned int numberOfnewControlPointsV = oldData.getControlPoints().size()[1];
    typename NurbsPatchData::ControlPointNetType newCPv({numberOfnewControlPointsU, numberOfnewControlPointsV});

    using ControlPoint = typename NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

    //    auto newCPs = line(newCPv);

    const int degreeU = oldData.getOrder()[0];
    const int a       = findUpperSpanIndex(degreeU, newKnots[0].front(), oldKnotVecU);
    const int b       = findUpperSpanIndex(degreeU, newKnots[0].back(), oldKnotVecU);

    for (unsigned int i = 0; i < a - degreeU + 1; ++i)
      std::ranges::transform(line(oldCPv,1, at(i)), line(newCPv,1, at(i)).begin(), scaleCPWithW);
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    for (unsigned int i = 0; i < newCPv.size()[0]; ++i)
      for (auto& cp : line(newCPv,1, at(i)))
        std::cout << cp.p << " " << cp.w << std::endl;
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "===========oldCPv.size()[0]=============="<<oldCPv.size()[0] << std::endl;
    std::cout << "===========b + degreeU - 1=============="<<b + degreeU - 1<< std::endl;
    for (unsigned int i = b ; i < oldCPv.size()[0]; ++i)
      std::ranges::transform(line(oldCPv,1, at(i)),
                             line(newCPv,1, at(i+newKSize)).begin(), scaleCPWithW);
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    for (unsigned int i = 0; i < newCPv.size()[0]; i++)
      for (auto&& cp : line(newCPv,1, at(i)))
        std::cout << cp.p << " " << cp.w << std::endl;
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "=========================" << std::endl;
    for (unsigned int i = 0; i < oldCPv.size()[0]; i++)
      for (auto&& cp : line(oldCPv,1, at(i)))
        std::cout << cp.p << " " << cp.w << std::endl;
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;

    int k = b +newKSize;
    int i = b;
    auto newVLineAt = [&newCPv](const int index) { return line(newCPv,1, at(index));};
    auto oldVLineAt = [&oldCPv](const int index) { return line(oldCPv,1, at(index));};
          auto newCpLine = newVLineAt(k)  ;
          auto oldCpLine = oldVLineAt(i) ;

        auto currentNewKnot = newKnotVecU.end();
        auto currentOldKnot = oldKnotVecU.end();
    //
        for (auto const& currentAdditionalKnot : reverse_view(newKnots[0])) {
          while (currentAdditionalKnot <= *(currentOldKnot - 1)
                 && std::distance(oldKnotVecU.begin(), currentOldKnot) > a + 1) {
            std::ranges::transform(oldCpLine,newCpLine.begin(),scaleCPWithW);
            --currentNewKnot;
            --currentOldKnot;
            newCpLine =  newVLineAt(--k);
            oldCpLine =  oldVLineAt(--i);
          }

//          *newCp = *(newCp + 1);
          std::ranges::copy(newVLineAt(k+1),newCpLine.begin());
          newCpLine =  newVLineAt(++k);
          currentOldKnot -= degreeU;
          for ([[maybe_unused]] int _ : iota_view{0, degreeU}) {
            newCpLine =  newVLineAt(++k);

            auto alpha = *currentNewKnot - currentAdditionalKnot;
            using std::abs;
            if (abs(alpha) < 1e-7)
              std::ranges::copy(newCpLine, newVLineAt(k-1).begin());
            else {
              alpha            = alpha / (*currentNewKnot - *(currentOldKnot));
              std::ranges::transform(newCpLine, newVLineAt(k-1),newVLineAt(k-1).begin(),[&alpha](auto& cp, auto& cpL){return  alpha * cpL + cp * (1.0 - alpha);});
//              leftOfcurCPValue = alpha * leftOfcurCPValue + curCPValue * (1.0 - alpha);
            }
            ++currentOldKnot;
            ++currentNewKnot;
          }

          currentNewKnot -= degreeU + 1;
          k-=degreeU+2;
          newCpLine =  newVLineAt(k);
        }
        std::cout << "=========================" << std::endl;
        for (unsigned int i = 0; i < newCPv.size()[0]; i++)
          for (auto&& cp : line(newCPv,1, at(i)))
            std::cout << cp.p << " " << cp.w << std::endl;
        std::cout << std::endl;
        std::cout << "=========================" << std::endl;
        std::ranges::transform(newCPv.directGetAll(), newCPv.directGetAll().begin(), [](auto& cp) -> ControlPoint { return {.p = cp.p / cp.w,.w=cp.w}; });
        std::ranges::for_each(newCPv.directGetAll(),[](auto& cp){std::cout<<"cp.p: "<<cp.p<<" "<<"cp.w: "<<cp.w<<std::endl;});

    //
    //    newKnotsArray[0].insert(newKnotsArray[0].begin(), newKnotsArray[0].front());
    //    newKnotsArray[0].push_back(newKnotsArray[0].back());

    return NURBSPatchData<2, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.getOrder());
  }

  //  template <int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  //  auto surveKnotRefinementu(const NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData,
  //                           const std::array<std::vector<double>, 1>& newKnots) {
  //    using NurbsPatchData = NURBSPatchData<1, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
  //    using namespace std::ranges;
  //    auto oldKnotVec = oldData.getKnots()[0];
  //    oldKnotVec.erase(oldKnotVec.begin());
  //    oldKnotVec.erase(oldKnotVec.end() - 1);
  //
  //    const auto newKSize = newKnots[0].size();
  ////    const unsigned int newKSize = static_cast<unsigned int>(knotsU.size());
  //
  //    const int p        = oldData.getOrder()[0];
  //    const int q        = oldData.getOrder()[1];newCpLine
  //
  //    const unsigned int nb_poles_u = geometry.nb_poles_u();
  //    const unsigned int nb_poles_v = geometry.nb_poles_v();
  //
  //    const unsigned int nb_knots_u = geometry.nb_knots_u();
  //    const unsigned int nb_knots_v = geometry.nb_knots_v();
  //
  //    const int a = findUpperSpanIndex(p, newKnots[0].front(), oldKnotVec);
  //    const int b = findUpperSpanIndex(p, newKnots[0].back(), oldKnotVec);
  //
  ////    const unsigned int numberOfnewControlPoints = oldData.getControlPoints().directSize() + newKSize;
  //
  //
  //    Pointer<SurfaceGeometry> refined = new_<SurfaceGeometry>(p, q, numberOfnewControlPointsInU, nb_poles_v, true);
  //    // FIXME: check is_rational
  //
  //    for (unsigned int i = 0; i < a + 1 - p + 1; i++) {
  //      for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //        refined->pole(i, m) = geometry.pole(i, m) * geometry.weight(i, m);
  //        refined->weight(i, m) = geometry.weight(i, m);
  //      }
  //    }
  //
  //    for (unsigned int i = b + 2 - 1; i < nb_poles_u; i++) {
  //      for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //        refined->pole(newKSize + i, m) = geometry.pole(i, m) * geometry.weight(i, m);
  //        refined->weight(newKSize + i, m) = geometry.weight(i, m);
  //      }
  //    }
  //
  //    for (unsigned int i = 0; i < a + 1; i++) {
  //      refined->knot_u(i) = geometry.knot_u(i);
  //    }
  //
  //    for (unsigned int i = b + 2 + p - 1; i < geometry.nb_knots_u(); i++) {
  //      refined->knot_u(i + newKSize) = geometry.knot_u(i);
  //    }
  //
  //    for (unsigned int i = 0; i < geometry.nb_knots_v(); i++) {
  //      refined->knot_v(i) = geometry.knot_v(i);
  //    }
  //
  //    const unsigned int n = nb_poles_u - 1;
  //    const unsigned int m = n + p + 1;
  //    const unsigned int r = newKSize - 1;
  //
  //    unsigned int i = b + 2 + p - 1;
  //    unsigned int i = b + 2 + p - 1;
  //    unsigned int k = b + 2 + p + r;
  //    unsigned int k = b  newKSize - 1;


  //    unsigned int j = r;
  //
  //    while (j >= 0) {
  //      while (knotsU[j] <= geometry.knot_u(-1 + i) && i > a + 1) {
  //        for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //          const auto pole = geometry.pole(i - p - 1, m);
  //          const auto weight = geometry.weight(i - p - 1, m);
  //          refined->pole(k - p -1, m) = pole * weight;
  //          refined->weight(k - p -1, m) = weight;
  //        }
  //
  //        refined->knot_u(-1+k) = geometry.knot_u(-1 + i);
  //
  //        k -= 1;
  //        i -= 1;
  //      }
  //
  //      for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //        refined->pole(k - p -1, m) = refined->pole(k - p, m);
  //        refined->weight(k - p -1, m) = refined->weight(k - p, m);
  //      }
  //
  //      for (unsigned int l = 1; l < p + 1; l++) {
  //        const unsigned int index = k - p + l;
  //        auto alpha = refined->knot_u(-1+k + l) - knotsU[j];
  //
  //        if (std::abs(alpha) < 1e-7) {
  //          for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //            refined->pole(index - 1, m) = refined->pole(index, m);
  //            refined->weight(index - 1, m) = refined->weight(index, m);
  //          }
  //        } else {
  //          alpha = alpha / (refined->knot_u(k + l - 1) - geometry.knot_u(i + l - p - 1));
  //          for (unsigned int m = 0; m < geometry.nb_poles_v(); m++) {
  //            refined->pole(index - 1, m) = refined->pole(index - 1, m) * alpha + refined->pole(index, m) * (1 -
  //            alpha); refined->weight(index - 1, m) = refined->weight(index - 1, m) * alpha + refined->weight(index,
  //            m) * (1 - alpha);
  //          }
  //        }
  //      }
  //
  //      refined->knot_u(-1 + k) = knotsU[j];
  //
  //      k -= 1;
  //      j -= 1;
  //    }
  //
  //    for (unsigned int i = 0; i < refined->nb_poles(); i++) {
  //      refined->pole(i) = refined->pole(i) / refined->weight(i);
  //    }
  //
  //    return refined;
  //  }
}  // namespace Dune::IGA