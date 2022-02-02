//
// Created by lex on 21.10.21.
//

#pragma once

#include <concepts>
#include <numbers>
#include <ranges>

#include <dune/iga/bsplinealgorithms.hh>
#include <dune/iga/linearalgebra.hh>
#include <dune/iga/multidimensionNet.hh>
#include <dune/iga/nurbspatchdata.hh>
#include <dune/iga/typetraits.hh>
namespace Dune::IGA {

  namespace Impl {
    template <std::integral auto dim>
    auto ordersFromDegrees(const std::array<int, dim>& degree) {
      std::array<int, dim> order;
      std::ranges::transform(degree, order.begin(), [](const auto& p) { return p + 1; });
      return order;
    }

    template <std::integral auto dim>
    auto createPartialSubDerivativPermutations(const FieldVector<int, dim>& v) {
      std::vector<FieldVector<int, dim>> perm;
      for (int i = 1; i < std::pow(2, v.size()); ++i) {
        FieldVector<int, dim> x{};
        for (int j = 0; j < v.size(); ++j)
          if ((i & (1 << j)) != 0) {
            x[j] = v[j];
            if (v[j] == 0) goto outer;
          }
        perm.push_back(x);
      outer:;
      }
      return perm;
    }

    template <Vector VectorType>
    int binom(const VectorType& n, const VectorType& k) {
      return std::inner_product(n.begin(), n.end(), k.begin(), 1, std::multiplies{},
                                [](auto& ni, auto& ki) { return Dune::binomial(ni, ki); });
    }

  }  // namespace Impl

  template <std::floating_point ScalarType, std::integral auto dim>
  auto generateRefinedKnots(const std::array<std::vector<ScalarType>, dim>& knotSpans, const int dir, const int refinementLevel) {
    const int newKnotsSizeForEachSpan = Dune::power(2, refinementLevel);
    auto unique_Knots                 = knotSpans;
    std::vector<double> additionalKnots;
    auto& unique_KnotPerDim = unique_Knots[dir];
    unique_KnotPerDim.erase(std::ranges::begin(std::ranges::unique(unique_KnotPerDim)), std::end(unique_KnotPerDim));

    for (int i = 0; i < unique_KnotPerDim.size() - 1; ++i) {
      const double spanLength = unique_KnotPerDim[i + 1] - unique_KnotPerDim[i];
      const double increment  = spanLength / newKnotsSizeForEachSpan;
      for (int j = 1; j < newKnotsSizeForEachSpan; ++j) {
        additionalKnots.emplace_back(unique_KnotPerDim[i] + increment * j);
      }
    }
    return additionalKnots;
  }




  /** \brief Function return a copy of the subNet of a knotspan, i.e. it is used to return the weights or controlpoints defined on one element
   *
   * @tparam dim dimension of the patch/net
   * @tparam NetValueType The type inside the net. It should be scalar for weights, Controlpoint<> or some vector type which satisfies Vector for the controlpoint coordinates
   * @param subNetStart The start of the controlpointnet in terms of knot span indices. From this the degree is substracted to get the correct Controlpoint index start
   * @param degree array of degrees per direction
   * @param net the full net itself
   * @return subNet with the size degree+1 per direction
   */
  template <std::integral auto dim, typename NetValueType> requires (std::floating_point<NetValueType> || Vector<NetValueType> || is_instantiation_of<ControlPoint,NetValueType>::value )
  auto netOfSpan(std::array<int, dim> subNetStart, const std::array<int, dim>& degree, const MultiDimensionNet<dim, NetValueType>& net) {
    std::array<int, dim> order = Impl::ordersFromDegrees(degree);
    for (std::size_t i = 0; i < dim; ++i)
      subNetStart[i] -= degree[i];
    return net.subNet(subNetStart, order);
  }

  /** \brief Same as netOfSpan above but the start is searched for using the knotvector value   */
  template <std::floating_point ScalarType, std::integral auto dim, typename NetValueType> requires (std::floating_point<NetValueType> || Vector<NetValueType> || is_instantiation_of<ControlPoint,NetValueType>::value )
  auto netOfSpan(const std::array<ScalarType, dim>& u, const std::array<std::vector<ScalarType>, dim>& knots,
                 const std::array<int, dim>& degree, const MultiDimensionNet<dim, NetValueType>& net) {
    auto subNetStart = findSpan(degree, u, knots);
    return netOfSpan(subNetStart, degree, net);
  }

  template <std::integral auto dim, ControlPointConcept ValueType>
  auto extractWeights(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverWeights = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.w; });
    return MultiDimensionNet<dim, typename ValueType::VectorType::value_type>(cpsandWeight.size(), viewOverWeights);
  }

  template <std::integral auto dim, ControlPointConcept ValueType>
  auto extractControlCoordinates(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverCps = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.p; });
    return MultiDimensionNet<dim, typename ValueType::VectorType>(cpsandWeight.size(), viewOverCps);
  }

  template <std::size_t dim, LinearAlgebra NurbsGridLinearAlgebraTraits = DuneLinearAlgebraTraits<double>>
  class Nurbs {
  public:
    Nurbs()                 = default;
    using ScalarType        = typename NurbsGridLinearAlgebraTraits::value_type;
    using DynamicVectorType = typename NurbsGridLinearAlgebraTraits::DynamicVectorType;
    using DynamicMatrixType = typename NurbsGridLinearAlgebraTraits::DynamicMatrixType;

    template <std::integral auto dimworld>
    Nurbs(const Dune::IGA::NURBSPatchData<dim, dimworld>& data, const std::optional<std::array<int, dim>>& spIndex = std::nullopt)
        : knots_{data.knotSpans}, degree_{data.degree}, weights_{extractWeights(data.controlPoints)}, spIndex_{spIndex} {}

    Nurbs(const std::array<std::vector<ScalarType>, dim>& knots, const std::array<int, dim>& degree, const std::vector<ScalarType>& weights)
        : knots_{knots}, degree_{degree}, weights_{weights} {}

    auto operator()(const std::array<ScalarType, dim>& u) { return basisFunctions(u, knots_, degree_, weights_, spIndex_).directGetAll(); }

    auto basisFunctionNet(const std::array<ScalarType, dim>& u) const { return basisFunctions(u, knots_, degree_, weights_, spIndex_); }

    static auto basisFunctions(std::array<ScalarType, dim> u, const std::array<std::vector<ScalarType>, dim>& knots,
                               const std::array<int, dim>& degree, const MultiDimensionNet<dim, ScalarType>& weights,
                               std::optional<std::array<int, dim>> spIndex = std::nullopt) {
      const std::array<int, dim> order = Impl::ordersFromDegrees(degree);
      std::array<DynamicVectorType, dim> bSplines;

      for (std::size_t i = 0; i < dim; ++i)
        bSplines[i] = BsplineBasis1D<ScalarType>::basisFunctions(u[i], knots[i], degree[i],
                                                                 (spIndex ? std::optional<int>(spIndex.value()[i]) : std::nullopt));

      auto Nnet        = MultiDimensionNet<dim, ScalarType>(bSplines);
      auto subNetStart = spIndex ? spIndex.value() : findSpan(degree, u, knots);

      const auto subNetWeights = netOfSpan(subNetStart, degree, weights);

      const ScalarType invSumWeight = dot(Nnet, subNetWeights);
      Nnet *= subNetWeights;
      Nnet /= invSumWeight;
      return Nnet;
    }

    /** \*brief This function return the basis function and the corresponding derivatives
     *
     * @param u dim-dimensional position in current span
     * @param knots dim-dimensional array of knot vectors
     * @param degree   dim-dimensional array of polynomial degree
     * @param weights dim-dimensional net of controlpoint weights
     * @param derivativeOrder  up to which order should derivatives be computed
     * @param triangleDerivatives
     * @param spIndex
     * @return
     */
    static auto basisFunctionDerivatives(std::array<ScalarType, dim> u, const std::array<std::vector<ScalarType>, dim>& knots,
                                         const std::array<int, dim>& degree, const MultiDimensionNet<dim, double>& weights,
                                         const int derivativeOrder, const bool triangleDerivatives = false,
                                         std::optional<std::array<int, dim>> spIndex = std::nullopt) {
      std::array<typename NurbsGridLinearAlgebraTraits::DynamicMatrixType, dim> bSplineDerivatives;
      for (int i = 0; i < dim; ++i)
        bSplineDerivatives[i] = BsplineBasis1D<ScalarType>::basisFunctionDerivatives(
            u[i], knots[i], degree[i], derivativeOrder, (spIndex ? std::optional<int>(spIndex.value()[i]) : std::nullopt));

      std::array<std::vector<ScalarType>, dim> dimArrayOfVectors;
      FieldVector<int, dim> dimSize(derivativeOrder + 1);
      MultiDimensionNet<dim, MultiDimensionNet<dim, ScalarType>> netOfDerivativeNets(dimSize);

      for (int j = 0; auto& derivNet : netOfDerivativeNets.directGetAll()) {
        auto multiIndex = netOfDerivativeNets.directToMultiIndex(j++);
        for (int i = 0; i < multiIndex.size(); ++i)
          dimArrayOfVectors[i] = bSplineDerivatives[i][multiIndex[i]].container();
        derivNet = MultiDimensionNet<dim, ScalarType>(dimArrayOfVectors);
      }

      MultiDimensionNet<dim, MultiDimensionNet<dim, ScalarType>> R = netOfDerivativeNets;
      MultiDimensionNet<dim, ScalarType> netsOfWeightfunctions(dimSize);
      auto subNetStart         = spIndex ? spIndex.value() : findSpan(degree, u, knots);
      const auto subNetWeights = netOfSpan(subNetStart, degree, weights);

      for (int j = 0; j < R.directSize(); ++j) {
        R.directGet(j) *= subNetWeights;
        netsOfWeightfunctions.directGet(j) = dot(netOfDerivativeNets.directGet(j), subNetWeights);
      }

      for (int j = 0; j < R.directSize(); ++j) {
        const auto derivOrders = R.template directToMultiIndex<FieldVector<int, dim>>(j);
        if (triangleDerivatives)
          if (std::accumulate(derivOrders.begin(), derivOrders.end(), derivativeOrder)) continue;
        const auto perms = Impl::createPartialSubDerivativPermutations(derivOrders);

        for (const auto& perm : perms) {
          MultiDimensionNet<dim, int> kNet(perm + FieldVector<int, dim>(1));
          auto startMultiIndex = perm;
          std::ranges::transform(startMultiIndex, startMultiIndex.begin(), [](auto& v) { return (v != 0); });

          for (int kk = kNet.index(startMultiIndex); kk < kNet.directSize(); ++kk) {
            const auto multik = kNet.template directToMultiIndex<FieldVector<int, dim>>(kk);
            R.directGet(j) -= Impl::binom(perm, multik)
                              * (R.get(derivOrders - multik) * netsOfWeightfunctions.get(multik));  // generalized Piegl Tiller (4.20)
          }
        }
        R.directGet(j) /= netsOfWeightfunctions.directGet(0);
      }
      return R;
    }

    auto basisFunctionDerivatives(const std::array<ScalarType, dim>& u, const int derivativeOrder,
                                  const bool triangleDerivatives = false) const {
      return basisFunctionDerivatives(u, knots_, degree_, weights_, derivativeOrder, triangleDerivatives, spIndex_);
    }

  private:
    std::array<std::vector<ScalarType>, dim> knots_;
    std::array<int, dim> degree_;
    MultiDimensionNet<dim, ScalarType> weights_;
    std::optional<std::array<int, dim>> spIndex_;
  };

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto degreeElevate(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData, const int refinementDirection,
                     const int elevationFactor) {
    using DynamicMatrixType = typename NurbsGridLinearAlgebraTraitsImpl::DynamicMatrixType;
    using ScalarType        = typename NurbsGridLinearAlgebraTraitsImpl::value_type;
    const int t             = elevationFactor;
    const int p             = oldData.degree[refinementDirection];
    const int m             = oldData.controlPoints.size()[refinementDirection] + p;
    const int ph            = p + t;
    const int ph2           = ph / 2;
    /* Compute Bezier degree elevation coefficients */
    DynamicMatrixType bezalfs(p + t + 1, p + 1);
    bezalfs[0][0] = bezalfs[ph][p] = 1.0;
    for (int i = 1; i <= ph2; ++i) {
      const ScalarType inv = 1.0 / Dune::binomial(ph, i);
      for (int j = std::max(0, i - t); j <= std::min(p, i); ++j)
        bezalfs[i][j] = inv * Dune::binomial(p, j) * Dune::binomial(t, i - j);
    }
    for (int i = ph2 + 1; i <= ph - 1; ++i)
      for (int j = std::max(0, i - t); j <= std::min(p, i); ++j)
        bezalfs[i][j] = bezalfs[ph - i][p - j];

    double ua = oldData.knotSpans[refinementDirection][0];
    typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType newCPv(oldData.controlPoints.size());
    using ControlPointType = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;

    auto oldCPs = oldData.controlPoints;

    auto scaleCPWithW = [](const auto& cp) -> ControlPointType { return {.p = cp.w * cp.p, .w = cp.w}; };
    std::ranges::transform(oldCPs.directGetAll(), oldCPs.directGetAll().begin(), scaleCPWithW);

    std::vector<ControlPointType> bpts(p + 1), ebpts(ph + 1), nextbpts(p - 1);

    NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl> newData;
    newData.degree = oldData.degree;
    newData.degree[refinementDirection] += t;
    const auto& U = oldData.knotSpans[refinementDirection];
    std::vector<ScalarType> Uh;
    for (int j = 0, i = 0; i < U.size(); ++i)  // insert knot t times for each unique knot
      if (j < U.size()) {
        do {
          Uh.push_back(U[j]);
          ++j;
        } while (j != U.size() && Dune::FloatCmp::eq(U[j - 1], U[j]));
        Uh.insert(Uh.end(), t, U[j - 1]);
        if (j == U.size()) break;
      }

    for (int i = 0; i < dim; ++i)
      if (i == refinementDirection)
        newData.knotSpans[i] = Uh;
      else
        newData.knotSpans[i] = oldData.knotSpans[i];

    auto dimSize                 = oldData.controlPoints.size();
    dimSize[refinementDirection] = Uh.size() - ph - 1;
    newData.controlPoints        = MultiDimensionNet<dim, ControlPointType>(dimSize);
    auto& newCPs                 = newData.controlPoints;

    int totalCurvesInPatch = 1;
    std::array<int, dim - 1> indicesOther;
    for (int counter = 0, i = 0; i < dim; ++i)
      if (i != refinementDirection) {
        totalCurvesInPatch *= dimSize[i];
        indicesOther[counter++] = dimSize[i];
      }

    auto newAndOldCurve = [&refinementDirection, indicesOther, &oldCPs, &newCPs](auto i) {
      const auto mI = MultiDimensionNet<dim - 1, int>::template directToMultiIndex<std::array<int, dim - 1>>(indicesOther, i);
      std::array<int, dim> multiIndex{};
      for (int c = 0, j = 0; j < dim; ++j)
        if (j != refinementDirection) multiIndex[j] = mI[c++];

      auto oldCurve = [&oldCPs, &refinementDirection, multiIndex](int i) mutable {
        multiIndex[refinementDirection] = i;
        return oldCPs.get(multiIndex);
      };

      auto newCurve = [&newCPs, &refinementDirection, multiIndex ](int i) mutable -> auto& {
        multiIndex[refinementDirection] = i;
        return newCPs.get(multiIndex);
      };
      return std::make_pair(oldCurve, newCurve);
    };

    for (auto [oldCurve, newCurve] : std::ranges::iota_view(0, totalCurvesInPatch) | std::views::transform(newAndOldCurve)) {
      int kind = ph + 1;
      int r    = -1;
      int a    = p;
      int b    = p + 1;
      int cind = 1;

      newCurve(0) = oldCurve(0);
      for (int i = 0; i <= p; ++i)
        bpts[i] = oldCurve(i);

      std::vector<ScalarType> alphas(p - 1);
      while (b < m) {
        const int bold = b;
        while (b < m && Dune::FloatCmp::eq(U[b], U[b + 1]))
          ++b;
        const int mul  = b - bold + 1;
        const auto ub  = U[b];
        const int oldr = r;
        r              = p - mul;
        const int lbz  = (oldr > 0) ? ((oldr + 2) / 2) : 1;  // Insert knot u(b) r times
        const int rbz  = (r > 0) ? (ph - (r + 1) / 2) : ph;
        if (r > 0)  // Insert knot to get Bezier segment
        {
          const auto numer = ub - ua;
          for (int k = p; k > mul; --k)
            alphas[k - mul - 1] = numer / (U[a + k] - ua);
          for (int j = 1; j <= r; ++j) {
            const int s = mul + j;
            for (int k = p; k >= s; --k)
              bpts[k] = alphas[k - s] * bpts[k] + (1.0 - alphas[k - s]) * bpts[k - 1];
            nextbpts[r - j] = bpts[p];
          }
        }  // End of insert knot

        for (int i = lbz; i <= ph; ++i) {  // Degree elevate Bezier
          ebpts[i].setZero();
          for (int j = std::max(0, i - t); j <= std::min(p, i); ++j)
            ebpts[i] += bezalfs[i][j] * bpts[j];
        }                // End of degree elevating Bezier
        if (oldr > 1) {  // Must remove knot u = U[a] oldr times
          int first            = kind - 2;
          int last             = kind;
          const ScalarType den = ub - ua;
          const ScalarType bet = (ub - Uh[kind - 1]) / den;
          for (int tr = 1; tr < oldr; ++tr) {  // Knot removal loop
            int i  = first;
            int j  = last;
            int kj = j - kind + 1;
            while (j - i > tr) {
              if (i < cind) {
                const ScalarType alf = (ub - Uh[i]) / (ua - Uh[i]);
                newCurve(i)          = alf * newCurve(i) + (1.0 - alf) * newCurve(i - 1);
              }
              if (j >= lbz) {
                const auto factor = (j - tr <= kind - ph + oldr) ? (ub - Uh[j - tr]) / den : bet;
                ebpts[kj]         = factor * ebpts[kj] + (1.0 - factor) * ebpts[kj + 1];
              }
              ++i;
              --j;
              --kj;
            }
            --first;
            ++last;
          }
        }  // End of removing knot, u= U[a]
        if (a != p) kind += ph - oldr;
        for (int j = lbz; j <= rbz; ++j, ++cind)  // Load ctrl pts into Qw
          newCurve(cind) = ebpts[j];
        if (b < m) {  // Set up for next pass thru loop
          std::ranges::copy(nextbpts | std::views::take(r), bpts.begin());
          for (int j = r; j <= p; ++j)
            bpts[j] = oldCurve(b - p + j);
          a = b;
          ++b;
          ua = ub;
        }
      }  // End of while-loop (b<m)
    }
    std::ranges::transform(newData.controlPoints.directGetAll(), newData.controlPoints.directGetAll().begin(),
                           [](auto& cp) -> ControlPointType { return {.p = cp.p / cp.w, .w = cp.w}; });
    return newData;
  }

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto knotRefinement(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& oldData, const std::vector<double>& newKnots,
                      const int refinementDirection) {
    assert(std::ranges::is_sorted(newKnots)
           && "You passed a non-sorted vector of new knots. Sort it first or check if you passed the correct vector.");
    using NurbsPatchData = NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using ControlPoint   = typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;
    using namespace std::ranges;
    std::array<std::vector<double>, dim> newKnotsArray;
    std::array<int, dim - 1> otherDirections;
    for (int counter = 0, i = 0; i < dim; ++i) {
      if (i == refinementDirection) continue;
      newKnotsArray[i]           = oldData.knotSpans[i];
      otherDirections[counter++] = i;
    }

    typename NurbsPatchData::ControlPointNetType newCPv(oldData.controlPoints.size());

    auto oldCPv     = oldData.controlPoints;
    auto oldKnotVec = oldData.knotSpans[refinementDirection];

    auto& newKnotVec = newKnotsArray[refinementDirection];

    newKnotVec.reserve(oldKnotVec.size() + newKnots.size());
    merge(oldKnotVec, newKnots, std::back_inserter(newKnotVec));

    const auto newKSize                   = newKnots.size();
    std::array<int, dim> numberOfCPperDir = oldCPv.size();
    numberOfCPperDir[refinementDirection] += newKSize;

    newCPv.resize(numberOfCPperDir);

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

    const int p = oldData.degree[refinementDirection];
    const int a = findSpan(p, newKnots.front(), oldKnotVec);
    const int b = findSpan(p, newKnots.back(), oldKnotVec);

    auto newSurfI = newCPv.hyperSurfBegin(otherDirections);
    auto oldSurfI = oldCPv.hyperSurfBegin(otherDirections);

    for (auto hSOld = oldSurfI, hSNew = newSurfI; hSOld != oldCPv.hyperSurfBegin(otherDirections) + (a - p + 1); ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);
    for (auto hSOld = oldSurfI + b, hSNew = newSurfI + b; hSOld != oldCPv.hyperSurfEnd(otherDirections); ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);

    newSurfI = newCPv.hyperSurfEnd(otherDirections) - 1;
    oldSurfI = oldCPv.hyperSurfEnd(otherDirections) - 1;

    auto currentNewKnot = newKnotVec.end();
    auto currentOldKnot = oldKnotVec.end();

    for (auto const& currentAdditionalKnot : reverse_view(newKnots)) {
      while (currentAdditionalKnot <= *(currentOldKnot - 1) && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1) {
        std::ranges::transform(*oldSurfI, (*newSurfI).begin(), scaleCPWithW);
        --currentNewKnot, --currentOldKnot;
        --newSurfI, --oldSurfI;
      }

      std::ranges::copy(*(newSurfI + 1), (*newSurfI).begin());
      ++newSurfI;
      currentOldKnot -= p;
      for ([[maybe_unused]] int _ : iota_view{0, p}) {
        ++newSurfI;

        auto alpha = *currentNewKnot - currentAdditionalKnot;
        using std::abs;
        if (abs(alpha) < 1e-7)
          std::ranges::copy(*newSurfI, (*(newSurfI - 1)).begin());
        else {
          alpha /= (*currentNewKnot - *(currentOldKnot));
          std::ranges::transform(*newSurfI, *(newSurfI - 1), (*(newSurfI - 1)).begin(),
                                 [&alpha](auto& cp, auto& cpL) { return alpha * cpL + cp * (1.0 - alpha); });
        }
        ++currentOldKnot, ++currentNewKnot;
      }
      currentNewKnot -= p + 1;
      newSurfI -= p + 2;
    }

    std::ranges::transform(newCPv.directGetAll(), newCPv.directGetAll().begin(),
                           [](auto& cp) -> ControlPoint { return {.p = cp.p / cp.w, .w = cp.w}; });

    return NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, newCPv, oldData.degree);
  }

  namespace Impl {
    template <Vector VectorType>
    auto projectPointOntoLine(const VectorType basepoint, const VectorType revolutionaxis, const VectorType point) {
      const VectorType e1 = basepoint + revolutionaxis - basepoint;
      const VectorType e2 = point - basepoint;
      using Dune::dot;
      const typename VectorType::value_type angle = dot(e1, e2);
      const typename VectorType::value_type len2  = two_norm(e1);

      VectorType p = basepoint + angle * e1 / len2;
      return p;
    }

    template <Vector VectorType>
    auto Intersect3DLines(const VectorType basepoint1, const VectorType direction1, const VectorType basepoint2,
                          const VectorType direction2) {
      using ScalarType                     = typename VectorType::value_type;
      const ScalarType tol                 = 1e-8;
      const VectorType basePointDifference = basepoint1 - basepoint2;
      using Dune::dot;
      const ScalarType dir1squaredLength = dot(direction1, direction1);
      const ScalarType angle             = dot(direction1, direction2);
      const ScalarType dir2squaredLength = dot(direction2, direction2);
      const ScalarType dir1BasePdiff     = dot(direction1, basePointDifference);
      const ScalarType dir2BasePdiff     = dot(direction2, basePointDifference);
      const ScalarType D                 = dir1squaredLength * dir2squaredLength - angle * angle;

      if (std::abs(D) < tol) throw std::logic_error("The two lines are almost parallel.");

      const ScalarType parameter1 = (angle * dir2BasePdiff - dir2squaredLength * dir1BasePdiff) / D;
      const ScalarType parameter2 = (dir1squaredLength * dir2BasePdiff - angle * dir1BasePdiff) / D;

      VectorType c1 = basepoint1 + parameter1 * direction1;
      assert(dot(c1 - (basepoint2 + parameter2 * direction2), c1 - (basepoint2 + parameter2 * direction2)) < tol
             && "Both calculated points do not coincide");
      return c1;
    }
  }  // namespace Impl

  // Piegl and Tiller Algo A7.1
  template <LinearAlgebra NurbsGridLinearAlgebraTraitsImpl = Dune::IGA::DuneLinearAlgebraTraits<double>>
  auto makeCircularArc(const typename NurbsGridLinearAlgebraTraitsImpl::value_type radius                  = 1.0,
                       const typename NurbsGridLinearAlgebraTraitsImpl::value_type startAngle              = 0.0,
                       typename NurbsGridLinearAlgebraTraitsImpl::value_type endAngle                      = 360.0,
                       const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3> origin = {0, 0, 0},
                       const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3> X      = {1, 0, 0},
                       const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3> Y      = {0, 1, 0}) {
    using ScalarType           = typename NurbsGridLinearAlgebraTraitsImpl::value_type;
    using GlobalCoordinateType = typename NURBSPatchData<1, 3, NurbsGridLinearAlgebraTraitsImpl>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<typename NurbsGridLinearAlgebraTraitsImpl::value_type>;

    if (endAngle < startAngle) endAngle += 360.0;
    const ScalarType theta = endAngle - startAngle;
    const int narcs        = std::ceil(theta / 90);

    typename NURBSPatchData<1, 3, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType circleCPs(2 * narcs + 1);
    const ScalarType dtheta  = theta / narcs * pi / 180;
    const int n              = 2 * narcs;
    const ScalarType w1      = cos(dtheta / 2.0);
    GlobalCoordinateType PO  = origin + radius * cos(startAngle) * X + radius * sin(startAngle) * Y;
    GlobalCoordinateType TO  = -sin(startAngle) * X + cos(startAngle) * Y;
    circleCPs.directGet(0).p = PO;
    ScalarType angle         = startAngle;
    for (int index = 0, i = 0; i < narcs; ++i) {
      angle += dtheta;
      const GlobalCoordinateType P2    = origin + radius * cos(angle) * X + radius * sin(angle) * Y;
      circleCPs.directGet(index + 2).p = P2;
      const GlobalCoordinateType T2    = -sin(angle) * X + cos(angle) * Y;
      const GlobalCoordinateType P1    = Impl::Intersect3DLines(PO, TO, P2, T2);
      circleCPs.directGet(index + 1)   = {.p = P1, .w = w1};
      index += 2;
      if (i < narcs - 1) {
        PO = P2;
        TO = T2;
      }
    }

    auto knotVec = std::array<std::vector<double>, 1>{};
    auto& U      = knotVec[0];
    U.resize(2 * (narcs + 2));
    if (narcs == 1) {
    } else if (narcs == 2) {
      U[3] = U[4] = 1.0 / 2.0;
    } else if (narcs == 3) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 1.0 / 4.0;
      U[5] = U[6] = 2.0 / 4.0;
      U[7] = U[8] = 3.0 / 4.0;
    }

    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);
    return NURBSPatchData<1, 3, NurbsGridLinearAlgebraTraitsImpl>(knotVec, circleCPs, {2});
  }

  // Algo 8.1
  template <LinearAlgebra NurbsGridLinearAlgebraTraitsImpl = Dune::IGA::DuneLinearAlgebraTraits<double>>
  auto makeSurfaceOfRevolution(const NURBSPatchData<1, 3, NurbsGridLinearAlgebraTraitsImpl>& generatrix,
                               const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3> point,
                               const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3> revolutionaxisI,
                               const typename NurbsGridLinearAlgebraTraitsImpl::template FixedVectorType<3>::value_type revolutionAngle
                               = 360.0) {
    const auto& genCP          = generatrix.controlPoints;
    using ScalarType           = typename NurbsGridLinearAlgebraTraitsImpl::value_type;
    using ControlPoint         = typename NURBSPatchData<2, 3, NurbsGridLinearAlgebraTraitsImpl>::ControlPointType;
    using GlobalCoordinateType = typename NURBSPatchData<2, 3, NurbsGridLinearAlgebraTraitsImpl>::GlobalCoordinateType;
    const auto pi              = std::numbers::pi_v<typename NurbsGridLinearAlgebraTraitsImpl::value_type>;

    const auto revolutionaxis = revolutionaxisI / revolutionaxisI.two_norm();
    auto newKnotsArray        = std::array<std::vector<double>, 2UL>();
    newKnotsArray[1]          = generatrix.knotSpans[0];
    auto& U                   = newKnotsArray[0];

    const int narcs = std::ceil(revolutionAngle / 90);
    U.resize(2 * (narcs + 2));
    if (revolutionAngle <= 90.0) {
    } else if (revolutionAngle <= 180.0) {
      U[3] = U[4] = 0.5;
    } else if (revolutionAngle <= 270.0) {
      U[3] = U[4] = 1.0 / 3.0;
      U[5] = U[6] = 2.0 / 3.0;
    } else {
      U[3] = U[4] = 0.25;
      U[5] = U[6] = 0.5;
      U[7] = U[8] = 0.75;
    }
    const ScalarType dtheta = revolutionAngle / narcs * pi / 180;
    std::ranges::fill_n(U.begin(), 3, 0.0);
    std::ranges::fill_n(std::ranges::reverse_view(U).begin(), 3, 1.0);

    typename NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType surfaceCP(2 * narcs + 1,
                                                                                                       generatrix.controlPoints.size()[0]);
    using std::cos;
    using std::sin;
    const ScalarType wm = cos(dtheta / 2.0);
    std::vector<ScalarType> cosines(narcs);
    std::vector<ScalarType> sines(narcs);
    ScalarType angle = 0.0;
    for (int i = 0; i < narcs; i++) {
      angle += dtheta;
      cosines[i] = cos(angle);
      sines[i]   = sin(angle);
    }
    ControlPoint PO = genCP.directGet(0);
    for (int j = 0; j < genCP.size()[0]; j++) {
      const GlobalCoordinateType Om = Impl::projectPointOntoLine(point, revolutionaxis, genCP.directGet(j).p);
      GlobalCoordinateType X        = genCP.directGet(j).p - Om;
      const ScalarType r            = X.two_norm();
      X /= r;
      const GlobalCoordinateType Y = cross(revolutionaxis, X);
      surfaceCP(0, j) = PO    = genCP.directGet(j);
      GlobalCoordinateType TO = Y;
      for (int index = 0, i = 0; i < narcs; ++i) {
        const GlobalCoordinateType P2 = Om + r * cosines[i] * X + r * sines[i] * Y;
        surfaceCP(index + 2, j)       = {.p = P2, .w = genCP.directGet(j).w};

        const GlobalCoordinateType T2 = -sines[i] * X + cosines[i] * Y;
        surfaceCP(index + 1, j).p     = Impl::Intersect3DLines(PO.p, TO, P2, T2);
        surfaceCP(index + 1, j).w     = wm * genCP.directGet(j).w;
        index += 2;
        if (i < narcs - 1) {
          PO.p = P2;
          TO   = T2;
        }
      }
    }
    return NURBSPatchData<2UL, 3UL, NurbsGridLinearAlgebraTraitsImpl>(newKnotsArray, surfaceCP, {2, generatrix.degree[0]});
  }
}  // namespace Dune::IGA