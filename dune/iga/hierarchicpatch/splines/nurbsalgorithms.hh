// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <concepts>
#include <numbers>
#include <ranges>

#include <dune/iga/hierarchicpatch/concepts.hh>
#include <dune/iga/hierarchicpatch/splines/bsplinealgorithms.hh>
#include <dune/iga/hierarchicpatch/splines/nurbspatchdata.hh>
#include <dune/iga/hierarchicpatch/utils/linearalgebra.hh>
#include <dune/iga/hierarchicpatch/utils/mdnet.hh>
#include <dune/iga/hierarchicpatch/utils/typetraits.hh>
namespace Dune::IGANEW {

  namespace Impl {
    template <std::integral auto dim>
    auto ordersFromDegrees(const std::array<int, dim>& degree) {
      std::array<int, dim> order;
      std::ranges::transform(degree, order.begin(), [](const auto& p) { return p + 1; });
      return order;
    }

    template <std::integral auto dim>
    void createPartialSubDerivativPermutations(const FieldVector<int, dim>& v,
                                               std::vector<FieldVector<int, dim>>& perm) {
      perm.resize(0);
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
    }

    template <Concept::Vector VectorType>
    int binom(const VectorType& n, const VectorType& k) {
      return std::inner_product(n.begin(), n.end(), k.begin(), 1, std::multiplies{},
                                [](auto& ni, auto& ki) { return Dune::binomial(ni, ki); });
    }

  }  // namespace Impl

  template <std::floating_point ScalarType, std::integral auto dim>
  auto generateRefinedKnots(const std::array<std::vector<ScalarType>, dim>& knotSpans, const int dir,
                            const int refinementLevel) {
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

  /** \brief Function return a copy of the subNet of a knotspan, i.e. it is used to return the weights or controlpoints
   * defined on one element
   *
   * @tparam dim dimension of the patch/net
   * @tparam NetValueType The type inside the net. It should be scalar for weights, Controlpoint<> or some vector type
   * which satisfies Vector for the controlpoint coordinates
   * @param subNetStart The start of the control point net in terms of knot span indices. From this the degree is
   * subtracted to get the correct control point index start
   * @param degree array of degrees per direction
   * @param net the full net itself
   * @return subNet with the strideSizes degree+1 per direction
   */
  template <std::integral auto dim, typename NetValueType>
  requires(std::floating_point<NetValueType> || Concept::Vector<NetValueType> || is_instantiation_of<ControlPoint, NetValueType>::value) auto netOfSpan(
      std::array<int, dim> subNetStart, const std::array<int, dim>& degree,
      const MultiDimensionNet<dim, NetValueType>& net) {
    std::array<int, dim> order = Impl::ordersFromDegrees(degree);
    for (std::size_t i = 0; i < dim; ++i)
      subNetStart[i] -= degree[i];
    return net.subNet(subNetStart, order);
  }

  /** \brief Same as netOfSpan above but the start is searched for using the knotvector value   */
  template <std::floating_point ScalarType, std::integral auto dim, std::integral auto dim2, typename NetValueType>
  requires(std::floating_point<NetValueType> || Concept::Vector<NetValueType> || is_instantiation_of<ControlPoint, NetValueType>::value) auto netOfSpan(
      const Dune::FieldVector<ScalarType, dim>& u, const std::array<std::vector<ScalarType>, dim2>& knots,
      const std::array<int, dim2>& degree, const MultiDimensionNet<dim2, NetValueType>& net) requires(dim == dim2) {
    auto subNetStart = findSpan(degree, u, knots);
    return netOfSpan(subNetStart, degree, net);
  }

  template <std::integral auto dim, Concept::ControlPoint ValueType>
  auto extractWeights(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverWeights = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.w; });
    return MultiDimensionNet<dim, typename ValueType::VectorType::value_type>(cpsandWeight.strideSizes(),
                                                                              viewOverWeights);
  }

  template <std::integral auto dim, Concept::ControlPoint ValueType>
  auto extractControlCoordinates(const MultiDimensionNet<dim, ValueType>& cpsandWeight) {
    auto viewOverCps = std::ranges::transform_view(cpsandWeight.directGetAll(), [](auto& cp) { return cp.p; });
    return MultiDimensionNet<dim, typename ValueType::VectorType>(cpsandWeight.strideSizes(), viewOverCps);
  }

  /** \brief A `dim` dimensional NURBS function
   *
   * /note This class explicitly does not relay on some controlpoints, but only on the weights
   * @tparam dim The dimension of the domain of the function
   * @tparam ScalarType_ The type for the functions values and arguments
   *
   */
  template <int dim, typename ScalarType_ = double>
  class Nurbs {
   public:
    Nurbs()                        = default;
    using ScalarType               = ScalarType_;
    using DynamicVectorType        = Dune::DynamicVector<ScalarType>;
    using DynamicMatrixType        = Dune::DynamicMatrix<ScalarType>;
    static constexpr int dimension = dim;

    /**
     *
     * @tparam dimworld
     * @param data
     * @param spIndex
     */
    template <std::integral auto dimworld>
    explicit Nurbs(const NURBSPatchData<dim, dimworld, ScalarType>& data)
        : knots_{data.knotSpans}, degree_{data.degree}, weights_{extractWeights(data.controlPoints)} {}

    auto localView() const { return LocalView(*this); }

    struct LocalView {
      static constexpr int dimension = dim;
      LocalView()                    = default;
      explicit LocalView(const Nurbs& nurbs) : nurbs_{&nurbs} {}

      void bind(const std::array<int, dim>& spIndex) { spIndex_ = spIndex; }

      auto basisFunctionDerivatives(const Dune::FieldVector<ScalarType, dim>& u, const int derivativeOrder) const {
        assert(spIndex_ && "Bind the local view first!");
        return nurbs_->basisFunctionDerivatives(u, derivativeOrder, spIndex_.value());
      }

      auto basisFunctions(const Dune::FieldVector<ScalarType, dim>& u) const {
        assert(spIndex_ && "Bind the local view first!");
        return nurbs_->basisFunctions(u, spIndex_.value());
      }

      std::optional<std::array<int, dim>> spIndex_;
      Nurbs const* nurbs_{nullptr};
    };

    Nurbs(const std::array<std::vector<ScalarType>, dim>& knots, const std::array<int, dim>& degree,
          const std::vector<ScalarType>& weights)
        : knots_{knots}, degree_{degree}, weights_{weights} {}

    // auto operator()(const Dune::FieldVector<ScalarType, dim>& u) {
    //   return basisFunctions(u, knots_, degree_, weights_).directGetAll();
    // }

    auto basisFunctions(const Dune::FieldVector<ScalarType, dim>& u) const {
      auto subNetStart = findSpanCorrecte(degree_, u, knots_);
      return basisFunctions(u, knots_, degree_, weights_, subNetStart);
    }

    auto basisFunctions(const Dune::FieldVector<ScalarType, dim>& u, const std::array<int, dim>& spIndex) const {
      return basisFunctions(u, knots_, degree_, weights_, spIndex);
    }

    static auto basisFunctions(const Dune::FieldVector<ScalarType, dim>& u,
                               const std::array<std::vector<ScalarType>, dim>& knots,
                               const std::array<int, dim>& degree, const MultiDimensionNet<dim, ScalarType>& weights,
                               std::array<int, dim> spIndex) {
      std::array<DynamicVectorType, (size_t)dim> bSplines;

      for (std::size_t i = 0; i < dim; ++i)
        bSplines[i] = BsplineBasis<ScalarType>::basisFunctions(u[i], knots[i], degree[i], spIndex[i]);

      auto Nnet = MultiDimensionNet<dim, ScalarType>(bSplines);
      // TODO Alex this should also be cached when called with a localView
      const auto subNetWeights = netOfSpan(spIndex, degree, weights);

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
    static auto basisFunctionDerivatives(const Dune::FieldVector<ScalarType, dim>& u,
                                         const std::array<std::vector<ScalarType>, dim>& knots,
                                         const std::array<int, dim>& degree,
                                         const MultiDimensionNet<dim, double>& weights, const int derivativeOrder,
                                         const bool triangleDerivatives              = false,
                                         std::optional<std::array<int, dim>> spIndex = std::nullopt) {
      std::array<DynamicMatrixType, dim> bSplineDerivatives;
      for (int i = 0; i < dim; ++i)
        bSplineDerivatives[i] = BsplineBasis<ScalarType>::basisFunctionDerivatives(
            u[i], knots[i], degree[i], derivativeOrder,
            (spIndex ? std::optional<int>(spIndex.value()[i]) : std::nullopt));

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

      for (int j = 0; j < R.size(); ++j) {
        R.directGet(j) *= subNetWeights;
        netsOfWeightfunctions.directGet(j) = dot(netOfDerivativeNets.directGet(j), subNetWeights);
      }

      std::vector<FieldVector<int, dim>> perms;
      for (int j = 0; j < R.size(); ++j) {
        const auto derivOrders = R.template directToMultiIndex<FieldVector<int, dim>>(j);
        if (triangleDerivatives)
          if (std::accumulate(derivOrders.begin(), derivOrders.end(), derivativeOrder)) continue;
        Impl::createPartialSubDerivativPermutations(derivOrders, perms);

        for (const auto& perm : perms) {
          const MultiDimensionNetIndex<dim> kNet(perm + FieldVector<int, dim>(1));
          auto startMultiIndex = perm;
          std::ranges::transform(startMultiIndex, startMultiIndex.begin(), [](const auto& v) { return (v != 0); });
          for (int kk = kNet.index(startMultiIndex); kk < kNet.size(); ++kk) {
            const auto multik     = kNet.template directToMultiIndex<FieldVector<int, dim>>(kk);
            const ScalarType fac  = (Impl::binom(perm, multik) * netsOfWeightfunctions.get(multik));
            const auto& Rdirect_d = R.get(derivOrders - multik);
            auto& Rdirect         = R.directGet(j);
            for (int i = 0; i < R.directGet(j).size(); ++i)
              Rdirect.directGet(i) -= fac * Rdirect_d.directGet(i);  // generalized Piegl Tiller (4.20)
          }
        }
        R.directGet(j) /= netsOfWeightfunctions.directGet(0);
      }
      return R;
    }

    auto basisFunctionDerivatives(const Dune::FieldVector<ScalarType, dim>& u, const int derivativeOrder) const {
      return basisFunctionDerivatives(u, knots_, degree_, weights_, derivativeOrder);
    }

    auto basisFunctionDerivatives(const Dune::FieldVector<ScalarType, dim>& u, const int derivativeOrder,
                                  const std::array<int, dim>& spIndex) const {
      return basisFunctionDerivatives(u, knots_, degree_, weights_, derivativeOrder, false, spIndex);
    }

   private:
    std::array<std::vector<ScalarType>, dim> knots_;
    std::array<int, dim> degree_;
    MultiDimensionNet<dim, ScalarType> weights_;
  };

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType_>
  auto degreeElevate(const NURBSPatchData<dim, dimworld, ScalarType_>& oldData, const int refinementDirection,
                     const int elevationFactor) {
    using ScalarType          = ScalarType_;
    using DynamicMatrixType   = Dune::DynamicMatrix<ScalarType>;
    using NURBSPatchData      = NURBSPatchData<dim, dimworld, ScalarType_>;
    using ControlPointNetType = typename NURBSPatchData::ControlPointNetType;
    using ControlPointType    = typename NURBSPatchData::ControlPointType;
    const int t               = elevationFactor;
    const int p               = oldData.degree[refinementDirection];
    const int m               = oldData.controlPoints.strideSizes()[refinementDirection] + p;
    const int ph              = p + t;
    const int ph2             = ph / 2;
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

    ControlPointNetType newCPv(oldData.controlPoints.strideSizes());

    auto oldCPs = oldData.controlPoints;

    auto scaleCPWithW = [](const auto& cp) -> ControlPointType { return {.p = cp.w * cp.p, .w = cp.w}; };
    std::ranges::transform(oldCPs.directGetAll(), oldCPs.directGetAll().begin(), scaleCPWithW);

    std::vector<ControlPointType> bpts(p + 1), ebpts(ph + 1), nextbpts(p - 1);

    NURBSPatchData newData;
    newData.degree = oldData.degree;
    newData.degree[refinementDirection] += t;
    const auto& U = oldData.knotSpans[refinementDirection];
    double ua     = U[0];
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

    auto dimSize                 = oldData.controlPoints.strideSizes();
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
      const auto mI
          = MultiDimensionNet<dim - 1, int>::template directToMultiIndex<std::array<int, dim - 1>>(indicesOther, i);
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

    for (auto [oldCurve, newCurve] :
         std::ranges::iota_view(0, totalCurvesInPatch) | std::views::transform(newAndOldCurve)) {
      ua       = U[0];
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
        if (b < m) {  // Set up for next pass through loop
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

  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType>
  auto knotRefinement(const NURBSPatchData<dim, dimworld, ScalarType>& oldData, const std::vector<double>& newKnots,
                      const int refinementDirection) {
    assert(std::ranges::is_sorted(newKnots)
           && "You passed a non-sorted vector of new knots. Sort it first or check if you passed the correct vector.");
    using NurbsPatchData = NURBSPatchData<dim, dimworld, ScalarType>;
    using ControlPoint   = typename NurbsPatchData::ControlPointType;
    using namespace std::ranges;
    std::array<std::vector<double>, dim> newKnotsArray;
    std::array<int, dim - 1> otherDirections;
    for (int counter = 0, i = 0; i < dim; ++i) {
      if (i == refinementDirection) continue;
      newKnotsArray[i]           = oldData.knotSpans[i];
      otherDirections[counter++] = i;
    }

    typename NurbsPatchData::ControlPointNetType newCPv(oldData.controlPoints.strideSizes());

    auto oldCPv     = oldData.controlPoints;
    auto oldKnotVec = oldData.knotSpans[refinementDirection];

    auto& newKnotVec = newKnotsArray[refinementDirection];

    newKnotVec.reserve(oldKnotVec.size() + newKnots.size());
    merge(oldKnotVec, newKnots, std::back_inserter(newKnotVec));

    const auto newKSize                   = newKnots.size();
    std::array<int, dim> numberOfCPperDir = oldCPv.strideSizes();
    numberOfCPperDir[refinementDirection] += newKSize;

    newCPv.resize(numberOfCPperDir);

    auto scaleCPWithW = [](const auto& cp) -> ControlPoint { return {.p = cp.w * cp.p, .w = cp.w}; };

    const int p = oldData.degree[refinementDirection];
    const int a = findSpan(p, newKnots.front(), oldKnotVec);
    const int b = findSpan(p, newKnots.back(), oldKnotVec);

    auto newSurfI = newCPv.hyperSurfBegin(otherDirections);
    auto oldSurfI = oldCPv.hyperSurfBegin(otherDirections);

    for (auto hSOld = oldSurfI, hSNew = newSurfI; hSOld != oldCPv.hyperSurfBegin(otherDirections) + (a - p + 1);
         ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);
    for (auto hSOld = oldSurfI + b, hSNew = newSurfI + b; hSOld != oldCPv.hyperSurfEnd(otherDirections);
         ++hSOld, ++hSNew)
      std::ranges::transform(*hSOld, (*hSNew).begin(), scaleCPWithW);

    newSurfI = newCPv.hyperSurfEnd(otherDirections) - 1;
    oldSurfI = oldCPv.hyperSurfEnd(otherDirections) - 1;

    auto currentNewKnot = newKnotVec.end();
    auto currentOldKnot = oldKnotVec.end();

    for (auto const& currentAdditionalKnot : reverse_view(newKnots)) {
      while (currentAdditionalKnot <= *(currentOldKnot - 1)
             && std::distance(oldKnotVec.begin(), currentOldKnot) > a + 1) {
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

    return NURBSPatchData<dim, dimworld, ScalarType>(newKnotsArray, newCPv, oldData.degree);
  }






}  // namespace Dune::IGANEW
