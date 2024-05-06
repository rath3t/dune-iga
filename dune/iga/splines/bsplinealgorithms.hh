// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>

#include <dune/common/dynmatrix.hh>
#include <dune/common/float_cmp.hh>

namespace Dune::IGANEW::Splines {

/** @brief Finds the spanIndex in range [u_0,...,u_0,...,u_a      ,u,...,u_n,...,u_n] which is first index lower than
 * u
 * @tparam Range knotvector range
 * @param p polynomial degree of the underlying spline
 * @param u span value which is searched for
 * @param U range where to search
 * @param offset adjust range where to start searching to improve efficiency
 * @return
 */
template <std::ranges::random_access_range Range>
auto findSpan(const int p, typename std::remove_cvref_t<Range>::value_type u, Range&& U,
              typename std::remove_cvref_t<Range>::difference_type offset = 0) {
  if (u <= U[0])
    return static_cast<long int>(p);
  if (u >= U.back())
    return static_cast<long int>(U.size() - p - 2); // if the coordinate is to big we return to the last non-end span
  auto it = std::upper_bound(U.begin() + p - 1 + offset, U.end(), u);
  return static_cast<long int>(std::distance(U.begin(), it) - 1);
}

/** @brief Same as findSpan but for dim  knotvectors @see findSpan */
template <int dim, size_t dim2, typename ValueType>
auto findSpan(const std::array<int, dim2>& p, const Dune::FieldVector<ValueType, dim>& u,
              const std::array<std::vector<ValueType>, dim2>& U)
requires(dim2 == dim)
{
  std::array<int, dim> res;
  for (auto i = 0; i < dim; ++i)
    res[i] = findSpan(p[i], u[i], U[i]);
  return res;
}

/**
 * @brief Finds the multiplicity in a B-spline knot vector within the specified range.
 *
 * This function determines the multiplicity of a given span value 'u' within the knot vector 'U' in the specified
 * range. The knot vector 'U' is assumed to be non-decreasing. The function returns a pair containing the multiplicity
 * of the knot at position 'u' and the corresponding span index, which is the findSpan result.
 *
 * @tparam Range Random access range type representing the knot vector.
 * @param U The knot vector range.
 * @param u The span value to search for.
 * @return A pair containing the multiplicity of the knot at position 'u' and the corresponding span index.
 *
 * @note The knot vector 'U' should be non-decreasing.
 *
 * @example
 * @code
 *
 * std::vector<double> v({0,0,0,1,1,1,2,2,3,3,3});
 * std::vector<double> vUnique;
 *
 * std::ranges::unique_copy(v, std::back_inserter(vUnique));
 *
 * for (double u: vUnique)
 * {
 *   auto [m, span] = multiplicity(v, u);
 *   std::cout << u << " is at pos (" << span << ") with multiplicity: " << m << "\n";
 * }
 * @endcode
 *
 * The above code results in the following console output:
 * @code
 * 0 is at pos (2) with multiplicity: 3
 * 1 is at pos (5) with multiplicity: 3
 * 2 is at pos (7) with multiplicity: 2
 * 3 is at pos (7) with multiplicity: 3
 * @endcode
 */
template <std::ranges::random_access_range Range>
auto multiplicity(Range&& U, typename std::remove_cvref_t<Range>::value_type u) {
  auto [lit, rit]        = std::ranges::equal_range(U, u);
  auto rightPos          = std::distance(U.begin(), rit);
  const int multiplicity = static_cast<long int>(std::distance(lit, rit));
  if (u >= U.back())
    return std::make_pair(multiplicity,
                          static_cast<long int>(rightPos - multiplicity -
                                                1)); // if the coordinate is to big we return to the last non-end span
  if (u <= U[0])
    return std::make_pair(multiplicity, static_cast<long int>(multiplicity - 1));
  return std::make_pair(multiplicity, static_cast<long int>(std::distance(U.begin(), rit) - 1));
}

/** @brief One dimensional b-spline basis
 *
 * @tparam ScalarType_ Scalar type of the coordinates
 */
template <typename ScalarType_>
class BsplineBasis
{
public:
  using ScalarType        = ScalarType_;
  using DynamicVectorType = Dune::DynamicVector<ScalarType>;
  using DynamicMatrixType = Dune::DynamicMatrix<ScalarType>;
  using RowFixedMatrix    = std::array<DynamicVectorType, 2>;

  BsplineBasis(const std::vector<ScalarType>& knots, const int degree)
      : knots_{knots},
        degree_{degree} {}

  auto operator()(ScalarType u) {
    return basisFunctions(u, knots_, degree_);
  }

  /** @brief The evaluation function modified version of The Nurbs Book Algorithm A2.2
   *
   * @tparam Range range of knotvector
   * @param u scalar point where to evaluate the spline
   * @param knots knotvector
   * @param p polynomial degree of the underlying b-spline basis
   * @param spIndex optional index in which range the evaluation point lies, if omitted it is searched for
   * @return Non-zero B-spline basis functions evaluated at u
   */
  template <std::ranges::random_access_range Range>
  static auto basisFunctions(typename std::remove_cvref_t<Range>::value_type u, Range&& knots, const int p,
                             std::optional<int> spIndex = std::nullopt) {
    assert(std::ranges::count(knots.begin(), knots.begin() + p + 1, knots.front()) == p + 1);
    assert(std::ranges::count(knots.end() - p - 1, knots.end(), knots.back()) == p + 1);
    assert(spIndex < knots.size() - p - 1);
    DynamicVectorType N;
    N.resize(p + 1, 0.0);
    u = std::clamp(u, knots.front(), knots.back()); // we clamp the value since sometimes we overshoot slightly
    // due to machine precision
    if (FloatCmp::eq(u, knots.back())) // early exit
    {
      N.back() = 1;
      return N;
    }
    if (FloatCmp::eq(u, knots.front())) {
      N.front() = 1;
      return N;
    }

    // if the index in the knot span is not given we search for it
    const int sp = spIndex ? spIndex.value() : findSpan(p, u, knots);
    using namespace std::ranges;
    auto lDiff = transform_view(reverse_view(std::views::counted(knots.begin() + sp + 1 - p, p)),
                                [&u](auto& kn) { return u - kn; });
    auto rDiff = transform_view(std::views::counted(knots.begin() + sp + 1, p), [&u](auto& kn) { return kn - u; });

    N[0] = ScalarType(1);

    int j = 0;
    for (auto lDp = lDiff.begin(); lDp != lDiff.end(); lDp += j + 2, ++j) {
      ScalarType saved{};
      for (int r = 0; auto&& rD : rDiff | std::views::take(j + 1)) {
        const auto& lD       = (*lDp);
        const ScalarType tmp = N[r] / (rD + lD);
        N[r++]               = saved + rD * tmp;
        saved                = lD * tmp;
        --lDp;
      }
      N[j + 1] = saved;
    }
    // The basis functions are always >=0!
    assert(std::ranges::all_of(N, [](const auto& Ni) { return not Dune::FloatCmp::lt(Ni, -1e-8); }));

    return N;
  }

  /** @brief The evaluation function modified version of The Nurbs Book Algorithm A2.3
   *
   * @param u scalar point where to evaluate the spline
   * @param knots knotvector
   * @param p polynomial p of the underlying b-spline basis
   * @param derivativeOrder maximum derivative degree that are constructed
   * @param spIndex optional index in which range the evaluation point lies, if omitted it is searched for
   * @return Non-zero B-spline derivatives  evaluated at u
   */
  static auto basisFunctionDerivatives(ScalarType u, const std::vector<ScalarType>& knots, const int p,
                                       const int derivativeOrder, std::optional<int> spIndex = std::nullopt) {
    assert(spIndex < knots.size() - p - 1);
    const int order = p + 1;
    const int sp    = spIndex ? spIndex.value() : findSpan(p, u, knots);
    using namespace std::ranges;

    DynamicMatrixType dN(derivativeOrder + 1, order);

    std::vector<ScalarType> left(order);
    std::vector<ScalarType> right(order);
    DynamicMatrixType ndu(order, order);

    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; j++) {
      left[j]      = u - knots[sp + 1 - j];
      right[j]     = knots[sp + j] - u;
      double saved = 0.0;
      for (int r = 0; r < j; r++) { /* Lower triangle */
        ndu[j][r]   = right[r + 1] + left[j - r];
        double temp = ndu[r][j - 1] / ndu[j][r];
        ndu[r][j]   = saved + right[r + 1] * temp;
        saved       = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j)
      dN[0][j] = ndu[j][p];

    // Compute the derivatives
    RowFixedMatrix a;
    a[0].resize(p + 1);
    a[1].resize(p + 1);
    auto& a1Row = a[0];
    auto& a2Row = a[1];
    for (int r = 0; r <= p; ++r) {
      a[0][0] = 1.0;

      // Compute the k-th derivative
      for (int k = 1; k <= derivativeOrder; ++k) {
        const int rk = r - k;
        const int pk = p - k;

        auto& nduRowpk1 = ndu[pk + 1];
        auto& dNcur     = dN[k][r];
        dNcur           = 0.0;
        if (r >= k) {
          a2Row[0] = a1Row[0] / nduRowpk1[rk];
          dNcur    = a2Row[0] * ndu[rk][pk];
        }

        const int j1 = (rk >= -1) ? 1 : -rk;
        const int j2 = (r - 1 <= pk) ? k - 1 : p - r;

        for (int j = j1; j <= j2; ++j) {
          a2Row[j] = (a1Row[j] - a1Row[j - 1]) / nduRowpk1[rk + j];
          dNcur += a2Row[j] * ndu[rk + j][pk];
        }

        if (r <= pk) {
          a2Row[k] = -a1Row[k - 1] / nduRowpk1[r];
          dNcur += a2Row[k] * ndu[r][pk];
        }
        std::swap(a2Row, a1Row); // Switch rows
      }
    }
    /* Multiply through by the correct factors */
    /* (Eq. [2.9])                             */
    for (int r = p, k = 1; k <= derivativeOrder; ++k) {
      for (int j = 0; j <= p; ++j)
        dN[k][j] *= r;
      r *= (p - k);
    }
    return dN;
  }

private:
  std::vector<ScalarType> knots_;
  int degree_;
};

template <std::size_t dim, typename ScalarType>
auto createUniqueKnotSpans(const std::array<std::vector<ScalarType>, dim>& knotSpans) {
  std::array<std::vector<double>, dim> uniqueKnotVector;
  for (int i = 0; i < dim; ++i) // create unique knotspan vectors
    std::ranges::unique_copy(knotSpans[i], std::back_inserter(uniqueKnotVector[i]),
                             [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

  return uniqueKnotVector;
}
} // namespace Dune::IGANEW::Splines
