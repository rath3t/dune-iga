//
// Created by lex on 07.11.21.
//

#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>

#include <dune/common/dynmatrix.hh>
#include <dune/common/float_cmp.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/concepts.hh>

namespace Dune::IGA {

  /** \brief Finds the spanIndex in range [u_0,...,u_0,...,u_a      ,u,...,u_n,...,u_n] which is first index lower than u
   *                                        -p+1-times-      returned       -p+1-times-
   * @tparam Range knotvector range
   * @param p polynomial degree of the underlying spline
   * @param u span value which is searched for
   * @param U range where to search (knotvector
   * @param offset adjust range where to start searching to improve efficiency
   * @return
   */
  template <std::ranges::random_access_range Range>
  auto findSpan(const int p, const typename std::remove_cvref_t<Range>::value_type u, Range&& U, int offset = 0) {
    if (u <= U[0]) return static_cast<long int>(p);
    auto it = std::upper_bound(U.begin() + p - 1 + offset, U.end(), u);
    return std::distance(U.begin(), it) - 1;
  }

  /** \brief Same as findSpan() but for dim - knotvectors  */
  template <std::size_t dim, typename ValueType>
  auto findSpan(const std::array<int, dim>& p, const std::array<ValueType, dim>& u, const std::array<std::vector<ValueType>, dim>& U) {
    std::array<int, dim> res;
    for (auto i = 0; i < dim; ++i)
      res[i] = findSpan(p[i], u[i], U[i]);
    return res;
  }

  /** \brief One dimensional b-spline basis
   *
   * @tparam T is either the scalar type pf the point where to evaluate or a specif non-default LinearAlgebraTraits where
   * DynamicMatrixType and DynamicVectorType,... is derived from
   */
  template <typename T = DuneLinearAlgebraTraits<double>>
  requires LinearAlgebra<T> || std::floating_point<T>
  class BsplineBasis1D;

  /** \brief One dimensional b-spline basis
   *
   * @tparam NurbsGridLinearAlgebraTraits Specialization where the Traits are directly given
   */
  template <LinearAlgebra NurbsGridLinearAlgebraTraits>
  class BsplineBasis1D<NurbsGridLinearAlgebraTraits> {
  public:
    using ScalarType        = typename NurbsGridLinearAlgebraTraits::value_type;
    using DynamicVectorType = typename NurbsGridLinearAlgebraTraits::DynamicVectorType;
    using DynamicMatrixType = typename NurbsGridLinearAlgebraTraits::DynamicMatrixType;
    template <int cols = 0>
    using RowFixedMatrix = typename NurbsGridLinearAlgebraTraits::template RowFixedMatrix<cols>;

    BsplineBasis1D(const std::vector<ScalarType>& knots, const int degree) : knots_{knots}, degree_{degree} {}

    auto operator()(ScalarType u) { return basisFunctions(u, knots_, degree_); }

    /** \brief The evaluation function modified version of The Nurbs Book Algorithm A2.2
     *
     * @tparam Range range of knotvector
     * @tparam ContainerType Container where the returned values should be stored, it needs to satisfy StdVectorLikeContainer
     * @param u scalar point where to evaluate the spline
     * @param knots knotvector
     * @param p polynomial degree of the underlying b-spline basis
     * @param spIndex optional index in which range the evaluationg point lies, if omited it is searched for
     * @return Non-zero B-spline basisfunctions evaluated at u
     */
    template <std::ranges::random_access_range Range>
    static auto basisFunctions(ScalarType u, Range&& knots, const int p, std::optional<int> spIndex = std::nullopt) {
      assert(std::ranges::count(knots.begin(), knots.begin() + p + 1, knots.front()) == p + 1);
      assert(std::ranges::count(knots.end() - p - 1, knots.end(), knots.back()) == p + 1);
      assert(spIndex < knots.size() - p - 1);
      DynamicVectorType N;
      N.resize(p + 1, 0.0);
      if (Dune::FloatCmp::eq(u, knots.back()))  // early exit
      {
        N.back() = 1;
        return N;
      } else if (Dune::FloatCmp::eq(u, knots.front())) {
        N.front() = 1;
        return N;
      }

      const int sp = spIndex ? spIndex.value() : findSpan(p, u, knots);
      using namespace std::ranges;
      auto lDiff = transform_view(reverse_view(std::views::counted(knots.begin() + sp + 1 - p, p)), [&u](auto& kn) { return u - kn; });
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
      for ([[maybe_unused]] auto& Ni : N)
        assert(Dune::FloatCmp::ge(Ni, 0.0)); // The basis functions are always >=0!

      return N;
    }

    /** \brief The evaluation function modified version of The Nurbs Book Algorithm A2.3
     *
     * @param u scalar point where to evaluate the spline
     * @param knots knotvector
     * @param p polynomial p of the underlying b-spline basis
     * @param derivativeOrder maximum derivative degree that are constructed
     * @param spIndex optional index in which range the evaluationg point lies, if omited it is searched for
     * @return Non-zero B-spline derivatives  evaluated at u
     */
    static auto basisFunctionDerivatives(ScalarType u, const std::vector<ScalarType>& knots, const int p, const int derivativeOrder,
                                         std::optional<int> spIndex = std::nullopt) {
      assert(spIndex < knots.size() - p - 1);
      const int order = p + 1;
      const int sp    = spIndex ? spIndex.value() : findSpan(p, u, knots);
      using namespace std::ranges;
      auto lDiff = transform_view(reverse_view(std::views::counted(begin(knots) + sp + 1 - p, p)), [&u](auto& kn) { return u - kn; });
      auto rDiff = transform_view(std::views::counted(begin(knots) + sp + 1, p), [&u](auto& kn) { return kn - u; });

      DynamicMatrixType dN(derivativeOrder + 1, order);

      std::vector<ScalarType> left(order);
      std::vector<ScalarType> right(order);
      DynamicMatrixType ndu(order, order);

      ndu[0][0] = 1.0;
      {
        int j = 0;
        for (auto lDp = lDiff.begin(); lDp != lDiff.end(); lDp += j + 2, ++j) {
          ScalarType saved{};
          for (int r = 0; auto&& rD : rDiff | std::views::take(j + 1)) {
            const auto& lD = (*lDp);
            /* Lower triangle */
            ndu[j + 1][r]  = rD + lD;
            ScalarType tmp = ndu[r][j] / ndu[j + 1][r];
            /* Upper triangle */
            ndu[r++][j + 1] = saved + rD * tmp;
            saved           = lD * tmp;
            --lDp;
          }
          ndu[j + 1][j + 1] = saved;
        }
      }
      for (int j = p; j >= 0; --j)
        dN[0][j] = ndu[j][p];

      // Compute the derivatives
      RowFixedMatrix<2> a(2, p + 1);
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
          std::swap(a2Row, a1Row);  // Switch rows
        }
      }
      /* Multiply through by the correct factors */
      /* (Eq. [2.9])                             */
      for (int r = p, k = 1; k <= derivativeOrder; ++k) {
        for (int j = p; j >= 0; --j)
          dN[k][j] *= r;
        r *= p - k;
      }
      return dN;
    }

  private:
    std::vector<ScalarType> knots_;
    int degree_;
  };

  /** \brief One dimensional b-spline basis
   *
   * @tparam ScalarType Specialization where the scalar type is given and the LinearAlgebraTraits are defaulted
   */
  template <std::floating_point ScalarType>
  class BsplineBasis1D<ScalarType> : public BsplineBasis1D<DuneLinearAlgebraTraits<ScalarType>> {
    using Base = BsplineBasis1D<DuneLinearAlgebraTraits<ScalarType>>;

  public:
    BsplineBasis1D(const std::vector<ScalarType>& knots, const int degree) : Base(knots, degree) {}
  };

}  // namespace Dune::IGA
