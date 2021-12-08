//
// Created by lex on 07.11.21.
//

#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>

#include <dune/common/dynmatrix.hh>
#include <dune/iga/dunelinearalgebratraits.hh>

namespace Dune::IGA {
  template <std::ranges::random_access_range Range>
  auto findSpan(const int p, const typename std::remove_cvref_t<Range>::value_type u, Range&& U, int offset = 0) {
    if (u <= U[0]) return static_cast<long int>(p);
    auto it = std::upper_bound(U.begin() + p - 1 + offset, U.end(), u);
    return std::distance(U.begin(), it) - 1;
  }

  template <std::size_t dim, typename ValueType>
  auto findSpan(const std::array<int, dim>& p, const std::array<ValueType, dim>& u, const std::array<std::vector<ValueType>, dim>& U) {
    std::array<int, dim> res;
    for (auto i = 0; i < dim; ++i)
      res[i] = findSpan(p[i], u[i], U[i]);
    return res;
  }

  template <typename ContainerType>
  void resize(ContainerType&& c, int size) {
    if constexpr (requires(ContainerType c2) { c2.resize(size); }) c.resize(size);
  }

  template <std::floating_point ScalarType>
  class Bspline {
  public:
    Bspline(const std::vector<ScalarType>& knots, const int degree) : knots_{knots}, degree_{degree} {}

    template <typename ContainerType = std::vector<ScalarType>>
    auto operator()(ScalarType u) {
      return basisFunctions<ContainerType>(u, knots_, degree_);
    }

    // The Nurbs Book Algorithm A2.2
    template <typename ContainerType = std::vector<ScalarType>>
    static auto basisFunctions(ScalarType u, const std::vector<ScalarType>& knots, const int degree,
                               std::optional<int> spIndex = std::nullopt) {
      assert(std::ranges::count(knots.begin(), knots.begin() + degree + 1, knots.front()) == degree + 1);
      assert(std::ranges::count(knots.end() - degree - 1, knots.end(), knots.back()) == degree + 1);
      assert(spIndex < knots.size() - degree - 1);
      ContainerType N;
      const int p = degree;
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

      return N;
    }

    // The Nurbs Book Algorithm A2.3
    static auto basisFunctionDerivatives(ScalarType u, const std::vector<ScalarType>& knots, const int degree, const int derivativeOrder,
                                         std::optional<int> spIndex = std::nullopt) {
      assert(spIndex < knots.size() - degree - 1);
      //      if(spIndex.has_value() and spIndex.value()==knots.size()) {
      ////        spIndex.value()-= degree-2;
      //        u -= 1e-4;
      ////        --spIndex.value();
      //      }
      const int order = degree + 1;
      const int p     = degree;
      const int sp    = spIndex ? spIndex.value() : findSpan(p, u, knots);
      using namespace std::ranges;
      auto lDiff = transform_view(reverse_view(std::views::counted(begin(knots) + sp + 1 - p, p)), [&u](auto& kn) { return u - kn; });
      auto rDiff = transform_view(std::views::counted(begin(knots) + sp + 1, p), [&u](auto& kn) { return kn - u; });

      DynamicMatrix<ScalarType> dN(derivativeOrder + 1, order);

      std::vector<ScalarType> left(order);
      std::vector<ScalarType> right(order);
      DynamicMatrix<ScalarType> ndu(order, order);

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
      DynamicMatrix<ScalarType> a(2, p + 1);
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
}  // namespace Dune::IGA
