//
// Created by lex on 07.11.21.
//

#pragma once
#include <dune/iga/traits.hh>
#include <concepts>
#include <ranges>
#include <algorithm>

namespace Dune::IGA
{
  template <std::ranges::random_access_range Range>
  auto findUpperSpanIndex(const int p, const typename std::remove_cvref_t<Range>::value_type u, Range&& U) {
    if (u <= U[0]) return static_cast<long int>(p);
    auto it = std::upper_bound(U.begin() + p-1, U.end() - p-1, u);
    return std::distance(U.begin(), it) - 1;
  }


  template<typename ContainerType>
  void resize(ContainerType&& c, int size)
  {
    if constexpr (requires (ContainerType c2 ){  c2.resize(size);})
        c.resize(size);
  }

  template <std::floating_point ScalarType, int degreeT= 0>
  struct Bspline{

//    static constexpr long unsigned int degree = (degreeT==-1 ) ? 0 : degreeT;

// The Nurbs Book Algorithm A2.2
template<typename ContainerType= std::array<double,degreeT+1>>
  static auto basisFunctions(     ScalarType u,
      int degree,  const std::span<ScalarType>& knots
  )
  {
    ContainerType N;
    const int p = degree;
    resize(N,p+1);
    const int i = findUpperSpanIndex( degree,u, knots);

    auto lDiff = std::ranges::transform_view(std::ranges::reverse_view(std::views::counted(knots.begin()+i+1-p,p)) ,[&u](auto& knot){return u-knot;}) ;
    auto rDiff = std::ranges::transform_view(std::views::counted(knots.begin()+i+1,p) ,[&u](auto& knot){return knot-u;}) ;

    N[0] = ScalarType(1);
    for (int j=1; j<=p; ++j)
    {
      ScalarType saved = ScalarType(0);
      for (int r=0; r<j; r++)
      {
        const ScalarType tmp = N[r]/(rDiff[r]+lDiff[j-r-1]);
        N[r] = saved + rDiff[r]*tmp;
        saved = lDiff[j-r-1]*tmp;
      }
      N[j] = saved;
    }
    return N;
  }
  };
}
