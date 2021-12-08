//
// Created by alex on 12/8/21.
//

#pragma once
#include <dune/iga/concepts.hh>
namespace Dune::IGA {

  /** \brief Free cross product function it calls the member function cross of VectorType if it exists and falls back to an implementation
   * by hand otherwise
   */
  template <Vector VectorType>
  requires(VectorType::dimension == 3) inline VectorType cross(const VectorType& a, const VectorType& b) {
    if constexpr (requires { a.cross(b); })
      return a.cross(b);
    else
      return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }
}  // namespace Dune::IGA