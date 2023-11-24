

#pragma once
#include "dune/iga/hierarchicpatch/concepts.hh"

namespace Dune::IGANEW {
  template <Concept::Vector VectorType>
  VectorType projectPointOntoLine(const VectorType basepoint, const VectorType revolutionaxis, const VectorType point) {
    const VectorType e1 = basepoint + revolutionaxis - basepoint;
    const VectorType e2 = point - basepoint;
    using Dune::dot;
    const typename VectorType::value_type angle = dot(e1, e2);
    const typename VectorType::value_type len2  = e1.two_norm();

    VectorType p = basepoint + angle * e1 / len2;
    return p;
  }

  template <Concept::Vector VectorType>
  VectorType intersect3DLines(const VectorType basepoint1, const VectorType direction1, const VectorType basepoint2,
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

    if (std::abs(D) < tol) DUNE_THROW(Dune::MathError, "The two lines are almost parallel.");

    const ScalarType parameter1 = (angle * dir2BasePdiff - dir2squaredLength * dir1BasePdiff) / D;
    const ScalarType parameter2 = (dir1squaredLength * dir2BasePdiff - angle * dir1BasePdiff) / D;

    VectorType c1 = basepoint1 + parameter1 * direction1;
    assert(dot(c1 - (basepoint2 + parameter2 * direction2), c1 - (basepoint2 + parameter2 * direction2)) < tol
           && "Both calculated points do not coincide");
    return c1;
  }

}  // namespace Dune::IGANEW
