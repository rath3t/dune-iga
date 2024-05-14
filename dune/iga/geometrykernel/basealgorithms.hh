

#pragma once
#include "dune/iga/hierarchicpatch/concepts.hh"

namespace Dune::IGA {
/**
 * @brief Projects a point onto a line defined by a basepoint and direction.
 *
 * This function calculates the projection of a given point onto a line
 * defined by a basepoint and a direction vector.
 *
 * @tparam VectorType A vector type satisfying the Concept::Vector concept.
 * @param basepoint The basepoint of the line.
 * @param direction The direction vector of the line.
 * @param point The point to be projected onto the line.
 * @return The projected of the point.
 *
 * Example:
 * @code
 * using Vec3D = Dune::FieldVector<double, 3>; // Assuming Dune::FieldVector for vector operations
 * Vec3D basepoint(0.0, 0.0, 0.0);
 * Vec3D direction(1.0, 0.0, 0.0);
 * Vec3D point(2.0, 1.0, 1.0);
 * Vec3D projectedPoint = projectPointOntoLine(basepoint, direction, point);
 * // projectedPoint = [2.0, 0.0, 0.0]
 * @endcode
 */
template <Concept::Vector VectorType>
VectorType projectPointOntoLine(const VectorType basepoint, const VectorType direction, const VectorType point) {
  const VectorType e1 = basepoint + direction - basepoint;
  const VectorType e2 = point - basepoint;
  using Dune::dot;
  const typename VectorType::value_type angle = dot(e1, e2);
  const typename VectorType::value_type len2  = e1.two_norm();

  VectorType p = basepoint + angle * e1 / len2;
  return p;
}

/**
 * @brief Finds the intersection point of two 3D lines.
 *
 * This function calculates the intersection point of two 3D lines defined
 * by their basepoints and direction vectors. The lines are represented by
 * the equations:
 *
 * Line 1: \f$P_1(t) = \text{{basepoint1}} + t \cdot \text{{direction1}}\f$
 *
 * Line 2: \f$P_2(s) = \text{{basepoint2}} + s \cdot \text{{direction2}}\f$
 *
 * where \f$t\f$ and \f$s\f$ are parameters. The function returns the point
 * where the two lines intersect.
 *
 * @tparam VectorType A vector type satisfying the Concept::Vector concept.
 * @param basepoint1 The basepoint of the first line.
 * @param direction1 The direction vector of the first line.
 * @param basepoint2 The basepoint of the second line.
 * @param direction2 The direction vector of the second line.
 * @return The intersection point of the two lines.
 *
 * @throws Dune::MathError if the lines are almost parallel (within tolerance).
 */
template <Concept::Vector VectorType>
VectorType intersect3DLines(const VectorType basepoint1, const VectorType direction1, const VectorType basepoint2,
                            const VectorType direction2) {
  using ScalarType     = typename VectorType::value_type;
  const ScalarType tol = 1e-8;

  const VectorType basePointDifference = basepoint1 - basepoint2;
  using Dune::dot;
  const ScalarType dir1squaredLength = dot(direction1, direction1);
  const ScalarType angle             = dot(direction1, direction2);
  const ScalarType dir2squaredLength = dot(direction2, direction2);
  const ScalarType dir1BasePdiff     = dot(direction1, basePointDifference);
  const ScalarType dir2BasePdiff     = dot(direction2, basePointDifference);
  const ScalarType D                 = dir1squaredLength * dir2squaredLength - angle * angle;

  if (std::abs(D) < tol)
    DUNE_THROW(Dune::MathError, "The two lines are almost parallel.");

  const ScalarType parameter1 = (angle * dir2BasePdiff - dir2squaredLength * dir1BasePdiff) / D;
  const ScalarType parameter2 = (dir1squaredLength * dir2BasePdiff - angle * dir1BasePdiff) / D;

  VectorType c1 = basepoint1 + parameter1 * direction1;
  assert(dot(c1 - (basepoint2 + parameter2 * direction2), c1 - (basepoint2 + parameter2 * direction2)) < tol &&
         "Both calculated points do not coincide");
  return c1;
}

} // namespace Dune::IGA
