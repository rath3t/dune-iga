//
// Created by lex on 16.11.21.
//

#pragma once
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>

namespace Dune::IGA {
  /** \brief a geometry implementation for NURBS*/
  template <std::integral auto mydim, std::integral auto dimworld, std::integral auto griddim,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits> class NURBSLocalGeometry {
  public:
    /** coordinate type */
    typedef double ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto coorddimension = griddim;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<mydim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<griddim>;

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposed = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<mydim, coorddimension>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<coorddimension, mydim>;

    using ControlPointType = typename NURBSPatchData<griddim, dimworld , NurbsGridLinearAlgebraTraits>::ControlPointType;
    using ControlPointNetType =
        typename NURBSPatchData<griddim, dimworld , NurbsGridLinearAlgebraTraits>::ControlPointNetType;

  private:
    /* Helper class to compute a matrix pseudo inverse */
    typedef MultiLinearGeometryTraits<ctype>::MatrixHelper MatrixHelper;

  public:
    /** \brief Constructor from NURBSPatchData and an iterator to a specific knot
     *
     *  \param[in] Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
     *  \param[in] corner Iterator (for each dimension) to the Knot span where the Geometry object is supposed to operate
     */
    NURBSLocalGeometry(int localSubEntityIndex)
        : localIndexInElement_{localSubEntityIndex}

    {
//      referenceElement_= Dune::Geo::ReferenceElements<ctype,mydim>::cube();
      assert((localIndexInElement_ < 2 && griddim == 1) || (localIndexInElement_ < 4 && griddim == 2)
             || (localIndexInElement_ < 6 && griddim == 3));
    }

    /** \brief Map the center of the element to the geometry */
    GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    [[nodiscard]] double volume() const { return 1; }

    /** \brief Return the number of corners of the element */
    [[nodiscard]] int corners() const { return 1 << mydimension; }

    /** \brief Return world coordinates of the k-th corner of the element */
    GlobalCoordinate corner(int k) const {
      LocalCoordinate localcorner;
      for (size_t i = 0; i < mydimension; i++) {
        localcorner[i] = (k & (1 << i)) ? 1 : 0;
      }
      return global(localcorner);
    }

    [[nodiscard]] bool affine() const { return false; }

    /** \brief Type of the element: a hypercube of the correct dimension */
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydim); }

    /** \brief evaluates the NURBS mapping
     *
     *  \param[in] local local coordinates for each dimension
     */
    GlobalCoordinate global(const LocalCoordinate& local) const {
      const double offset = 1e-15;
      if constexpr (mydim == 0)
        return (localIndexInElement_ == 0) ? GlobalCoordinate(0) : GlobalCoordinate(1-offset);
      else if constexpr (mydim == 1)
        switch (localIndexInElement_) {
          case 0:
            return GlobalCoordinate({0, local[0]});
          case 1:
            return GlobalCoordinate({1-offset, local[0]});
          case 2:
            return GlobalCoordinate({local[0], 0});
          case 3:
            return GlobalCoordinate({local[0], 1-offset});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydim == 3)
        switch (localIndexInElement_) {
          case 0:
            return GlobalCoordinate({0, local[0], local[1]});
          case 1:
            return GlobalCoordinate({1-offset, local[0], local[1]});
          case 2:
            return GlobalCoordinate({local[0], 0, local[1]});
          case 3:
            return GlobalCoordinate({local[0], 1-offset, local[1]});
          case 4:
            return GlobalCoordinate({local[0], local[1], 0});
          case 5:
            return GlobalCoordinate({local[0], local[1], 1-offset});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();
    }

    LocalCoordinate local(const GlobalCoordinate& global) const {
      if constexpr (mydim == 0)
        return (localIndexInElement_ == 0) ? GlobalCoordinate(0) : GlobalCoordinate(1);
      else if constexpr (mydim == 1)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate({global[1]});
          case 2:
          case 3:
            return LocalCoordinate({global[0]});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydim == 3)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate({global[1], global[2]});
          case 2:
          case 3:
            return LocalCoordinate({global[0], global[2]});
          case 4:
          case 5:
            return LocalCoordinate({global[0], global[1]});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param[in] local local coordinates for each dimension
     */
    JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      if constexpr (mydim == 0)
        return 0;
      else if constexpr (mydim == 1)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return JacobianTransposed({FieldVector<ctype,mydim>({0,1})});
          case 2:
          case 3:
            return JacobianTransposed({FieldVector<ctype,mydim>({1,0})});
          default:
            __builtin_unreachable();
        }
      else if constexpr (mydim == 3)
        switch (localIndexInElement_) {
          case 0:
          case 1:
            return LocalCoordinate({FieldVector<ctype,mydim>({0,1,0}), FieldVector<ctype,mydim>({0,0,1})});
          case 2:
          case 3:
            return LocalCoordinate({FieldVector<ctype,mydim>({1,0,0}), FieldVector<ctype,mydim>({0,0,1})});
          case 4:
          case 5:
            return LocalCoordinate({FieldVector<ctype,mydim>({1,0,0}), FieldVector<ctype,mydim>({0,1,0})});
          default:
            __builtin_unreachable();
        }
      __builtin_unreachable();

    }

    ctype integrationElement(const LocalCoordinate& local) const {
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(jacobianTransposed(local));
    }

    JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianInverseTransposed jacobianInverseTransposed1;
      MatrixHelper::template rightInvA<mydimension, coorddimension>(jacobianTransposed(local), jacobianInverseTransposed1);
      return jacobianInverseTransposed1;
    }

    GlobalCoordinate unitNormal(const LocalCoordinate& local) const requires(mydimension == 2) {
      auto J = jacobianTransposed(local);
      auto N = cross(J[0], J[1]);
      return N / N.two_norm();
    }

  private:
      Dune::Geo::ReferenceElements<ctype,mydim> referenceElement_;
    int localIndexInElement_;
  };

}  // namespace Dune::IGA