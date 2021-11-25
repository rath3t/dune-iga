//
// Created by lex on 16.11.21.
//

#pragma once
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>

namespace Dune::IGA {
  namespace Impl {
    enum class FixedOrFree { fixed, free };
  }

  /** \brief a geometry implementation for NURBS*/
  template <std::integral auto mydim, std::integral auto dimworld, std::integral auto griddim,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
  class NURBSGeometry {
  public:
    /** coordinate type */
    typedef double ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto coorddimension = dimworld;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<mydim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = typename NurbsGridLinearAlgebraTraits::template FixedVectorType<coorddimension>;

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposed = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<mydim, coorddimension>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = typename NurbsGridLinearAlgebraTraits::template FixedMatrixType<coorddimension, mydim>;

    using ControlPointType    = typename NURBSPatchData<griddim, dimworld, NurbsGridLinearAlgebraTraits>::ControlPointType;
    using ControlPointNetType = typename NURBSPatchData<griddim, dimworld, NurbsGridLinearAlgebraTraits>::ControlPointNetType;

  private:
    /* Helper class to compute a matrix pseudo inverse */
    typedef MultiLinearGeometryTraits<ctype>::MatrixHelper MatrixHelper;

  public:
    /** \brief Constructor from NURBSPatchData and an iterator to a specific knot
     *
     *  \param[in] Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
     *  \param[in] corner Iterator (for each dimension) to the Knot span where the Geometry object is supposed to operate
     */
    NURBSGeometry(std::shared_ptr<NURBSPatchData<griddim, dimworld, NurbsGridLinearAlgebraTraits>> patchData,
                  const std::array<std::vector<double>::const_iterator, griddim>& varyingSpans,
                  //                  const std::array<double, griddim - mydimension>& fixedSpans,
                  const std::array<Impl::FixedOrFree, griddim>& fixedOrVaryingDirections)
        : patchData_(patchData),
          varyingSpans_(varyingSpans),
          fixedOrVaryingDirections_{fixedOrVaryingDirections},
          nurbs_{*patchData}
    //          cpCoordinateNet_{netOfSpan(transform01ToKnotSpan(std::array<double, griddim>{}, varyingSpans_), patchData->knotSpans,
    //                                     patchData->order, extractControlpoints(patchData->controlPoints))}
    {}

    /** \brief Map the center of the element to the geometry */
    GlobalCoordinate center() const {
      LocalCoordinate localcenter(0.5);
      return global(localcenter);
    }

    [[nodiscard]] double volume() const {
      const auto rule = Dune::QuadratureRules<ctype, mydimension>::rule(this->type(), (*std::ranges::max_element(patchData_->order)));
      ctype vol       = 0.0;
      for (auto& gp : rule)
        vol += integrationElement(gp.position()) * gp.weight();
      return vol;
    }

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

    template <typename ReturnType = std::array<typename LocalCoordinate::value_type, griddim>>
    auto transformLocalToSpan(const LocalCoordinate& local) const {
      ReturnType localInSpan;
      if constexpr (local.size() != 0)
        for (int loci = 0, i = 0; i < griddim; ++i) {
          localInSpan[i] = (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::free)
                               ? local[loci++] * (*(varyingSpans_[i] + 1) - *varyingSpans_[i]) + *varyingSpans_[i]
                               : *varyingSpans_[i];
        }
      else
        for (int i = 0; i < griddim; ++i)
          localInSpan[i] = *(varyingSpans_[i]);
      return localInSpan;
    }
    /** \brief evaluates the NURBS mapping
     *
     *  \param[in] local local coordinates for each dimension
     */
    GlobalCoordinate global(const LocalCoordinate& local) const {
      const auto localInSpan = transformLocalToSpan(local);
      auto basis             = nurbs_.basisFunctionNet(localInSpan);
      auto cpNetofSpan = netOfSpan(localInSpan, patchData_->knotSpans, patchData_->order, extractControlpoints(patchData_->controlPoints));
      return dot(basis, cpNetofSpan);
    }

    LocalCoordinate local(const GlobalCoordinate& global) const {
      const ctype tolerance = ctype(16) * std::numeric_limits<ctype>::epsilon();
      LocalCoordinate x     = LocalCoordinate(0.5);
      LocalCoordinate dx{};
      do {  // from multilinearGeometry
        // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
        const GlobalCoordinate dglobal = (*this).global(x) - global;
        MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);
        const bool invertible = MatrixHelper::template xTRightInvA<mydimension, coorddimension>(jacobianTransposed(x), dglobal, dx);

        if (!invertible) return LocalCoordinate(std::numeric_limits<ctype>::max());
        x -= dx;
      } while (dx.two_norm2() > tolerance);
      return x;
    }

    /** \brief compute the Jacobian transposed matrix
     *
     *  \param[in] local local coordinates for each dimension
     */
    JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
      JacobianTransposed result;
      std::array<unsigned int, mydim> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }
      const auto localInSpan              = transformLocalToSpan(local);
      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(localInSpan, 1);
      auto cpNetofSpan = netOfSpan(localInSpan, patchData_->knotSpans, patchData_->order, extractControlpoints(patchData_->controlPoints));
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpNetofSpan);
        result[dir] *= *(varyingSpans_[subDirs[dir]] + 1) - *varyingSpans_[subDirs[dir]];  // transform back to 0..1 domain
      }
      return result;
    }

    ctype integrationElement(const LocalCoordinate& local) const {
      return MatrixHelper::template sqrtDetAAT<mydimension, coorddimension>(jacobianTransposed(local));
    }

    JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
      JacobianInverseTransposed jacobianInverseTransposed1;
      MatrixHelper::template rightInvA<mydimension, coorddimension>(jacobianTransposed(local), jacobianInverseTransposed1);
      return jacobianInverseTransposed1;
    }

    ctype gaussianCurvature(const LocalCoordinate& local) const requires(mydimension == 2) {
      auto metricDet = metric(local).determinant();
      auto secondF   = secondFundamentalForm(local).determinant();
      return secondF / metricDet;
    }

    auto metric(const LocalCoordinate& local) const {
      const auto J = jacobianTransposed(local);
      FieldMatrix<ctype, mydimension, mydimension> metric;
      MatrixHelper::AAT(J, metric);
      return metric;
    }

    GlobalCoordinate unitNormal(const LocalCoordinate& local) const requires(mydimension == 2 && coorddimension=3) {
      auto N = normal(local);
      return N / N.two_norm();
    }

    GlobalCoordinate normal(const LocalCoordinate& local) const requires(mydimension == 2 && coorddimension=3) {
      auto J = jacobianTransposed(local);
      return cross(J[0], J[1]);
    }

    auto secondFundamentalForm(const LocalCoordinate& local) const requires(mydimension == 2) {
      const auto secDerivatives = secondDerivativeOfPosition(local);
      const auto unitnormal     = unitNormal(local);
      FieldMatrix<ctype, mydimension, mydimension> b;
      b[0][0] = secDerivatives[0] * unitnormal;
      b[1][1] = secDerivatives[1] * unitnormal;
      b[0][1] = b[1][0] = secDerivatives[2] * unitnormal;
      return b;
    }

    auto secondDerivativeOfPosition(const LocalCoordinate& local) const {
      FieldMatrix<ctype, coorddimension, mydimension*(mydimension + 1) / 2> result;
      std::array<unsigned int, mydim> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if (fixedOrVaryingDirections_[i] == Dune::IGA::Impl::FixedOrFree::fixed) continue;
        subDirs[subI++] = i;
      }
      const auto localInSpan              = transformLocalToSpan(local);
      const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(localInSpan, 2);
      auto cpNetofSpan = netOfSpan(localInSpan, patchData_->knotSpans, patchData_->order, extractControlpoints(patchData_->controlPoints));
      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 2;  // second derivative in dir direction
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpNetofSpan);
        result[dir] *= Dune::power(*(varyingSpans_[subDirs[dir]] + 1) - *varyingSpans_[subDirs[dir]], 2);  // transform back to 0..1
      }
      std::array<int, mydimension> mixeDerivs;
      std::ranges::fill_n(mixeDerivs.begin(), mydimension, 1);  // first mixed derivatives
      int mixedDireCounter = mydimension;
      do {
        result[mixedDireCounter++] = dot(basisFunctionDerivatives.get(mixeDerivs), cpNetofSpan);
        for (int dir = 0; dir < mixeDerivs.size(); ++dir) {
          if (mixeDerivs[dir] == 0) continue;
          result[mixedDireCounter - 1] *= *(varyingSpans_[subDirs[dir]] + 1) - *varyingSpans_[subDirs[dir]];
        }
      } while (std::ranges::next_permutation(mixeDerivs, std::greater()).found);

      return result;
    }

  private:
    std::shared_ptr<NURBSPatchData<griddim, dimworld, NurbsGridLinearAlgebraTraits>> patchData_;
    std::array<std::vector<double>::const_iterator, griddim> varyingSpans_;
    std::array<Impl::FixedOrFree, griddim> fixedOrVaryingDirections_{free};
    Dune::IGA::Nurbs<double, griddim> nurbs_;
    //    MultiDimensionNet<griddim, typename ControlPointType::VectorType> cpCoordinateNet_;
  };

  template <std::integral auto mydim, std::integral auto dimworld, std::integral auto griddim,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits>
  auto referenceElement(const NURBSGeometry<mydim, dimworld, griddim, NurbsGridLinearAlgebraTraits>& geo) {
    return Dune::ReferenceElements<typename NurbsGridLinearAlgebraTraits::value_type, mydim>::cube();
  };
}  // namespace Dune::IGA