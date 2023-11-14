// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDGEOMETRY_HH
#define DUNE_IDENTITYGRIDGEOMETRY_HH

/** \file
 * \brief The PatchGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/iga/hierarchicpatch/nurbspatchgeometry.hh>

namespace Dune::IGA {
  namespace Impl {
    enum class FixedOrFree2 { fixed, free };
  }



  template<int mydim, int coorddim, class GridImp>
  class PatchGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, PatchGridGeometry>
  {
  public:

    static constexpr std::integral auto mydimension = mydim;

    static constexpr std::integral auto coorddimension = coorddim;
    static constexpr std::integral auto griddim        = GridImp::dimension;
    using ctype                     = typename GridImp::ctype;
    using PatchGeometry =NURBSPatchGeometry<GridImp::dimension,coorddimension,ctype>;
    using LocalCoordinateInPatch           = typename PatchGeometry::LocalCoordinate;
    using LocalCoordinate           = FieldVector<ctype, mydimension>;
    using GlobalCoordinate          = FieldVector<ctype, coorddimension>;
    using JacobianTransposed        = FieldMatrix<ctype, mydimension, coorddimension>;
    using Hessian                   = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
    using Jacobian                  = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;
    using JacobianInverse           = FieldMatrix<ctype, mydimension, coorddimension>;
    using Volume                    = ctype;

    // The codimension of this entitypointer wrt the host grid
    constexpr static int CodimInHostGrid = GridImp::HostGridType::dimension - mydim;
    constexpr static int DimensionWorld = GridImp::HostGridType::dimensionworld;

    // select appropriate hostgrid geometry via typeswitch
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;

    typedef typename std::conditional<coorddim==DimensionWorld, HostGridGeometryType, HostGridLocalGeometryType>::type HostGridGeometry;
    auto getDirectionsOfSubEntityInParameterSpace() const {
      std::array<unsigned int, mydimension> subDirs;
      for (int subI = 0, i = 0; i < griddim; ++i) {
        if constexpr (mydimension != griddim)
          if (fixedOrFree[i] == Impl::FixedOrFree2::fixed) continue;
        subDirs[subI++] = i;
      }
      return subDirs;
    }

    //! type of the LocalView of the patch geometry
    using ElementGeometryLocalView= typename NURBSPatchGeometry<GridImp::dimension,coorddimension,ctype>::LocalView;

    /** constructor from host geometry
     */
    PatchGridGeometry(const HostGridGeometry& hostGeometry, ElementGeometryLocalView&& geometryLocalView)
      : hostGeometry_(hostGeometry),
    geometryLocalView_(std::forward<ElementGeometryLocalView>(geometryLocalView)){
      geometryLocalView_.bind(hostGeometry_.center());

      std::vector<LocalCoordinateInPatch> cornersVector;
      std::cout<<"Number of corners:"<<hostGeometry.corners()<<std::endl;
      for (int i = 0; i < hostGeometry.corners(); ++i)
        cornersVector.push_back(hostGeometry.corner(i));
      for (int i = 0; i < hostGeometry.corners(); ++i)
        std::cout<<cornersVector[i]<<std::endl;
      std::cout<<std::endl;

      const ctype tolerance = std::sqrt( std::numeric_limits< ctype >::epsilon() );
      LocalCoordinateInPatch average(0.0);
      for (int i = 0; i < hostGeometry.corners(); ++i)
        average+=hostGeometry.corner(i);
      average/=hostGeometry.corners();

      for (int i = 0; i < griddim; ++i)
        if(abs(average[i]-cornersVector[0][i])<tolerance)
          fixedOrFree[i]=Impl::FixedOrFree2::fixed;
        else
          fixedOrFree[i]=Impl::FixedOrFree2::free;
      //
      // }
      // hostGeometry_.corner(0);
      const auto& thisSpanIndices =geometryLocalView_.spanIndices();
      const auto& patchData =geometryLocalView_.patchData();

      for (int i = 0; i < griddim; ++i) {
        if (thisSpanIndices[i] + 1 < patchData.knotSpans[i].size())
          scaling_[i] = patchData.knotSpans[i][thisSpanIndices[i] + 1] - patchData.knotSpans[i][thisSpanIndices[i]];
        offset_[i] = patchData.knotSpans[i][thisSpanIndices[i]];
      }
    }

    template <typename ReturnType = LocalCoordinate>
LocalCoordinate spanToLocal(const LocalCoordinateInPatch& inSpan) const {
      LocalCoordinate localL;
      // for (int i = 0; i < griddim; ++i) {
      //   localL[i] = (inSpan[i] - offset_[i]) / scaling_[i];
      //   localL[i] = clampToDomain(localL[i], domain());
      // }

      return localL;
    }


decltype(auto) localToSpan(const LocalCoordinate& local) const {
      if constexpr (griddim==mydim)
        return local;
      else {
        LocalCoordinateInPatch localInSpan;
        // if constexpr (LocalCoordinate::dimension != 0) {
        //   for (int loci = 0, i = 0; i < griddim; ++i) {
        //     localInSpan[i] = (fixedOrVaryingDirections_[i] == Impl::FixedOrFree::free)
        //                          ? local[loci++] * scaling_[i] + offset_[i]
        //                          : offset_[i];
        //   }
        // } else
        //   for (int i = 0; i < griddim; ++i)
        //     localInSpan[i] = offset_[i];
        return localInSpan;
      }
    }


    /** \brief Return the element type identifier
     */
    [[nodiscard]] GeometryType type () const {
      return geometryLocalView_.type();
    }

    // return whether we have an affine mapping
    [[nodiscard]] bool affine() const {
      return geometryLocalView_.affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    [[nodiscard]] int corners () const {
      return hostGeometry_.corners();
    }


    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner (int i) const {
      return geometryLocalView_.corner(i);
    }


    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    GlobalCoordinate global (const LocalCoordinate& local) const {
      auto localInGrid = hostGeometry_.global(local);
      return geometryLocalView_.global(localInGrid);
    }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed
    jacobianTransposed ( const LocalCoordinate& local ) const {
      JacobianTransposed result;
      // auto localInGrid = hostGeometry_.global(local);
      // std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();
      //
      // const auto basisFunctionDerivatives = nurbs_.basisFunctionDerivatives(localInSpan, 1);
      //
      // for (int dir = 0; dir < mydimension; ++dir) {
      //   std::array<unsigned int, griddim> ithVec{};
      //   ithVec[subDirs[dir]] = 1;
      //   result[dir]          = dot(basisFunctionDerivatives.get(ithVec), cpCoordinateNet_);
      //   result[dir] *= scaling_[subDirs[dir]];  // transform back to 0..1 domain
      // }
      // return result;
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    LocalCoordinate local (const GlobalCoordinate& global) const {

      // return spanToLocal(geometryLocalView_.local(global));
    }


    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<ctype, mydim> &local) const {
      return geometryLocalView_.checkInside(local);
    }


    /**
     */
    [[nodiscard]] ctype integrationElement (const LocalCoordinate& local) const {
      // return geometryLocalView_.integrationElement(local);
    }


    //! The Jacobian matrix of the mapping from the reference element to this element
    [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      JacobianTransposed result;
      std::array<unsigned int, mydimension> subDirs = getDirectionsOfSubEntityInParameterSpace();

      const auto localInSpan              = localToSpan(local);
      const auto basisFunctionDerivatives = geometryLocalView_.nurbs().basisFunctionDerivatives(localInSpan, 1);

      for (int dir = 0; dir < mydimension; ++dir) {
        std::array<unsigned int, griddim> ithVec{};
        ithVec[subDirs[dir]] = 1;
        result[dir]          = dot(basisFunctionDerivatives.get(ithVec), geometryLocalView_.controlPointCoordinates());
        result[dir] *= scaling_[subDirs[dir]];  // transform back to 0..1 domain
      }
      return result;
    }
      // return geometryLocalView_.jacobianInverseTransposed(local);


  private:
    auto getDirectionsOfSubEntityInParameterSpace() const {
      // std::array<unsigned int, mydimension> subDirs;
      // for (int subI = 0, i = 0; i < griddim; ++i) {
      //   if constexpr (mydimension != griddim)
      //     if (fixedOrVaryingDirections_[i] == Dune::IGA::Impl::FixedOrFree2::fixed) continue;
      //   subDirs[subI++] = i;
      // }
      // return subDirs;
    }

    HostGridGeometry hostGeometry_;
    ElementGeometryLocalView geometryLocalView_{};
    std::array<ctype, griddim> offset_;
    std::array<ctype, griddim> scaling_;
    std::array<Impl::FixedOrFree2,griddim> fixedOrFree;

  };

}  // namespace Dune

#endif
