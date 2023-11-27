
#pragma once

#include "elementtrimdata.hh"
namespace Dune {
  namespace IGANEW {
    namespace Trim {

      template <int dim, typename ScalarType = double>
      struct DefaultTrimmer;

      template <int mydim_, int coorddim, typename ScalarType>
      class DefaultTrimmedPatchLocalGeometry {
       public:
        using ctype = ScalarType;

        static constexpr int mydimension = mydim_;
        using TrimmerType                = DefaultTrimmer<mydimension, ctype>;
        using TrimDataType               = DefaultElementTrimData<mydimension, ctype>;

        static constexpr int coorddimension = coorddim;
        static constexpr int codim          = coorddim - mydimension;
        using PatchGeometry                 = GeometryKernel::NURBSPatch<mydimension, coorddimension, ctype>;
        using LocalCoordinateInPatch        = typename PatchGeometry::LocalCoordinate;
        using LocalCoordinate               = FieldVector<ctype, mydimension>;
        using GlobalCoordinate              = FieldVector<ctype, coorddimension>;
        using JacobianTransposed            = FieldMatrix<ctype, mydimension, coorddimension>;
        using Hessian                       = FieldMatrix<ctype, mydimension*(mydimension + 1) / 2, coorddimension>;
        using Jacobian                      = FieldMatrix<ctype, coorddimension, mydimension>;
        using JacobianInverseTransposed     = FieldMatrix<ctype, coorddimension, mydimension>;
        using JacobianInverse               = FieldMatrix<ctype, mydimension, coorddimension>;
        using Volume                        = ctype;

        //! type of the LocalView of the patch geometry
        using GeometryLocalView =
            typename GeometryKernel::NURBSPatch<mydimension, coorddimension,
                                                ctype>::template GeometryLocalView<codim, TrimmerType>;

        /** constructor from host geometry
         */
        DefaultTrimmedPatchLocalGeometry(const TrimDataType& trimData) : trimData_{&trimData} {}

        /** \brief Return the element type identifier
         */
        [[nodiscard]] GeometryType type() const { return trimData_->type(); }

        // return whether we have an affine mapping
        [[nodiscard]] bool affine() const {
          // handy if the trims are straight lines for other algorithms
          if constexpr (codim == 0)
            return true;
          else
            return false;  // TODO for straight lines this should return true
        }

        //! return the number of corners of this element. Corners are numbered 0...n-1
        [[nodiscard]] int corners() const {
          // TODO
          //  return hostGeometry_.corners();
          return {};
        }

        //! access to coordinates of corners. Index is the number of the corner
        GlobalCoordinate corner(int i) const {
          if constexpr (codim == 0)
            return cubeGeometry.corner(i);
          else
            return GlobalCoordinate{};  // get corner at trimdata intersections
        }

        /** \brief Maps a local coordinate within reference element to
         * global coordinate in element  */
        GlobalCoordinate global(const LocalCoordinate& local) const {
          if constexpr (codim == 0) return cubeGeometry.global(local);
        }

        /** \brief Return the transposed of the Jacobian
         */
        JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const {
          if constexpr (codim == 0) return cubeGeometry.jacobianTransposed(local);
        }

        /** \brief Maps a global coordinate within the element to a
         * local coordinate in its reference element */
        LocalCoordinate local(const GlobalCoordinate& global) const {
          if constexpr (codim == 0) return cubeGeometry.local(global);
          ;
        }

        //! Returns true if the point is in the current element
        bool checkInside(const LocalCoordinate& local) const { return trimData_->checkInside(local); }

        [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const {
          if constexpr (codim == 0) return cubeGeometry.integrationElement(local);
        }

        //! The Jacobian matrix of the mapping from the reference element to this element
        [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
          if constexpr (codim == 0) return cubeGeometry.jacobianInverseTransposed(local);
        }

       private:
        const typename ReferenceElements<ctype, mydimension>::ReferenceElement& cubeGeometry{
            ReferenceElements<ctype, mydimension>::cube()};
        const TrimDataType* trimData_{nullptr};
      };
    }  // namespace Trim
  }    // namespace IGANEW
}  // namespace Dune
