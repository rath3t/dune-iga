
#pragma once

#include "elementtrimdata.hh"
namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      enum class LocalGeometryTag {
        InParameterSpace,
        InReferenceElement
      };

      template <int dim, typename ScalarType = double>
      struct Trimmer;

      template <int mydim_, int coorddim, typename ScalarType,LocalGeometryTag localGeometryTag>
      class TrimmedLocalGeometry {
       public:
        using ctype = ScalarType;

        static constexpr int mydimension = mydim_;
        using TrimmerType                = Trimmer<mydimension, ctype>;
        using TrimDataType               = ElementTrimData<mydimension, ctype>;

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
        explicit TrimmedLocalGeometry(const TrimDataType& trimData) : trimData_{&trimData} {}

        /** @brief Return the element type identifier
         */
        [[nodiscard]] GeometryType type() const { return GeometryTypes::none(mydimension); }

        // return whether we have an affine mapping
        [[nodiscard]] bool affine() const {
          // handy if the trims are straight lines for other algorithms
          if constexpr (codim == 0)
            return true;
          else
            return false;  // TODO for straight lines this should return true
        }

        //! return the number of corners of this element. Corners are numbered 0...n-1
        [[nodiscard]] int corners() const { return 0; }

        GlobalCoordinate center() const { return {}; }

        //! access to coordinates of corners. Index is the number of the corner
        GlobalCoordinate corner(int i) const { return GlobalCoordinate{}; }

        /** @brief Maps a local coordinate within reference element to
         * global coordinate in element  */
        GlobalCoordinate global(const LocalCoordinate& local) const { return GlobalCoordinate{}; }

        /** @brief Return the transposed of the Jacobian
         */
        JacobianTransposed jacobianTransposed(const LocalCoordinate& local) const { return JacobianTransposed{}; }

        /** @brief Maps a global coordinate within the element to a
         * local coordinate in its reference element */
        LocalCoordinate local(const GlobalCoordinate& global) const { return LocalCoordinate{}; }

        //! Returns true if the point is in the current element
        bool checkInside(const LocalCoordinate& local) const { return true; }

        [[nodiscard]] Volume integrationElement(const LocalCoordinate& local) const { return Volume{}; }

        //! The Jacobian matrix of the mapping from the reference element to this element
        [[nodiscard]] JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& local) const {
          return JacobianInverseTransposed{};
        }

       private:
        const TrimDataType* trimData_{nullptr};
      };
    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune
