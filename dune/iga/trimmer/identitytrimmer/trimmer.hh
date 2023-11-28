
#pragma once
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"

namespace Dune::IGANEW {
  namespace GeometryKernel {
    template <int dim_, int dimworld_, typename ScalarType>
    class NURBSPatch;
  }
  namespace Trim {

    struct IdentityTrimParameter {};

    template <int mydim_, typename ScalarType>
    struct IdentityElementTrimData {
      // outer loops
      // inner loops ....
    };

    template <typename ParameterSpaceGrid>
    struct IdentityElementTrimDataContainer {};
    template <int mydim_, typename ScalarType>
    struct IdentityPatchTrimData {
      // outer loops
      // inner loops ....
    };

    template <int dim, typename ScalarType = double>
    struct IdentityTrimmer {
      static constexpr int mydimension = dim;         ///< Dimension of the patch.
      using ctype                      = ScalarType;  ///< Scalar type for the coordinates.

      template <int codim>
      static constexpr bool isLocalGeometryLinear
          = true;  ///< boolean for the linearity of the local geometry, for the untrimmed case this is always true
      static constexpr bool isAlwaysTrivial = true;  ///< Boolean indicating if the trimming is always trivial, no
      ///< trimming or simple deletion of element.

      IdentityTrimmer() = default;

      // using Grid= PatchGrid<dim, dimworld, IdentityTrimmer,ScalarType>;

      using ParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

      using UntrimmedParameterSpaceGrid = Empty;

      using ReferenceElementType = typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement;

      template <int codim>
      using LocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;

      using ParameterType = IdentityTrimParameter;  ///< Type for trimming parameters.

      template </* Dune::Concept::EntityExtended */ typename EntityType>
      static auto referenceElement(const EntityType& entity) {
        return Dune::referenceElement<ctype, EntityType::mydimension>(entity.type());
      }

      using ElementTrimData
          = IdentityElementTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

      using PatchTrimData = IdentityPatchTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

      using ElementTrimDataContainer = IdentityElementTrimDataContainer<ParameterSpaceGrid>;

      template <int dimworld>
      void createParameterSpaceGrid(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch,
                                    const std::optional<PatchTrimData>&) {
        parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(patch.uniqueKnotVector());
      }

      template <int dimworld>
      IdentityTrimmer(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch,
                      const std::optional<PatchTrimData>& trimData) {
        createParameterSpaceGrid(patch, trimData);
      }

      const ParameterSpaceGrid& parameterSpaceGrid() const { return *parameterSpaceGrid_; }
      ParameterSpaceGrid& parameterSpaceGrid() { return *parameterSpaceGrid_; }

      std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;
    };
  }  // namespace Trim
}  // namespace Dune::IGANEW
