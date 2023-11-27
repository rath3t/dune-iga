
#pragma once
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
namespace Dune {
  namespace IGANEW {

    namespace Trim {
      template <int mydim_, typename ScalarType>
      struct NoOpElementTrimData {
        // outer loops
        // inner loops ....
      };

      template <typename ParameterSpaceGrid>
      struct NoOpElementTrimDataContainer {};
      template <int mydim_, typename ScalarType>
      struct NoOpPatchTrimData {
        // outer loops
        // inner loops ....
      };

      template <int dim, int dw, typename ScalarType = double>
      struct NoOpTrimmer {
        static constexpr int mydimension = dim;
        static constexpr int dimworld    = dw;
        using ctype                      = ScalarType;

        // using Grid= PatchGrid<dim, dimworld, NoOpTrimmer,ScalarType>;

        // boolean for the linearity of the local geometry, for the untrimmed case this is always true
        template <int codim>
        static constexpr bool isLocalGeometryLinear = true;
        static constexpr bool isAlwaysTrivial       = true;

        using ParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

        using ReferenceElementType = typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement;

        template <Dune::Concept::EntityExtended EntityType>
        static auto referenceElement(const EntityType& entity) {
          return Dune::referenceElement<ctype, mydimension>(entity.type());
        }

        template <int codim>
        using LocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;

        using ElementTrimData = NoOpElementTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

        using PatchTrimData = NoOpPatchTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

        using ElementTrimDataContainer = NoOpElementTrimDataContainer<ParameterSpaceGrid>;
      };
    }  // namespace Trim
  }    // namespace IGANEW
}  // namespace Dune
