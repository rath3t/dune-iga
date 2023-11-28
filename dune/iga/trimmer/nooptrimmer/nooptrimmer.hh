
#pragma once
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"

  namespace Dune::IGANEW {
    namespace GeometryKernel {
      template <int dim_, int dimworld_, typename ScalarType>
  class NURBSPatch ;
    }
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

      template <int dim, typename ScalarType = double>
      struct NoOpTrimmer {
        static constexpr int mydimension = dim;
        using ctype                      = ScalarType;
        NoOpTrimmer()=default;


        // using Grid= PatchGrid<dim, dimworld, NoOpTrimmer,ScalarType>;

        // boolean for the linearity of the local geometry, for the untrimmed case this is always true
        template <int codim>
        static constexpr bool isLocalGeometryLinear = true;
        static constexpr bool isAlwaysTrivial       = true;

        using ParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

        using UntrimmedParameterSpaceGrid = Empty;

        using ReferenceElementType = typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement;

        template <Dune::Concept::EntityExtended EntityType>
        static auto referenceElement(const EntityType& entity) {
          return Dune::referenceElement<ctype, mydimension>(entity.type());
        }

        template <int codim>
        using LocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;

        // template <int codim>
        // using LocalHostGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;


        using ElementTrimData = NoOpElementTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

        using PatchTrimData = NoOpPatchTrimData<ParameterSpaceGrid::dimension, typename ParameterSpaceGrid::ctype>;

        using ElementTrimDataContainer = NoOpElementTrimDataContainer<ParameterSpaceGrid>;

        template<int dimworld>
        void createParameterSpaceGrid(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch, const std::optional<PatchTrimData>& ) {
          parameterSpaceGrid_= std::make_unique<ParameterSpaceGrid>(patch.uniqueKnotVector());
        }

        template<int dimworld>
NoOpTrimmer(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch, const std::optional<PatchTrimData>& trimData) {
          createParameterSpaceGrid(patch,trimData);
        }

        const ParameterSpaceGrid& parameterSpaceGrid()const {return *parameterSpaceGrid_;}
         ParameterSpaceGrid& parameterSpaceGrid() {return *parameterSpaceGrid_;}

        std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;

      };
    }  // namespace Trim
  } // namespace Dune::IGANEW

