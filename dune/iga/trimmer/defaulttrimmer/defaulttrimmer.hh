#pragma once

#include "elementtrimdata.hh"
#include "patchtrimdata.hh"
#include "referenceelement.hh"
#include "trimmedlocalgeometry.hh"

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"

#include <dune/subgrid/subgrid.hh>

namespace Dune {
  namespace IGANEW {

    namespace Trim {

      /**
       * @brief DefaultTrimParameter struct representing parameters for the trimming operation.
       */
      struct DefaultTrimParameter {
        int dummy            = 7;      ///< Dummy variable.
        double trimPrecision = 1e-10;  ///< Precision for trimming.
      };

      /**
       * @brief DefaultTrimmer Trimmer with resonable defaults.
       * @tparam dim Dimension of the patch.
       * @tparam ScalarType Scalar type for geometric calculations.
       */
      template <int dim, typename ScalarType>
      struct DefaultTrimmer {
        static constexpr int mydimension = dim;  ///< Dimension of the patch.

        using ctype = ScalarType;  ///< Scalar type for the coordinates.

        /**
         * @brief Boolean for the linearity of the local geometry.
         * For codim==0, the parameter geometry is linear.
         */
        template <int codim>
        static constexpr bool isLocalGeometryLinear = codim == 0;

        static constexpr bool isAlwaysTrivial = false;  ///< Boolean indicating if the trimming is always trivial, no
                                                        ///< trimming or simple deletion of element.

        using UntrimmedParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

        using ParameterSpaceGrid
            = SubGrid<mydimension,   UntrimmedParameterSpaceGrid>;  ///< Type of the Parametric
                                                                                             ///< grid

        using ReferenceElementType = DefaultTrimmedReferenceElement<mydimension, ctype>;  ///< Reference element type.
        using ParameterType        = DefaultTrimParameter;  ///< Type for trimming parameters.

        /**
         * @brief Get the reference element for a given entity.
         * @tparam EntityType Type of the entity.
         * @param entity The entity for which the reference element is requested.
         * @return Reference element for the entity.
         */
        template <Dune::Concept::EntityExtended EntityType>
        static auto referenceElement(const EntityType& entity) {
          return ReferenceElementType(entity.impl().trimData());
        }

        /**
         * @brief Type alias for local geometry of a specified codimension.
         * @tparam codim Codimension of the local geometry.
         */
        template <int codim>
        using LocalGeometry = typename ReferenceElementType::template Codim<codim>::Geometry;

        using ElementTrimData
            = DefaultElementTrimData<ParameterSpaceGrid::dimension,
                                     typename ParameterSpaceGrid::ctype>;  ///< Element trim data type.
        using PatchTrimData = PatchTrimData<ParameterSpaceGrid::dimension,
                                            typename ParameterSpaceGrid::ctype>;  ///< Patch trim data type.

        using ElementTrimDataContainer = std::map<typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
                                                  ElementTrimData>;  ///< Container for element trim data.

        /**
         * @brief Trim elements based on patch data and trim data.
         * @tparam dimworld Dimension of the world.
         * @param patchData NURBS patch data.
         * @param patchTrimData Patch trim data.
         */
        template <int dimworld>
        auto trimElements(const NURBSPatchData<dim, dimworld, ctype>& patchData, const PatchTrimData& patchTrimData) {
          // fill up container
          // patchTrimData,trimDatas_;
          ;
        }

        ElementTrimDataContainer trimDatas_;  ///< Container for element trim data.
        PatchTrimData patchTrimData;          ///< Patch trim data.
        DefaultTrimParameter parameter;       ///< Trimming parameters.
      };

    }  // namespace Trim
  }    // namespace IGANEW
}  // namespace Dune
