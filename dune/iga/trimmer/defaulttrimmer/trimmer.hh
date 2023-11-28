#pragma once

#include "elementtrimdata.hh"
#include "patchtrimdata.hh"
#include "referenceelement.hh"
#include "trimmedlocalgeometry.hh"

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include "dune/iga/trimmer/localGeometryVariant.hh"

#include <dune/subgrid/subgrid.hh>
#include <dune/subgrid/test/common.hh>

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
        static constexpr int mydimension = dim;         ///< Dimension of the patch.
        using ctype                      = ScalarType;  ///< Scalar type for the coordinates.

        /**
         * @brief Boolean for the linearity of the local geometry.
         * For codim==0, the parameter geometry is linear.
         */
        template <int codim>
        static constexpr bool isLocalGeometryLinear = codim == 0;
        static constexpr bool isAlwaysTrivial = false;  ///< Boolean indicating if the trimming is always trivial, no
        ///< trimming or simple deletion of element.

        DefaultTrimmer() = default;

        template <int dimworld>
        DefaultTrimmer(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patchData,
                       const std::optional<PatchTrimData>& trimData) {
          createParameterSpaceGrid(patchData, trimData);
        }

       private:
        using UntrimmedParameterSpaceGrid = YaspGrid<mydimension, TensorProductCoordinates<ctype, mydimension>>;

       public:
        using ParameterSpaceGrid = SubGrid<mydimension, UntrimmedParameterSpaceGrid>;  ///< Type of the Parametric
                                                                                       ///< grid
        template <int mydim>
        using ReferenceElementType = DefaultTrimmedReferenceElement<mydim, ctype>;  ///< Reference element type.

       private:
        template <int codim>
        using MyLocalGeometry = DefaultTrimmedPatchLocalGeometry<mydimension - codim, mydimension, ctype>;
        template <int codim>
        using UntrimmedLocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;

       public:
        /**
         * @brief Type alias for local geometry of a specified codimension.
         * @tparam codim Codimension of the local geometry.
         */
        template <int codim>
        using LocalGeometry
            = LocalGeometryVariant<DefaultTrimmer, UntrimmedLocalGeometry<codim>, MyLocalGeometry<codim>>;

        using ParameterType = DefaultTrimParameter;  ///< Type for trimming parameters.

        /**
         * @brief Get the reference element for a given entity.
         * @tparam EntityType Type of the entity.
         * @param entity The entity for which the reference element is requested.
         * @return Reference element for the entity.
         */
        template </* Dune::Concept::EntityExtended */ typename EntityType>
        static auto referenceElement(const EntityType& entity) {
          return ReferenceElementType<EntityType::mydimension>(entity.trimData());
        }

        using ElementTrimData
            = DefaultElementTrimData<ParameterSpaceGrid::dimension,
                                     typename ParameterSpaceGrid::ctype>;  ///< Element trim data type.
        using PatchTrimData = PatchTrimData<ParameterSpaceGrid::dimension,
                                            typename ParameterSpaceGrid::ctype>;  ///< Patch trim data type.

        using ElementTrimDataContainer = std::map<typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
                                                  ElementTrimData>;  ///< Container for element trim data.

        template </* Dune::Concept::EntityExtended */ typename EntityType, typename GlobalIdSet>
        std::optional<std::reference_wrapper<const ElementTrimData>> trimData(const EntityType& element,
                                                                              const GlobalIdSet& globalIdSet) const {
          auto iter = trimDatas_.find(globalIdSet.template id<0>(element));
          if (iter != trimDatas_.end())
            return std::make_optional<std::reference_wrapper<const ElementTrimData>>(std::cref(iter->second));

          return std::nullopt;
        }

        /**
         * @brief Refine grid and create new trim
         * @tparam dimworld Dimension of the world.
         * @param patchData NURBS patch data.
         * @param patchTrimData Patch trim data.
         */
        auto globalRefine(int ref) {
          // fill up container
          // patchTrimData,trimDatas_;
          ;
        }

        const ParameterSpaceGrid& parameterSpaceGrid() const { return *parameterSpaceGrid_; }
        ParameterSpaceGrid& parameterSpaceGrid() { return *parameterSpaceGrid_; }
        const UntrimmedParameterSpaceGrid& unTrimmedParameterSpaceGrid() const { return *untrimmedParameterSpaceGrid_; }
        UntrimmedParameterSpaceGrid& unTrimmedParameterSpaceGrid() { return *untrimmedParameterSpaceGrid_; }

       private:
        template <int dimworld>
        void createParameterSpaceGrid(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch,
                                      const std::optional<PatchTrimData>&) {
          untrimmedParameterSpaceGrid_ = std::make_unique<UntrimmedParameterSpaceGrid>(patch.uniqueKnotVector());

          parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);
          parameterSpaceGrid_->createBegin();
          for (auto hostEntity : elements(untrimmedParameterSpaceGrid_->leafGridView())) {
            // if decide which elements are full or trim and add them to the subgrid
            // subGrid.insert(hostEntity);
          }
          parameterSpaceGrid_->insertLeaf();
          parameterSpaceGrid_->createEnd();
        }

        std::unique_ptr<UntrimmedParameterSpaceGrid> untrimmedParameterSpaceGrid_;
        std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;

        ElementTrimDataContainer trimDatas_;  ///< Container for element trim data.
        PatchTrimData patchTrimData;          ///< Patch trim data.
        DefaultTrimParameter parameter;       ///< Trimming parameters.
      };

    }  // namespace Trim
  }    // namespace IGANEW
}  // namespace Dune
