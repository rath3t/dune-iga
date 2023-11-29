// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file trimmer.hh
 * @brief Definition of the default trimmer class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include "trimmedentity.hh"
#include "elementtrimdata.hh"
#include "patchtrimdata.hh"
#include "referenceelement.hh"
#include "trimmedlocalgeometry.hh"

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include "dune/iga/trimmer/localGeometryVariant.hh"
#include "dune/iga/trimmer/entityvariant.hh"

#include <dune/subgrid/subgrid.hh>
#include <dune/subgrid/test/common.hh>

namespace Dune {
  namespace IGANEW {
    namespace DefaultTrim {

      /**
       * @brief DefaultTrimParameter struct representing parameters for the trimming operation.
       */
      struct Parameter {
        int dummy            = 7;      ///< Dummy variable.
        double trimPrecision = 1e-10;  ///< Precision for trimming.
      };

      /**
       * @brief DefaultTrimmer class representing a trimmer with reasonable defaults.
       * @ingroup Trimmer
       * @tparam dim Dimension of the patch.
       * @tparam ScalarType Scalar type for geometric calculations.
       */
      template <int dim, typename ScalarType>
      class Trimmer {
       public:
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

       private:
        using UntrimmedParameterSpaceGrid
            = YaspGrid<mydimension,
                       TensorProductCoordinates<ctype, mydimension>>;  ///< Type of the untrimmed Parametric grid

       public:
        using ParameterSpaceGrid = SubGrid<mydimension, UntrimmedParameterSpaceGrid>;  ///< Type of the Parametric grid
        template <int mydim>
        using ReferenceElementType = TrimmedReferenceElement<mydim, ctype>;  ///< Reference element type.

        using ElementTrimData = ElementTrimData<ParameterSpaceGrid::dimension,
                                                typename ParameterSpaceGrid::ctype>;  ///< Element trim data type.
        using PatchTrimData   = PatchTrimData<ParameterSpaceGrid::dimension,
                                            typename ParameterSpaceGrid::ctype>;  ///< Patch trim data type.

        using ElementTrimDataContainer = std::map<typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
                                                  ElementTrimData>;  ///< Container for element trim data.

        /**
         * @brief Default constructor for Trimmer.
         */
        Trimmer() = default;

        /**
         * @brief Constructor for Trimmer with patch and trim data.
         * @tparam dimworld Dimension of the world.
         * @param patchData NURBS patch data.
         * @param trimData Optional patch trim data.
         */
        template <int dimworld>
        Trimmer(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patchData,
                const std::optional<PatchTrimData>& trimData) {
          createParameterSpaceGrid(patchData, trimData);
        }

       private:
        template <int codim_, int dim_, class GridImp_>
friend class TrimmedParameterSpaceGridEntity;
        template <int codim>
        using TrimmedLocalParameterSpaceGeometry
            = TrimmedLocalGeometry<mydimension - codim, mydimension, ctype, LocalGeometryTag::InParameterSpace>;
        template <int codim>
        using TrimmedLocalGeometry
            = TrimmedLocalGeometry<mydimension - codim, mydimension, ctype, LocalGeometryTag::InReferenceElement>;
        template <int codim>
        using UntrimmedLocalParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
        template <int codim>
        using UntrimmedLocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;
        template <int codim, class GridImp>
        using TrimmedParameterSpaceGridEntity
            = TrimmedParameterSpaceGridEntity<codim, mydimension,GridImp>;
        template <int codim>
        using UntrimmedParameterSpaceGridEntity = typename ParameterSpaceGrid::template Codim<codim>::Entity;


       public:
        template <int codim, class GridImp>
using ParameterSpaceGridEntity=  TrimmedParameterSpaceGridEntity<codim,GridImp>;
        /**
         * @brief Type alias for local geometry of a specified codimension.
         * @tparam codim Codimension of the local geometry.
         */
        template <int codim>
        using LocalParameterSpaceGeometry
            = Trim::LocalGeometryVariant<Trimmer, UntrimmedLocalParameterSpaceGeometry<codim>,
                                         TrimmedLocalParameterSpaceGeometry<codim>>;

        template <int codim>
        using LocalGeometry
            = Trim::LocalGeometryVariant<Trimmer, UntrimmedLocalGeometry<codim>, TrimmedLocalGeometry<codim>>;

        using ParameterType = Parameter;  ///< Type for trimming parameters.

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

        /**
         * @brief Get the trim data for a given element and global ID set.
         * @tparam EntityType Type of the entity.
         * @tparam GlobalIdSet Type of the global ID set.
         * @param element The entity for which trim data is requested.
         * @param globalIdSet The global ID set associated with the element.
         * @return Trim data for the element if available, std::nullopt otherwise.
         */
        template </* Dune::Concept::Entity */ typename EntityType, typename GlobalIdSet>
        std::optional<std::reference_wrapper<const ElementTrimData>> trimData(const EntityType& element,
                                                                              const GlobalIdSet& globalIdSet) const {
          auto iter = trimDatas_.find(globalIdSet.template id<0>(element));
          if (iter != trimDatas_.end())
            return std::make_optional<std::reference_wrapper<const ElementTrimData>>(std::cref(iter->second));

          return std::nullopt;
        }

        /**
         * @brief Refine the grid globally.
         * @param ref Number of refinement levels.
         */
        auto globalRefine(int ref) {
          // fill up container
          // patchTrimData,trimDatas_;
          ;
        }

        /**
         * @brief Get a const reference to the parameter space grid.
         * @return Const reference to the parameter space grid.
         */
        const ParameterSpaceGrid& parameterSpaceGrid() const { return *parameterSpaceGrid_; }

        /**
         * @brief Get a reference to the parameter space grid.
         * @return Reference to the parameter space grid.
         */
        ParameterSpaceGrid& parameterSpaceGrid() { return *parameterSpaceGrid_; }

        /**
         * @brief Get a const reference to the untrimmed parameter space grid.
         * @return Const reference to the untrimmed parameter space grid.
         */
        // const UntrimmedParameterSpaceGrid& unTrimmedParameterSpaceGrid() const { return
        // *untrimmedParameterSpaceGrid_; }

        /**
         * @brief Get a reference to the untrimmed parameter space grid.
         * @return Reference to the untrimmed parameter space grid.
         */
        // UntrimmedParameterSpaceGrid& unTrimmedParameterSpaceGrid() { return *untrimmedParameterSpaceGrid_; }

        /**
         * @brief Pass parameters to the trimmer.
         * @param par The parameters.
         */
        void setup(const ParameterType& par) { parameter = par; }

       private:
        /**
         * @brief Create the parameter space grid based on the patch and trim data.
         * @tparam dimworld Dimension of the world.
         * @param patch NURBS patch.
         * @param patchTrimData Optional patch trim data.
         */
        template <int dimworld>
        void createParameterSpaceGrid(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch,
                                      const std::optional<PatchTrimData>& patchTrimData) {
          untrimmedParameterSpaceGrid_ = std::make_unique<UntrimmedParameterSpaceGrid>(patch.uniqueKnotVector());

          parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);
          parameterSpaceGrid_->createBegin();
          for (auto hostEntity : elements(untrimmedParameterSpaceGrid_->leafGridView())) {
            // if decide which elements are full or trim and add them to the subgrid
             // parameterSpaceGrid_->insert(hostEntity);
          }
          parameterSpaceGrid_->insertLeaf();
          parameterSpaceGrid_->createEnd();
        }

        std::unique_ptr<UntrimmedParameterSpaceGrid>
            untrimmedParameterSpaceGrid_;                         ///< The untrimmed parameter space grid.
        std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;  ///< The trimmed parameter space grid.

        ElementTrimDataContainer trimDatas_;  ///< Container for element trim data.
        PatchTrimData patchTrimData;          ///< Patch trim data.
        ParameterType parameter;              ///< Trimming parameters.
      };

    }  // namespace DefaultTrim
  }    // namespace IGANEW
}  // namespace Dune

