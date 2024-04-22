// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file trimmer.hh
 * @brief Definition of the default trimmer class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include "elementtrimdata.hh"
#include "entitycontainer.hh"
#include "idset.hh"
#include "patchtrimdata.hh"
#include "referenceelement.hh"
#include "trimmedentity.hh"
#include "trimmedlocalgeometry.hh"

#include <clipper2/clipper.h>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include "dune/iga/trimmer/entityvariant.hh"
#include "dune/iga/trimmer/localGeometryVariant.hh"
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/subgrid/subgrid.hh>
#include <dune/subgrid/test/common.hh>

namespace Dune {
namespace IGANEW {
  namespace DefaultTrim {
    // template<typename  ParameterSpaceGrid,typename ElementTrimData>
    // struct ElementTrimDataContainer {
    //
    //   // trimData()
    //
    //   using ElementTrimDataContainer = std::map<typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
    //                                     ElementTrimData>;  ///< Container for element trim data.
    // };

    /**
     * @brief DefaultTrimParameter struct representing parameters for the trimming operation.
     */
    struct Parameter
    {
      int dummy            = 7;     ///< Dummy variable.
      double trimPrecision = 1e-10; ///< Precision for trimming.
    };

    template <typename HostIdType>
    struct IdType
    {
      enum class HostOrTrimmed
      {
        host,
        trimmed
      };
      HostOrTrimmed host_or_trimmed;
      HostIdType id;

      friend std::ostream& operator<<(std::ostream& stream, const IdType& id) {
        stream << "Type: " << std::string(id.host_or_trimmed == HostOrTrimmed::host ? "Host" : "Trimmed")
               << ", Key: " << id.id << "\n";
        return stream;
      }
    };

    template <typename HostIdType>
    bool operator<(const IdType<HostIdType>& lhs, const IdType<HostIdType>& rhs) {
      if (lhs.host_or_trimmed == rhs.host_or_trimmed)
        return lhs.id < rhs.id;
      else if (lhs.host_or_trimmed < rhs.host_or_trimmed)
        return true;
      else
        return false;
    }

    /**
     * @brief DefaultTrimmer class representing a trimmer with reasonable defaults.
     * @ingroup Trimmer
     * @tparam dim Dimension of the patch.
     * @tparam ScalarType Scalar type for geometric calculations.
     */
    template <int dim, int dimworld, typename ScalarType>
    class Trimmer
    {
    public:
      static constexpr int mydimension = dim;        ///< Dimension of the patch.
      using ctype                      = ScalarType; ///< Scalar type for the coordinates.

      /**
       * @brief Boolean for the linearity of the local geometry.
       * For codim==0, the parameter geometry is linear.
       */
      template <int codim>
      static constexpr bool isLocalGeometryLinear = codim == 0;
      static constexpr bool isAlwaysTrivial       = false; ///< Boolean indicating if the trimming is always trivial, no
                                                           ///< trimming or simple deletion of element.

    private:
      using UntrimmedParameterSpaceGrid =
          YaspGrid<mydimension,
                   TensorProductCoordinates<ctype, mydimension>>; ///< Type of the untrimmed Parametric grid

    public:
      using TrimmingCurve      = GeometryKernel::NURBSPatch<mydimension - 1, mydimension, ctype>;
      using ParameterSpaceGrid = SubGrid<mydimension, UntrimmedParameterSpaceGrid>; ///< Type of the Parametric grid

      template <int mydim>
      using ReferenceElementType = TrimmedReferenceElement<mydim, Trimmer>; ///< Reference element type.

      using ElementTrimData = ElementTrimData<ParameterSpaceGrid::dimension,
                                              typename ParameterSpaceGrid::ctype>; ///< Element trim data type.
      using PatchTrimData   = PatchTrimData<Trimmer>;                              ///< Patch trim data type.

      using ElementTrimDataContainer = std::map<typename ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
                                                ElementTrimData>; ///< Container for element trim data.
      template <typename GridImp>
      using EntityContainerType = EntityContainer<GridImp>;

      /**
       * @brief Default constructor for Trimmer.
       */
      Trimmer() = default;

      /**
       * @brief Constructor for Trimmer with patch and trim data.
       * @param patchData NURBS patch data.
       * @param trimData Optional patch trim data.
       */
      Trimmer(const GeometryKernel::NURBSPatch<dim, dimworld, ctype>& patch,
              const std::optional<PatchTrimData>& trimData)
          : patch_{patch},
            patchTrimData_{trimData} {
      }

    private:
      template <int codim_, int dim_, class GridImp_>
      friend class TrimmedParameterSpaceGridEntity;
      template <class GridImp_>
      friend class EntityContainer;
      template <int codim>
      using TrimmedLocalParameterSpaceGeometry =
          TrimmedLocalGeometry<mydimension - codim, Trimmer, LocalGeometryTag::InParameterSpace>;
      template <int codim>
      using TrimmedLocalGeometry =
          TrimmedLocalGeometry<mydimension - codim, Trimmer, LocalGeometryTag::InReferenceElement>;
      template <int codim>
      using UntrimmedLocalParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;
      template <int codim>
      using UntrimmedLocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;
      template <int codim, class GridImp>
      using TrimmedParameterSpaceGridEntity = TrimmedParameterSpaceGridEntity<codim, mydimension, GridImp>;
      template <int codim>
      // using UntrimmedParameterSpaceGridEntity = typename UntrimmedParameterSpaceGrid::template
      // Codim<codim>::Entity;
      using HostParameterSpaceGridEntity = typename ParameterSpaceGrid::template Codim<codim>::Entity;

    public:
      using HostIdType = typename UntrimmedParameterSpaceGrid::Traits::GlobalIdSet::IdType;

      template <class GridImp>
      using GlobalIdSet = DefaultTrim::PatchGridGlobalIdSet<GridImp>;

      template <class GridImp>
      friend class DefaultTrim::PatchGridGlobalIdSet;
      template <int codim, class GridImp>
      using ParameterSpaceGridEntity = TrimmedParameterSpaceGridEntity<codim, GridImp>;
      /**
       * @brief Type alias for local geometry of a specified codimension.
       * @tparam codim Codimension of the local geometry.
       */
      template <int codim>
      using LocalParameterSpaceGeometry =
          Trim::LocalGeometryVariant<Trimmer, UntrimmedLocalParameterSpaceGeometry<codim>,
                                     TrimmedLocalParameterSpaceGeometry<codim>>;

      template <int codim>
      using LocalGeometry =
          Trim::LocalGeometryVariant<Trimmer, UntrimmedLocalGeometry<codim>, TrimmedLocalGeometry<codim>>;

      using ParameterType = Parameter; ///< Type for trimming parameters.

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
      const ParameterSpaceGrid& parameterSpaceGrid() const {
        return *parameterSpaceGrid_;
      }

      /**
       * @brief Get a reference to the parameter space grid.
       * @return Reference to the parameter space grid.
       */
      ParameterSpaceGrid& parameterSpaceGrid() {
        return *parameterSpaceGrid_;
      }

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
      void setup(const ParameterType& par) {
        parameter = par;
      }

      /**
       * @brief Create the parameter space grid based on the patch and trim data.
       * @tparam dimworld Dimension of the world.
       * @param patch NURBS patch.
       * @param patchTrimData Optional patch trim data.
       */
      template <typename GridImp>
      void createIdSetAndParameterGrid(GridImp& grid) {
        untrimmedParameterSpaceGrid_ = std::make_unique<UntrimmedParameterSpaceGrid>(patch_.uniqueKnotVector());

        if (not patchTrimData_) {
          parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);
          createSubGridLevel(grid);
          createLevel(grid);
          return;
        } else {
          parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(*untrimmedParameterSpaceGrid_);
          createSubGridLevel(grid);
          createLevel(grid);
        }
      }

      template <typename GridImp>

      void createSubGridLevel(GridImp& grid) {
        if (not patchTrimData_) {
          parameterSpaceGrid_->createBegin();
          parameterSpaceGrid_->insertLeaf();
          parameterSpaceGrid_->createEnd();
        } else {
          parameterSpaceGrid_->createBegin();
          parameterSpaceGrid_->insertLeaf();
          parameterSpaceGrid_->createEnd();
        }
      }

      template <typename GridImp>
      void createLevel(GridImp& grid) {
        using IdType                 = typename GridImp::Traits::GlobalIdSet::IdType;
        using EdgeHostType           = typename UntrimmedParameterSpaceGrid::template Codim<1>::Entity;
        using EdgeGridType           = typename GridImp::template Codim<1>::Entity;
        using EleGridType            = typename GridImp::template Codim<0>::Entity;
        using EdgeParameterSpaceType = TrimmedParameterSpaceGridEntity<1, const GridImp>;
        using EleParameterSpaceType  = TrimmedParameterSpaceGridEntity<0, const GridImp>;

        auto& entityContainer           = grid.entityContainer_;
        auto& globalIdSet               = grid.globalIdSet_;
        auto gv                         = parameterSpaceGrid_->leafGridView();
        auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

        for (const auto& ele : elements(gv)) {
          // auto trimedEdges= trimCubeReferenceElement(ele,clipCurve);
          // std::cout<<"T"<<trimedEdges.size()<<std::endl;
          // for(auto v: trimedEdges)
          //   std::cout<<v<<std::endl;
          std::cout << "Element" << std::endl;
          // auto elementId = globalIdSet->getStableIndex();
          IdType elementId = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = globalIdSetParameterSpace.id(ele)};

          auto& elementEdgeIndices   = entityContainer.globalEdgesIdOfElementsMap_[elementId];
          auto& elementVertexIndices = entityContainer.globalVerticesIdOfElementsMap[elementId];
          std::set<IdType> elementVertexSet;
          elementVertexSet.clear();
          for (const auto& intersection : intersections(gv, ele)) {
            auto localEdgeIndex = intersection.indexInInside();
            std::cout << "Edge " << localEdgeIndex << std::endl;

            auto edge    = intersection.inside().template subEntity<1>(localEdgeIndex);
            auto vertex0 = intersection.inside().template subEntity<2>(0);
            auto vertex1 = intersection.inside().template subEntity<2>(1);
            auto refEle  = Dune::ReferenceElements<ctype, mydimension>::cube();

            if (true /* untrimmed */) {
              auto eleParameterSpace = EleParameterSpaceType(&grid, ele, elementId);
              // auto realEntity = typename EleGridType::Implementation(&grid,eleParameterSpace);
              entityContainer.globalElementMap.insert({elementId, eleParameterSpace});
              auto edgeId = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
              // IndexType id =globalIdSet->g(edgeId);
              IdType id = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = edgeId};
              elementEdgeIndices.push_back(id);
              auto localVertexid0 = refEle.subEntity(localEdgeIndex, 1, 0, 2);
              auto localVertexid1 = refEle.subEntity(localEdgeIndex, 1, 1, 2);
              auto vertexId0      = globalIdSetParameterSpace.subId(ele, localVertexid0, 2);
              auto vertexId1      = globalIdSetParameterSpace.subId(ele, localVertexid1, 2);

              IdType idVert0 = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = vertexId0};
              IdType idVert1 = {.host_or_trimmed = IdType::HostOrTrimmed::host, .id = vertexId1};

              elementVertexSet.insert(idVert0);
              elementVertexSet.insert(idVert1);
              // globalEdgeMap.insert({id,Edge(edge)});
              auto paraEnt = EdgeParameterSpaceType(&grid, edge, id);
              // auto realEdge = typename EdgeGridType::Implementation(&grid,paraEnt);
              entityContainer.globalEdgeMap.insert({id, paraEnt});
            } else /* trimmed */
            {
              if (true /* edge is trimmed but has a rest on the original intersection*/)
              // if edge is trimmed but contains a part of the original intersection,
              //  then this part of the edge that belongs to the original intersection will be inserted using the id
              //  of the host the rest of the edge will be created using a new index
              {
                // auto edgeId = globalIdSet.subId(ele,localEdgeIndex,1);
                // IndexType id =getStableIndex(edgeId);
                // elementEdgeIndices.push_back(id);
                // globalEdgeMap.insert({id,"Part of trimmed edge is part of untrimmed edge" +
                // std::to_string(id.id.touint())}); IndexType id2 = getStableIndex(myTrimmedEdgeIndex);
                // globalEdgeMap.insert({id2,"Part of Trimmed edge with ID, that is not completely new but this part
                // is"+std::to_string(myTrimmedEdgeIndex.id.touint())}); elementEdgeIndices.push_back(id);
                //
                // ++myTrimmedEdgeIndex.id;

              } else if (false /* the edge is completly gone?*/) {
                //
              }
              //
              // IndexType id =getStableIndex(myTrimmedEdgeIndex);
              // elementEdgeIndices.push_back(id);
              // //  globalEdgeMap.insert({id,Edge(edge)});
              //
              // globalEdgeMap.insert({getStableIndex(myTrimmedEdgeIndex),"Trimmed edge with ID
              // "+std::to_string(myTrimmedEdgeIndex.id.touint())});
              // ++myTrimmedEdgeIndex.id;
            }
            // auto edgeGeo = intersection.geometry();
            // auto firstCorner = edgeGeo.corner(0);
            // auto secondCorner = edgeGeo.corner(1);
            // edges[0][indexMapper(i)](edgeGeo.corner(0), pos0[1]);
          }
          elementVertexIndices.resize(4);
          std::ranges::copy(elementVertexSet, elementVertexIndices.begin());
        }

        for (const auto& [key, value] : entityContainer.globalEdgeMap)
          if (key.host_or_trimmed == IdType::HostOrTrimmed::host)
            std::cout << "Host: Index value: " << key.id << std::endl;
          else
            std::cout << "Trimmed: Index value: " << key.id << std::endl;

        for (const auto& [key, value_vec] : entityContainer.globalEdgesIdOfElementsMap_) {
          std::cout << "Host: Index value: " << key.id << std::endl;
          for (auto value : value_vec)
            std::cout << " Edge:  id: " << value.id << std::endl;
        }

        for (const auto& [key, value_vec] : entityContainer.globalVerticesIdOfElementsMap) {
          std::cout << "Host: Index value: " << key.id << std::endl;
          for (auto value : value_vec)
            std::cout << " Vertex: " << value.id << std::endl;
        }
      }

      template <typename LevelIndexSets, typename LeafIndexSet>
      void updateIndices(LevelIndexSets& lvlIndexSets, LeafIndexSet& leafIndexSet) {
      }

    private:
      std::unique_ptr<UntrimmedParameterSpaceGrid>
          untrimmedParameterSpaceGrid_;                        ///< The untrimmed parameter space grid.
      std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_; ///< The trimmed parameter space grid.

      ElementTrimDataContainer trimDatas_; ///< Container for element trim data.
      GeometryKernel::NURBSPatch<dim, dimworld, ctype> patch_;
      std::optional<PatchTrimData> patchTrimData_; ///< Patch trim data.
      ParameterType parameter;                     ///< Trimming parameters.
    };

  } // namespace DefaultTrim
} // namespace IGANEW
} // namespace Dune
