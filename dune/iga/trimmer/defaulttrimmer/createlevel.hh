// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
namespace Dune::IGANEW::DefaultTrim {
  template <int dim, int dimworld, typename ScalarType>
  void Trimmer<dim,dimworld,ScalarType>::createLevel( GridImp& grid, int lvl) {
    using IdType = typename GridImp::Traits::GlobalIdSet::IdType;
    using EdgeHostType = typename UntrimmedParameterSpaceGrid::template Codim<1>::Entity;
    using EdgeGridType = typename GridImp::template Codim<1>::Entity;
    using EleGridType = typename GridImp::template Codim<0>::Entity;
    using EdgeParameterSpaceType = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
    using EleParameterSpaceType = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;

    entityContainer_.push_back();

    auto& entityContainer = entityContainer_;
    auto& globalIdSet=grid.globalIdSet_;
    auto gv= parameterSpaceGrid_->levelGridView(lvl);
    auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

    for (const auto& ele: elements(gv)) {
      // auto trimedEdges= trimCubeReferenceElement(ele,clipCurve);
      // std::cout<<"T"<<trimedEdges.size()<<std::endl;
      // for(auto v: trimedEdges)
      //   std::cout<<v<<std::endl;
      std::cout<<"Element"<<std::endl;
      // auto elementId = globalIdSet->getStableIndex();
      IdType elementId = {.host_or_trimmed=IdType::HostOrTrimmed::host,.id=globalIdSetParameterSpace.id(ele)};

      auto& elementEdgeIndices= entityContainer.globalEdgesIdOfElementsMap_[elementId];
      auto& elementVertexIndices= entityContainer.globalVerticesIdOfElementsMap[elementId];
      std::set<IdType> elementVertexSet;
      elementVertexSet.clear();
      for(const auto& intersection: intersections(gv,ele)){
        auto localEdgeIndex= intersection.indexInInside();
        std::cout<<"Edge "<<localEdgeIndex<<std::endl;

        auto edge= intersection.inside().template subEntity<1>(localEdgeIndex);
        auto vertex0= intersection.inside().template subEntity<2>(0);
        auto vertex1= intersection.inside().template subEntity<2>(1);
        auto refEle =Dune::ReferenceElements<ctype,mydimension>::cube();


        if (true /* untrimmed */) {
          auto eleParameterSpace = EleParameterSpaceType(&grid,ele,elementId);
          // auto realEntity = typename EleGridType::Implementation(&grid,eleParameterSpace);
          entityContainer.globalElementMap.insert({elementId,eleParameterSpace});
          auto edgeId = globalIdSetParameterSpace.subId(ele,localEdgeIndex,1);
          // IndexType id =globalIdSet->g(edgeId);
          IdType id = {.host_or_trimmed=IdType::HostOrTrimmed::host,.id=edgeId};
          elementEdgeIndices.push_back(id);
          auto localVertexid0 =      refEle.subEntity(localEdgeIndex,1,0,2);
          auto localVertexid1 =      refEle.subEntity(localEdgeIndex,1,1,2);
          auto vertexId0 = globalIdSetParameterSpace.subId(ele,localVertexid0,2);
          auto vertexId1 = globalIdSetParameterSpace.subId(ele,localVertexid1,2);

          IdType idVert0 = {.host_or_trimmed=IdType::HostOrTrimmed::host,.id=vertexId0};
          IdType idVert1 = {.host_or_trimmed=IdType::HostOrTrimmed::host,.id=vertexId1};

          elementVertexSet.insert(idVert0);
          elementVertexSet.insert(idVert1);
          // globalEdgeMap.insert({id,Edge(edge)});
          auto paraEnt= EdgeParameterSpaceType(&grid,edge,id);
          // auto realEdge = typename EdgeGridType::Implementation(&grid,paraEnt);
          entityContainer.globalEdgeMap.insert({id,paraEnt});
        }
        else /* trimmed */
        {
          if(true/* edge is trimmed but has a rest on the original intersection*/)
            //if edge is trimmed but contains a part of the original intersection,
              // then this part of the edge that belongs to the original intersection will be inserted using the id of the host
                // the rest of the edge will be created using a new index
          {
            // auto edgeId = globalIdSet.subId(ele,localEdgeIndex,1);
            // IndexType id =getStableIndex(edgeId);
            // elementEdgeIndices.push_back(id);
            // globalEdgeMap.insert({id,"Part of trimmed edge is part of untrimmed edge" + std::to_string(id.id.touint())});
            // IndexType id2 = getStableIndex(myTrimmedEdgeIndex);
            // globalEdgeMap.insert({id2,"Part of Trimmed edge with ID, that is not completely new but this part is"+std::to_string(myTrimmedEdgeIndex.id.touint())});
            // elementEdgeIndices.push_back(id);
            //
            // ++myTrimmedEdgeIndex.id;

          }
          else if (false /* the edge is completly gone?*/) {
            //
          }
          //
          // IndexType id =getStableIndex(myTrimmedEdgeIndex);
          // elementEdgeIndices.push_back(id);
          // //  globalEdgeMap.insert({id,Edge(edge)});
          //
          // globalEdgeMap.insert({getStableIndex(myTrimmedEdgeIndex),"Trimmed edge with ID "+std::to_string(myTrimmedEdgeIndex.id.touint())});
          // ++myTrimmedEdgeIndex.id;
        }
        //  auto edgeGeo = intersection.geometry();
        //  auto firstCorner = edgeGeo.corner(0);
        // auto secondCorner = edgeGeo.corner(1);
        // edges[0][indexMapper(i)](edgeGeo.corner(0), pos0[1]);

      }
      elementVertexIndices.resize(4);
      std::ranges::copy(elementVertexSet,elementVertexIndices.begin());

    }


    for (const auto& [key, value] : entityContainer.globalEdgeMap)
      if(key.host_or_trimmed==IdType::HostOrTrimmed::host)
        std::cout << "Host: Index value: " << key.id<< std::endl;
      else
        std::cout << "Trimmed: Index value: " << key.id<< std::endl;

    for (const auto& [key, value_vec] : entityContainer.globalEdgesIdOfElementsMap_) {
      std::cout << "Host: Index value: " << key.id << std::endl;
      for(auto value:value_vec)
        std::cout << " Edge:  id: "<<value.id << std::endl;
    }

    for (const auto& [key, value_vec] : entityContainer.globalVerticesIdOfElementsMap) {
      std::cout << "Host: Index value: " << key.id << std::endl;
      for(auto value:value_vec)
        std::cout << " Vertex: "<<value.id<< std::endl;
    }
  }

  /**
   * \brief Create the paramter grid levels
   * \param grid
   */
  template <int dim, int dimworld, typename ScalarType>
  void Trimmer<dim,dimworld,ScalarType>::refineParameterSpaceGrid(int refCount) {
          using IdType = typename GridImp::Traits::GlobalIdSet::IdType;
          using EdgeHostType = typename UntrimmedParameterSpaceGrid::template Codim<1>::Entity;
          using EdgeGridType = typename GridImp::template Codim<1>::Entity;
          using EleGridType = typename GridImp::template Codim<0>::Entity;
          using EdgeParameterSpaceType = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
          using EleParameterSpaceType = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;

          entityContainer_.entityImps_.emplace_back();
            auto& entityContainer = entityContainer_;
            // auto& globalIdSet=grid.globalIdSet_;
            untrimmedParameterSpaceGrid_->globalRefine(refCount);
            auto gv= untrimmedParameterSpaceGrid_->leafGridView();
            auto& globalIdSetParameterSpace = untrimmedParameterSpaceGrid_->globalIdSet();
                parameterSpaceGrid_->createBegin();
            for (const auto& ele: elements(gv)) {

              parameterSpaceGrid_->insert(ele);
            }
    parameterSpaceGrid_->createEnd();
        }
}