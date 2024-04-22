// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>

namespace Dune::IGANEW::DefaultTrim {

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::refineParameterSpaceGrid(int refCount, bool initFlag) {
  const int oldLevel = untrimmedParameterSpaceGrid_->maxLevel();
  untrimmedParameterSpaceGrid_->globalRefine(refCount);
  // auto gvu = untrimmedParameterSpaceGrid_->leafGridView();
  parameterSpaceGrid_->createBegin();
  parameterSpaceGrid_->insertLeaf();
  parameterSpaceGrid_->createEnd();

  assert((initFlag and oldLevel == 0 and refCount == 0) or
         !initFlag && "If we initialize the grid, the untrimmedParameterSpaceGrid_ should only have one level");
  // if we init we start at 0 otherwise at 1
  // if the grid is refined once later this would yield int i = 1; i < 2;
  for (int i = !initFlag; i < refCount + 1; ++i) {
    const int newLevel = oldLevel + i;
    entityContainer_.entityImps_.emplace_back();

    auto gv                            = parameterSpaceGrid_->levelGridView(newLevel);
    auto& globalIdSetParameterSpace    = parameterSpaceGrid_->globalIdSet();
    unsigned int unTrimmedElementIndex = 0;
    unsigned int trimmedElementIndex   = 0;

    auto indices = [&]() { return std::make_tuple(unTrimmedElementIndex, trimmedElementIndex, newLevel); };

    auto& indexSet              = gv.indexSet();
    const auto elementTrimDatas = trimElements(newLevel);

    entityContainer_.idToVertexInfoMap.emplace_back();
    entityContainer_.trimmedVertexIds_.emplace_back();
    entityContainer_.edgeCount.emplace_back(gv.size(1));
    entityContainer_.vertexCount.emplace_back(gv.size(2));

    for (const auto& ele : elements(gv)) {
      const ElementTrimData& eleTrimData = elementTrimDatas[indexSet.index(ele)];
      const ElementTrimFlag eleTrimFlag  = eleTrimData.flag();

      // For now we are exiting if ele is empty
      if (eleTrimFlag == ElementTrimFlag::empty)
        continue;

      // The elementID stays persistent, it doenst matter if ele is trimmed or not (for now)
      GlobalIdType elementId = makeElementID(ele);

      // Father detection is done in createAndSaveElementInfo
      if (eleTrimFlag == ElementTrimFlag::full) {
        createAndSaveElementInfo(indices(), ele, false);
        ++unTrimmedElementIndex;
      } else if (eleTrimFlag == ElementTrimFlag::trimmed) {
        createAndSaveElementInfo(indices(), ele, true);
        ++trimmedElementIndex;
      }

      // ******************* Now sub entitities, collect all ids of subentities  ************************

      collectElementEdges(newLevel, ele, eleTrimData);
      collectElementVertices(newLevel, ele, eleTrimData);
    }

    // save numbers of untrimmed and trimmed elements per level
    entityContainer_.numberOfTrimmedElements.push_back(trimmedElementIndex);
    entityContainer_.numberOfUnTrimmedElements.push_back(unTrimmedElementIndex);

    createElements(newLevel, elementTrimDatas);
    createSubEntities(newLevel);
  }
}
} // namespace Dune::IGANEW::DefaultTrim
