// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>
#include <clipper2/clipper.rectclip.h>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/trimrectangle.hh>

namespace Dune::IGANEW::DefaultTrim {





  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
    const PatchTrimData& trimmingCurves) {
    std::cout<<"START "<<std::endl;
    using namespace Clipper2Lib;
  auto geo  = element.geometry();
  std::array cornerIndices= {0,1,3,2};// see dune book page 127 Figure 5.12
  std::array<FieldVector<ctype,2>,4> corners;// see dune book page 127 Figure 5.12
  corners[0] = geo.corner(0);
  corners[1] = geo.corner(1);
  corners[2] = geo.corner(3);
  corners[3]  = geo.corner(2);

    PathsD elementPath;
    elementPath.push_back({});
    for (int i=0; i<4 ; ++i)
      elementPath[0].push_back({corners[i][0],corners[i][1],i+1});


    PathsD trimmedSampledCurve;

    for (auto trimmingCurve:trimmingCurves.curves()) {
      trimmedSampledCurve.emplace({});
      for (int i = 0;auto v: std::ranges::reverse_view(Utilities::linspace(trimmingCurve.domain()[0],2))) {
        auto fieldVectorPoint = trimmingCurve.global({v});
        trimmedSampledCurve.back().push_back({fieldVectorPoint[0],fieldVectorPoint[1],500+i++});
      }
    }

    trimRectangle(elementPath,trimmedSampledCurve);
  

}
}
