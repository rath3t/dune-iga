// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>

namespace Dune::IGANEW::DefaultTrim {
  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
    const PatchTrimData& trimmingCurves) {
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
      elementPath[0].push_back({corners[i][0],corners[i][1],cornerIndices[i]});
    // elementPath[0].push_back({corners[0][0],corners[0][1],cornerIndices[0]});
    ClipperD c(5);
    int intersectionCounter= 4;
    c.SetZCallback([&](const PointD& e1bot, const PointD& e1top,
  const PointD& e2bot, const PointD& e2top, PointD& pt){pt.z=intersectionCounter++;});
    Clipper2Lib::RectD elementRect{corners[0][0], corners[0][1], corners[2][0], corners[2][1]};
    //
    //  Clipper2Lib::Intersect(ClipperPaths{elementEdges}, clip_, Clipper2Lib::FillRule::NonZero);
    // c.AddSubject(elementPath);

    PathsD trimmedSampledCurve;

    for (auto trimmingCurve:trimmingCurves.curves()) {
      trimmedSampledCurve.emplace({});
      for (auto v: Utilities::linspace(trimmingCurve.domain()[0],100)) {
        auto fieldVectorPoint = trimmingCurve.global({v});
        trimmedSampledCurve.back().push_back({fieldVectorPoint[0],fieldVectorPoint[1]});
      }
    }
    PathsD clippedEdges = RectClip(elementRect,trimmedSampledCurve);
    std::cout<<clippedEdges.size()<<std::endl;

}
}
