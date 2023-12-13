// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>
#include <clipper2/clipper.rectclip.h>

namespace Dune::IGANEW::DefaultTrim {

  template<typename ctype>
struct TrimmingTracker {
  std::map<int,std::tuple<int,int,FieldVector<ctype,2>,bool>> betweenMap;
};


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

    // auto tmp =  elementPath[0];
    // std::ranges::reverse_copy(tmp, elementPath[0].begin());
    // elementPath[0]=tmp;
     // elementPath[0]=std::ranges::reverse( elementPath[0]);
    // elementPath[0].push_back({corners[0][0],corners[0][1],cornerIndices[0]});
    ClipperD c(5);
    int intersectionCounter= 5;

    TrimmingTracker<ctype> tracker;
    c.SetZCallback([&](const PointD& e1bot, const PointD& e1top,
  const PointD& e2bot, const PointD& e2top, PointD& pt) {
    std::cout<<"e1bot: "<<e1bot.z<<std::endl;
    std::cout<<"e1top: "<<e1top.z<<std::endl;
    std::cout<<"e2bot: "<<e2bot.z<<std::endl;
    std::cout<<"e2top: "<<e2top.z<<std::endl;
    std::cout<<"pt: "<<pt<<std::endl;
    std::cout<<"intersectionCounter: "<<intersectionCounter<<std::endl;
    bool in{false};
    auto areWeAtEdge = [&](int i1,int i2){return ((e2top.z== i1 and e2bot.z ==i2) or (e2top.z== i2 and e2bot.z ==i1)); };
    const bool areWeAtTopEdge =areWeAtEdge(3,4);
     const bool areWeAtBottomEdge =areWeAtEdge(1,2);
          const bool areWeAtLeftEdge =areWeAtEdge(4,1);
                 const bool areWeAtRightEdge =areWeAtEdge(2,3);
    if (areWeAtTopEdge  )
     in=(e1bot.z<e1top.z);
    else if (areWeAtLeftEdge)
        in=(e1bot.z>e1top.z);
  else if (areWeAtBottomEdge)
      in=(e1bot.z<e1top.z);
   else if (areWeAtRightEdge)
      in=(e1bot.z<e1top.z);
      tracker.betweenMap.insert({{intersectionCounter},{e2top.z,e2bot.z,FieldVector<ctype,2>({pt.x,pt.y}),in}});
      pt.z=intersectionCounter++;
    });
    const auto offsetX= corners[0][0];
    const auto offsetY= corners[0][1];
    const auto right= corners[2][0];
    const auto bottom= corners[2][1];
    Clipper2Lib::RectD elementRect{offsetX, offsetY, right, bottom};
    //
   //   c.AddSubject(elementPath);
std::cout<<"Lower left corner: "<<corners[0]<<std::endl;
    PathsD trimmedSampledCurve;

    for (auto trimmingCurve:trimmingCurves.curves()) {
      trimmedSampledCurve.emplace({});
      for (int i = 0;auto v: Utilities::linspace(trimmingCurve.domain()[0],2)) {
        auto fieldVectorPoint = trimmingCurve.global({v});
        trimmedSampledCurve.back().push_back({fieldVectorPoint[0],fieldVectorPoint[1],500+i++});
      }
    }
    std::cout<<"Trimming Curve"<<std::endl;

    for( auto point: trimmedSampledCurve[0])
      std::cout<<point<<std::endl;
    std::cout<<"Trimming Curve End"<<std::endl;

    c.AddOpenSubject(trimmedSampledCurve );
    c.AddClip(elementPath );
    PathsD clippedOpenEdges;
    PathsD clippedClosedEdges;
    PolyTreeD tree;
    //  c.Execute(ClipType::Intersection,FillRule::NonZero,tree);
    //  std::cout<<"Tree size"<<tree<<std::endl;
    //       c.Execute(ClipType::Intersection,FillRule::EvenOdd,tree);
    //  std::cout<<"Tree size"<<tree<<std::endl;

    if(not c.Execute(ClipType::Intersection,FillRule::EvenOdd,clippedClosedEdges,clippedOpenEdges))
      DUNE_THROW(InvalidStateException,"Trimming failed of element with lower left corner "<<corners[0]);
    // = Clipper2Lib::Intersect(elementPath, trimmedSampledCurve, Clipper2Lib::FillRule::NonZero);


    std::cout<<"Result Open"<<std::endl;
    for( auto& path: clippedOpenEdges)
      for( auto point: path)
        std::cout<<point<<std::endl;

    std::cout<<"Result Closed"<<std::endl;
    for( auto& path: clippedClosedEdges)
      for( auto point: path)
        std::cout<<point<<std::endl;

    std::cout<<"Tracked intersections: "<<std::endl;
    for(const auto& [key,value]: tracker.betweenMap)
      std::cout<<key<<" v:"<<std::get<0>(value)<<" " <<std::get<1>(value)<<" "<<std::get<2>(value)<<" "<<(std::get<3>(value)? " in": " out")<<std::endl;
     // PathsD clippedEdges =Clipper2Lib::Intersect(elementPath, trimmedSampledCurve, Clipper2Lib::FillRule::NonZero);
    // c.AddClip(trimmedSampledCurve);
    // c.Execute(ClipType::Intersection,FillRule::NonZero);
    // PathsD clippedEdges = RectClipLines(elementRect,trimmedSampledCurve,5);
    // std::cout<<clippedEdges.size()<<std::endl;
    std::cout<<"END "<<std::endl;

    //Compute final path
    std::vector<std::pair<int,int>> edges;
    edges.push_back({1,2});
    edges.push_back({2,3});
    edges.push_back({3,4});
    edges.push_back({4,1});



    std::cout<<"Final path "<<std::endl;

    // for(const auto& [key,value]: tracker.betweenMap) {
    //   trimmedPath[std::get<1>(value)]=std::get<2>(value);
    // }

    // for( auto point: trimmedPath)
    //   std::cout<<point<<std::endl;

}
}
