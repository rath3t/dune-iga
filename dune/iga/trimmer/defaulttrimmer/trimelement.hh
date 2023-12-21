// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>
#include <clipper2/clipper.rectclip.h>

#include <dune/iga/trimmer/defaulttrimmer/elementtrimdata.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>

namespace Dune::IGANEW::DefaultTrim {

  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(
      const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
      const PatchTrimData& patchTrimData) {
    std::cout << "START " << std::endl;
    using namespace Clipper2Lib;
    auto geo = element.geometry();
    std::array<FieldVector<ctype, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    PathD elementPath;
    for (const auto i : std::views::iota(0, 4))
      elementPath.push_back({corners[i][0], corners[i][1], i});

    PathsD clipPaths;
    PathD tempPath;
    constexpr int N = 5;

    for (auto loop : patchTrimData.loops()) {
      tempPath.clear();
      for (int i = 100; auto& curve : loop.curves()) {
        for (auto v : Utilities::linspace(curve.domain()[0], N)) {
          auto fV = curve.global({v});
          tempPath.emplace_back(fV[0], fV[1], i++);
        }
        // Add an additional point just outside of the element
        auto localLastPoint = curve.domain()[0][1];
        auto lastPoint      = curve.global(localLastPoint);

        auto scale = element.geometry().volume() / 10;

        auto dx = curve.jacobian({localLastPoint}) * scale;
        tempPath.emplace_back(lastPoint[0] + dx[0], lastPoint[1] + dx[1], i++);

        i += 99 - N;
      }
      clipPaths.push_back(tempPath);
    }

    auto [flag, result] = Impl::clipElementRectangle(elementPath, clipPaths);

    // Create ElementTrimData with exact intersection Points and a geometry representation of the element edges
    ElementTrimData trimData(flag);

    if (flag != ElementTrimFlag::trimmed) return trimData;

    auto nextEntity  = [&](const int i) { return (i + 1) % result.vertices_.size(); };
    auto isNewVertex = [](const auto& vV) { return std::holds_alternative<Impl::ClippingResult::NewVertex>(vV); };
    auto getPt       = [](auto&& vV) { return vV.pt; };
    auto getHostIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::HostVertex>(vV).hostIdx; };
    auto getEdgeIdx  = [](const auto& vV) { return std::get<Impl::ClippingResult::NewVertex>(vV).pt.z; };

    for (const auto i : std::views::iota(0u, result.vertices_.size())) {
      auto vV1 = result.vertices_[i];
      auto vV2 = result.vertices_[nextEntity(i)];

      // First case edge is completly untrimmed
      if (!isNewVertex(vV1) and !isNewVertex(vV2)) {
        trimData.addEdge(Impl::giveEdgeIdx(getHostIdx(vV1), getHostIdx(vV2)));
        continue;
      }

      auto pt1 = std::visit(getPt, vV1);
      auto pt2 = std::visit(getPt, vV2);

      // Second case edge begins on a hostVertes and ends on a newVertex
      if (!isNewVertex(vV1) and isNewVertex(vV2)) {
        // Find intersection point on edge(vV2)
        auto edge = getEdgeIdx(vV2);
      }

      // Third case edge begins on a newVertex and ends in a HostVertex
      if (isNewVertex(vV1), !isNewVertex(vV2)) {
        // Intersection Point should already be found in case 2
      }
    }

    return trimData;
  }
}  // namespace Dune::IGANEW::DefaultTrim
