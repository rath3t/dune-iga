// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <clipper2/clipper.h>
#include <clipper2/clipper.rectclip.h>

#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/clipelementrectangle.hh>

namespace Dune::IGANEW::DefaultTrim {


  template <int dim, int dimworld, typename ScalarType>
  auto TrimmerImpl<dim, dimworld, ScalarType>::trimElement(
      const typename GridFamily::TrimmerTraits::template Codim<0>::UnTrimmedHostParameterSpaceGridEntity& element,
      const PatchTrimData& trimData) {
    std::cout << "START " << std::endl;
    using namespace Clipper2Lib;
    auto geo                 = element.geometry();
    std::array<FieldVector<ctype, 2>, 4> corners;  // see dune book page 127 Figure 5.12
    corners[0] = geo.corner(0);
    corners[1] = geo.corner(1);
    corners[2] = geo.corner(3);
    corners[3] = geo.corner(2);

    PathD elementPath;
    for (int i = 0; i < 4; ++i)
      elementPath.push_back({corners[i][0], corners[i][1], i});

    PathsD clipPaths;
    PathD tempPath;

    for (auto loop : trimData.loops()) {
      tempPath.clear();
      for (auto& curve: loop.curves()) {
        for (auto v : Utilities::linspace(curve.domain()[0], 2)) {
          auto fV = curve.global({v});
          tempPath.emplace_back(fV[0], fV[1]);
        }
        // Add an additional point just outside of the element
        auto localLastPoint = curve.domain()[0][1];
        auto lastPoint = curve.global(localLastPoint);
        auto scale = element.geometry().volume() / 10;
        auto dx = curve.jacobian({localLastPoint}) * scale;
        tempPath.emplace_back(dx[0] + lastPoint[0], dx[1] + lastPoint[1]);
      }
      // For curve other than the first one, delete first node
      if (!clipPaths.empty())
        tempPath.erase(tempPath.begin());
      clipPaths.push_back(tempPath);
    }


    auto result = clipElementRectangle(elementPath, clipPaths);
  }
}  // namespace Dune::IGANEW::DefaultTrim
