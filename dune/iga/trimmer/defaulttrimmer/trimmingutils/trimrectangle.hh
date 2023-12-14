// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <clipper2/clipper.h>
#include <variant>


namespace Dune::IGANEW::DefaultTrim {

  struct TrimmingResult {
    std::size_t nVertices;
    std::size_t newVertices = 4;
    std::array<bool, 4> verticesVisited{false, false, false, false};
    std::array<bool, 4> edgesVisited{false, false, false, false};

    struct HostVertex {
      std::size_t idx{};
      std::size_t originalIdx{};
    };
    struct NewVertex {
      std::size_t idx{};
      int onEdgeIdx{};
      std::size_t getsNewID{};

      Clipper2Lib::PointD pt;
    };

    using VertexVariant = std::variant<HostVertex, NewVertex>;
    std::vector<VertexVariant> vertices;

    std::size_t addOriginalVertex(const size_t originalIdx) {
      verticesVisited[originalIdx] = true;
      vertices.emplace_back(HostVertex(nVertices++, originalIdx));
      return nVertices;
    }
    std::size_t addNewVertex(const int edgeIdx, const Clipper2Lib::PointD& pt) {
      edgesVisited[edgeIdx] = true;
      vertices.emplace_back(NewVertex(nVertices++, edgeIdx, newVertices++, pt));
      return nVertices;
    }
    void finish() {
      // We have to somehow close the loop
      for (size_t i = 0; const auto didIvisit : verticesVisited) {
        if (not didIvisit) {
          if (not (edgesVisited[i] and edgesVisited[i-1])) {
            // @todo unlikely ? but edge 0 is not visited than -1
            addOriginalVertex(i);
          }
        }
        i++;
      }
    }
  };

  inline int giveEdgeIdx(int e1, int e2) {
    if ((e1 == 0 and e2 == 1) or (e1 == 1 and e2 == 0))
      return 0;
    if ((e1 == 1 and e2 == 2) or (e1 == 2 and e2 == 1))
      return 1;
    if ((e1 == 2 and e2 == 3) or (e1 == 3 and e2 == 2))
      return 2;
    if ((e1 == 3 and e2 == 0) or (e1 == 0 and e2 == 3))
      return 3;
    __builtin_unreachable();
  }

  inline auto trimElementImpl(Clipper2Lib::PathD& eleRect, Clipper2Lib::PathsD& trimmingCurves) {
    using namespace Clipper2Lib;

    ClipperD clipper(8);
    // todo just add cuve ID to z value of curve points, than we can just same this numer in NewVertex
    // Counters
    std::size_t counter     = 0;
    std::size_t loopCounter = 0;
    std::vector<size_t> counterBreaks{};

    // Give points in eleRect the z-Val from 0 to 4 counter-clockwise
    std::ranges::for_each(eleRect, [&](auto& point) { point.z = counter++; });
    counterBreaks.push_back(counter);

    for (auto& trimmingCurve : trimmingCurves) {
      std::ranges::for_each(trimmingCurve, [&](auto& point) { point.z = counter++; });
      counterBreaks.push_back(counter);
      loopCounter++;
    }

    // std::cout << "Counter: " << counter << std::endl;
    // std::cout << "loopCounter: " << loopCounter << std::endl;
    // std::cout << "Counter Breaks:\n";
    // for (const auto& c : counterBreaks)
    //   std::cout << c << std::endl;

    TrimmingResult result{};
    std::size_t vertexCounter = 0;

    clipper.SetZCallback(
        [&](const PointD& e1bot, const PointD& e1top, const PointD& e2bot, const PointD& e2top, const PointD& pt) {
          std::cout << "New Intersection x: " << pt.x << " y: " << pt.y << std::endl;
          std::cout << e1bot.z << " " << e1top.z << std::endl;
          std::cout << e2bot.z << " " << e2top.z << std::endl;

          if (e2bot.z > e2top.z and (int) e2top.z != 0) {
            // Trimming path comes in ... add Vertices before
            for (const auto i : std::views::iota(vertexCounter, static_cast<std::size_t>(e2bot.z)))
              vertexCounter = result.addOriginalVertex(i);

            // add new Intersection
            vertexCounter = result.addNewVertex(giveEdgeIdx((int) e2bot.z, (int) e2top.z), pt);
          }
          else {
            // Trimming path goes out
            // @todo there is an edge Case where trimming path goes out on an host vertex
            vertexCounter = result.addNewVertex(giveEdgeIdx((int) e2bot.z, (int) e2top.z), pt);
          }
        });

    clipper.AddClip({eleRect});
    clipper.AddOpenSubject(trimmingCurves);

    PathsD solution{};

    clipper.Execute(ClipType::Intersection, FillRule::EvenOdd, solution);
    result.finish();

    std::cout << "Vertices found\n";

    struct Visitor {
      void operator()(const TrimmingResult::HostVertex& v) {
        std::cout << v.idx << std::endl;
      }
      void operator()(const TrimmingResult::NewVertex& v) {
        std::cout << v.idx << " Edge: " << v.onEdgeIdx << ", gets Idx "  << v.getsNewID << std::endl;
      }
    };

    for (auto& vV : result.vertices)
      std::visit(Visitor{}, vV);
  }

}  // namespace Dune::IGANEW::DefaultTrim
