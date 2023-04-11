//
// Created by Henri on 10.02.2023.
//

#pragma once

#include <clipper2/clipper.h>
#include <matplot/matplot.h>
#include <string>
#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/nurbspatchgeometry.h>
#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/nurbstrimfunctionality.hh>
#include <dune/grid/io/file/vtk.hh>

namespace Plot {

namespace Draw {

template <typename GridView>
void drawElements(const GridView& gridView, double lineWidth = 1, std::string&& lineColor = "black",
                  std::string&& nodeColor = "red") {
    using namespace matplot;

    auto ax = gca();
    hold(ax, true);

    // Code from Ikarus
    constexpr int edgeCodim = GridView::dimension - 1;

    for (auto& element : elements(gridView)) {
        std::array<std::array<double, 2>, GridView::dimensionworld> edgeCoords{};
        for (size_t edgeIndex = 0; edgeIndex < element.subEntities(edgeCodim); ++edgeIndex) {
            auto edge = element.template subEntity<edgeCodim>(edgeIndex);
            for (int i = 0; i < 2; ++i) {
                const auto vertCoords = edge.geometry().corner(i);
                for (int j = 0; j < GridView::dimensionworld; ++j)
                    edgeCoords[j][i] = vertCoords[j];
            }
            auto l = ax->plot(edgeCoords[0], edgeCoords[1], "-o");
            l->line_width(lineWidth);
            l->color(lineColor);
            l->marker_size(lineWidth * 2);
            l->marker_face_color(nodeColor);
        }
    }
}


void drawPoints(std::vector<Clipper2Lib::PointD>& points) {
    using namespace matplot;

    // auto figure{gcf()};
    auto ax{gca()};
    hold(ax, true);

    // prepare data
    std::vector<double> x;
    std::vector<double> y;
    for (const auto& point : points) {
        x.push_back(point.x);
        y.push_back(point.y);
    }
    ax->scatter(x, y);
}

template <typename T>
void drawPaths(Clipper2Lib::Paths<T> paths, std::string&& lineColor, auto lineWidth, bool isLoop = true,
               bool fill = true) {
    std::vector<double> x;
    std::vector<double> y;

    using namespace matplot;

    auto ax{gca()};
    hold(ax, true);

    // Loop über paths
    for (auto& path : paths) {
        x.clear();
        y.clear();

        // Loop over points in path
        for (auto& point : path) {
            x.push_back(point.x);
            y.push_back(point.y);
        }
        // Loop add first point
        if (isLoop) {
            x.push_back(path.front().x);
            y.push_back(path.front().y);
        }

        // Draw path
        auto l{ax->plot(x, y)};
        l->color(lineColor);
        l->line_width(lineWidth);

        if (fill) {
            auto fi = matplot::fill(x, y);
            auto c{l->color()};
            c[0] = 0.8;
            fi->color(c);
        }
    }
}

}  // namespace Draw

void plotGridView(auto& gridView, std::string&& file_name) {
    auto figure = matplot::figure(true);
    Draw::drawElements(gridView);

    axis(matplot::gca(), matplot::equal);
    figure->size(1000, 800);

    figure->save(file_name + ".svg");
    figure->save(file_name + ".jpg");
}

void plotPaths(auto& paths, std::string&& file_name) {
    auto figure = matplot::figure(true);
    Draw::drawPaths(paths, "black", 1);

    axis(matplot::gca(), matplot::equal);
    figure->size(1000, 800);

    //figure->save(file_name + ".png");
    figure->save(file_name + ".jpg");
}

void plotGridViewAndPaths(auto& gridView, Clipper2Lib::PathsD& paths, std::string&& file_name, bool markEndPoints = true) {
    auto figure = matplot::figure(true);
    Draw::drawElements(gridView, 0.5);
    Draw::drawPaths(paths, "red", 1, false, false);

    if (markEndPoints) {
      std::vector<Clipper2Lib::PointD> endPoints;
      for (auto& path : paths) {
        endPoints.push_back(path.front());
      }

      Draw::drawPoints(endPoints);
    }

    axis(matplot::gca(), matplot::equal);
    figure->size(1000, 800);

    //figure->save(file_name + ".png");
    figure->save(file_name + ".jpg");
}

void plotParametricGridAndPhysicalGrid(const std::shared_ptr<Dune::IGA::NURBSGrid<2, 2>>& grid, std::string&& postfix = "") {
    if (!(grid->trimData_.has_value()))
      return;

    auto geometry = Dune::IGA::NURBSPatchGeometry<2, 2>(std::make_shared<Dune::IGA::NURBSPatchData<2, 2>>(grid->currentPatchRepresentation_));
    auto boundarieLoops = grid->trimData_.value();

    Clipper2Lib::PathsD transformedPaths;
    Clipper2Lib::PathsD parametricPaths;

    for (auto& loop : boundarieLoops->boundaryLoops) {
      for (auto& trimInfo : loop.boundaries) {
        auto path = trimInfo.path(200, false);
        Clipper2Lib::PathD transformedPath;
        for (auto& point : path) {
            auto transformedPoint = geometry.global({point.x, point.y});
            transformedPath.emplace_back(transformedPoint[0], transformedPoint[1]);
        }
        parametricPaths.push_back(path);
        transformedPaths.push_back(transformedPath);
      }
    }

    auto paraGrid = grid->getPatch().parameterSpaceGrid();

    auto paraGridView = paraGrid->leafGridView();
    plotGridViewAndPaths(paraGridView, parametricPaths, "plot" + postfix + "/parametricGrid");

    // Print Brep in physical Space
    auto patchGridView = grid->leafGridView();
    Plot::plotGridViewAndPaths(patchGridView, transformedPaths, "plot"+ postfix +"/grid");

}


void plotEveryReconstructedGrid(const std::shared_ptr<Dune::IGA::NURBSGrid<2, 2>>& grid, std::string&& postfix = "") {
  for (int i = 0; auto& ele : elements(grid->leafGridView())) {
      if (ele.impl().getTrimFlag() == Dune::IGA::ElementTrimFlag::trimmed) {
        auto gV = grid->getPatch().getTrimmedElementRepresentation(i).value()->getGridView();
        Plot::plotGridView(gV, "plot" + postfix + "/reconstruction/grid_" + std::to_string(i));
      }
      ++i;
  }
}

void saveEveryReconstructedGrid(const std::shared_ptr<Dune::IGA::NURBSGrid<2, 2>>& grid, std::string&& postfix = "") {
  for (int i = 0; auto& ele : elements(grid->leafGridView())) {
      if (ele.impl().getTrimFlag() == Dune::IGA::ElementTrimFlag::trimmed) {
        auto gV = grid->getPatch().getTrimmedElementRepresentation(i).value()->getGridView();
        Dune::VTKWriter vtkWriter(gV);
        vtkWriter.write("plot" + postfix + "/reconstruction/grid_" + std::to_string(i));
      }
      ++i;
  }
}


}  // namespace Plot
