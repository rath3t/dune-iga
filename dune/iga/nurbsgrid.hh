// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once
#include <clipper2/clipper.h>
#include <utility>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgridindexsets.hh>
#include <dune/iga/nurbsgridtraits.hh>
#include <dune/iga/nurbsidset.hh>
#include <dune/iga/nurbsintersection.hh>
#include <dune/iga/nurbsleafgridview.hh>
#include <dune/iga/nurbslocalgeometry.hh>
#include <dune/iga/nurbspatch.hh>
#include <dune/iga/nurbstrimmedpatch.hh>
#include <dune/iga/nurbstrimboundary.hh>
#include <dune/iga/nurbstrimfunctionality.hh>
#include <dune/iga/reconstructedgridhandler.hh>
#include <dune/iga/nurbspatchgeometry.h>

namespace Dune::IGA {
  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  class NURBSGrid;

  template <int cd, typename GridImpl>
  class EntitySeedStruct {
   public:
    [[nodiscard]] bool isValid() const { return valid_; }
    static constexpr int codimension = cd;

   private:
    bool valid_{false};
    template <int codim, int dim, typename GridImpl1>
    friend class NURBSGridEntity;
    template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
    friend class NURBSGrid;

    int index_{-1};
  };

  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily;

  /** \brief NURBS grid manager */
  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  class NURBSGrid : public Dune::Grid<dim, dimworld, typename NurbsGridLinearAlgebraTraitsImpl::value_type,
                                      NurbsGridFamily<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>> {
   public:
    using LinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;

    static constexpr std::integral auto dimension      = dim;
    static constexpr std::integral auto dimensionworld = dimworld;
    using ctype                                        = typename LinearAlgebraTraits::value_type;

    using ControlPointNetType =
        typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>::ControlPointNetType;

    using GridFamily = NurbsGridFamily<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;

    using Traits = typename GridFamily::Traits;
    template <int cd>
    using Codim = typename GridFamily::Traits::template Codim<cd>;

    using LeafIndexSet  = typename Traits::LeafIndexSet;
    using LevelIndexSet = typename Traits::LevelIndexSet;
    using GridView      = typename Traits::LeafGridView;
    using ElementEntity = typename Traits::template Codim<0>::Entity;
    using GlobalIdSet   = typename Traits::GlobalIdSet;
    NURBSGrid()         = default;

    // Boundaries can be given as an optional
    using BoundariesOpt = std::optional<std::vector<Boundary>>;

    explicit NURBSGrid(const NURBSPatchData<dim, dimworld, LinearAlgebraTraits>& nurbsPatchData,
                       BoundariesOpt _boundaries = std::nullopt)
        : coarsestPatchRepresentation_{nurbsPatchData},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()},
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this))},
          idSet_{std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView())},
          boundaries(std::move(_boundaries)),
          nurbsSurface{Dune::IGA::Nurbs<dim>(currentPatchRepresentation_)} {
      static_assert(dim <= 3, "Higher grid dimensions are unsupported");
      assert(nurbsPatchData.knotSpans[0].size() - nurbsPatchData.degree[0] - 1 == nurbsPatchData.controlPoints.size()[0]
             && "The size of the controlpoints and the knotvector size do not match in the first direction");
      if constexpr (dim > 1)
        assert(nurbsPatchData.knotSpans[1].size() - nurbsPatchData.degree[1] - 1
                   == nurbsPatchData.controlPoints.size()[1]
               && "The size of the controlpoints and the knotvector size do not match in the second direction");
      if constexpr (dim > 2)
        assert(nurbsPatchData.knotSpans[2].size() - nurbsPatchData.degree[2] - 1
                   == nurbsPatchData.controlPoints.size()[2]
               && "The size of the controlpoints and the knotvector size do not match in the third direction");
      // FIXME check sanity of knotvector and degree
      finestPatches_->emplace_back(currentPatchRepresentation_);
      leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      idSet_        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());
      indexdSet_    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());

      trimFlags = std::vector<ElementTrimFlag>(size(0));
      std::fill(trimFlags.begin(), trimFlags.end(), ElementTrimFlag::full);

      if (boundaries.has_value()) trimElement();
    }

    /** \brief  constructor
     *
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order degree of the B-Spline structure for each dimension
     */
    NURBSGrid(const std::array<std::vector<double>, dim>& knotSpans, const ControlPointNetType& controlPoints,
              const std::array<int, dim>& order)
        : coarsestPatchRepresentation_{NURBSPatchData<dim, dimworld, LinearAlgebraTraits>(knotSpans, controlPoints,
                                                                                          order)},
          currentPatchRepresentation_{coarsestPatchRepresentation_},
          finestPatches_{std::make_shared<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>>()},
          leafGridView_{std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(currentPatchRepresentation_, *this))},
          idSet_{std::make_unique<Dune::IGA::IgaIdSet<NURBSGrid>>(*this)},
          indexdSet_{std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl())},
          nurbsSurface{Dune::IGA::Nurbs<dim>(currentPatchRepresentation_)} {}

    void globalRefine(int refinementLevel) {
      if (refinementLevel == 0) return;
      for (int refDirection = 0; refDirection < dim; ++refDirection) {
        auto additionalKnots
            = generateRefinedKnots(currentPatchRepresentation_.knotSpans, refDirection, refinementLevel);
        currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, refDirection);
      }
      updateStateAfterRefinement();
    }

    void globalRefineInDirection(const int dir, const int refinementLevel) {
      if (refinementLevel == 0) return;
      auto additionalKnots        = generateRefinedKnots(currentPatchRepresentation_.knotSpans, dir, refinementLevel);
      currentPatchRepresentation_ = knotRefinement<dim>(currentPatchRepresentation_, additionalKnots, dir);
      updateStateAfterRefinement();
    }

    [[nodiscard]] int size(int codim) const { return finestPatches_.get()->front().size(codim); }

    /** \brief returns the number of boundary segments within the macro grid */
    [[nodiscard]] int numBoundarySegments() const {
      if constexpr (dimension == 1)
        return 2;
      else if constexpr (dimension == 2)
        return (finestPatches_.get()->front().validKnotSize()[0] + finestPatches_.get()->front().validKnotSize()[1])
               * 2;
      else if constexpr (dimension == 3)
        return 2
               * (finestPatches_.get()->front().validKnotSize()[0] * finestPatches_.get()->front().validKnotSize()[1]
                  + finestPatches_.get()->front().validKnotSize()[1] * finestPatches_.get()->front().validKnotSize()[2]
                  + finestPatches_.get()->front().validKnotSize()[0]
                        * finestPatches_.get()->front().validKnotSize()[2]);
      __builtin_unreachable();
    }
    [[nodiscard]] int size(int level, int codim) const { return this->size(codim); }

    const GridView& leafGridView() const { return *leafGridView_; }
    const GridView& levelGridView([[maybe_unused]] int level) const { return *leafGridView_; }
    int getMark(const ElementEntity& element) const { return 0; }
    bool mark(int refCount, const ElementEntity& element) { return false; }

    template <typename Seed>
    typename Codim<Seed::codimension>::Entity entity(const Seed& seed) const {
      return leafGridView_->impl().template getEntity<Seed::codimension>(seed.impl().index_);
    }

    [[nodiscard]] int size(const GeometryType& type) const {
      if (type == Dune::GeometryTypes::vertex || type == Dune::GeometryTypes::cube(1)
          || type == Dune::GeometryTypes::cube(2) || type == Dune::GeometryTypes::cube(3))
        return this->leafGridView().size(dimension - type.dim());
      else
        return 0;
    }
    int size(int lvl, const GeometryType& type) const { return this->size(type); }

    const GlobalIdSet& globalIdSet() const { return *idSet_; }
    const LevelIndexSet& levelIndexSet(int lvl) const { return *indexdSet_; }
    const LeafIndexSet& leafIndexSet() const { return *indexdSet_; }

    [[nodiscard]] int maxLevel() const { return 0; }

    const auto& localIdSet() const { return this->globalIdSet(); }

    [[nodiscard]] const typename Traits::CollectiveCommunication& comm() const { return ccobj; }

    // The following 3 functions are an implementation to get global Coordinates for the whole patch (not on element
    // lvl)
    [[nodiscard]] ElementEntity::Geometry::Implementation geometry(std::array<int, dim> spanIndex) const {
      using Impl = ElementEntity::Geometry::Implementation;
      return Impl(std::make_shared<NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits>>(
                      currentPatchRepresentation_),
                  {Dune::IGA::Impl::FixedOrFree::free}, spanIndex);
    }

    [[nodiscard]] std::array<int, dim> findSpanIndex(std::array<double, dim> u) const {
      std::array<int, dim> spans{};

      for (int i = 0; i < dim; ++i)
        spans[i] = static_cast<int>(Dune::IGA::findSpanUncorrected(currentPatchRepresentation_.degree[i], u[i],
                                                                   currentPatchRepresentation_.knotSpans[i]));

      return spans;
    }

    [[nodiscard]] FieldVector<ctype, dimworld> global(Dune::FieldVector<double, dim> u) const {
      // Bad workaround -> into array<double, dim>
      std::array<double, dim> pointArray{};
      for (int i = 0; i < dim; ++i)
        pointArray[i] = u[i];

      auto spanIndex = findSpanIndex(pointArray);
      auto geo       = geometry(spanIndex);
      auto basis     = nurbsSurface.basisFunctionNet(pointArray);

      return Dune::IGA::dot(basis, geo.cpCoordinateNet_);
    }

    // For dim != 2
    std::optional<typename ReconstructedGridHandler<dimworld>::GridView> getReconstructedGridViewForTrimmedElement(
        int index) {
      return std::nullopt;
    }

    // Get if available the reconstructedGridView for a particular Element (only valid for dim == 2 (and also only for
    // dimwold == 2))
    std::optional<typename ReconstructedGridHandler<dimworld>::GridView> getReconstructedGridViewForTrimmedElement(
        int index)
      requires(dimworld == 2 && dim == 2)
    {
      if ((boundaries.has_value()) && (trimFlags[index] == ElementTrimFlag::trimmed))
        return std::make_optional<typename ReconstructedGridHandler<dimworld>::GridView>(
            trimResultMap[index]->grid->leafGridView());
      else
        return std::nullopt;
    }

   private:
    void updateStateAfterRefinement() {
      finestPatches_.get()->front() = NURBSPatch<dim, dimworld, LinearAlgebraTraits>(currentPatchRepresentation_);
      leafGridView_                 = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this));
      indexdSet_                    = std::make_unique<NURBSGridLeafIndexSet<NURBSGrid>>(this->leafGridView().impl());
      idSet_                        = std::make_unique<IgaIdSet<NURBSGrid>>(this->leafGridView());

      if (boundaries.has_value()) {
        nurbsSurface = Dune::IGA::Nurbs<dim>(currentPatchRepresentation_);
        trimResultMap.clear();
        trimElement();
      }
    }

    typename Traits::CollectiveCommunication ccobj;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> coarsestPatchRepresentation_;
    NURBSPatchData<(size_t)dim, (size_t)dimworld, LinearAlgebraTraits> currentPatchRepresentation_;
    std::shared_ptr<std::vector<NURBSPatch<dim, dimworld, LinearAlgebraTraits>>> finestPatches_;
    std::shared_ptr<GridView> leafGridView_;
    std::unique_ptr<LeafIndexSet> indexdSet_;
    std::unique_ptr<IgaIdSet<NURBSGrid>> idSet_;

    Dune::IGA::Nurbs<dim> nurbsSurface;
    std::vector<ElementTrimFlag> trimFlags;
    BoundariesOpt boundaries;
    std::map<int, std::unique_ptr<ReconstructedGridHandler<dimworld>>> trimResultMap;

    auto parameterSpaceGrid() {
      auto gV        = leafGridView();
      auto patchData = gV.impl().getPatchData();

      using TensorSpaceCoordinates = TensorProductCoordinates<double, dimension>;
      using ParametricGrid         = YaspGrid<dimension, TensorSpaceCoordinates>;

      std::array<std::vector<double>, dimension> tensorProductCoordinates;
      for (int i = 0; i < dimension; ++i) {
        auto unique_it              = std::unique(patchData.knotSpans[i].begin(), patchData.knotSpans[i].end());
        tensorProductCoordinates[i] = std::vector<double>(patchData.knotSpans[i].begin(), unique_it);
      }
      // assert(gV.size(0) == ParametricGrid(tensorProductCoordinates).leafGridView().size(0));

      return std::make_unique<ParametricGrid>(tensorProductCoordinates);
    }

    void trimElement() {
      // Nothing to do
      std::cout << "No trim functionality for grids that arent dim == 2 and worldDim = 2" << std::endl;
    }

    void trimElement()
      requires(dim == 2 && dimworld == 2)
    {
      // Set up trimFlags
      trimFlags = std::vector<ElementTrimFlag>(size(0));
      std::fill(trimFlags.begin(), trimFlags.end(), ElementTrimFlag::full);

      if (!boundaries.has_value()) return;

      // Local variable for the result (idk if this is really necessary)
      std::map<int, std::unique_ptr<ReconstructedGridHandler<dimworld>>> _trimResultMap;

      // Parameters and Preferences
      using namespace Impl::Trim;
      const int clipperPrecision = 8;
      const int pathSamples      = 200;

      auto paraGrid               = parameterSpaceGrid();
      auto parameterSpaceGridView = paraGrid->leafGridView();

      // Some Lambdas // (move into functionality)
      constexpr double tolerance  = 1e-8;
      auto approximatelyEqual = [tolerance](auto a, auto b) {
        // const double tolerance = double(16) * std::numeric_limits<double>::epsilon();
        return fabs(a - b) < tolerance;
      };

      auto distance = [](ClipperPoint p1, ClipperPoint p2) {
        return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
      };

      auto pointOnLine = [&approximatelyEqual, distance](ClipperPoint p, ClipperPoint a, ClipperPoint b) {
        // https://stackoverflow.com/a/17693146
        return (approximatelyEqual(distance(a, p) + distance(b, p), distance(a, b)));
      };

      auto approximatelySamePoint = [&approximatelyEqual](ClipperPoint a, ClipperPoint b) {
        return (approximatelyEqual(a.x, b.x) && approximatelyEqual(a.y, b.y));
      };

      // Get Clip as ClipperPath
      Clipper2Lib::PathD pathParametric;
      for (auto& boundary : boundaries.value())
        for (auto& point : boundary.path(pathSamples))
          pathParametric.push_back(point);

      Clipper2Lib::PathsD clip = {pathParametric};

      // Prepping Clip Path
      // This is imported later when we search for intersection Points with the Grid, so we don't get duplicate
      // intersection Points
      assert(clip.size() == 1);
      clip.front() = Clipper2Lib::TrimCollinear(clip.front(), clipperPrecision);

      const auto& indexSet = parameterSpaceGridView.indexSet();
      for (auto& element : elements(parameterSpaceGridView)) {
        auto index{indexSet.index(element)};
        auto geo = element.geometry();
        std::vector<FieldVector<double, 2>> corners;
        corners.push_back({geo.corner(0)[0], geo.corner(0)[1]});
        corners.push_back({geo.corner(1)[0], geo.corner(1)[1]});
        corners.push_back({geo.corner(3)[0], geo.corner(3)[1]});  // see dune book page 127 Figure 5.12
        corners.push_back({geo.corner(2)[0], geo.corner(2)[1]});

        Clipper2Lib::PathD elementEdges;
        for (int i = 0; i < 4; ++i)
          elementEdges.emplace_back(corners[i][0], corners[i][1]);

        Clipper2Lib::RectD elementRect{corners[0][0], corners[1][1], corners[1][0], corners[3][1]};

        // Clip with RectClip (more efficient than normal Intersect Clip, works because our ParameterSpace Grid
        // has always rectangular elements
        Clipper2Lib::PathsD clippedEdges = Clipper2Lib::RectClip(elementRect, clip, clipperPrecision);

        // If there are no clipped edges, the element is not impacted by the trimming (empty Element)
        // Also here the FloatCompare sometimes has to a lot more liberal, especially with 3D shell-like structures

        if (clippedEdges.empty())
          trimFlags[index] = ElementTrimFlag::empty;
        else {
          if (FloatCmp::eq(Clipper2Lib::Area(clippedEdges), Clipper2Lib::Area(elementEdges)))
            trimFlags[index] = ElementTrimFlag::full;
          else
            trimFlags[index] = ElementTrimFlag::trimmed;
        }
        // Update the GridView with the set trimFlags
        leafGridView_ = std::make_shared<GridView>(NURBSLeafGridView<NURBSGrid>(finestPatches_, *this, trimFlags));

        if (!(trimFlags[index] == ElementTrimFlag::trimmed)) continue;

        // Now we only look at the trimmed elements

        // Wenn Code hier → Elemente sind betroffen
        assert(clippedEdges.size() == 1);
        clippedEdges.front() = Clipper2Lib::TrimCollinear(clippedEdges.front(), clipperPrecision);

        // The following can be made a lot more concise. E.g. the lambda pointOnAnyEdge could report the edge it has
        // found an intersection Point, so we don't have to loop over the points again
        // Also accuracy is a huge problem here. There has to be a balance that all Intersection Points are all found
        // but also not inexact enough that there are more than one IP found next to each other (next Point on
        // ClipperPath)

        // Find the intersection Points and determine
        int intersectionPointCount = 0;
        std::vector<Clipper2Lib::PointD> intersectionPoints;
        for (auto& point : clippedEdges[0])
          // Todo: etwas in Parameter einbauen, dass hier wechselt evtl
          // if (Clipper2Lib::PointInPolygon(point, elementEdges) == Clipper2Lib::PointInPolygonResult::IsOn) {
          if (pointOnAnyEdge(point, elementEdges)) {
            intersectionPointCount++;
            intersectionPoints.push_back(point);
          }
        int n_edges = 4;
        // Determine which edges are cut
        // Result: -1 not on edge (should definitely not happen), 10-13 is on corner, else 0-3 for edges
        auto whichEdgeOrNode
            = [&pointOnLine, &elementEdges, &approximatelySamePoint, n_edges](const Clipper2Lib::PointD& point) -> int {
          // Check if point lies in corner (node)
          for (int i = 0; const auto& p : elementEdges) {
            if (approximatelySamePoint(p, point))
              return 10 + i;
            else
              ++i;
          }
          for (int i = 0; i < n_edges; ++i)
            if (pointOnLine(point, elementEdges[i], elementEdges[iPlusX(i, 1)])) return i;

          return -1;
        };

        // Create a local pointMap to store the intersection points in regard to their corresponding edges or nodes
        IntersectionPointMap pointMap;

        // The following numbering is used for edges and nodes:
        /*
         * 3 -- 2 -- 2
         * |         |
         * 3         1
         * |         |
         * 0 -- 0 -- 1
         */

        // This keeps track of how many intersection Points there are at a given edge
        std::array<int, 4> edgeCounter{0, 0, 0, 0};
        std::array<int, 4> nodeCounter{0, 0, 0, 0};

        for (int i = 0; const auto& point : intersectionPoints) {
          int edgeNodeNr{whichEdgeOrNode(point)};
          pointMap[edgeNodeNr].emplace_back(intersectionPoints[i]);

          if (edgeNodeNr >= 0 && edgeNodeNr < 10)
            edgeCounter[edgeNodeNr]++;
          else if (edgeNodeNr >= 10)
            nodeCounter[edgeNodeNr - 10]++;

          ++i;
        }

        IntersectionPointMapPtr pointMapPtr = std::make_shared<IntersectionPointMap>(pointMap);

        // Sort IntersectionPoints
        for (int i = 0; i < 4; ++i)
          sortIntersectionPoints(pointMapPtr, i);

        auto hasIntersectionPointOnEdgeNr = [&edgeCounter](int nr) -> bool { return edgeCounter[nr] > 0; };
        auto hasIntersectionPointOnNodeNr = [&nodeCounter](int nr) -> bool { return nodeCounter[nr] > 0; };

        int edgePointCount = std::accumulate(edgeCounter.begin(), edgeCounter.end(), 0);
        int nodePointCount = std::accumulate(nodeCounter.begin(), nodeCounter.end(), 0);

        // Jetzt gehen wir systematisch vor

        std::vector<Boundary> elementBoundaries;

        struct LoopState {
          bool isCurrentlyOnNode            = true;
          bool moreIntersectionPointsOnEdge = false;
          int edgeTheLoopIsOn               = -1;     // This variable is only filed when !loopIsCurrentlyOnNode
          bool onlyOneSide                  = false;  // Currently not used
        };

        // Find the first IntersectionPointOnNode
        auto it_nodes = std::find_if(nodeCounter.begin(), nodeCounter.end(), [](auto x) { return x > 0; });

        // Find the first IntersectionPointOnEdge (currently not used)
        auto it_edges = std::find_if(edgeCounter.begin(), edgeCounter.end(), [](auto x) { return x > 0; });

        int nodeToBegin;
        LoopState state{};
        if (it_nodes == nodeCounter.end()) {
          std::cout << "Not implemented" << std::endl;
          // TODO: For now those elements are skipped and marked as empty
          trimFlags[index] = ElementTrimFlag::empty;
          continue;
        } else
          nodeToBegin = std::distance(nodeCounter.begin(), it_nodes);

        int node                = nodeToBegin;
        Point current_NodePoint = corners[node];
        while (true) {
          int edge = node;  // node 0 korrespondiert zu edge 0, die von Node 0 zu Node 1 läuft

          // Check if edge has no intersections and next node has one
          if (state.isCurrentlyOnNode && !hasIntersectionPointOnEdgeNr(edge)
              && hasIntersectionPointOnNodeNr(iPlusX(node, 1))) {
            elementBoundaries.push_back(boundaryForEdge(corners, edge));

            // Update Node and Point
            node              = iPlusX(node, 1);
            current_NodePoint = corners[node];
            if (node == nodeToBegin)
              break;
            else
              continue;
          }

          if (state.isCurrentlyOnNode && !hasIntersectionPointOnEdgeNr(edge)
              && !hasIntersectionPointOnNodeNr(iPlusX(node, 1))) {
            std::cerr << "Nicht implementiert, höchstwahrscheinlich unendlicher Loop ab hier" << std::endl;

            // Todo: So lange durch die Nodes loopen, bis ein Intersection Point gefunden wurde, aufpassen dass dieser
            // auch der StartNode sein kann
          }
          if (state.isCurrentlyOnNode) {
            // If we are here, then there is an intersectionPoint on the current edge, we have to find the curve
            // that is intersecting there and then trace it_nodes to find the next intersection Point on the next edge
            // (or node)
            assert(hasIntersectionPointOnEdgeNr(edge));

            int current_node  = node;
            current_NodePoint = corners[current_node];

            ////////
            // Find the intersection Point on the current Edge
            ////////

            auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult, found]
                = findBoundaryThatHasIntersectionPoint(pointMapPtr, edge, corners, *boundaries);

            // There has to be an intersectionPoint
            if (!found) throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

            // Make a parametrisation of the section before the intersectionPoint
            PointVector newEdgeParametrisation{current_NodePoint, edgeIntersectResult.globalResult};
            elementBoundaries.emplace_back(newEdgeParametrisation);

            ////////
            // Curve Tracing
            ////////

            // Then we need the parametrisation for the curve to come
            Boundary boundaryToTrace       = (*boundaries)[boundaryIndexThatHasIntersectionPoint];
            CurveGeometry curveToTrace = boundaryToTrace.nurbsGeometry;

            // Call traceCurve
            // Todo: as for now we are only checking against edge Intersection Points, next intersection could be on a
            // node (highly unlikely but with bad luck due to precision in Intersection Point, an intersection could be
            // registered on node, instead of an edge)
            double beginU = edgeIntersectResult.localResult;
            auto [intersectU, foundInterSectionPoint, foundOnEdgeNr]
                = traceCurve(pointMapPtr, corners, curveToTrace, beginU);

            // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
            Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{beginU, intersectU});
            elementBoundaries.push_back(newBoundary);

            // Now we assume we are on an edge (has to be atm)
            state.isCurrentlyOnNode = false;
            current_NodePoint       = foundInterSectionPoint;
            state.edgeTheLoopIsOn   = foundOnEdgeNr;

            state.moreIntersectionPointsOnEdge
                = getNextIntersectionPointOnEdge(pointMapPtr, state.edgeTheLoopIsOn).second;
          }
          // Wir sollten jetzt nach gefundenem IntersectionPoint auf der Edge in einem der folgenden If-Satement
          // befinden

          // Todo: There is some code repetition here, but its not exactly the same behaviour than in the above if
          if (!state.isCurrentlyOnNode && state.moreIntersectionPointsOnEdge && !state.onlyOneSide) {
            // Call traceCurve
            edge = state.edgeTheLoopIsOn;

            ////////
            // Find the intersection Point on the current Edge
            ////////

            auto [boundaryIndexThatHasIntersectionPoint, edgeIntersectResult, found]
                = findBoundaryThatHasIntersectionPoint(pointMapPtr, edge, corners, *boundaries);

            // There has to be an intersectionPoint
            if (!found) throw std::runtime_error("No Intersection Point found on a edge where there has to be one.");

            // Make a parametrisation of the section before the intersectionPoint
            PointVector newEdgeParametrisation{current_NodePoint, edgeIntersectResult.globalResult};
            elementBoundaries.emplace_back(newEdgeParametrisation);

            ////////
            // Curve Tracing
            ////////

            // Then we need the parametrisation for the curve to come
            Boundary boundaryToTrace       = (*boundaries)[boundaryIndexThatHasIntersectionPoint];
            CurveGeometry curveToTrace = boundaryToTrace.nurbsGeometry;

            // Call traceCurve
            double beginU = edgeIntersectResult.localResult;
            auto [intersectU, foundInterSectionPoint, foundOnEdgeNr]
                = traceCurve(pointMapPtr, corners, curveToTrace, beginU);

            // After we found the next IntersectionPoint on this curve, add the so found part of the elementBoundary
            Boundary newBoundary = Boundary(curveToTrace, std::array<double, 2>{beginU, intersectU});
            elementBoundaries.push_back(newBoundary);

            // Now for now we assume we are still on an edge (Todo: again check against node Intersection Points)
            state.isCurrentlyOnNode = false;
            current_NodePoint       = foundInterSectionPoint;
            state.edgeTheLoopIsOn   = foundOnEdgeNr;

            state.moreIntersectionPointsOnEdge
                = getNextIntersectionPointOnEdge(pointMapPtr, state.edgeTheLoopIsOn).second;
          }

          if (!state.isCurrentlyOnNode && !state.moreIntersectionPointsOnEdge && !state.onlyOneSide) {
            PointVector v1{current_NodePoint, corners[iPlusX(state.edgeTheLoopIsOn, 1)]};
            elementBoundaries.emplace_back(v1);

            state.isCurrentlyOnNode = true;

            // Update Node and Point
            node                  = iPlusX(state.edgeTheLoopIsOn, 1);
            current_NodePoint     = corners[node];
            state.edgeTheLoopIsOn = -1;
            if (node == nodeToBegin)
              break;
            else
              continue;
          }
        }

        // ReconstructGrid and save in a GridHandler
        auto transferToGlobal
            = std::make_shared<std::function<FieldVector<double, dimworld>(FieldVector<double, dim>)>>(
                [this](FieldVector<double, dim> u) -> FieldVector<double, dimworld> { return this->global(u); });

        auto recoPtr = std::make_unique<ReconstructedGridHandler<dimworld>>(elementBoundaries, transferToGlobal,
                                                                            currentPatchRepresentation_.degree);
        _trimResultMap[index] = std::move(recoPtr);

      }  // Element Loop End

      // Pass the local trimResultMap to the member variable, idk why it wouldnt work directly (segfaults)
      trimResultMap = std::move(_trimResultMap);
    }
  };

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto levelGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid, int level) {
    return grid.levelGridView(level);
  }

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  auto leafGridView(const NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>& grid) {
    return grid.leafGridView();
  }

  template <int dim, int dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl>
  struct NurbsGridFamily {
    using GridImpl = Dune::IGA::NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraitsImpl>;
    using Traits   = NurbsGridTraits<dim, dimworld, GridImpl, NURBSGeometry, NURBSGridEntity, NURBSGridLeafIterator,
                                   NURBSintersection, NURBSintersection, NURBSGridInterSectionIterator,
                                   NURBSGridInterSectionIterator, NurbsHierarchicIterator, NURBSGridLeafIterator,
                                   NURBSGridLeafIndexSet<GridImpl>, NURBSGridLeafIndexSet<GridImpl>, IgaIdSet<GridImpl>,
                                   int, IgaIdSet<GridImpl>, int, Communication<No_Comm>, NurbsLeafGridViewTraits,
                                   NurbsLeafGridViewTraits, EntitySeedStruct, NURBSLocalGeometry>;
  };

}  // namespace Dune::IGA
