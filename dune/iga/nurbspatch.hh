// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "nurbsalgorithms.hh"

#include "dune/iga/controlpoint.hh"
#include "dune/iga/nurbsgeometry.hh"
#include "dune/iga/trim/nurbstrimboundary.hh"
#include "dune/iga/trim/nurbstrimmer.hh"
#include "dune/iga/trim/subgrid.hh"
#include "dune/iga/utils/concepts.hh"
#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
namespace Dune::IGA {

  template <int dim, int dimworld, typename ScalarType>
  class NURBSGrid;

  template <typename GridImpl>
  class NURBSGridLeafIndexSet;

  template <typename Grid>
  class NURBSLeafGridView;

  /** \brief Class where the NURBS geometry can work on */
  template <std::integral auto dim, std::integral auto dimworld, typename ScalarType = double>
  class NURBSPatch {
   public:
    using GridImpl = NURBSGrid<dim, dimworld, ScalarType>;
    using Types    = typename GridImpl::Traits::LeafIndexSet::Types;
    friend class NURBSLeafGridView<const GridImpl>;
    template <int codim, int dim1, typename GridImpl1>
    friend class NURBSGridEntity;

    using DirectIndex = unsigned int;
    using RealIndex   = unsigned int;

    using SubGridType = GridImpl::Traits::SubGridType;

    /** \brief Constructor of NURBS from knots, control points, weights and degree
     *  \param knotSpans vector of knotSpans for each dimension
     *  \param controlPoints a dim-dimensional net of control points
     *  \param degree degree of the spline basis for each dimension
     */
    NURBSPatch(const std::array<std::vector<double>, dim>& knotSpans,
               const typename NURBSPatchData<dim, dimworld, ScalarType>::ControlPointNetType controlPoints,
               const std::array<int, dim> degree)
        : NURBSPatch(NURBSPatchData<dim, dimworld, ScalarType>(knotSpans, controlPoints, degree)) {}

    explicit NURBSPatch(const NURBSPatchData<dim, dimworld, ScalarType>& patchData,
                        std::optional<std::shared_ptr<TrimData>> trimData = std::nullopt)
        : patchData_{std::make_shared<NURBSPatchData<dim, dimworld, ScalarType>>(patchData)},
          patchGeometry_{std::make_shared<NURBSPatchGeometry<dim, dimworld>>(patchData_)},
          trimData_(std::move(trimData)) {
      for (int i = 0; i < dim; ++i)  // create unique knotspan vectors
        std::ranges::unique_copy(patchData_->knotSpans[i], std::back_inserter(uniqueKnotVector_[i]),
                                 [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

      for (int i = 0; i < dim; ++i)
        uniqueSpanSize_[i] = uniqueKnotVector_[i].size() - 1;

      elementNet_ = std::make_shared<MultiDimensionNet<dim, double>>(uniqueSpanSize_);
      vertexNet_  = std::make_shared<MultiDimensionNet<dim, double>>(uniqueKnotVector_);

      computeElementTrimInfo();
    }

    template <unsigned int codim, std::integral T = int>
    requires(codim == 0 || codim == dim) [[nodiscard]] std::array<unsigned int, dim> originalGridSize() const {
      if constexpr (codim == 0)
        return elementNet_->template sizeAsT<T>();
      else
        return vertexNet_->template sizeAsT<T>();
    }

    /** \brief Returns the number of entities of a given codimension in the untrimmed case */
    [[nodiscard]] int originalSize(int codim) const {
      assert((codim <= dim) and (codim <= 3) and (codim >= 0));

      if (codim == 0)
        return elementNet_->size();
      else if (codim == dim)
        return vertexNet_->size();
      else if (dim - codim == 1)  // edge case
      {
        int edgeSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? uniqueSpanSize_[i] : uniqueSpanSize_[i] + 1;
          edgeSize += subs;
        }
        return edgeSize;
      } else if (dim - codim == 2 && dim == 3)  // surface case
      {
        int surfSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? uniqueSpanSize_[i] + 1 : uniqueSpanSize_[i];
          surfSize += subs;
        }
        return surfSize;
      }
      __builtin_unreachable();
    }

    /** \brief Returns the number of entities of a given codimension (considers trimmed elements) */
    [[nodiscard]] int size(int codim) const {
      assert((codim <= dim) and (codim <= 3) and (codim >= 0));
      if constexpr (dim != 2 && dimworld != 2) return originalSize(codim);

      if (!(trimData_.has_value())) return originalSize(codim);

      assert(!(std::isnan(n_fullElement)));
      if (codim == 0)
        return (n_fullElement + n_trimmedElement);
      else if (codim == 1)
        return edgeIndexMap.size();
      else
        return vertexIndexMap.size();
    }

    /** \brief Returns the patch for this patch, i.e. the control points and knotvectors */
    const auto& getPatchData() { return patchData_; }

    /** \brief Checks if the element given by it's id is on the patch boundary */
    [[nodiscard]] bool isPatchBoundary(const int& id) const {
      throw std::logic_error("Not Implemented: ask the element if it has a neighbor");
    }

    [[nodiscard]] int numBoundarySegments() const {
      if constexpr (dim == 1)
        return 2;
      else if constexpr (dim == 2) {
        if (trimData_.has_value()) return boundarySegmentList().size();
        return (validKnotSize()[0] + validKnotSize()[1]) * 2;
      } else if constexpr (dim == 3)
        return 2
               * (validKnotSize()[0] * validKnotSize()[1] + validKnotSize()[1] * validKnotSize()[2]
                  + validKnotSize()[0] * validKnotSize()[2]);
      __builtin_unreachable();
    }

    [[nodiscard]] int patchBoundaryIndex(const RealIndex intersectionRealIndex) const requires(dim == 2) {
      auto realIndexOfBoundaryIntersections = boundarySegmentList();

      auto it = std::ranges::find(realIndexOfBoundaryIntersections, intersectionRealIndex);
      assert(it != realIndexOfBoundaryIntersections.end());

      return static_cast<int>(std::ranges::distance(realIndexOfBoundaryIntersections.begin(), it));
    }

    /** \brief Returns the boundary index of the passed codim 1 entity index
     * If the   passed codim 1 entity does not lie on the boundary this method returns random indices!
     * @param ocdim1Id
     * @return Increasing boundary index
     */

    [[nodiscard]] int patchBoundaryIndex(const int ocdim1Id) const {
      int index = ocdim1Id;

      if constexpr (dim == 1) {  // early exit curves only have a left or right boundary point
        if (ocdim1Id == 0) return 0;
        if (ocdim1Id == uniqueSpanSize_[0]) return 1;
      } else {
        const auto codimSizes = this->sizeOfCodim1PerDirection<1>();

        int surfFixedDir = 0;
        int edgesPerDir  = codimSizes[0];
        for (int j = 0; j < dim; ++j)
          if (ocdim1Id < edgesPerDir) {  // find out to which edge direction this edge index belongs
            surfFixedDir = j;
            break;
          } else
            edgesPerDir += codimSizes.at(j + 1);

        // move index to correct subset, e.g. edges in y-dir or surfaces with fixed z-dir
        index -= std::accumulate(codimSizes.begin(), codimSizes.begin() + surfFixedDir, 0);

        if constexpr (dim == 2) {
          if (surfFixedDir == 1) {
            const double fac = std::ceil((double)index / (uniqueSpanSize_[0] + 1));
            index -= fac * (uniqueSpanSize_[0] - 1);
          }
          if (surfFixedDir == 0) {
            const int fac2 = std::floor((double)(index) / ((uniqueSpanSize_[1]) * uniqueSpanSize_[0]));
            index -= fac2 * (uniqueSpanSize_[0] * (uniqueSpanSize_[1] - 1));
          }

          std::array<int, dim> boundaries;
          std::ranges::transform(uniqueSpanSize_, boundaries.begin(), [](auto& i) { return i * 2; });

          // move back to correct index, e.g. the first y - boundary index comes after all x-boundary indices
          index += std::accumulate(boundaries.begin(), boundaries.begin() + surfFixedDir, 0);
        } else if constexpr (dim == 3) {
          if (surfFixedDir == 0) {
            const int fac = (index + 1) / (uniqueSpanSize_[0] + 1);
            index -= fac * (uniqueSpanSize_[0] - 1);
          } else if (surfFixedDir == 1) {
            const int fac = (index + uniqueSpanSize_[0]) / (uniqueSpanSize_[0] * (uniqueSpanSize_[1] + 1));
            index -= fac * (uniqueSpanSize_[0] * (uniqueSpanSize_[1] - 1));
          } else if (surfFixedDir == 2) {
            const int fac = index / ((uniqueSpanSize_[0] * uniqueSpanSize_[1]) * (uniqueSpanSize_[2]));
            index -= fac * (uniqueSpanSize_[1] * uniqueSpanSize_[0] * (uniqueSpanSize_[2] - 1));
          }

          std::array<int, dim> boundaries;
          std::ranges::fill(boundaries, 2);
          for (int i = 0; i < dim; ++i)
            for (int k = 0; k < dim; ++k) {
              if (i == k) continue;
              boundaries[i] *= uniqueSpanSize_[k];
            }

          // move back to correct index, e.g. the first y - boundary index comes after all x-boundary indices
          index += std::accumulate(boundaries.begin(), boundaries.begin() + surfFixedDir, 0);
        }
        assert(index >= 0);
        return index;
      }
      __builtin_unreachable();
    }

    /** \brief Returns span indices  in the knot vectors and the fixed or free directions of the entity from the
     * directindex
     *
     * @tparam codim codim of the entity
     * @param directIndex index of the entity This index starts from 0 for each codim.
     * @return pair of spanindices and free or fixed direction array
     */
    template <std::integral auto codim>
    std::pair<std::array<int, dim>, std::array<Impl::FixedOrFree, dim>> spanAndDirectionFromDirectIndex(
        const int directIndex) const {
      std::array<int, dim> currentKnotSpan;
      std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
      if constexpr (codim == 0) {  // elements
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::free);
        const auto multiIndex = elementNet_->directToMultiIndex(directIndex);
        for (size_t i = 0; i < dim; i++)
          currentKnotSpan[i] = Dune::IGA::findSpanCorrected(patchData_->degree[i], uniqueKnotVector_[i][multiIndex[i]],
                                                            patchData_->knotSpans[i], multiIndex[i]);
      } else if constexpr (codim == dim) {  // vertex
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        const auto multiIndex = vertexNet_->directToMultiIndex(directIndex);

        for (size_t i = 0; i < dim; ++i) {
          if (Dune::FloatCmp::eq(uniqueKnotVector_[i][multiIndex[i]], patchData_->knotSpans[i].back())) {
            // If we are the vertex on the rightmost end of the knotspan, we set the knot span index by hand to the last
            // entry This is needed since findSpan return to use the end - degree -2 index, which is wrong for the
            // rightmost end
            currentKnotSpan[i] = patchData_->knotSpans[i].size() - 1;
          } else
            currentKnotSpan[i] = Dune::IGA::findSpanCorrected(
                patchData_->degree[i], uniqueKnotVector_[i][multiIndex[i]], patchData_->knotSpans[i], multiIndex[i]);
        }

      } else if constexpr (dim - codim == 1 && dim > 1)  // edge case
      {
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        std::array<int, dim> edgeSizes = sizeOfCodim1PerDirection<codim>();

        int edgeDir = 0, edgesPerDir = edgeSizes[0];
        for (int j = 0; j < dim; ++j)
          if (directIndex < edgesPerDir) {  // find out to which edge direction this edge index belongs
            edgeDir                 = j;
            fixedOrFreeDirection[j] = Impl::FixedOrFree::free;
            break;
          } else
            edgesPerDir += edgeSizes.at(j + 1);

        assert(((edgeDir < 2 && dim == 2) || (edgeDir < 3 && dim == 3)) && edgeDir >= 0);
        const int indexInDir
            = (directIndex - std::accumulate(edgeSizes.begin(), edgeSizes.begin() + edgeDir, 0)) % edgeSizes[edgeDir];
        std::array<int, dim> spanIndices{};
        spanIndices[edgeDir]
            = std::fmod(std::floor(indexInDir
                                   / (std::accumulate(uniqueSpanSize_.begin(), uniqueSpanSize_.begin() + edgeDir, 1,
                                                      [](auto res, auto r) { return res * (r + 1); }))),
                        uniqueSpanSize_[edgeDir]);

        if (edgeDir == 0) {
          spanIndices[1] = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[0])), uniqueSpanSize_[1] + 1);
          if constexpr (dim == 3)
            spanIndices[2] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0]) * (uniqueSpanSize_[1] + 1))),
                                       uniqueSpanSize_[2] + 1);
        } else {
          spanIndices[0] = std::fmod(indexInDir, uniqueSpanSize_[0] + 1);
          if (edgeDir == 1)
            if constexpr (dim == 3) {
              spanIndices[2] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[1]) * uniqueSpanSize_[0] + 1)),
                                         uniqueSpanSize_[2] + 1);
            }
          if constexpr (dim == 3)
            if (edgeDir == 2)
              spanIndices[1] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0] + 1))), uniqueSpanSize_[1] + 1);
        }
        for (size_t i = 0; i < dim; i++)
          if (fixedOrFreeDirection[i] == Impl::FixedOrFree::free) {
            currentKnotSpan[i] = Dune::IGA::findSpanCorrected(
                patchData_->degree[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
          } else {
            currentKnotSpan[i] = Dune::IGA::findSpanUncorrected(
                patchData_->degree[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
          }
      } else if constexpr (dim - codim == 2 && dim > 2)  // surface case
      {
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::free);
        std::array<int, dim> surfaceSizes = sizeOfCodim1PerDirection<codim>();

        int surfFixedDir = 0;
        int edgesPerDir  = surfaceSizes[0];
        for (int j = 0; j < dim; ++j)
          if (directIndex < edgesPerDir) {  // find out to which edge direction this edge index belongs
            surfFixedDir            = j;
            fixedOrFreeDirection[j] = Impl::FixedOrFree::fixed;
            break;
          } else
            edgesPerDir += surfaceSizes.at(j + 1);
        assert((surfFixedDir < 3 && dim == 3) && surfFixedDir >= 0);
        const int indexInDir
            = (directIndex - std::accumulate(surfaceSizes.begin(), surfaceSizes.begin() + surfFixedDir, 0))
              % surfaceSizes[surfFixedDir];
        std::array<int, dim> spanIndices{};
        spanIndices[surfFixedDir]
            = std::fmod(std::floor(indexInDir
                                   / (std::accumulate(uniqueSpanSize_.begin(), uniqueSpanSize_.begin() + surfFixedDir,
                                                      1, std::multiplies{}))),
                        uniqueSpanSize_[surfFixedDir] + 1);

        if (surfFixedDir == 0) {
          spanIndices[1] = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[0] + 1)), uniqueSpanSize_[1]);
          spanIndices[2]
              = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0] + 1) * uniqueSpanSize_[1])), uniqueSpanSize_[2]);
        } else {
          spanIndices[0] = std::fmod(indexInDir, uniqueSpanSize_[0]);
          if (surfFixedDir == 1)
            spanIndices[2] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[1] + 1) * uniqueSpanSize_[0])),
                                       uniqueSpanSize_[2]);
          else if (surfFixedDir == 2)
            spanIndices[1] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0]))), uniqueSpanSize_[1]);
        }
        for (size_t i = 0; i < dim; i++)
          currentKnotSpan[i] = Dune::IGA::findSpanCorrected(patchData_->degree[i], uniqueKnotVector_[i][spanIndices[i]],
                                                            patchData_->knotSpans[i], spanIndices[i]);
        for (size_t i = 0; i < dim; i++)
          if (fixedOrFreeDirection[i] == Impl::FixedOrFree::free) {
            currentKnotSpan[i] = Dune::IGA::findSpanCorrected(
                patchData_->degree[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
          } else {
            currentKnotSpan[i] = Dune::IGA::findSpanUncorrected(
                patchData_->degree[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
          }
      }
      return std::make_pair(currentKnotSpan, fixedOrFreeDirection);
    };

    /** \brief creates a NURBSGeometry object
     *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
     *  and generates a Geometry object
     *
     *  \param[in] ijk array of indices for each dimension
     */
    template <std::integral auto codim>
    typename GridImpl::template Codim<codim>::Geometry geometry(const RealIndex realIndex) const {
      DirectIndex directIndex = getDirectIndex<codim>(realIndex);

      auto [currentKnotSpan, fixedOrFreeDirection] = spanAndDirectionFromDirectIndex<codim>(directIndex);

      auto geo = (codim == 0 && trimData_ && isTrimmed(realIndex))
                     ? NURBSGeometry<dim - codim, dimworld, const GridImpl>(
                         patchData_, fixedOrFreeDirection, currentKnotSpan, subGridInfoMap.at(directIndex).subgrid)
                     : NURBSGeometry<dim - codim, dimworld, const GridImpl>(patchData_, fixedOrFreeDirection,
                                                                            currentKnotSpan);
      return typename GridImpl::template Codim<codim>::Geometry(geo);
    }

    std::pair<std::array<double, dim>, std::array<double, dim>> scalingAndOffset(const DirectIndex directIndex) {
      // This is a mirroring of nurbsGeometryCode
      std::array<double, dim> scaling;
      std::array<double, dim> offset;

      auto [thisSpanIndices, fixedOrFreeDirection] = spanAndDirectionFromDirectIndex<0>(directIndex);
      for (int i = 0; i < dim; ++i) {
        if (thisSpanIndices[i] + 1 < patchData_->knotSpans[i].size())
          scaling[i] = patchData_->knotSpans[i][thisSpanIndices[i] + 1] - patchData_->knotSpans[i][thisSpanIndices[i]];
        offset[i] = patchData_->knotSpans[i][thisSpanIndices[i]];
      }
      return std::make_pair(scaling, offset);
    }

    /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
    std::array<int, dim> validKnotSize() const { return uniqueSpanSize_; }

    Types typesInCodim(const int codim) const {
      assert(!std::isnan(n_fullElement));  // This would mean that the trim Code hasn't been run
      if (codim != 2)
        return {GeometryTypes::cube(codim)};
      else {
        if (n_fullElement == 0)
          return {GeometryTypes::none(codim)};
        else if (n_trimmedElement > 0)
          return {GeometryTypes::cube(codim), GeometryTypes::none(codim)};
        else
          return {GeometryTypes::cube(codim)};
      }
    }

    /// \brief returns false if an error happened during processing the trim
    [[nodiscard]] bool reportTrimError() const { return not trimErrorFlag_; }

   private:
    [[nodiscard]] int getGlobalVertexIndexFromElementIndex(const RealIndex realIndex, const int localVertexIndex,
                                                           const bool returnOriginal = false) const {
      DirectIndex elementDirectIndex = getDirectIndex<0>(realIndex);
      auto multiIndexOfVertex        = elementNet_->directToMultiIndex(elementDirectIndex);
      const auto bs                  = std::bitset<dim>(localVertexIndex);
      for (int i = 0; i < bs.size(); ++i)
        multiIndexOfVertex[i] += bs[i];
      if (returnOriginal) return vertexNet_->index(multiIndexOfVertex);
      return getRealIndex<2>(vertexNet_->index(multiIndexOfVertex));
    }

    template <int codim>
    std::array<int, dim> sizeOfCodim1PerDirection() const {
      std::array<int, dim> codim1Sizes;
      std::ranges::fill(codim1Sizes, 1);
      if constexpr (dim > 1) {
        for (int i = 0; i < dim; ++i)
          for (int k = 0; k < dim; ++k)
            if constexpr (dim - codim == 2)
              codim1Sizes[i] *= (i == k) ? uniqueSpanSize_[k] + 1 : uniqueSpanSize_[k];  // surfaces
            else
              codim1Sizes[i] *= (i == k) ? uniqueSpanSize_[k] : uniqueSpanSize_[k] + 1;  // edges
      } else
        codim1Sizes[0] = uniqueSpanSize_[0] + 1;  // vertices

      return codim1Sizes;
    }

    [[nodiscard]] int getGlobalEdgeIndexFromElementIndex(const int realIndex, const int eI,
                                                         const bool returnOriginal = false) const {
      DirectIndex elementDirectIndex = getDirectIndex<0>(realIndex);
      const auto eleI                = elementNet_->directToMultiIndex(elementDirectIndex);
      std::array<int, dim> edgeSizes;
      std::ranges::fill(edgeSizes, 1);
      for (int i = 0; i < dim; ++i)
        for (int k = 0; k < dim; ++k)
          edgeSizes[i] *= (i == k) ? uniqueSpanSize_[k] : uniqueSpanSize_[k] + 1;
      assert((dim == 2 && eI < 4) || (dim == 3 && eI < 12));
      int edgeDir;
      if constexpr (dim == 2)
        edgeDir = (eI == 0 || eI == 1) ? 1 : 0;  // 2,3 edgeDir =0
      else
        edgeDir
            = (eI >= 0 && eI < 4) ? 2 : ((eI == 4 || eI == 5 || eI == 8 || eI == 9) ? 1 : 0);  // 6,7,10,11 edgeDir=0

      assert((edgeDir < 2 && dim == 2) || (edgeDir < 3 && dim == 3));
      int dIndex
          = std::accumulate(edgeSizes.begin(), edgeSizes.begin() + edgeDir, 0);  // move index to correct subset
                                                                                 //    directIndex += eleI[edgeDir];
      if constexpr (dim == 3) {
        if (edgeDir == 0)
          dIndex += eleI[0] + (uniqueSpanSize_[0]) * (eleI[1])
                    + (uniqueSpanSize_[0]) * (uniqueSpanSize_[1] + 1) * (eleI[2]);
        else if (edgeDir == 1)
          dIndex += eleI[0] + ((uniqueSpanSize_[0] + 1) * eleI[1])
                    + ((uniqueSpanSize_[1]) * (uniqueSpanSize_[0] + 1)) * eleI[2];
        else if (edgeDir == 2)
          dIndex += eleI[0] + ((uniqueSpanSize_[0] + 1) * eleI[1])
                    + ((uniqueSpanSize_[1] + 1) * (uniqueSpanSize_[0] + 1)) * eleI[2];
      } else if constexpr (dim == 2) {
        if (edgeDir == 0)
          dIndex += eleI[0] + (uniqueSpanSize_[0]) * (eleI[1]);
        else
          dIndex += eleI[0] + ((uniqueSpanSize_[0] + 1) * eleI[1]);
      }

      switch (dim) {
        case 2:
          dIndex += (eI == 1) ? 1 : (eI == 3) ? uniqueSpanSize_[0] : 0;
          break;
        case 3:
          assert(uniqueSpanSize_.size() == dim);
          switch (eI) {
            case 1:
              dIndex += 1;
              break;
            case 2:
              dIndex += uniqueSpanSize_[0] + 1;
              break;
            case 3:
              dIndex += uniqueSpanSize_[0] + 2;
              break;
            case 5:
              dIndex += 1;
              break;
            case 8:
              dIndex += (uniqueSpanSize_[0] + 1) * uniqueSpanSize_[1];
              break;
            case 9:
              dIndex += (uniqueSpanSize_[0] + 1) * uniqueSpanSize_[1] + 1;
              break;
            case 7:
              dIndex += uniqueSpanSize_[0];
              break;
            case 10:
              dIndex += uniqueSpanSize_[0] * (uniqueSpanSize_[1] + 1);
              break;
            case 11:
              dIndex += uniqueSpanSize_[0] * (uniqueSpanSize_[1] + 2);
              break;
            default:  // edges 0,4,6
              break;
          }
      }
      if (returnOriginal) return dIndex;
      return getRealIndex<1>(dIndex);
    }

    [[nodiscard]] int getGlobalSurfaceIndexFromElementIndex(const int realIndex, const int eI) const {
      DirectIndex elementDirectIndex = getDirectIndex<0>(realIndex);
      auto eleI                      = elementNet_->directToMultiIndex(elementDirectIndex);

      std::array<int, dim> edgeSizes;
      std::ranges::fill(edgeSizes, 1);
      for (int i = 0; i < dim; ++i)
        for (int k = 0; k < dim; ++k)
          edgeSizes[i] *= (i == k) ? uniqueSpanSize_[k] + 1 : uniqueSpanSize_[k];
      assert((dim == 2 && eI < 4) || (dim == 3 && eI < 12));
      const int surfaceFixedDir = std::floor(eI / 2);
      int dIndex
          = std::accumulate(edgeSizes.begin(), edgeSizes.begin() + surfaceFixedDir, 0);  // move index to correct subset

      assert(surfaceFixedDir < 3 && dim == 3);

      if (surfaceFixedDir == 0)
        dIndex += eleI[0] + (uniqueSpanSize_[0] + 1) * (eleI[1])
                  + (uniqueSpanSize_[0] + 1) * (uniqueSpanSize_[1]) * (eleI[2]);
      else if (surfaceFixedDir == 1)
        dIndex += eleI[0] + (uniqueSpanSize_[0] * eleI[1]) + ((uniqueSpanSize_[1] + 1) * uniqueSpanSize_[0]) * eleI[2];
      else if (surfaceFixedDir == 2)
        dIndex += eleI[0] + (uniqueSpanSize_[0] * eleI[1]) + (uniqueSpanSize_[1] * uniqueSpanSize_[0]) * eleI[2];

      switch (eI) {
        case 1:
          ++dIndex;
          break;
        case 3:
          dIndex += uniqueSpanSize_[0];
          break;
        case 5:
          dIndex += uniqueSpanSize_[1] * uniqueSpanSize_[0];
          break;
        default:  // surface local indices 0,2,4
          break;
      }
      return dIndex;
    }

   public:
    template <typename ctype = double>
    auto parameterSpaceGrid(int scale = 0) const {
      using TensorSpaceCoordinates = TensorProductCoordinates<ctype, dim>;
      using ParameterSpaceGrid     = YaspGrid<dim, TensorSpaceCoordinates>;

      std::array<std::vector<double>, dim> uniqueKnots;
      for (int i = 0; i < dim; ++i)
        std::ranges::unique_copy(patchData_->knotSpans[i], std::back_inserter(uniqueKnots[i]));

      if constexpr (std::is_same<double, ctype>::value) return std::make_unique<ParameterSpaceGrid>(uniqueKnots);

      std::array<std::vector<ctype>, dim> uniqueKnotsTransformed;
      for (int i = 0; i < dim; ++i)
        std::ranges::transform(uniqueKnots[i], std::back_inserter(uniqueKnotsTransformed[i]),
                               [scale](const auto& knot) { return knot * static_cast<ctype>(std::pow(10, scale)); });

      return std::make_unique<ParameterSpaceGrid>(uniqueKnotsTransformed);
    }

    // In the following a DirectIndex is the flat index of the knot span in the tensor-product net of the untrimmed
    // patch The RealIndex refers to a flat index that runs from 0 ... j with j = n_T + n_F

    [[nodiscard]] bool hasEmptyElements() const { return (n_fullElement + n_trimmedElement < originalSize(0)); }
    [[nodiscard]] int emptyElements() const { return originalSize(0) - (n_fullElement + n_trimmedElement); }

    template <unsigned int codim>
    [[nodiscard]] auto getDirectIndex(RealIndex idx) const -> DirectIndex {
      return idx;
    }
    template <unsigned int codim>
    [[nodiscard]] auto getRealIndex(DirectIndex idx) const -> RealIndex {
      return idx;
    }
    template <unsigned int codim>
    requires(dim == 2) [[nodiscard]] auto getDirectIndex(RealIndex idx) const -> DirectIndex {
      return getEntityMap<codim>()[idx];
    }

    template <unsigned int codim>
    requires(dim == 2) [[nodiscard]] auto getRealIndex(DirectIndex idx) const -> RealIndex {
      constexpr auto orVal = std::numeric_limits<DirectIndex>::max();
      auto res             = getRealIndexOr<codim>(idx, orVal);
      if (res != orVal)
        return res;
      else
        throw std::runtime_error("No corresponding realIndex");
    }
    template <unsigned int codim>
    requires(dim == 2) [[nodiscard]] auto getRealIndexOr(auto idx, auto orValue) const noexcept
        requires(std::same_as<decltype(idx), decltype(orValue)>) {
      auto& map = this->getEntityMap<codim>();
      auto it   = std::ranges::find(map, idx);
      return it != map.end() ? std::ranges::distance(map.begin(), it) : orValue;
    }

    [[nodiscard]] auto getTrimFlagForDirectIndex(DirectIndex directIndex) const -> ElementTrimFlag {
      return trimFlags_[directIndex];
    }
    [[nodiscard]] auto getTrimFlag(RealIndex realIndex) const -> ElementTrimFlag {
      int nurbsIdx = getDirectIndex<0>(realIndex);
      return getTrimFlagForDirectIndex(nurbsIdx);
    }
    [[nodiscard]] bool isTrimmed(RealIndex realIndex) const {
      return getTrimFlag(realIndex) == ElementTrimFlag::trimmed;
    }

    [[nodiscard]] auto getElementSubGrid(RealIndex realIndex) const -> std::shared_ptr<SubGridType> {
      assert(isTrimmed(realIndex) && "You can only obtain trimmedElementRepresentations for trimmed Elements");

      auto directIndex = getDirectIndex<0>(realIndex);
      return subGridInfoMap.at(directIndex).subgrid;
    }

    void prepareForNoTrim() {
      auto n_ele = originalSize(0);
      trimFlags_ = std::vector<ElementTrimFlag>(n_ele);
      std::ranges::fill(trimFlags_, ElementTrimFlag::full);

      if constexpr (dim == 2) {
        fill1to1Maps<0>();
        fill1to1Maps<1>();
        fill1to1Maps<2>();
      }

      n_fullElement = n_ele;
    }

    template <unsigned int codim>
    requires(dim == 2) void fill1to1Maps() {
      int n_entity = originalSize(codim);
      auto& map    = getEntityMap<codim>();

      auto iView = std::views::iota(0, n_entity);
      map.reserve(n_entity);
      map.insert(map.end(), iView.begin(), iView.end());
    }

    // This is for grid dims where no trimming functionality is enabled
    void computeElementTrimInfo() { prepareForNoTrim(); }

    // Trim Info gets computed for dim == 2 && dimworld == 2
    void computeElementTrimInfo() requires(dim == 2) {
      if (!trimData_.has_value()) {
        prepareForNoTrim();
        return;
      }
      // Prepare Results
      n_fullElement    = 0;
      n_trimmedElement = 0;
      trimFlags_.resize(originalSize(0));

      using Trimmer = Trim::NURBSPatchTrimmer<>;

      auto paraGrid               = parameterSpaceGrid<Trimmer::intType>(Trimmer::scale);
      auto parameterSpaceGridView = paraGrid->leafGridView();

      Trimmer trimmer(parameterSpaceGridView, trimData_.value().get());

      const auto& indexSet = parameterSpaceGridView.indexSet();
      for (auto& element : elements(parameterSpaceGridView)) {
        DirectIndex directIndex = indexSet.index(element);
        RealIndex realIndex     = n_fullElement + n_trimmedElement;

        auto [trimFlag, boundaries, errorFlag] = trimmer.trimElement(directIndex);
        trimFlags_[directIndex]                = trimFlag;
        if (errorFlag) trimErrorFlag_ = true;

        if (trimFlag == ElementTrimFlag::trimmed) {
          elementIndexMap.push_back(directIndex);
          ++n_trimmedElement;
          subGridInfoMap.emplace(
              directIndex, ElementSubGridInfo{.realIndex = realIndex,
                                              .subgrid   = std::make_unique<SubGridType>(boundaries.value(),
                                                                                       scalingAndOffset(directIndex))});
        } else if (trimFlag == ElementTrimFlag::full) {
          elementIndexMap.push_back(directIndex);
          ++n_fullElement;
        }
      }

      constructSubEntityMaps<2>();
      constructSubEntityMaps<1>();
    }
    template <unsigned int codim>
    requires(codim == 2 || codim == 1) && (dim == 2) void constructSubEntityMaps() {
      int n_entities_original = originalSize(codim);
      std::set<DirectIndex> indicesOfEntityInTrim;

      for (auto realIndex : std::views::iota(0u, n_fullElement + n_trimmedElement)) {
        for (int i = 0; i < 4; ++i) {
          if constexpr (codim == 2)
            indicesOfEntityInTrim.insert(getGlobalVertexIndexFromElementIndex(realIndex, i, true));
          else
            indicesOfEntityInTrim.insert(getGlobalEdgeIndexFromElementIndex(realIndex, i, true));
        }
      }

      unsigned int realIndexCounter = 0;
      auto& map                     = getEntityMap<codim>();

      for (auto i : std::views::iota(0, n_entities_original)) {
        auto it = std::ranges::find(indicesOfEntityInTrim, i);
        if (it != indicesOfEntityInTrim.end()) {
          map.push_back(i);
          ++realIndexCounter;
        }
      }
    }

    // TODO Cache this
    [[nodiscard]] std::set<RealIndex> boundarySegmentList() const requires(dim == 2) {
      // This is the same functionality as in entity<0> where the intersection are made, maybe cache this and use it
      // when the elements are created as well
      constexpr int noNeighbor = -1;

      auto getRealIndexForOuterIndex = [&](int outerIndex) -> int {
        if (outerIndex == noNeighbor) return noNeighbor;
        return getRealIndexOr<0>(outerIndex, noNeighbor);
      };

      std::set<RealIndex> realIndexOfBoundaryIntersections;

      for (int i : std::views::iota(0, size(0))) {
        auto nurbsDirectIndex = getDirectIndex<0>(i);
        for (int innerLocalIndex = 0; innerLocalIndex < 4; ++innerLocalIndex) {
          auto multiIndex = elementNet_->directToMultiIndex(nurbsDirectIndex);
          multiIndex[static_cast<int>(std::floor(innerLocalIndex / 2))]
              += ((innerLocalIndex % 2)
                      ? 1
                      : noNeighbor);  // increase the multiIndex depending on where the outer element should lie
          auto directOuterIndex = (elementNet_->isValid(multiIndex)) ? elementNet_->index(multiIndex) : noNeighbor;
          directOuterIndex      = getRealIndexForOuterIndex(directOuterIndex);
          if (directOuterIndex == noNeighbor)
            realIndexOfBoundaryIntersections.insert(getGlobalEdgeIndexFromElementIndex(i, innerLocalIndex));
        }
      }
      return realIndexOfBoundaryIntersections;
    }

   private : template <typename GridImpl>
             friend class NURBSGridLeafIndexSet;
    template <typename GridImpl>
    friend class NURBSintersection;
    std::shared_ptr<NURBSPatchData<dim, dimworld, ScalarType>> patchData_;
    std::shared_ptr<NURBSPatchGeometry<dim, dimworld>> patchGeometry_;
    std::array<int, dim> uniqueSpanSize_;
    std::array<std::vector<double>, dim> uniqueKnotVector_;

    std::shared_ptr<MultiDimensionNet<dim, double>> elementNet_;
    std::shared_ptr<MultiDimensionNet<dim, double>> vertexNet_;

    std::optional<std::shared_ptr<TrimData>> trimData_;
    std::vector<ElementTrimFlag> trimFlags_;
    bool trimErrorFlag_{false};

    std::vector<DirectIndex> elementIndexMap;
    std::vector<DirectIndex> vertexIndexMap;
    std::vector<DirectIndex> edgeIndexMap;

    unsigned int n_fullElement{std::numeric_limits<unsigned int>::signaling_NaN()};
    unsigned int n_trimmedElement{0};

    struct ElementSubGridInfo {
      RealIndex realIndex{};
      std::shared_ptr<SubGridType> subgrid;
    };

    std::map<DirectIndex, ElementSubGridInfo> subGridInfoMap;

    template <unsigned int codim>
    requires(dim == 2) auto& getEntityMap() const {
      if constexpr (codim == 0)
        return elementIndexMap;
      else if constexpr (codim == 1)
        return edgeIndexMap;
      else
        return vertexIndexMap;
    }
    template <unsigned int codim>
    requires(dim == 2) auto& getEntityMap() {
      if constexpr (codim == 0)
        return elementIndexMap;
      else if constexpr (codim == 1)
        return edgeIndexMap;
      else
        return vertexIndexMap;
    }
  };

}  // namespace Dune::IGA
