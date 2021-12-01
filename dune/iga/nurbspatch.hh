#pragma once

#include "igaalgorithms.hh"

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/nurbsgeometry.hh>

namespace Dune::IGA {

  template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double>>
  class NURBSGrid;

  template <typename GridImpl>
  class NURBSGridLeafIndexSet;

  template <typename Grid>
  class NURBSLeafGridView;

  /** \brief Class where the NURBS geometry can work on */
  template <std::integral auto dim, std::integral auto dimworld,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double>>
  class NURBSPatch {
  public:
    using GridImpl = NURBSGrid<dim, dimworld, NurbsGridLinearAlgebraTraits>;
    friend class NURBSLeafGridView<GridImpl>;
    template <int codim, int dim1, typename GridImpl1>
    friend class NURBSGridEntity;

    /** \brief Constructor of NURBS from knots, control points, weights and order
     *  \param[in] knotSpans vector of knotSpans for each dimension
     *  \param[in] controlPoints a n-dimensional net of control points
     *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
     *  \param[in] order order of the NURBS structure for each dimension
     */
    NURBSPatch(const std::array<std::vector<double>, dim>& knotSpans,
               const typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>::ControlPointNetType controlPoints,
               const std::array<int, dim> order)
        : NURBSPatch(NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>(knotSpans, controlPoints, order)) {}

    explicit NURBSPatch(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>& patchData)
        : patchData_{std::make_shared<NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>>(patchData)} {
      for (int i = 0; i < dim; ++i)
        std::ranges::unique_copy(patchData_->knotSpans[i], std::back_inserter(uniqueKnotVector_[i]),
                                 [](auto& l, auto& r) { return Dune::FloatCmp::eq(l, r); });

      for (int i = 0; i < dim; ++i)
        uniqueSpanSize_[i] = uniqueKnotVector_[i].size() - 1;

      elementNet_ = std::make_shared<MultiDimensionNet<dim, double>>(uniqueSpanSize_);

      std::array<int, dim> vertexSize;
      for (int i = 0; i < dim; ++i)
        vertexSize[i] = uniqueSpanSize_[i] + 1;

      vertexNet_ = std::make_shared<MultiDimensionNet<dim, double>>(vertexSize);
    }

    [[nodiscard]] int size(int codim) const {
      assert(codim <= dim);
      assert(codim <= 3);
      assert(codim >= 0);

      if (codim == 0)
        return elementNet_->directSize();
      else if (codim == dim)
        return vertexNet_->directSize();
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

    const auto& getPatchData() { return patchData_; }

    bool isPatchBoundary(const int& id) {
      auto const& knotElementNet = this->elementNet_;
      auto const& multiIndex     = knotElementNet->directToMultiIndex(id);

      for (int i = 0; i < dim; ++i)
        if ((multiIndex[i] == knotElementNet->size()[i] - 1))
          return true;
        else if ((multiIndex[i] == 0))
          return true;

      return false;
    }

    int patchBoundaryIndex(const int& id) {
      auto const& knotElementNet = this->elementNet_;
      auto const& multiIndex     = knotElementNet->directToMultiIndex(id);
      int index=0;
      if constexpr (dim==1) {
        if ((multiIndex[0] == 0))
          return 0;
        else
          return 1;
      } else if constexpr (dim==2) {
        if ((multiIndex[0] == 0))
          return multiIndex[1];
        else if ((multiIndex[0] == knotElementNet->size()[0] - 1))
          return multiIndex[0]*multiIndex[1]+multiIndex[1];
        else if ((multiIndex[0] == 0))
          return multiIndex[1];
        else if ((multiIndex[0] == knotElementNet->size()[0] - 1))
          return multiIndex[0]*multiIndex[1]+multiIndex[1];
      }

      for (int dimCounter = 0, i = 0; i < dim; ++i) {
        if ((multiIndex[0] == knotElementNet->size()[0] - 1))
          return dimCounter + 1;
        else if ((multiIndex[0] == 0))
          return dimCounter;
        dimCounter += 2;
      }
      return -1;
    }

    /** \brief creates a NURBSGeometry object
     *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
     *  and generates a Geometry object
     *
     *  \param[in] ijk array of indices for each dimension
     */
    template <std::integral auto codim>
    NURBSGeometry<dim - codim, dimworld, GridImpl> geometry(const int directIndex) const {
      if constexpr (codim == 0) {  // elements
        const auto multiIndex = elementNet_->directToMultiIndex(directIndex);
        std::array<int, dim> currentKnotSpan;
        for (size_t i = 0; i < dim; i++)
          currentKnotSpan[i]
              = Dune::IGA::findSpan(patchData_->order[i], uniqueKnotVector_[i][multiIndex[i]], patchData_->knotSpans[i], multiIndex[i]);

        std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::free);
        return NURBSGeometry<dim - codim, dimworld, GridImpl>(patchData_, fixedOrFreeDirection, currentKnotSpan);
      } else if constexpr (codim == dim) {  // vertex
        const auto vertexSpanIndex = vertexNet_->directToMultiIndex(directIndex);
        std::array<int, dim> currentKnotSpan;
        for (size_t i = 0; i < dim; ++i)
          currentKnotSpan[i] = Dune::IGA::findSpan(patchData_->order[i], uniqueKnotVector_[i][vertexSpanIndex[i]], patchData_->knotSpans[i],
                                                   vertexSpanIndex[i]);

        std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        return NURBSGeometry<0, dimworld, GridImpl>(patchData_, fixedOrFreeDirection, currentKnotSpan);
      } else if constexpr (dim - codim == 1 && dim > 1)  // edge case
      {
        std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        std::array<int, dim> edgeSizes;
        std::ranges::fill(edgeSizes, 1);
        for (int i = 0; i < dim; ++i)
          for (int k = 0; k < dim; ++k)
            edgeSizes[i] *= (i == k) ? uniqueSpanSize_[k] : uniqueSpanSize_[k] + 1;

        int edgeDir = 0, edgesPerDir = edgeSizes[0];
        for (int j = 0; j < dim; ++j)
          if (directIndex < edgesPerDir) {  // find out to which edge direction this edge index belongs
            edgeDir                 = j;
            fixedOrFreeDirection[j] = Impl::FixedOrFree::free;
            break;
          } else
            edgesPerDir += edgeSizes.at(j + 1);

        assert(((edgeDir < 2 && dim == 2) || (edgeDir < 3 && dim == 3)) && edgeDir >= 0);
        const int indexInDir = (directIndex - std::accumulate(edgeSizes.begin(), edgeSizes.begin() + edgeDir, 0)) % edgeSizes[edgeDir];

        std::array<int, dim> spanIndices{};
        spanIndices[edgeDir] = indexInDir % (uniqueSpanSize_[edgeDir]);
        if constexpr (dim == 3) {
          const bool dirFac = !(edgeDir == 0 || edgeDir == 2);
          spanIndices[(edgeDir + 1 + dirFac) % dim]
              = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[edgeDir])), uniqueSpanSize_[(edgeDir + 1 + dirFac) % dim] + 1);
          spanIndices[(edgeDir + 1 + !dirFac) % dim]
              = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[edgeDir] * (uniqueSpanSize_[(edgeDir + 1 + dirFac) % dim] + 1))),
                          uniqueSpanSize_[(edgeDir + 1 + !dirFac) % dim] + 1);
        } else if constexpr (dim == 2)
          spanIndices[(edgeDir + 1) % dim]
              = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[edgeDir])), uniqueSpanSize_[(edgeDir + 1) % dim] + 1);
        std::array<int, dim> currentKnotSpan;
        for (size_t i = 0; i < dim; i++)
          currentKnotSpan[i]
              = Dune::IGA::findSpan(patchData_->order[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
        return NURBSGeometry<1, dimworld, GridImpl>(patchData_, fixedOrFreeDirection, currentKnotSpan);
      } else if constexpr (dim - codim == 2 && dim > 2)  // surface case
      {
          std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
          std::ranges::fill(fixedOrFreeDirection,Impl::FixedOrFree::free);
          std::array<int, dim> edgeSizes;
          std::ranges::fill(edgeSizes, 1);
          for (int i = 0; i < dim; ++i)
            for (int k = 0; k < dim; ++k)
              edgeSizes[i] *= (i == k) ? uniqueSpanSize_[k] +1: uniqueSpanSize_[k];

          int  surfFixedDir=0;
          int   edgesPerDir = edgeSizes[0];
          for (int j = 0; j < dim; ++j)
            if (directIndex < edgesPerDir) {  // find out to which edge direction this edge index belongs
              surfFixedDir                 = j;
              fixedOrFreeDirection[j] = Impl::FixedOrFree::fixed;
              break;
            } else
              edgesPerDir += edgeSizes.at(j + 1);
          std::cout<<directIndex<<" "<<surfFixedDir <<" ";
          assert((surfFixedDir < 3 && dim == 3) && surfFixedDir >= 0);
          const int indexInDir = (directIndex - std::accumulate(edgeSizes.begin(), edgeSizes.begin()+surfFixedDir , 0)) % edgeSizes[surfFixedDir];
          std::cout<<indexInDir<<" || ";
          std::array<int, dim> spanIndices{};
          spanIndices[surfFixedDir] = std::fmod(std::floor(indexInDir / (std::accumulate(uniqueSpanSize_.begin(), uniqueSpanSize_.begin()+surfFixedDir , 1,std::multiplies{}))),uniqueSpanSize_[surfFixedDir]+1);

          if (surfFixedDir == 0)
          {spanIndices[1] = std::fmod(std::floor(indexInDir / (uniqueSpanSize_[0]+1)),uniqueSpanSize_[1]);
            spanIndices[2] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0]+1)*uniqueSpanSize_[1])),uniqueSpanSize_[2]);
          }
          else if (surfFixedDir !=0)
          { spanIndices[0] = std::fmod(indexInDir ,uniqueSpanSize_[0]);
            if (surfFixedDir == 1)
              spanIndices[2] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[1]+1)*uniqueSpanSize_[0])),uniqueSpanSize_[2]);
            else if (surfFixedDir== 2)
              spanIndices[1] = std::fmod(std::floor(indexInDir / ((uniqueSpanSize_[0]))),uniqueSpanSize_[1]);
          }
          std::array<int, dim> currentKnotSpan;
          for (size_t i = 0; i < dim; i++)
            currentKnotSpan[i]
                = Dune::IGA::findSpan(patchData_->order[i], uniqueKnotVector_[i][spanIndices[i]], patchData_->knotSpans[i], spanIndices[i]);
          return NURBSGeometry<2, dimworld, GridImpl>(patchData_, fixedOrFreeDirection, currentKnotSpan);
      }
    }

    /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
    std::array<int, dim> validKnotSize() const { return uniqueSpanSize_; }

  private:
    auto getGlobalVertexIndexFromElementIndex(const int elementDirectIndex, const int localVertexIndex) const {
      auto multiIndexOfVertex = elementNet_->directToMultiIndex(elementDirectIndex);
      const auto bs           = std::bitset<dim>(localVertexIndex);
      for (int i = 0; i < bs.size(); ++i)
        multiIndexOfVertex[i] += bs[i];
      return vertexNet_->index(multiIndexOfVertex);
    }

    int getGlobalEdgeIndexFromElementIndex(const int elementDirectIndex, const int eI) {
      const auto eleI = elementNet_->directToMultiIndex(elementDirectIndex);
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
        edgeDir = (eI >= 0 && eI < 4) ? 2 : ((eI == 4 || eI == 5 || eI == 8 || eI == 9) ? 1 : 0);  // 6,7,10,11 edgeDir=0

      assert((edgeDir < 2 && dim == 2) || (edgeDir < 3 && dim == 3));
      int dIndex = std::accumulate(edgeSizes.begin(), edgeSizes.begin() + edgeDir, 0);  // move index to correct subset
      dIndex += eleI[edgeDir];
      if constexpr (dim == 3) {
        const bool df = (edgeDir == 2 || edgeDir == 0);
        dIndex += uniqueSpanSize_[edgeDir]
                  * (eleI[(edgeDir + df + 1) % dim] * (uniqueSpanSize_[(edgeDir + !df + 1) % dim] + 1) + eleI[(edgeDir + !df + 1) % dim]);
      } else if constexpr (dim == 2)
        dIndex += uniqueSpanSize_[edgeDir] * eleI[((edgeDir + 1) % dim)];

      if constexpr (dim == 2)
        dIndex += (eI == 1) ? uniqueSpanSize_[1] : (eI == 3) ? uniqueSpanSize_[0] : 0;
      else if constexpr (dim == 3)
        switch (eI) {
          case 1: dIndex += uniqueSpanSize_[2]; break;
          case 2:
          case 3: dIndex += (uniqueSpanSize_[0] + (eI - 1) % 6) * uniqueSpanSize_[2]; break;
          case 5: dIndex += uniqueSpanSize_[1]; break;
          case 8:
          case 9: dIndex += (uniqueSpanSize_[0] + (eI - 1) % 6) * uniqueSpanSize_[1]; break;
          case 7: dIndex += uniqueSpanSize_[0]; break;
          case 10:
          case 11: dIndex += uniqueSpanSize_[0] * (uniqueSpanSize_[1] + 1 + eI % 2); break;
          default:  // edge local indices 0,4,6
            break;
        }
      return dIndex;
    }

    int  getGlobalSurfaceIndexFromElementIndex(const int elementDirectIndex, const int eI) {
      auto eleI = elementNet_->directToMultiIndex(elementDirectIndex);

      std::array<int, dim> edgeSizes;
      std::ranges::fill(edgeSizes, 1);
      for (int i = 0; i < dim; ++i)
        for (int k = 0; k < dim; ++k)
          edgeSizes[i] *= (i == k) ? uniqueSpanSize_[k]+1 : uniqueSpanSize_[k];
      assert((dim == 2 && eI < 4) || (dim == 3 && eI < 12));
      const int  surfaceFixedDir = std::floor(eI/2);
      int dIndex = std::accumulate(edgeSizes.begin(), edgeSizes.begin() + surfaceFixedDir, 0); // move index to correct subset

      assert( surfaceFixedDir < 3 && dim == 3);

      if(surfaceFixedDir==0)
        dIndex += eleI[0]+(uniqueSpanSize_[0]+1)*(eleI[1])+(uniqueSpanSize_[0]+1)*(uniqueSpanSize_[1])*(eleI[2]);
      else if(surfaceFixedDir==1)
        dIndex += eleI[0]+(uniqueSpanSize_[0]*eleI[1]) +((uniqueSpanSize_[1]+1)*uniqueSpanSize_[0])*eleI[2];//(uniqueSpanSize_[0]+1)*(eleI[0]);
      else if(surfaceFixedDir==2)
        dIndex += eleI[0]+(uniqueSpanSize_[0]*eleI[1]) +(uniqueSpanSize_[1]*uniqueSpanSize_[0])*eleI[2];//(uniqueSpanSize_[0]+1)*(eleI[0]);

      switch (eI) {
        case 1: ++dIndex; break;
        case 3: dIndex += uniqueSpanSize_[0]; break;
        case 5: dIndex += uniqueSpanSize_[1]*uniqueSpanSize_[0]; break;
        default:  // edge local indices 0,2,4
          break;
      }
      return dIndex;
    }

    template <typename GridImpl>
    friend class NURBSGridLeafIndexSet;
    std::shared_ptr<NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>> patchData_;
    std::array<int, dim> uniqueSpanSize_;
    std::array<std::vector<double>, dim> uniqueKnotVector_;
    std::shared_ptr<MultiDimensionNet<dim, double>> elementNet_;
    std::shared_ptr<MultiDimensionNet<dim, double>> vertexNet_;
  };

}  // namespace Dune::IGA
