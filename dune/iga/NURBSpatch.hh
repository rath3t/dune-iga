#pragma once

#include "igaalgorithms.hh"

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/nurbsgeometry.hh>
#include <dune/iga/traits.hh>


namespace Dune::IGA {

  template <std::integral auto dim, std::integral auto dimworld,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double, dim, dimworld>>
  class NURBSGrid;

  template<typename GridImpl>
  class NURBSGridLeafIndexSet;

  template <typename Grid>
  class NURBSLeafGridView;

  /** \brief Class where the NURBS geometry can work on */
  template <std::integral auto dim, std::integral auto dimworld,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double, dim, dimworld>>
  class NURBSPatch {
  public:
    friend class NURBSLeafGridView<NURBSGrid<dim, dimworld>>;
    template <std::integral auto codim, class GridViewImp>
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

      // Build a knot net to make iterator operations easier
      // Here each "point" of the net is a element(knot span)
      elementNet_ = std::make_shared<MultiDimensionNet<dim, double>>(uniqueSpanSize_);

      std::array<int, dim> vertexSize;
      for (int i = 0; i < dim; ++i)
        vertexSize[i] = uniqueSpanSize_[i] + 1;

      vertexNet_ = std::make_shared<MultiDimensionNet<dim, double>>(vertexSize);
    }

    int size(int codim) {
      assert(codim <= dim);
      assert(codim <= 3);

      if (codim == 0)
        return elementNet_->directSize();
      else if (codim == dim)
        return vertexNet_->directSize();
      else if (codim == 1 && dim == 2)  // edge case
      {
        int edgeSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? uniqueSpanSize_[i]  : uniqueSpanSize_[i]+1;
          edgeSize += subs;
        }
        return edgeSize;
      } else if (codim == 1 && dim == 3)  // surface case
      {
        int surfSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? uniqueSpanSize_[i] +1: uniqueSpanSize_[i] ;
          surfSize += subs;
        }
        return surfSize;
      }

    }

    const auto& getPatchData() { return patchData_; }

    // this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
//    auto findSpanIndex(const std::array<int, dim>& ijk) const {
//      const auto& knotSpans = patchData_->knotSpans;
//
//      std::array<int, dim> index;
//      std::fill(index.begin(), index.end(), 0);
//
//      /*finds the working geometry object ijk
//       *(working geometry objects are defined between 2 knots, where knot[i]<knot[i+1])*/
//      for (int count, j = 0; j < dim; ++j) {
//        count = 0;
//        while (count <= ijk[j]) {
//          if (index[j] == knotSpans[j].size()) break;
//
//          if (knotSpans[j][index[j] + 1] > knotSpans[j][index[j]]) count++;
//
//          ++index[j];
//        }
//      }
//      return index;
//    }

    bool isBorderElement(const int& id) {
      auto const& knotElementNet = this->elementNet_;
      auto const& multiIndex     = knotElementNet->directToMultiIndex(id);

      for (int i = 0; i < dim; ++i)
        if ((multiIndex[i] == knotElementNet->size()[i] - 1) || (multiIndex[i] == 0)) return true;

      return false;
    }



//    template <std::integral auto codim>
//    NURBSGeometry<dim - codim, dimworld, dim, NurbsGridLinearAlgebraTraits> geometry(const int directIndex,
//                                                                                     [[maybe_unused]] const int subIndex) const {
//      throw std::logic_error("The specialzations should be called");
//      //      if constexpr (codim == dim)  // vertices
//      //      {
//      //        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
//      //        switch (subIndex) {
//      //          case 1:
//      //            ++freeSpans[0];
//      //          case 2:
//      //            ++freeSpans[1];
//      //          case 3:
//      //            ++freeSpans[1];
//      //            ++freeSpans[0];
//      //          case 4:
//      //            ++freeSpans[2];
//      //          case 5:
//      //            ++freeSpans[2];
//      //            ++freeSpans[0];
//      //          case 6:
//      //            ++freeSpans[2];
//      //            ++freeSpans[1];
//      //          case 7:
//      //            ++freeSpans[2];
//      //            ++freeSpans[1];
//      //            ++freeSpans[0];
//      //          default:
//      //            assert(subIndex == 0);
//      //            break;
//      //        }
//      //        std::ranges::transform(freeSpans, fixedSpans.begin(), [](auto& vp) { return *vp; });
//      //      } else if constexpr (dim - codim == 1 && dim > 1)  // edges in a higher dim grid
//      //      {
//      //        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
//      //        // each direction has dim+1 edges
//      //        const int freeDirection             = std::floor(subIndex / dim);
//      //        fixedOrFreeDirection[freeDirection] = Impl::FixedOrFree::free;
//      //        ++freeSpans[freeDirection];
//      //        for (int fcounter = 0, i = 0; i < dim; ++i) {
//      //          if (i == freeDirection) continue;
//      //          fixedSpans[fcounter++] = freeSpans[i];
//      //        }
//      //
//      //        //        }
//      //        //        switch (subIndex) {
//      //        //          case 0:
//      //        //
//      //        //          case 1:subIndex
//      //        //            ++freeSpans[freeDirection];
//      //        //          case 2:
//      //        //
//      //        //          case 3:
//      //        //            ++freeSpans[freeDirection];
//      //        //          case 4:
//      //        //
//      //        //          case 5:
//      //        //            ++freeSpans[freeDirection];
//      //        //          case 6:
//      //        //
//      //        //          case 7:
//      //        //            ++freeSpans[freeDirection];
//      //        //          case 8:
//      //        //
//      //        //          case 9:
//      //        //            ++freeSpans[freeDirection];
//      //        //          case 10:
//      //        //
//      //        //          case 11:
//      //        //            ++freeSpans[freeDirection];
//      //        //          default:
//      //        //            assert(subIndex == 0);
//      //        //            break;
//      //        //        }
//      //      }
//    }
    /** \brief creates a NURBSGeometry object
     *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
     *  and generates a Geometry object
     *
     *  \param[in] ijk array of indices for each dimension
     */
    template <std::integral auto codim>
    NURBSGeometry<dim-codim, dimworld, dim, NurbsGridLinearAlgebraTraits> geometry(const int directIndex) const {
      if constexpr (codim==0) {
        const auto multiIndex = elementNet_->directToMultiIndex(directIndex);
        const auto& knotSpans = patchData_->knotSpans;

//        std::array<int, dim> index = findSpanIndex(multiIndex);

        /*the iterator on each dim-knotspan for geometry ijk is stored in an array named corners*/
        std::array<std::vector<double>::const_iterator, dim> freeSpans;
        for (int i = 0; i < dim; ++i)
          freeSpans[i] = uniqueKnotVector_[i].begin() + multiIndex[i] ;

        std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::free);
        return NURBSGeometry<dim-codim, dimworld, dim, NurbsGridLinearAlgebraTraits>(patchData_, freeSpans, fixedOrFreeDirection);
      }
      else if (codim==dim)
      {
        const auto vertexSpanIndex = vertexNet_->directToMultiIndex(directIndex);
        std::array<std::vector<double>::const_iterator, dim> freeSpans;
        for (int i = 0; i < dim; ++i)
          freeSpans[i] = uniqueKnotVector_[i].begin() + vertexSpanIndex[i] ;
        std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        return NURBSGeometry<0UL, dimworld, dim, NurbsGridLinearAlgebraTraits>(patchData_, freeSpans, fixedOrFreeDirection);

      }


    }

    /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
    std::array<int, dim> validKnotSize() const { return uniqueSpanSize_; }

  private:
    auto getGlobalVertexIndexFromElementIndex(const int elementDirectIndex, const int localVertexIndex) const
    {
      auto multiIndexOfVertex = elementNet_->directToMultiIndex(elementDirectIndex);
      const auto bs = std::bitset<dim>(localVertexIndex);
      for(int i =0; i < bs.size();++i)
        multiIndexOfVertex[i]+= bs[i];

      return vertexNet_->index(multiIndexOfVertex);
    }
    template<typename GridImpl>
    friend class NURBSGridLeafIndexSet;
    std::shared_ptr<NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>> patchData_;
    std::array<int, dim> uniqueSpanSize_;
    std::array<std::vector<double>, dim> uniqueKnotVector_;
    std::shared_ptr<MultiDimensionNet<dim, double>> elementNet_;
    std::shared_ptr<MultiDimensionNet<dim, double>> vertexNet_;
  };

}  // namespace Dune::IGA
