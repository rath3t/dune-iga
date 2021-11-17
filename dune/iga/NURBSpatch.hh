#pragma once

#include "igaalgorithms.hh"

#include <memory>

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/controlpoint.hh>
#include <dune/iga/nurbsgeometry.hh>
#include <dune/iga/traits.hh>

namespace Dune::IGA {

  template <std::integral auto dim, std::integral auto dimworld,
            NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits = LinearAlgebraTraits<double, dim, dimworld>>
  class NURBSGrid;

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
    explicit NURBSPatch(const NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>& patchData)
        : NURBSPatch(patchData.getKnots(), patchData.getControlPoints(), patchData.getWeights(), patchData.getOrder()) {}

    NURBSPatch(const std::array<std::vector<double>, dim>& knotSpans,
               const typename NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>::ControlPointNetType controlPoints,
               const std::array<int, dim> order)
        : patchData_(std::make_shared<NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>>(knotSpans, controlPoints, order)) {
      validKnotSize_ = this->validKnotSize();
      // Build a knot net to make iterator operations easier
      // Here each "point" of the net is a element(knot span)
      knotElementNet_ = std::make_shared<MultiDimensionNet<dim, double>>(validKnotSize_);
    }

    int size(int codim) {
      assert(codim <= dim);
      assert(codim <= 3);

      if (codim == 0)
        return knotElementNet_->directSize();
      else if (codim == dim)
        return patchData_->controlPoints.directSize();
      else if (codim == 1 && dim == 2)  // edge case
      {
        int edgeSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? validKnotSize_[i] - 1 : validKnotSize_[i];
          edgeSize += subs;
        }
        return edgeSize;
      } else if (codim == 1 && dim == 3)  // surface case
      {
        int surfSize = 0;
        for (int j = 0; j < dim; ++j) {
          int subs = 1;
          for (int i = 0; i < dim; ++i)
            subs *= (i == j) ? validKnotSize_[i] : validKnotSize_[i] - 1;
          surfSize += subs;
        }
        return surfSize;
      }
    }

    const auto& getPatchData() { return patchData_; }

    // this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
    auto findSpanIndex(const std::array<int, dim>& ijk) const {
      const auto& knotSpans = patchData_->knotSpans;

      std::array<int, dim> index;
      std::fill(index.begin(), index.end(), 0);

      /*finds the working geometry object ijk
       *(working geometry objects are defined between 2 knots, where knot[i]<knot[i+1])*/
      for (int count, j = 0; j < dim; ++j) {
        count = 0;
        while (count <= ijk[j]) {
          if (index[j] == knotSpans[j].size()) break;

          if (knotSpans[j][index[j] + 1] > knotSpans[j][index[j]]) count++;

          ++index[j];
        }
      }
      return index;
    }

    bool isBorderElement(const int& id) {
      auto const& knotElementNet = this->knotElementNet_;
      auto const& multiIndex     = knotElementNet->directToMultiIndex(id);

      for (int i = 0; i < dim; ++i)
        if ((multiIndex[i] == knotElementNet->size()[i] - 1) || (multiIndex[i] == 0)) return true;

      return false;
    }

    /** \brief creates a NURBSGeometry object
     *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
     *  and generates a Geometry object
     *
     *  \param[in] ijk array of indices for each dimension
     */
    template <std::integral auto entDim = dim>
    auto geometry(const std::array<int, dim>& ijk, [[maybe_unused]] const int subIndex = 0) const {
      const auto& knotSpans = patchData_->knotSpans;

      std::array<int, dim> index = findSpanIndex(ijk);

      /*the iterator on each dim-knotspan for geometry ijk is stored in an array named corners*/
      std::array<std::vector<double>::const_iterator, dim> freeSpans;
      for (int i = 0; i < dim; ++i)
        freeSpans[i] = (knotSpans[i]).begin() + index[i] - 1;

      std::array<Impl::FixedOrFree, dim> fixedOrFreeDirection;
      std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::free);
      std::array<double, dim - entDim> fixedSpans;
      if constexpr (entDim == 0) {
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        switch (subIndex) {
          case 1:
            ++freeSpans[0];
          case 2:
            ++freeSpans[1];
          case 3:
            ++freeSpans[1];
            ++freeSpans[0];
          case 4:
            ++freeSpans[2];
          case 5:
            ++freeSpans[2];
            ++freeSpans[0];
          case 6:
            ++freeSpans[2];
            ++freeSpans[1];
          case 7:
            ++freeSpans[2];
            ++freeSpans[1];
            ++freeSpans[0];
          default:
            assert(subIndex == 0);
            break;
        }
        std::ranges::transform(freeSpans, fixedSpans.begin(), [](auto& vp) { return *vp; });
      } else if constexpr (entDim == 1 && dim > 1) {
        std::ranges::fill(fixedOrFreeDirection, Impl::FixedOrFree::fixed);
        // each direction has dim+1 edges
        const int freeDirection             = std::floor(subIndex / dim);
        fixedOrFreeDirection[freeDirection] = Impl::FixedOrFree::free;
        ++freeSpans[freeDirection];
        for (int fcounter=0, i = 0; i < dim; ++i) {
          if (i == freeDirection) continue;
          fixedSpans[fcounter++] = freeSpans[i];
        }

        //        }
        //        switch (subIndex) {
//          case 0:
//
//          case 1:subIndex
//            ++freeSpans[freeDirection];
//          case 2:
//
//          case 3:
//            ++freeSpans[freeDirection];
//          case 4:
//
//          case 5:
//            ++freeSpans[freeDirection];
//          case 6:
//
//          case 7:
//            ++freeSpans[freeDirection];
//          case 8:
//
//          case 9:
//            ++freeSpans[freeDirection];
//          case 10:
//
//          case 11:
//            ++freeSpans[freeDirection];
//          default:
//            assert(subIndex == 0);
//            break;
//        }
      }

      return NURBSGeometry<entDim, dimworld, dim, NurbsGridLinearAlgebraTraits>(patchData_, freeSpans, fixedSpans, fixedOrFreeDirection);
    }

    /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
    std::array<int, dim> validKnotSize() const {
      const auto& knotSpans = patchData_->knotSpans;
      std::array<int, dim> validknotsize;
      std::fill(validknotsize.begin(), validknotsize.end(), 0);

      for (int j = 0; j < dim; ++j) {
        for (int i = 0; i < knotSpans[j].size() - 1; ++i) {
          if (knotSpans[j][i + 1] > knotSpans[j][i]) ++validknotsize[j];
        }
      }

      return validknotsize;
    }

  private:
    std::shared_ptr<NURBSPatchData<dim, dimworld, NurbsGridLinearAlgebraTraits>> patchData_;
    std::array<int, dim> validKnotSize_;
    std::shared_ptr<MultiDimensionNet<dim, double>> knotElementNet_;
  };

}  // namespace Dune::IGA
