//
// Created by Henri on 17.03.2023.
//

#pragma once

namespace Dune::IGA {

  template <std::integral auto dim, std::integral auto dimworld, LinearAlgebra NurbsGridLinearAlgebraTraitsImpl =  DuneLinearAlgebraTraits<double>>
  class NURBSPatchGeometry {
   public:
    static constexpr std::integral auto patchDim = dim;
    static constexpr std::integral auto coorddimension = dimworld;

    using LinearAlgebraTraits = NurbsGridLinearAlgebraTraitsImpl;
    using ctype               = typename LinearAlgebraTraits::value_type;
    using LocalCoordinate     = typename LinearAlgebraTraits::template FixedVectorType<patchDim>;
    using GlobalCoordinate    = typename LinearAlgebraTraits::template FixedVectorType<coorddimension>;
    using JacobianTransposed  = typename LinearAlgebraTraits::template FixedMatrixType<patchDim, coorddimension>;
    using JacobianInverseTransposed =
        typename LinearAlgebraTraits::template FixedMatrixType<coorddimension, patchDim>;

    using ControlPointType = typename NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>::ControlPointType;
    using MultiDimNet = MultiDimensionNet<patchDim, ControlPointType>;

    NURBSPatchGeometry() = default;

    explicit NURBSPatchGeometry(std::shared_ptr<NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>> patchData) : patchData_(patchData) {
      nurbs_ = Nurbs<patchDim, LinearAlgebraTraits>(*patchData);
    }


    [[nodiscard]] FieldVector<ctype, dimworld> global(const LocalCoordinate& local) const {
      // Bad workaround -> into array<double, dim>
      std::array<double, dim> u{};
      for (int i = 0; i < dim; ++i)
        u[i] = local[i];

      auto spanIndex = findSpanUncorrected(patchData_->degree, u, patchData_->knotSpans);

      auto cpCoordinateNet = netOfSpan(spanIndex, patchData_->degree, extractControlCoordinates(patchData_->controlPoints));
      auto basis     = nurbs_.basisFunctionNet(u);

      return Dune::IGA::dot(basis, cpCoordinateNet);

    }

   private:
    std::shared_ptr<NURBSPatchData<patchDim, dimworld, LinearAlgebraTraits>> patchData_;
    Dune::IGA::Nurbs<patchDim, LinearAlgebraTraits> nurbs_;

  };

}


