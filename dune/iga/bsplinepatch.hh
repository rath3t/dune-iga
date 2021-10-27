// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <memory>
#include <dune/geometry/type.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/iga/multidimensionNet.hh>
#include <dune/iga/concepts.hh>
#include <dune/iga/traits.hh>
#include <dune/iga/controlpoint.hh>


  namespace Dune::IGA
  {

    /** \brief class that holds all data regarding the B-Spline stucture */
    template<int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits=LinearAlgebraTraits<double,dim,dimworld>>
    class BsplinePatchData
    {
    public:
      using GlobalCoordinateType = typename NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
      using LocalCoordinateType = typename NurbsGridLinearAlgebraTraits::LocalCoordinateType;
      using JacobianTransposedType = typename NurbsGridLinearAlgebraTraits::JacobianTransposedType;
      using JacobianInverseTransposed = typename NurbsGridLinearAlgebraTraits::JacobianInverseTransposed;

      using ControlPointType = GlobalCoordinateType;


      using ControlPointNetType = MultiDimensionNet<dim,ControlPointType>;


      /** \brief  constructor for a BSplinePatchData from knots, control points and order
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      BsplinePatchData(const std::array<std::vector<double>,dim>& knotSpans,
                       const ControlPointNetType& controlPoints,
                       const std::array<int,dim> order)
      : knotSpans_(knotSpans)
      , controlPoints_(controlPoints)
      , order_(order)
      {
      }

      /** \brief returns the Knot Span*/
      const std::array<std::vector<double>, dim> & getKnots() const
      {
        return knotSpans_;
      }

      /** \brief returns the Control Point*/
      const ControlPointNetType & getControlPoints () const
      {
        return controlPoints_;
      }

      /** \brief returns the order*/
      const std::array<int,dim> & getOrder() const
      {
        return order_;
      }
    private:

      const std::array<std::vector<double>,dim>& knotSpans_;

      const ControlPointNetType  controlPoints_;

      const std::array<int,dim> order_;

    };

    /** \brief a geometry implementation for B-Splines */
    template <int dim, int dimworld>
    class BSplineElementGeometry {
    public:

      /** coordinate type */
      typedef double ctype;

      /** \brief Dimension of the cube element */
      enum {mydimension = dim};

      /** \brief Dimension of the world space that the cube element is embedded in*/
      enum {coorddimension = dimworld};

      /** \brief Type used for a vector of element coordinates */
      typedef FieldVector<ctype,dim> LocalCoordinate;

      /** \brief Type used for a vector of world coordinates */
      typedef FieldVector<ctype,dimworld> GlobalCoordinate;

      /** \brief Type for the transposed Jacobian matrix */
      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

      /** \brief Type for the transposed inverse Jacobian matrix */
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

    private:
      /* Helper class to compute a matrix pseudo inverse */
      typedef MultiLinearGeometryTraits<ctype>::MatrixHelper MatrixHelper;

    public:
      /** \brief Constructor of BSplineElementGeometry from BSplinePatchData and an iterator to a specific knot
       *
       *  \param[in] Patchdata shared pointer to an object where the all the data of the BSplinePatch is stored
       *  \param[in] corner Pointer (for each dimension) to the Knot span where the Geometry object is supposed to operate
       */
      BSplineElementGeometry(std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata,
                      std::array<std::vector<double>::const_iterator,dim> corner)
      : patchData_(Patchdata)
      , corner_(corner)
      {
        // note it`s maybe not such a good a idea to store all of that!
        std::array<std::vector<std::vector<double>>,dim> _basis;
      }

      /** \brief Map the center of the element to the geometry */
      GlobalCoordinate center() const
      {
        LocalCoordinate localcenter;
        std::fill(localcenter.begin(), localcenter.end(), 0.5);
        return global(localcenter);
      }

      /** \brief Return the number of corners of the element */
      int corners() const
      {
        return 1<<mydimension;
      }

      /** \brief Return world coordinates of the k-th corner of the element */
      GlobalCoordinate corner(int k) const
      {
        LocalCoordinate localcorner;
        for (size_t i=0; i<mydimension; i++) {
            localcorner[i] = (k & (1<<i)) ? 0 : 1;
        }
        return global(localcorner);
      }

      /** \brief I think it is not an affine mapping */
      bool affine () const { return false; }

      /** \brief Type of the element: a hypercube of the correct dimension */
      GeometryType type() const
      {
        return {GeometryType::cube,dim};
      }

      /** \brief evaluates the B-Spline mapping
       *
       *  \param[in] local local coordinates for each dimension
       */
      GlobalCoordinate global(const LocalCoordinate& local) const
      {
        const auto& knotSpans = patchData_->getKnots();
        const auto& controlPoints = patchData_->getControlPoints();
        const auto& order = patchData_->getOrder();
        const auto& basis  = this->getBasis(local);
        GlobalCoordinate glob;
        std::fill(glob.begin(), glob.end(), 0.0);
        std::array<unsigned int,dim> multiIndexBasisfunction, multiIndexControlNet;
        std::array<unsigned int,dim> dimsize;
        std::array<unsigned int,dim> cornerIdx;
        for (unsigned int d=0; d<dim; ++d)
        {
          dimsize[d] = order[d]+1;
          cornerIdx[d] = corner_[d]-knotSpans[d].begin();
        }
        /*Index net for valid basis functions*/
        auto basisFunctionNet = MultiDimensionNet<dim,FieldVector<double,dimworld>>(dimsize);
        for (unsigned int i=0; i<basisFunctionNet.directSize(); ++i)
        {
            multiIndexBasisfunction =  basisFunctionNet.directToMultiIndex(i);
          double temp = 1;
          for (unsigned int d=0; d<dim; ++d)
          {
              temp *= basis[d][multiIndexBasisfunction[d]][order[d]];
              multiIndexControlNet[d] = multiIndexBasisfunction[d] + cornerIdx[d] - order[d];
          }
          auto cp = controlPoints.get(multiIndexControlNet);
          for (unsigned int k=0; k<dimworld; ++k)
            glob[k] += cp[k]*temp;
        }

        return glob;
      }

      /** \brief compute the Jacobian transposed matrix
       *
       *  \param[in] local local coordinates for each dimension
       */
      JacobianTransposed jacobianTransposed(const LocalCoordinate &local)
      {
        for (int i=0; i<dim;++i)
          if (local[i]<0 || local[i]>1)
            DUNE_THROW(RangeError, "Local coordinates have to be in [0,1]^dim!");

        JacobianTransposed result;
        const auto& knotSpans = patchData_->getKnots();
        const auto& controlPoints = patchData_->getControlPoints();
        const auto& order = patchData_->getOrder();
        const auto& basis = this->getBasis(local);

        //const auto& basis = (basis_.empty())?this->getBasis(local):basis_;
        const auto& derBasis = this->getDerBasis();
        std::array<unsigned int,mydimension> multiIndexBasisfunction, multiIndexControlNet;
        std::array<unsigned int,mydimension> dimsize;
        std::array<unsigned int,mydimension> cornerIdx;
        for (unsigned int d=0; d<mydimension; ++d)
        {
          dimsize[d] = order[d]+1;
          cornerIdx[d] = corner_[d]-knotSpans[d].begin();
        }
        /*Index net for valid basis functions*/
        auto basisFunctionNet = MultiDimensionNetFVd<mydimension,coorddimension>(dimsize);
        /*Iterate the matrix row*/
        for (unsigned int r=0; r<mydimension; ++r)
        {
          std::fill(result[r].begin(),result[r].end(),0);
          /*Iterate the non-vanishing basis function grid*/
          for (unsigned int i=0; i<basisFunctionNet.directSize(); ++i)
          {
              multiIndexBasisfunction =  basisFunctionNet.directToMultiIndex(i);
            double product = 1;
            for (unsigned int d=0; d<mydimension; ++d)
            {
              double temp = (d == r) ? derBasis[d][multiIndexBasisfunction[d]] : basis[d][multiIndexBasisfunction[d]][order[d]];
              product *= temp;
              multiIndexControlNet[d] = multiIndexBasisfunction[d] + cornerIdx[d] - order[d];
            }
            auto cp = controlPoints.get(multiIndexControlNet);
            /* Fill the matrix row */
            for (unsigned int k=0; k<coorddimension; ++k)
            {
              result[r][k] += cp[k]*product;
            }

          }
        }
        return result;

      }

      /** \brief obtain the integration element
       *
       *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
       *  integration element \f$\mu(x)\f$ is given by
       *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
       *
       *  \param[in]  local  local coordinate to evaluate the integration element in
       *
       *  \returns the integration element \f$\mu(x)\f$.
       *
       */
      ctype integrationElement(const LocalCoordinate &local) const
      {
        return MatrixHelper::template sqrtDetAAT< mydimension,coorddimension>(jacobianTransposed(local));
      }

//      /** \brief computes the first-order derivative of the B-Spline with respect to the nth variable
//       *
//       *  \param[in] local local coordinates for each dimension
//       */
//       const std:std::vector<double>,mydimension> getDerivative(int n)
//       {
//
//
//       }

      /** \brief evaluates the non-zero derivatives of basis functions, order+1 for each dimension
       *
       *  \param[in] local local coordinates for each dimension
       */
      std::array<std::vector<double>,mydimension> getDerBasis()
      {
        const auto & order = patchData_->getOrder();
        const auto & knotSpans = patchData_->getKnots();
        std::array<std::vector<double>, mydimension> der;

        /*redo basis calculation or not? */

        for (unsigned int d=0; d<mydimension; ++d)
        {
          der[d].resize(order[d]+1);
          /*The first one , what if it doesn't exist*/
          der[d][0] = -order[d]*basis_[d][0][order[d]-1]/basis_[d][order[d]][0];
          for (unsigned int j=1; j<order[d]; ++j)
          {
            der[d][j] = order[d]*(basis_[d][j-1][order[d]-1]/basis_[d][order[d]][j-1]-basis_[d][j][order[d]-1]/basis_[d][order[d]][j]);
          }
          /*The last one, what if it doesn't exist */
          der[d][order[d]] = order[d]*basis_[d][order[d]-1][order[d]-1]/basis_[d][order[d]][order[d]-1];
        }
        return der;
      }


      /** \brief evaluates the non-zero basis Functions at the given local coordinates, order+1 for each dimension
       *
       *  \param[in] local local coordinates for each dimension
       */
      const std::array<std::vector<std::vector<double>>,mydimension > & getBasis (const LocalCoordinate& local) const
      {
        for (int i=0; i<dim; ++i)
          if (local[i]<0 || local[i]>1)
            DUNE_THROW(RangeError, "Local coordinates have to be in [0,1]^dim!");

        const auto & order = patchData_->getOrder();
        const auto & knotSpans = patchData_->getKnots();

        /*Compute the non-vanishing basis function. See "The NURBS book"*/
        for (int d=0; d<dim; ++d)
        {
          double loc;
          loc = local[d]*(*(corner_[d]+1)-*corner_[d]) + *corner_[d];
          basis_[d].resize(order[d]+1);
          std::vector<double> left, right;
          left.resize(order[d]+1);
          right.resize(order[d]+1);
          basis_[d][0].resize(order[d]+1);
          basis_[d][0][0] = 1;
          /*Build a squared table to store the knot differences as well to make derivative computation easier */
          for (int o=1; o<=order[d]; ++o)
          {
            basis_[d][o].resize(order[d]+1);
            left[o] = loc-*(corner_[d]+1-o);
            right[o] = *(corner_[d]+o)-loc;
            double saved = 0;
            for (int r=0; r<o; ++r)
            {
              /*Lower triangle - knot differences */
              basis_[d][o][r] = right[r+1]+left[o-r];
              double temp = basis_[d][r][o-1]/basis_[d][o][r];

              /*Upper triangle - basis function*/
              basis_[d][r][o] = saved+right[r+1]*temp;
              saved = left[o-r]*temp;
            }
            basis_[d][o][o] = saved;
          }
        }

        return basis_;
      }

    private:

        std::shared_ptr <BsplinePatchData<dim,dimworld>> patchData_;
        std::array<std::vector<double>::const_iterator,dim> corner_;
        mutable std::array<std::vector<std::vector<double>>,dim> basis_;
    };


    template<int dim, int dimworld>
    class BSplineLeafGridView;
    /** \brief Class where the B-Spline geometry can work on  */
    template <int dim, int dimworld, NurbsGridLinearAlgebra NurbsGridLinearAlgebraTraits=LinearAlgebraTraits<double,dim,dimworld>>
    class BSplinePatch
    {
    public:
      using GlobalCoordinateType = typename NurbsGridLinearAlgebraTraits::GlobalCoordinateType;
      using LocalCoordinateType = typename NurbsGridLinearAlgebraTraits::LocalCoordinateType;
      using JacobianTransposedType = typename NurbsGridLinearAlgebraTraits::JacobianTransposedType;
      using JacobianInverseTransposed = typename NurbsGridLinearAlgebraTraits::JacobianInverseTransposed;

      using ControlPointType = GlobalCoordinateType;
      using ControlPointNetType = MultiDimensionNet<dim,ControlPointType>;

      friend class BSplineLeafGridView<dim, dimworld>;

      /** \brief Constructor of BSplinePatch from knots, control points and order
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] order order of the B-Spline structure for each dimension
       */
      BSplinePatch(const std::array<std::vector<double>,dim>& knotSpans,
                   const ControlPointNetType controlPoints,
                   const std::array<int,dim> order)
      : patchData_(std::make_shared<BsplinePatchData<dim,dimworld>>(knotSpans, controlPoints, order))
      {
        validKnotSize_ = this -> validKnotSize();
        //Build a knot net to make iterator operations easier
        //Here each "point" of the net is a element(knot span)
        knotElementNet_ = std::make_shared<MultiDimensionNet<dim,double>>(validKnotSize_);
      }

      /** \brief creates a BSplineElementGeometry object
       *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
       *  and generates a Geometry object
       *
       * \param[in] ijk array of indices for each dimension
       */
      BSplineElementGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk ) const
      {
        const auto & knotSpans = patchData_->getKnots();

        std::array<unsigned int,dim> index;
        std::fill(index.begin(), index.end(), 0);

        /*finds the working geometry object ijk
         *(working geometry objects are defined between 2 knots, where knot[i]<knot[i+1])*/
        for (int count, j=0; j<dim; ++j)
        {
          count = 0;
          while(count <= ijk[j])
          {
            if (index[j] == knotSpans[j].size())
                break;
            if (knotSpans[j][index[j]+1] > knotSpans[j][index[j]])
                count++;

            ++index[j];
          }
        }

        /*the pointer on each dim-knotspan for geometry ijk is stored in an array named corners*/
        std::array<std::vector<double>::const_iterator,dim> corners;
        for(int i=0; i<dim; ++i)
          corners[i] = (knotSpans[i]).begin()+index[i]-1; //get starting span for each dim of element

        return BSplineElementGeometry<dim,dimworld>(patchData_,corners);
      }

      /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
      std::array<unsigned int,dim> validKnotSize() const
      {
        const auto & knotSpans = patchData_->getKnots();
        std::array<unsigned int,dim> validknotsize;
        validknotsize.fill(0);

        for (int j=0; j<dim; ++j)
          for (int i=0; i<knotSpans[j].size()-1; ++i)
            if (knotSpans[j][i+1] > knotSpans[j][i])
              ++validknotsize[j];

        return validknotsize;
      }

    private:

      std::shared_ptr <BsplinePatchData<dim,dimworld>> patchData_;
      std::array<unsigned int,dim> validKnotSize_;
      std::shared_ptr <MultiDimensionNet<dim,double>> knotElementNet_;
    };

}

#endif  // DUNE_IGA_BSPLINEPATCH_HH
