#ifndef DUNE_IGA_NURBSPATCH_HH
#define DUNE_IGA_NURBSPATCH_HH

#include <memory>
#include <dune/iga/bsplinepatch.hh>

namespace Dune
{
  namespace IGA
  {

    /** \brief class that holds all data regarding the NURBS structure */
    template<int dim, int dimworld>
    class NURBSPatchData
    {
    public:
      /** \brief constructor for a NURBSPatchData from knots, control points, weights and order
       *
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
       *  \param[in] order order of the NURBS structure for each dimension
       */
      NURBSPatchData(const std::array<std::vector<double>,dim>& knotSpans,
                     const MultiDimensionNet<dim,dimworld> controlPoints,
                     const MultiDimensionNet<dim,1> weights,
                     const std::array<int,dim> order)
      : knotSpans_(knotSpans)
      , controlPoints_(controlPoints)
      , weights_(weights)
      , order_(order)
      {
      }

      /** \brief returns the Knot Span*/
      const std::array<std::vector<double>, dim> & getKnots() const
      {
        return knotSpans_;
      }

      /** \brief returns the Control Point*/
      const MultiDimensionNet<dim,dimworld> & getControlPoints () const
      {
        return controlPoints_;
      }

      /** \brief returns the weights*/
      const MultiDimensionNet<dim,1> & getWeights() const
      {
        return weights_;
      }

      /** \brief returns the order*/
      const std::array<int,dim> & getOrder() const
      {
        return order_;
      }
    private:

      const std::array<std::vector<double>,dim>& knotSpans_;

      const MultiDimensionNet<dim,dimworld>  controlPoints_;

      const MultiDimensionNet<dim,1>  weights_;

      const std::array<int,dim> order_;

    };

    /** \brief a geometry implementation for NURBS*/
    template <int dim, int dimworld>
    class NURBSGeometry
    {
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
      typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< double > > MatrixHelper;

    public:

      /** \brief Constructor from NURBSPatchData and an iterator to a specific knot
       *
       *  \param[in] Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
       *  \param[in] corner Iterator (for each dimension) to the Knot span where the Geometry object is supposed to operate
       */
      NURBSGeometry(std::shared_ptr <NURBSPatchData<dim,dimworld>> patchData,
                    std::array<std::vector<double>::const_iterator,dim> corner)
      : patchData_(patchData)
      , corner_(corner)
      {
        const auto& weightedControlPoints = computeWeightedControlPoints();
        const auto& knotSpans = patchData_->getKnots();
        const auto& order = patchData_->getOrder();

        const auto& BSplinepatchData = std::make_shared<BsplinePatchData<dim,dimworld+1>>(knotSpans,weightedControlPoints,order);
        weightedBSpline = std::make_shared<BSplineGeometry<dim,dimworld+1>>(BSplinepatchData,corner_);

        /* only order+1 basis functions are needed for each dimension*/
        std::array<std::vector<double>,dim> _basis;
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

      /** \brief I think it is not an affine mapping, but not sure if should be a permanent false here*/
      bool affine () const { return false; }

      /** \brief Type of the element: a hypercube of the correct dimension */
      GeometryType type() const
      {
        return GeometryType(GeometryType::cube,dim);
      }

      /** \brief evaluates the NURBS mapping
       *
       *  \param[in] local local coordinates for each dimension
       */
      FieldVector<double,dimworld> global(const FieldVector<double,dim>& local)
      {
        const auto& weightedGlobal = weightedBSpline->global(local);
        return reduceDimension(weightedGlobal);

      }

//      /** \brief evaluates the NURBS mapping
//       *
//       *  \param[in] local local coordinates for each dimension
//       */
//      FieldVector<double,dimworld> global(const FieldVector<double,dim>& local)
//      {
//
//        const auto& knotSpans = patchData_->getKnots();
//        const auto& controlPoints = patchData_->getControlPoints();
//        const auto& weights = patchData_->getWeights();
//        const auto& order = patchData_->getOrder();
//        const auto& basis  = this->getBasis(local);
//
//        FieldVector<double,dimworld> denominator, numerator;
//        FieldVector<double,dimworld> glob;
//        std::fill(glob.begin(), glob.end(), 0.0);
//        std::fill(denominator.begin(), denominator.end(), 0.0);
//        std::fill(numerator.begin(), numerator.end(), 0.0);
//
//        std::array<unsigned int,dim> multiIndexBasisfucntion, multiIndexControlNet;
//        double temp;
//        std::array<unsigned int,dim> dimsize;
//        std::array<unsigned int,dim> cornerIdx;
//        for (unsigned int d=0; d<dim; ++d)
//        {
//          dimsize[d] = order[d]+1;
//          cornerIdx[d] = corner_[d]-knotSpans[d].begin();
//        }
//        /*Index net for valid basis functions*/
//        auto basisFunctionNet = MultiDimensionNet<dim,dimworld>(dimsize);
//        for (unsigned int i=0; i<basisFunctionNet.directSize(); ++i)
//        {
//          multiIndexBasisfucntion =  basisFunctionNet.directToMultiIndex(i);
//          temp = 1;
//          for (unsigned int d=0; d<dim; ++d)
//          {
//            temp *= basis[d][multiIndexBasisfucntion[d]];
//            multiIndexControlNet[d] = multiIndexBasisfucntion[d]+cornerIdx[d]-order[d];
//          }
//          auto cp = controlPoints.get(multiIndexControlNet);
//          auto w = weights.get(multiIndexControlNet);
//          for (unsigned int k=0; k<dimworld; ++k)
//          {
//            numerator[k] += cp[k]*w[0]*temp;
//            denominator[k] += w[0]*temp;
//          }
//        }
//        for (unsigned int k=0; k<dimworld; ++k)
//          glob[k] = numerator[k] / denominator[k];
//
//        return glob;
//      }

//      /** \brief evaluates the basis Functions at the given local coordinates
//       *
//       *  \param[in] local local coordinates for each dimension
//       */
//      const std::array<std::vector<double>,dim > & getBasis (const FieldVector<double,dim>& local)
//      {
//        for (int i=0; i<dim;++i)
//          if (local[i]<0 || local[i]>1)
//            DUNE_THROW(RangeError, "Local coordinates have to be in [0,1]^dim!");
//
//        const auto & order = patchData_->getOrder();
//        const auto & knotSpans = patchData_->getKnots();
//
//        /*Compute the non-vanishing basis function. See "The NURBS book"*/
//        for (int d=0; d<dim; ++d)
//        {
//          double loc;
//          loc = local[d]*(*(corner_[d]+1)-*corner_[d]) + *corner_[d];
//          basis_[d].resize(order[d]+1);
//          std::vector<double> left, right;
//          left.resize(order[d]+1);
//          right.resize(order[d]+1);
//          basis_[d][0] = 1;
//          for (int o=1; o<=order[d]; ++o)
//          {
//            left[o] = loc-*(corner_[d]+1-o);
//            right[o] = *(corner_[d]+o)-loc;
//            double saved = 0;
//            for (int r=0; r<o; ++r)
//            {
//              double temp = basis_[d][r]/(right[r+1]+left[o-r]);
//              basis_[d][r] = saved+right[r+1]*temp;
//              saved = left[o-r]*temp;
//            }
//            basis_[d][o] = saved;
//          }
//        }
//
//        return basis_;
//      }

      /** \brief compute the Jacobian transposed matrix
       *
       *  \param[in] local local coordinates for each dimension
       */
      const JacobianTransposed jacobianTransposed(const LocalCoordinate &local)
      {
        JacobianTransposed result;
        const auto& BSplineJacobian = weightedBSpline->jacobianTransposed(local);
        const auto& weightedGlobal = weightedBSpline->global(local);
        const auto& dividedGlobal = reduceDimension(weightedGlobal);
        const ctype weight = weightedGlobal[coorddimension];

        //Each row is a partial derivative with respect to one variable
        for (int i=0; i<mydimension; ++i)
        {
           //Iterate the columns except the last one
           const ctype derWeight = BSplineJacobian[i][coorddimension];
           for (int j=0; j<coorddimension; ++j)
           {
             result[i][j] = (BSplineJacobian[i][j]-derWeight*dividedGlobal[j])/weight;
           }
        }
        return result;
      }

      /** \brief constructs the homogeneous coordinates
       *
       *  \param[in] inPoint input point
       *  \param[in] weight weight
       */
      FieldVector<double,coorddimension+1> liftDimension(const FieldVector<double,coorddimension>& inPoint, ctype weight)
      {
        FieldVector<ctype,coorddimension+1> homoCoordinate;
        for(unsigned int i=0; i<coorddimension; i++)
        {
          homoCoordinate[i] = inPoint[i]*weight;
        }
        homoCoordinate[coorddimension] = weight;
        return homoCoordinate;
      }

      /** \brief get the projected non-homogeneous coordinates in a 1 degree lower dimensional space
       *
       *  \param[in] inPoint input point, will be divided by the last element
       */
      FieldVector<double,coorddimension> reduceDimension(const FieldVector<double,coorddimension+1>& inPoint)
      {
        FieldVector<ctype,coorddimension> oriCoordinate;
        for(unsigned int i=0; i<coorddimension; i++)
        {
          oriCoordinate[i] = inPoint[i]/inPoint[coorddimension];
        }
        return oriCoordinate;
      }

      /** \brief gets weighted control points
       *
       */
      MultiDimensionNet<mydimension,coorddimension+1> computeWeightedControlPoints()
      {
        const auto& controlPoints = patchData_->getControlPoints();
        const auto& weights = patchData_->getWeights();
        const auto& netSize = controlPoints.size();
        MultiDimensionNet<mydimension,coorddimension+1> weightedControlPoints(netSize);
        const auto& pointSize = controlPoints.directSize();
        for(unsigned int i=0; i<pointSize; ++i)
        {
          //auto const& liftedControlPoint = ;
          weightedControlPoints.directSet(i,liftDimension(controlPoints.directGet(i), weights.directGet(i)));
        }
        return weightedControlPoints;

      }


    private:

      std::shared_ptr <NURBSPatchData<dim,dimworld>> patchData_;

      std::array<std::vector<double>::const_iterator,dim> corner_;

      std::shared_ptr <BSplineGeometry<dim, dimworld+1>> weightedBSpline;

    };

    template<int dim, int dimworld>
    class NURBSLeafGridView;

    /** \brief Class where the NURBS geometry can work on */
    template <int dim, int dimworld>
    class NURBSPatch
    {
    public:

      friend class NURBSLeafGridView<dim, dimworld>;
      template<int codim, class GridViewImp>
      friend class NURBSGridEntity;

      /** \brief Constructor of NURBS from knots, control points, weights and order
       *  \param[in] knotSpans vector of knotSpans for each dimension
       *  \param[in] controlPoints a n-dimensional net of control points
       *  \param[in] weights vector a n-dimensional net of weights for each corresponding control points
       *  \param[in] order order of the NURBS structure for each dimension
       */
      NURBSPatch(const std::array<std::vector<double>,dim>& knotSpans,
                   const MultiDimensionNet<dim,dimworld> controlPoints,
                   const MultiDimensionNet<dim,1> weights,
                   const std::array<int,dim> order)
      : patchData_(std::make_shared<NURBSPatchData<dim,dimworld>>(knotSpans, controlPoints, weights, order))
      {
        validKnotSize_ = this -> validKnotSize();
        //Build a knot net to make iterator operations easier
        //Here each "point" of the net is a element(knot span)
        knotElementNet_ = std::make_shared<MultiDimensionNet<dim,1>>(validKnotSize_);
      }

      /** \brief creates a NURBSGeometry object
       *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimension
       *  and generates a Geometry object
       *
       *  \param[in] ijk array of indices for each dimension
       */
      NURBSGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk ) const
      {
        const auto & knotSpans = patchData_->getKnots();

        unsigned int count = 0;
        std::array<unsigned int,dim> index;
        std::fill(index.begin(), index.end(), 0);

        /*finds the working geometry object ijk
         *(working geometry objects are defined between 2 knots, where knot[i]<knot[i+1])*/
        for (int j=0; j<dim; ++j)
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

        /*the iterator on each dim-knotspan for geometry ijk is stored in an array named corners*/
        std::array<std::vector<double>::const_iterator,dim> corners;
        for(int i=0; i<dim; ++i)
        {
          corners[i] = (knotSpans[i]).begin()+index[i]-1;
        }

        return NURBSGeometry<dim,dimworld>(patchData_,corners);
      }

      /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
      std::array<unsigned int,dim> validKnotSize() const
      {
        const auto & knotSpans = patchData_->getKnots();
        std::array<unsigned int,dim> validknotsize;
        std::fill(validknotsize.begin(), validknotsize.end(), 0);

        for (int j=0; j<dim; ++j)
        {
          for (int i=0; i<knotSpans[j].size()-1; ++i)
          {
            if (knotSpans[j][i+1] > knotSpans[j][i])
              ++validknotsize[j];
          }
        }

        return validknotsize;
      }

    private:

      std::shared_ptr <NURBSPatchData<dim,dimworld>> patchData_;
      std::array<unsigned int,dim> validKnotSize_;
      std::shared_ptr <MultiDimensionNet<dim,1>> knotElementNet_;

    };
  }
}

#endif  // DUNE_IGA_NURBSPATCH_HH
