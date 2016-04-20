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
        /* only order+1 basis functions are needed for each dimension*/
        std::array<std::vector<double>,dim> _basis;
      }

      /** \brief evaluates the NURBS mapping
       *
       *  \param[in] local local coordinates for each dimension
       */
      FieldVector<double,dimworld> global(const FieldVector<double,dim>& local)
      {

        const auto& knotSpans = patchData_->getKnots();
        const auto& controlPoints = patchData_->getControlPoints();
        const auto& weights = patchData_->getWeights();
        const auto& order = patchData_->getOrder();
        const auto& basis  = this->getBasis(local);

        FieldVector<double,dimworld> denominator, numerator;
        FieldVector<double,dimworld> glob;
        std::fill(glob.begin(), glob.end(), 0.0);
        std::fill(denominator.begin(), denominator.end(), 0.0);
        std::fill(numerator.begin(), numerator.end(), 0.0);

        std::array<unsigned int,dim> multiIndexBasisfucntion, multiIndexControlNet;
        double temp;
        std::array<unsigned int,dim> dimsize;
        std::array<unsigned int,dim> cornerIdx;
        for (unsigned int d=0; d<dim; ++d)
        {
          dimsize[d] = order[d]+1;
          cornerIdx[d] = corner_[d]-knotSpans[d].begin();
        }
        /*Index net for valid basis functions*/
        auto basisFunctionNet = MultiDimensionNet<dim,dimworld>(dimsize);
        for (unsigned int i=0; i<basisFunctionNet.directSize(); ++i)
        {
          multiIndexBasisfucntion =  basisFunctionNet.directToMultiIndex(i);
          temp = 1;
          for (unsigned int d=0; d<dim; ++d)
          {
            temp *= basis[d][multiIndexBasisfucntion[d]];
            multiIndexControlNet[d] = multiIndexBasisfucntion[d]+cornerIdx[d]-order[d];
          }
          auto cp = controlPoints.get(multiIndexControlNet);
          auto w = weights.get(multiIndexControlNet);
          for (unsigned int k=0; k<dimworld; ++k)
          {
            numerator[k] += cp[k]*w[0]*temp;
            denominator[k] += w[0]*temp;
          }
        }
        for (unsigned int k=0; k<dimworld; ++k)
          glob[k] = numerator[k] / denominator[k];

        return glob;
      }

      /** \brief evaluates the basis Functions at the given local coordinates
       *
       *  \param[in] local local coordinates for each dimension
       */
      const std::array<std::vector<double>,dim > & getBasis (const FieldVector<double,dim>& local)
      {
        for (int i=0; i<dim;++i)
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
          basis_[d][0] = 1;
          for (int o=1; o<=order[d]; ++o)
          {
            left[o] = loc-*(corner_[d]+1-o);
            right[o] = *(corner_[d]+o)-loc;
            double saved = 0;
            for (int r=0; r<o; ++r)
            {
              double temp = basis_[d][r]/(right[r+1]+left[o-r]);
              basis_[d][r] = saved+right[r+1]*temp;
              saved = left[o-r]*temp;
            }
            basis_[d][o] = saved;
          }
        }

        return basis_;
      }

    private:

      std::shared_ptr <NURBSPatchData<dim,dimworld>> patchData_;

      std::array<std::vector<double>::const_iterator,dim> corner_;

      std::array<std::vector<double>,dim> basis_;

    };

    /** \brief Class where the NURBS geometry can work on */
    template <int dim, int dimworld>
    class NURBSPatch
    {
    public:
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

    };
  }
}

#endif  // DUNE_IGA_NURBSPATCH_HH
