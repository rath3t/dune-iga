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
    /** \brief default Constructor
     *
     * \param[in] knotSpans vector of knotSpans for each dimension
     * \param[in] controlPoints a n-dimensional net of Control controlPoints
     * \param[in] order order of the NURBS structure for each dimension
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
    const std::array<std::vector<double>, dim> & getknots() const
    {
        return knotSpans_;
    }

    /** \brief returns the Control Point*/
    const MultiDimensionNet<dim,dimworld> & getcontrolPoints () const
    {
        return controlPoints_;
    }

    /** \brief returns the weights*/
    const MultiDimensionNet<dim,1> & getweights() const
    {
        return weights_;
    }

    /** \brief returns the order*/
    const std::array<int,dim> & getorder() const
    {
        return order_;
    }
private:

  const std::array<std::vector<double>,dim>& knotSpans_;

  const MultiDimensionNet<dim,dimworld>  controlPoints_;

  const MultiDimensionNet<dim,1>  weights_;

  const std::array<int,dim> order_;

};

/** \brief */
template <int dim, int dimworld>
class NURBSGeometry
{
public:
  /** \brief default Constructor
   *
   * \param[in] Patchdata shared pointer to an object where the all the data of the NURBSPatch is stored
   * \param[in] corner Iterator (for each dimension) to the Knot span where the Geometry object is supposed to operate
   */
  NURBSGeometry(std::shared_ptr <NURBSPatchData<dim,dimworld>> Patchdata, std::array<std::vector<double>::const_iterator,dim> corner)
  : Patchdata_(Patchdata)
  , corner_(corner)
  {
      //* only order+1 basis functions are needed for each dimension*/
      std::array<std::vector<double>,dim> _basis;
  }

  /** \brief evaluates the NURBS mapping
   *
   * \param[in] local local coordinates for each dimension
   */
    FieldVector<double,dimworld> global(const FieldVector<double,dim>& local)
  {

    const auto& knotSpans = Patchdata_->getknots();
    const auto& controlPoints = Patchdata_->getcontrolPoints();
    const auto& weights = Patchdata_->getweights();
    const auto& order = Patchdata_->getorder();
    const auto& basis  = this->getbasis(local);

    FieldVector<double,dimworld> denominator, numerator;
    FieldVector<double,dimworld> glob;
    std::fill(glob.begin(), glob.end(), 0.0);
    std::fill(denominator.begin(), denominator.end(), 0.0);
    std::fill(numerator.begin(), numerator.end(), 0.0);

    std::array<unsigned int,dim> multiIndex_Basisfucntion, multiIndex_ControlNet;
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
        multiIndex_Basisfucntion =  basisFunctionNet.directToMultiIndex(i);
        temp = 1;
        for (unsigned int d=0; d<dim; ++d)
        {
            temp *= basis[d][multiIndex_Basisfucntion[d]];
            multiIndex_ControlNet[d] = multiIndex_Basisfucntion[d]+cornerIdx[d]-order[d];
        }
        auto cp = controlPoints.get(multiIndex_ControlNet);
        auto w = weights.get(multiIndex_ControlNet);
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
   * \param[in] local loacal coordinates for each dimension
   */
  const std::array<std::vector<double>,dim > & getbasis (const FieldVector<double,dim>& local)
  {
    for (int i=0; i<dim;++i)
      if (local[i]<0 || local[i]>1)
        DUNE_THROW(RangeError, "Local coordinates have to be in [0,1]^dim!");

      /*note: define lower and upperbounds so a loop over all knots is not needed.... or fnd a better way to generate the right basis functions*/
        const auto & order = Patchdata_->getorder();
        const auto & knotSpans = Patchdata_->getknots();
        double loc;
        double saved;

        /*Compute the nonvanishing basis function. See "The NURBS book"*/
        for (int d=0; d<dim; ++d)
        {
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
                saved = 0;
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

    std::shared_ptr <NURBSPatchData<dim,dimworld>> Patchdata_;

    std::array<std::vector<double>::const_iterator,dim> corner_;

    std::array<std::vector<double>,dim> basis_;

};

/** \brief class */
template <int dim, int dimworld>
class NURBSPatch
{
public:
    /** \brief default Constructor
     * \param[in] knotSpans vector of knotSpans for each dimension
     * \param[in] controlPoints a n-dimensional net of Control controlPoints
     * \param[in] order orde of the NURBS structure for each dimension
     */
  NURBSPatch(const std::array<std::vector<double>,dim>& knotSpans,
               const MultiDimensionNet<dim,dimworld> controlPoints,
               const MultiDimensionNet<dim,1> weights,
               const std::array<int,dim> order)
  : Patchdata_(new NURBSPatchData<dim,dimworld>(knotSpans, controlPoints, weights, order))
  {
  }

  /** \brief creates a NURBSGeometry object
   *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimesion
   *  and generates a Geometry object
   *
   * \param[in] ijk array of indecies for each dimension
   */
  NURBSGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk ) const
  {

        const auto & knotSpans = Patchdata_->getknots();

        unsigned int count = 0;
        std::array<unsigned int,dim> index;
        std::fill(index.begin(), index.end(), 0);

        /*finds the working geometry object ijk
         *(working geometry objects are difined between 2 knots, were knot[i]<kont[i+1])*/
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


    return NURBSGeometry<dim,dimworld>(Patchdata_,corners);
  }

  /** \brief returns the size of knot spans where knot[i] < knot[i+1] of each dimension */
  std::array<unsigned int,dim> validknotsize() const
  {
      const auto & knotSpans = Patchdata_->getknots();
      std::array<unsigned int,dim> validknotsize;
      std::fill(validknotsize.begin(), validknotsize.end(), 0);

      for (int j=0; j<dim; ++j){
        for (int i=0; i<knotSpans[j].size()-1; ++i)
        {
                if (knotSpans[j][i+1] > knotSpans[j][i])
                    ++validknotsize[j];

            }
        }

        return validknotsize;
  }

private:

  std::shared_ptr <NURBSPatchData<dim,dimworld>> Patchdata_;

};


}

}

#endif  // DUNE_IGA_NURBSPATCH_HH
