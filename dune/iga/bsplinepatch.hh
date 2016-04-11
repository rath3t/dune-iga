// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <memory>
#include <dune/geometry/type.hh>

namespace Dune
{
	/** \brief class holds a n-dim net */
template <int netdim, int dimworld>
class MultiDimensionNet
{
public:
	/** \brief default constructor
	 *
	 * \param[in] dimSize array of the size of each dimension
	 */
	MultiDimensionNet (std::array<unsigned int, netdim> dimSize)
	: dimSize_(dimSize)
	{
		int size = 1;
		for (int i=0; i<netdim; ++i)
			size*=dimSize_[i];

		values_.resize(size);
	}

	/** \brief constructor intended for the 1-D if the values are already in a vector
	 * \note can also be used if the values are already mapped
	 *
	 *  \param[in] dimSize array of the size of each dimension
	 *  \param[in] values vector with values
	 */
	MultiDimensionNet (std::array<unsigned int, netdim> dimSize, std::vector<FieldVector<double, dimworld>> values)
	: values_(values)
	, dimSize_(dimSize)
	{}

	/** \brief constructor intended for the 2-D if the values are already in a matrix
	 *
	 *  \param[in] dimSize array of the size of each dimension
	 *  \param[in] values matrix with values
	 */
	MultiDimensionNet (std::array<unsigned int, netdim> dimSize, std::vector<std::vector<FieldVector<double, dimworld> > > values)
	: dimSize_(dimSize)
	{
		values_.resize(values.size()*values[0].size());

		for (unsigned int i=0; i<values.size(); ++i){
			for(unsigned int j=0; j<values[0].size(); ++j){
				std::array<unsigned int,netdim> multiIndex = {i,j};
				this->set(multiIndex, values[i][j]);
			}
		}
	}

    /** \brief constructor for a grid of the same value
     *
     *  \param[in] dimSize array of the size of each dimension
     *  \param[in] value a common value to fill the grid
     */
    MultiDimensionNet (std::array<unsigned int, netdim> dimSize, FieldVector<double, dimworld> value)
    : dimSize_(dimSize)
    {
        int size = 1;
        for (int i=0; i<netdim; ++i)
             size*=dimSize_[i];
        values_.resize(size);
        std::fill(values_.begin(), values_.end(), value);
    }

	/** \brief sets a value at the multiindex */
	void set (std::array<unsigned int, netdim> multiIndex, FieldVector<double,dimworld> value)
	{
		int index = this->index(multiIndex);
		values_[index] = value;
	}

	/** \brief returns the value at the multiindex */
	FieldVector<double, dimworld> get (std::array<unsigned int, netdim> multiIndex) const
	{
		int index = this->index(multiIndex);
		return values_[index];
	}

	/** \brief returns a value at an index (unmapped)
	  * \note only to be used when the mapping is known
	 */
	FieldVector<double, dimworld> directGet (int index) const
	{
		return values_[index];
	}

	/** \brief returns a multiindex for a scalar index */
	std::array<unsigned int,netdim> directToMultiIndex (unsigned int index) const
	{
		std::array<unsigned int, netdim> multiIndex;

		unsigned int help = index ;
		int temp;
		for (int i=0; i<netdim; ++i)
		{
				temp = help%(dimSize_[netdim-(i+1)]);
				multiIndex[netdim-(i+1)] = temp;
				help -= temp;
				help = help/dimSize_[netdim-(i+1)];
		}


		return multiIndex;
	}

	/** \brief returns an array with the size of each dimension */
	std::array<unsigned int,netdim> size() const
	{
		return dimSize_;
	}

	unsigned int directSize () const
	{
		return values_.size();
	}
private:
	int index (std::array<unsigned int, netdim> multiIndex)const
	{
		int index,help ;
		index = 0;
		for (int i=0; i<netdim; ++i)
		{
			help = 1;
			for (int j=(i+1); j<netdim ;++j)
				help *= dimSize_[j];

			index += help*multiIndex[i];
		}
		return index;
	}

private:

	std::array<unsigned int,netdim> dimSize_;
	std::vector<FieldVector<double,dimworld> > values_;
};

namespace IGA
{

/** \brief class that holds all data regarding the B-Spline stucture */
template<int dim, int dimworld>
class BsplinePatchData
{
public:
	/** \brief default Constructor
	 *
	 * \param[in] knotSpans vector of knotSpans for each dimension
	 * \param[in] controlPoints a n-dimensional net of Control controlPoints
	 * \param[in] order orde of the B-Spline structure for each dimension
	 */
    BsplinePatchData(const std::array<std::vector<double>,dim>& knotSpans,
               const MultiDimensionNet<dim,dimworld> controlPoints,
               const std::array<int,dim> order)
    : knotSpans_(knotSpans)
    , controlPoints_(controlPoints)
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

    /** \brief returns the order*/
    const std::array<int,dim> & getorder() const
    {
        return order_;
    }
private:

  const std::array<std::vector<double>,dim>& knotSpans_;

  const MultiDimensionNet<dim,dimworld>  controlPoints_;

  const std::array<int,dim> order_;

};

/** \brief */
template <int dim, int dimworld>
class BSplineGeometry
{
public:

  //! coordinate type
  typedef double ctype;

  /** \brief Dimension of the cube element */
  enum {mydimension = dim};

  /** \brief Dimension of the world space that the cube element is embedded in*/
  enum {coorddimension = dimworld};

  /** \brief Type used for a vector of element coordinates */
  typedef FieldVector<ctype,dim> LocalCoordinate;

  /** \brief Type used for a vector of world coordinates */
  typedef FieldVector<ctype,dimworld> GlobalCoordinate;

  /** \brief default Constructor
   *
   * \param[in] Patchdata shared pointer to an object where the all the data of the BSplinePatch is stored
   * \param[in] corner Pointer (for each dimension) to the Knot span where the Geometry object is supposed to operate
   */
  BSplineGeometry(std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata, std::array<std::vector<double>::const_iterator,dim> corner)
  : Patchdata_(Patchdata)
  , corner_(corner)
  {
	  // note it`s maybe not such a good a idea to store all of that!
      std::array<std::vector<double>,dim> _basis;
  }

  /** \brief Map the center of the element to the geometry */
  GlobalCoordinate center() const
  {
      LocalCoordinate l_center;
      std::fill(l_center.begin(), l_center.end(), 0.5);
      return global(l_center);
  }

  /** \brief Return the number of corners of the element */
  int corners() const
  {
      return 1<<dim;
  }

  /** \brief Return world coordinates of the k-th corner of the element */
  GlobalCoordinate corner(int k) const
  {
      LocalCoordinate l_corner;
      for (size_t i=0; i<dim; i++) {
          l_corner[i] = (k & (1<<i)) ? 0 : 1;
      }
      return global(l_corner);
  }

  /** \brief I think it is not an affine mapping */
  bool affine () const { return false; }

  /** \brief Type of the element: a hypercube of the correct dimension */
  GeometryType type() const
  {
      return GeometryType(GeometryType::cube,dim);
  }

  /** \brief evaluates the B-Spline mapping
   *
   * \param[in] local local coordinates for each dimension
   */
  GlobalCoordinate global(const LocalCoordinate& local)
  {

    const auto& knotSpans = Patchdata_->getknots();
    const auto& controlPoints = Patchdata_->getcontrolPoints();
    const auto& order = Patchdata_->getorder();
    const auto& basis  = this->getbasis(local);
    GlobalCoordinate glob;
    std::fill(glob.begin(), glob.end(), 0.0);
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
		for (unsigned int k=0; k<dimworld; ++k)
			glob[k] += cp[k]*temp;
	}

    return glob;


  }

  /** \brief evaluates the basis Functions at the given local coordinates
   *
   * \param[in] local loacal coordinates for each dimension
   */
  const std::array<std::vector<double>,dim > & getbasis (const LocalCoordinate& local)
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

    std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata_;

    std::array<std::vector<double>::const_iterator,dim> corner_;

    std::array<std::vector<double>,dim> basis_;

};

/** \brief class */
template <int dim, int dimworld>
class BSplinePatch
{
public:
	/** \brief default Constructor
	 *
	 * \param[in] knotSpans vector of knotSpans for each dimension
	 * \param[in] controlPoints a n-dimensional net of Control controlPoints
	 * \param[in] order orde of the B-Spline structure for each dimension
	 */
  BSplinePatch(const std::array<std::vector<double>,dim>& knotSpans,
               const MultiDimensionNet<dim,dimworld> controlPoints,
               const std::array<int,dim> order)
  : Patchdata_(new BsplinePatchData<dim,dimworld>(knotSpans, controlPoints, order))
  {
  }

  /** \brief creates a BSplineGeometry object
   *  this function finds the i-th knot span where knot[i] < knot[i+1] for each dimesion
   *  and generates a Geometry object
   *
   * \param[in] ijk array of indecies for each dimension
   */
  BSplineGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk ) const
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

        /*the pointer on each dim-knotspan for geometry ijk is stored in an array named corners*/
        std::array<std::vector<double>::const_iterator,dim> corners;
        for(int i=0; i<dim; ++i)
        {
          corners[i] = (knotSpans[i]).begin()+index[i]-1;
        }


    return BSplineGeometry<dim,dimworld>(Patchdata_,corners);
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

  std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata_;

};


}

}

#endif  // DUNE_IGA_BSPLINEPATCH_HH
