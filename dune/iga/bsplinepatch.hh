// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <memory>


namespace Dune
{
template <int netdim, int dimworld>
class MultiDimensionNet
{
public:
	MultiDimensionNet (std::array<unsigned int, netdim> dimSize)
	: dimSize_(dimSize)
	{
		int size = 1;
		for (int i=0; i<netdim; ++i)
			size*=dimSize_[i];

		values_.resize(size);
	}
	MultiDimensionNet (std::array<unsigned int, netdim> dimSize, std::vector<FieldVector<double, dimworld>> values)
	: values_(values)
	, dimSize_(dimSize)
	{}

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

	void set (std::array<unsigned int, netdim> multiIndex, FieldVector<double,dimworld> value)
	{
		int index = this->index(multiIndex);
		values_[index] = value;
	}

	FieldVector<double, dimworld> get (std::array<unsigned int, netdim> multiIndex) const
	{
		int index = this->index(multiIndex);
		return values_[index];
	}

	FieldVector<double, dimworld> directGet (int index) const
	{
		return values_[index];
	}

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
// used for testing
// 	void disp() const
// 	{
// 		for (int i=0; i<values_.size(); ++i)
// 			std::cout<<"i: "<<i<<" value= "<<values_[i][0]<<" "<<values_[i][1]<<" "<<values_[i][2]<<std::endl;
// 	}

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

template<int dim, int dimworld>
class BsplinePatchData
{
public:
    BsplinePatchData(const std::array<std::vector<double>,dim>& knotSpans,
               const MultiDimensionNet<dim,dimworld> controlPoints,
               const std::array<int,dim> order)
    : knotSpans_(knotSpans)
    , controlPoints_(controlPoints)
    , order_(order)
    {
    }

    const std::array<std::vector<double>, dim> & getknots() const
    {
        return knotSpans_;
    }
    const MultiDimensionNet<dim,dimworld> & getcontrolPoints () const
    {
        return controlPoints_;
    }
    const std::array<int,dim> & getorder() const
    {
        return order_;
    }
private:

  const std::array<std::vector<double>,dim>& knotSpans_;

  const MultiDimensionNet<dim,dimworld>  controlPoints_;

  const std::array<int,dim> order_;

};

template <int dim, int dimworld>
class BSplineGeometry
{
public:

  BSplineGeometry(std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata, std::array<const double*,dim> corner)
  : Patchdata_(Patchdata)
  , corner_(corner)
  {
      std::array<std::vector<std::vector<double>>,dim> basis;
      auto & knotSpans = Patchdata_->getknots();
      auto & order = Patchdata_->getorder();

        for (int d=0; d<dim; ++d)
        {
            basis[d].resize(order[d]+1);
            for (int j=0; j<knotSpans[d].size()-1; ++j)
            {
                if(&(knotSpans[d][j]) == (corner_[d]))
                    basis[d][0].push_back(1);
                else
                    basis[d][0].push_back(0);
            }
        }
      basis_ = basis;

}

  FieldVector<double,dimworld> global(const FieldVector<double,dim>& local)
  {

    const auto& knotSpans = Patchdata_->getknots();
    const auto& controlPoints = Patchdata_->getcontrolPoints();
    const auto& order = Patchdata_->getorder();
    const auto& basis  = this->getbasis(local);

    FieldVector<double,dimworld> glob;
    std::fill(glob.begin(), glob.end(), 0.0);

	std::array<unsigned int,dim> multiIndex;
	double temp;

	for (unsigned int i=0; i<controlPoints.directSize(); ++i)
	{
		auto cp = controlPoints.directGet(i);
		multiIndex =  controlPoints.directToMultiIndex(i);
		temp = 1;
		for (unsigned int d=0; d<dim; ++d)
			temp *= basis[d][order[d]][multiIndex[d]];

		for (unsigned int k=0; k<dimworld; ++k)
			glob[k] += cp[k]*temp;

	}

    return glob;


  }

  const std::array<std::vector<std::vector<double> >,dim > & getbasis (const FieldVector<double,dim>& local)
  {
    for (int i=0; i<dim;++i)
      if (local[i]<0 || local[i]>1)
        DUNE_THROW(RangeError, "Local coordinates have to be in [0,1]^dim!");

      /*note: define lower and upperbounds so a loop over all knots is not needed.... or fnd a better way to generate the right basis functions*/
        const auto & order = Patchdata_->getorder();
        const auto & knotSpans = Patchdata_->getknots();
        double loc, A, B;

        /*generate the basis functions using the Cox-de Boor recursion formula*/
        for (int d=0; d<dim; ++d)
        {
            loc = local[d]*(*(corner_[d]+1)-*corner_[d]) + *corner_[d];
            for (int o=1; o<=order[d]; ++o)
            {
                basis_[d][o].resize(knotSpans[d].size()-o);
                for (int k=0; k<(knotSpans[d].size()-(o+1)); ++k)
                {
                    if ((knotSpans[d][k+o]-knotSpans[d][k]) == 0)
                        A = 0;
                    else
                        A = ((loc)-knotSpans[d][k])/ (knotSpans[d][k+o]-knotSpans[d][k]);

                    if ((knotSpans[d][k+o+1]-knotSpans[d][k+1]) == 0)
                        B = 0;
                    else
                        B = (knotSpans[d][k+o+1]-(loc))/ (knotSpans[d][k+o+1]-knotSpans[d][k+1]);

                    basis_[d][o][k] =  A*basis_[d][o-1][k] + B*basis_[d][o-1][k+1];

                }
            }
        }

    return basis_;
  }

private:

    std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata_;

    std::array<const double*,dim> corner_;

    std::array<std::vector<std::vector<double>>,dim> basis_;

};

template <int dim, int dimworld>
class BSplinePatch
{
public:

  BSplinePatch(const std::array<std::vector<double>,dim>& knotSpans,
               const MultiDimensionNet<dim,dimworld> controlPoints,
               const std::array<int,dim> order)
  : Patchdata_(new BsplinePatchData<dim,dimworld>(knotSpans, controlPoints, order))
  {
  }

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
        std::array<const double* ,dim> corners;
        for(int i=0; i<dim; ++i)
            corners[i] = &knotSpans[i][index[i]-1];


    return BSplineGeometry<dim,dimworld>(Patchdata_,corners);
  }

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