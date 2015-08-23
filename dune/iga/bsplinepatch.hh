// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <memory>

namespace Dune
{
namespace IGA
{

template<int dim, int dimworld>
class BsplinePatchData
{
public:
    BsplinePatchData(const std::array<std::vector<double>,dim>& knotSpans,
               const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints,
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
    const std::vector<std::vector<FieldVector<double,dimworld> > > & getcontrolPoints () const
    {
        return controlPoints_;
    }
    const std::array<int,dim> & getorder() const
    {
        return order_;
    }
private:

  const std::array<std::vector<double>,dim>& knotSpans_;

  const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints_;

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

    switch (dim)
    {
      case 1:
      {
        for (int j=0; j<controlPoints[0].size(); ++j)
            for (int k=0; k<dimworld; ++k)
                glob[k] += controlPoints[0][j][k]*basis[0][order[0]][j];

        break;
      }
      case 2:
      {
        for (int i=0; i<controlPoints[1].size(); ++i)
            for (int j=0; j<controlPoints[0].size(); ++j)
                for (int k=0; k<dimworld; ++k)
                    glob[k] += controlPoints[i][j][k]*basis[0][order[0]][i]*basis[1][order[1]][j];

        break;
      }
      default:
        DUNE_THROW(Dune::NotImplemented, "General case not implemented yet!");
        /* find a way to handle n-D nets (maybe a new stuct can help a n-D net mapped into a 2-D vector )*/
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

                for (int k=0; k<(knotSpans[d].size()-o); ++k)
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
               const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints,
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