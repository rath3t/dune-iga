#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <dune/geometry/multilineargeometry.hh>
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
               const std::vector<FieldVector<double,dimworld> > controlPoints,
               const int order)
    : knotSpans_(knotSpans)
    , controlPoints_(controlPoints)
    , order_(order)
    {
    }

    const std::array<std::vector<double>, dim> & getknots()
    {
        return knotSpans_;
    }
    const std::vector<FieldVector<double,dimworld>> & getcontrolPoints ()
    {
        return controlPoints_;
    }
    const int getorder()
    {
        return order_;
    }
private:

  const std::array<std::vector<double>,dim>& knotSpans_;

  std::vector<FieldVector<double,dimworld> > controlPoints_;

  const int order_;

};

template <int dim, int dimworld>
class BSplineGeometry
{
public:

  BSplineGeometry(std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata, const double* corner)
  : Patchdata_(Patchdata)
  , corner_(corner)
  {}

  FieldVector<double,dimworld> global(const FieldVector<double,dim>& local) const
  {
    if (local[0]<0 || local[0]>1)
    {
        //DUNE_THROW(OutofRange, "local coordinates havae to be in [0,1]!");
        std::cout<<"fail"<<std::endl;
    }
    else
    {
        double loc = local[0]*(*(corner_+1)-*corner_) + *corner_;

        const auto& knotSpans = Patchdata_->getknots();
        const auto& controlPoints = Patchdata_->getcontrolPoints();
        auto order = Patchdata_->getorder();

        std::vector<std::vector<double>> basis;
        basis.resize(order+1);

        for (int j=0; j<knotSpans[0].size()-1; ++j)
        {
            if(&(knotSpans[0][j]) == (corner_))
                basis[0].push_back(1);
            else
            {
                basis[0].push_back(0);
            }
        }

         double A,B;
        for (int o=1; o<=order; ++o)
        {
            for (int k=0; k<(knotSpans[0].size()-o); ++k)
            {
                if ((knotSpans[0][k+o]-knotSpans[0][k]) == 0)
                    A = 0;
                else
                    A = ((loc)-knotSpans[0][k])/ (knotSpans[0][k+o]-knotSpans[0][k]);

                if ((knotSpans[0][k+o+1]-knotSpans[0][k+1]) == 0)
                    B = 0;
                else
                    B = (knotSpans[0][k+o+1]-(loc))/ (knotSpans[0][k+o+1]-knotSpans[0][k+1]);

                 basis[o].push_back(A*basis[o-1][k] + B*basis[o-1][k+1]);
            }
        }

        FieldVector<double,dimworld> glob = {0,0,0};

        for (int i=0; i<controlPoints.size(); ++i)
        {
            for (int j=0; j<dimworld; ++j)
                glob[j] += controlPoints[i][j]*basis[order][i];

        }

    return glob;

    }
  }

private:

    std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata_;

    const double* corner_;

};

template <int dim, int dimworld>
class BSplinePatch
{
public:

  BSplinePatch(/*const FieldVector<double,dim>& lower,
               const FieldVector<double,dim>& upper,*/
               const std::array<std::vector<double>,dim>& knotSpans,
               const std::vector<FieldVector<double,dimworld> > controlPoints,
               const int order)
  : Patchdata_(new BsplinePatchData<dim,dimworld>(knotSpans, controlPoints, order))
  {
  }

  BSplineGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk ) const
  {
//     std::vector<FieldVector<double,dimworld> > corners(1<<dim);
//         corners[0] = controlPoints_[ijk[0]];
//         corners[1] = controlPoints_[ijk[0]+1];

//     switch (dim)
//     {
//       case 1:
//       {
//         corners[0] = controlPoints_[ijk[0]];
//         corners[1] = controlPoints_[ijk[0]+1];
//
//         break;
//       }
//       case 2:
//       {
//         corners[0] = controlPoints_[ijk[1]*(knotSpans_[0]+1) + ijk[0]];
//         corners[1] = controlPoints_[ijk[1]*(knotSpans_[0]+1) + ijk[0]+1];
//         corners[2] = controlPoints_[(ijk[1]+1)*(knotSpans_[0]+1) + ijk[0]];
//         corners[3] = controlPoints_[(ijk[1]+1)*(knotSpans_[0]+1) + ijk[0]+1];
//
//         break;
//       }
//       default:
//         DUNE_THROW(Dune::NotImplemented, "General case not implemented yet!");
//     }


        unsigned int i = 0;
        unsigned int count = 0;

        auto& knotSpans = Patchdata_->getknots();


        for (int j=0; j<dim; ++j)
        {
            while(count <= ijk[0])
            {
                if (i == knotSpans[0].size())
                    break;

                if (knotSpans[j][i+1] > knotSpans[j][i])
                    count++;

                ++i;
            }

        }

    return BSplineGeometry<dim,dimworld>(Patchdata_,&knotSpans[0][i-1]);
  }

private:

  std::shared_ptr <BsplinePatchData<dim,dimworld>> Patchdata_;

};


}

}

#endif  // DUNE_IGA_BSPLINEPATCH_HH