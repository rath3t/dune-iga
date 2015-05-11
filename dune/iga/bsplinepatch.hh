#ifndef DUNE_IGA_BSPLINEPATCH_HH
#define DUNE_IGA_BSPLINEPATCH_HH

#include <dune/geometry/multilineargeometry.hh>

namespace Dune
{
namespace IGA
{

template <int dim, int dimworld>
class BSplineGeometry
{
public:

  BSplineGeometry(const std::vector<FieldVector<double,dimworld> >& corners)
  : multiLinearGeometry_(GeometryType(GeometryType::BasicType::cube,dim), corners)
  {}

  FieldVector<double,dimworld> global(const FieldVector<double,dim>& local) const
  {
    return multiLinearGeometry_.global(local);
  }

  MultiLinearGeometry<double,dim,dimworld> multiLinearGeometry_;
};

template <int dim, int dimworld>
class BSplinePatch
{
public:

  BSplinePatch(const FieldVector<double,dim>& lower,
               const FieldVector<double,dim>& upper,
               const std::array<unsigned int,dim>& knotSpans,
               const std::vector<FieldVector<double,dimworld> > controlPoints)
  : knotSpans_(knotSpans),
    controlPoints_(controlPoints)
  {

  }

  BSplineGeometry<dim,dimworld> geometry(const std::array<unsigned int,dim>& ijk) const
  {
    std::vector<FieldVector<double,dimworld> > corners(1<<dim);

    switch (dim)
    {
      case 1:
      {
        corners[0] = controlPoints_[ijk[0]];
        corners[1] = controlPoints_[ijk[0]+1];

        break;
      }
      case 2:
      {
        corners[0] = controlPoints_[ijk[1]*(knotSpans_[0]+1) + ijk[0]];
        corners[1] = controlPoints_[ijk[1]*(knotSpans_[0]+1) + ijk[0]+1];
        corners[2] = controlPoints_[(ijk[1]+1)*(knotSpans_[0]+1) + ijk[0]];
        corners[3] = controlPoints_[(ijk[1]+1)*(knotSpans_[0]+1) + ijk[0]+1];

        break;
      }
      default:
        DUNE_THROW(Dune::NotImplemented, "General case not implemented yet!");
    }
    return BSplineGeometry<dim,dimworld>(corners);
  }

private:

  const std::array<unsigned int,dim>& knotSpans_;

  std::vector<FieldVector<double,dimworld> > controlPoints_;

};

}

}

#endif  // DUNE_IGA_BSPLINEPATCH_HH