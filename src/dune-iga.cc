#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/vtkfile.hh>

using namespace Dune;

int main(int argc, char** argv) try
{
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

  // parameters
  unsigned int subSampling = 2;

  ////////////////////////////////////////////////////////////////
  //  Create a 2d B-spline surface in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 2;
  const int dimworld = 3;

  FieldVector<double,dim> lower = {0.0, 0.0};
  FieldVector<double,dim> upper = {1.0, 1.0};

  std::array<unsigned int,dim> knotSpans = {2, 2};

  // very simple geometry: a square
  std::vector<FieldVector<double,dimworld> > controlPoints = {{0,0,0},    {1,0,0.25}, {2,0,0},
                                                              {0,1,0.25}, {1,1,0.5},  {2,1,0.25},
                                                              {0,2,0},    {1,2,0.25}, {2,2,0}};

  IGA::BSplinePatch<dim,dimworld> patch(lower, upper, knotSpans, controlPoints);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  std::size_t offset = 0;

  for (unsigned int i=0; i<knotSpans[1]; i++)
  {
    for (unsigned int j=0; j<knotSpans[0]; j++)
    {
      auto geometry = patch.geometry({j,i});

      // Add vertex coordinates to the VTK file
      for (int iy=0; iy<=(1<<subSampling); iy++)
      {
        FieldVector<double,dim> localPos;
        localPos[1] = ((double)iy)/(1<<subSampling);

        // Add vertex coordinates to the VTK file
        for (int ix=0; ix<=(1<<subSampling); ix++)
        {
          localPos[0] = ((double)ix)/(1<<subSampling);

          vtkFile.points_.push_back(geometry.global(localPos));
        }
      }

      // Add elements to the VTK file
      for (int k=0; k<(1<<subSampling); k++)
      {
        for (int l=0; l<(1<<subSampling); l++)
        {
          vtkFile.cellConnectivity_.push_back(offset + k    *((1<<subSampling)+1) + l);
          vtkFile.cellConnectivity_.push_back(offset + k    *((1<<subSampling)+1) + l+1);
          vtkFile.cellConnectivity_.push_back(offset + (k+1)*((1<<subSampling)+1) + l+1);
          vtkFile.cellConnectivity_.push_back(offset + (k+1)*((1<<subSampling)+1) + l);

          // 4 corners per element
          if (vtkFile.cellOffsets_.size()==0)
            vtkFile.cellOffsets_.push_back(4);
          else
            vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+4);

          // Element type: a 4-node quadrilateral
          vtkFile.cellTypes_.push_back(9);
        }
      }

      offset += ((1<<subSampling)+1) * ((1<<subSampling)+1);
    }

    // Actually write the VTK file
    vtkFile.write("bsplinesurface");
  }

  return 0;
}
catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
}
