// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>


#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/iga/bsplinegrid.hh>
#include <dune/iga/bsplinepatch.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/NURBSgrid.hh>
#include <dune/iga/vtkfile.hh>


#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

using namespace Dune;

void testBSplineGridSurface()
{
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  //parameters
  unsigned int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  const int dim      = 2;
  const int dimworld = 3;
  const std::array<int,dim> order = {2,2};
  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,1},{0,0,0,1,1,1}}};
  const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints = {{{{0,0,1},{1,0,1},{2,0,2}}, {{0,1,0},{1,1,0},{2,1,0}}, {{0,2,1},{1,2,2},{2,2,2}} }};

  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size()),static_cast<unsigned int>(controlPoints[0].size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);
  //controlNet.disp();

  IGA::BSplineGrid<dim,dimworld> grid(knotSpans, controlNet, order);
  const auto gridView = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();
  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  //  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView,subSampling);
  //  vtkWriter.write("foo");
  ////////////////////////////////////////////////////////////////

  //! Original manual VTK writer, it works fine. Uncomment to have a working version
//  IGA::VTKFile vtkFile;
//
//  //  The number of vertices that have been inserted so far
//  std::size_t offset = 0;
//
//  //Range-based for loop to get each element and its corresponding geometry
//  for (auto const &element:gridView)
//  {
//    auto geometry = element.geometry();
//    std::cout<<"Element index: "<<indexSet.index(element)<<std::endl;
//    //Add vertex coordinates to the VTK file
//    for (int iy=0; iy<=(1<<subSampling); iy++)
//    {
//      FieldVector<double,dim> localPos;
//      localPos[1] = ((double)iy)/(1<<subSampling);
//
//      // Add vertex coordinates to the VTK file
//      for (int ix=0; ix<=(1<<subSampling); ix++)
//      {
//        localPos[0] = ((double)ix)/(1<<subSampling);
//
//        vtkFile.points_.push_back(geometry.global(localPos));
//      }
//    }
//
//    // Add elements to the VTK file
//    for (int k=0; k<(1<<subSampling); k++)
//    {
//      for (int l=0; l<(1<<subSampling); l++)
//      {
//        vtkFile.cellConnectivity_.push_back(offset + k    *((1<<subSampling)+1) + l);
//        vtkFile.cellConnectivity_.push_back(offset + k    *((1<<subSampling)+1) + l+1);
//        vtkFile.cellConnectivity_.push_back(offset + (k+1)*((1<<subSampling)+1) + l+1);
//        vtkFile.cellConnectivity_.push_back(offset + (k+1)*((1<<subSampling)+1) + l);
//
//        // 4 corners per element
//        if (vtkFile.cellOffsets_.size()==0)
//          vtkFile.cellOffsets_.push_back(4);
//        else
//          vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+4);
//
//        // Element type: a 4-node quadrilateral
//        vtkFile.cellTypes_.push_back(9);
//      }
//    }
//
//    offset += ((1<<subSampling)+1) * ((1<<subSampling)+1);
//  }
//  vtkFile.write("BSplineGridTest-Surface");

}

void testBSplineGridCurve()
{
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  //parameters
  unsigned int subSampling = 5;

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int,dim> order = {2};
  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,2,3,4,4,5,5,5}}};

  const std::vector<FieldVector<double,dimworld> > controlPoints = {{1,3,4}, {2,2,2}, {3,4,5}, {5,1,7}, {4,7,2}, {8,6,2}, {2,9,9}, {1,4,3},{1,7,1}};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);

  IGA::BSplineGrid<dim,dimworld> grid(knotSpans, controlNet, order);
  const auto& gridView = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  int offset = 0;

  //Range-based for loop to get each element and its corresponding geometry
  for (auto const &element:gridView)
  {
    std::cout<<"Element index: "<<indexSet.index(element)<<std::endl;
    auto geometry = element.geometry();
    FieldVector<double,dim> localPos;

    // Add vertex coordinates to the VTK file
    for (int ix=0; ix<=(1<<subSampling); ix++)
    {
      localPos[0] = ((double)ix)/(1<<subSampling);
      vtkFile.points_.push_back(geometry.global(localPos));
      //std::cout << "Jacob transposed:" << geometry.jacobianTransposed(localPos);
    }

    // Add elements to the VTK file
    for (int l=0; l<(1<<subSampling); l++)
    {
      vtkFile.cellConnectivity_.push_back(offset + l);
      vtkFile.cellConnectivity_.push_back(offset + l+1);

      // 2 corners per element
      if (vtkFile.cellOffsets_.empty())
        vtkFile.cellOffsets_.push_back(2);
      else
        vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+2);

      // Element type: a line segment
      vtkFile.cellTypes_.push_back(3);
    }

    offset += ((1<<subSampling)+1);
  }

  // Actually write the VTK file
  vtkFile.write("BSplineGridTest-Curve");



}

void testNURBSGridCurve()
{
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  //parameters
  unsigned int subSampling = 5;

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int,dim> order = {2};
  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,2,3,4,4,5,5,5}}};

  const std::vector<FieldVector<double,dimworld> > controlPoints = {{1,3,4}, {2,2,2}, {3,4,5}, {5,1,7}, {4,7,2}, {8,6,2}, {2,9,9}, {1,4,3},{1,7,1}};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);

  const std::vector<FieldVector<double,1> > weight = {{2},{2},{1},{1},{4},{2},{1},{2},{4}};
  //auto weightNet = MultiDimensionNet<dim,1>(dimsize,FieldVector<double,1>(1));
  auto weightNet = MultiDimensionNet<dim,1>(dimsize,weight);

  IGA::NURBSGrid<dim,dimworld> grid(knotSpans, controlNet, weightNet, order);
  const auto& gridView = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  int offset = 0;

  //Range-based for loop to get each element and its corresponding geometry
  for (auto const &element:gridView)
  {
    std::cout<<"Element index: "<<indexSet.index(element)<<std::endl;
    auto geometry = element.geometry();
    FieldVector<double,dim> localPos;

    // Add vertex coordinates to the VTK file
    for (int ix=0; ix<=(1<<subSampling); ix++)
    {
      localPos[0] = ((double)ix)/(1<<subSampling);
      vtkFile.points_.push_back(geometry.global(localPos));
      //std::cout << "Jacob transposed:" << geometry.jacobianTransposed(localPos);
    }

    // Add elements to the VTK file
    for (int l=0; l<(1<<subSampling); l++)
    {
      vtkFile.cellConnectivity_.push_back(offset + l);
      vtkFile.cellConnectivity_.push_back(offset + l+1);

      // 2 corners per element
      if (vtkFile.cellOffsets_.empty())
        vtkFile.cellOffsets_.push_back(2);
      else
        vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+2);

      // Element type: a line segment
      vtkFile.cellTypes_.push_back(3);
    }

    offset += ((1<<subSampling)+1);
  }

  // Actually write the VTK file
  vtkFile.write("NURBSGridTest-Curve");



}

void testNURBSGridSurface()
{
  unsigned int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim      = 2;
  const int dimworld = 3;
  const std::array<int,dim> order = {2,2};

  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,1},{0,0,0,1,1,1}}};
  const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints = {{{{0,0,1},{1,0,1},{2,0,2}}, {{0,1,0},{1,1,0},{2,1,0}}, {{0,2,1},{1,2,2},{2,2,2}} }};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size()),static_cast<unsigned int>(controlPoints[0].size())};
  //
  const std::vector<std::vector<FieldVector<double,1> > > weight = {{{2},{2},{1}}, {{1},{4},{1}}, {{1},{2},{4}}};
  //auto weightNet = MultiDimensionNet<dim,1>(dimsize,FieldVector<double,1>(1));
  auto weightNet = MultiDimensionNet<dim,1>(dimsize,weight);
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);

  IGA::NURBSGrid<dim,dimworld> grid(knotSpans, controlNet, weightNet, order);
  const auto& gridView = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

 ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////


  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  std::size_t offset = 0;

  //Range-based for loop to get each element and its corresponding geometry
  for (auto const &element:gridView)
  {
    auto geometry = element.geometry();
    std::cout<<"Element index: "<<indexSet.index(element)<<std::endl;
    //Add vertex coordinates to the VTK file
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
        if (vtkFile.cellOffsets_.empty())
          vtkFile.cellOffsets_.push_back(4);
        else
          vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+4);

        // Element type: a 4-node quadrilateral
        vtkFile.cellTypes_.push_back(9);
      }
    }

    offset += ((1<<subSampling)+1) * ((1<<subSampling)+1);
  }
  vtkFile.write("NURBSGridTest-Surface");
}



void testNURBSSurface()
{
  //parameters
  unsigned int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim      = 2;
  const int dimworld = 3;
  const std::array<int,dim> order = {2,2};

  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,1},{0,0,0,1,1,1}}};
  const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints = {{{{0,0,1},{1,0,1},{2,0,2}}, {{0,1,0},{1,1,0},{2,1,0}}, {{0,2,1},{1,2,2},{2,2,2}} }};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size()),static_cast<unsigned int>(controlPoints[0].size())};
  //
  const std::vector<std::vector<FieldVector<double,1> > > weight = {{{2},{2},{1}}, {{1},{4},{1}}, {{1},{2},{4}}};
  //auto weightNet = MultiDimensionNet<dim,1>(dimsize,FieldVector<double,1>(1));
  auto weightNet = MultiDimensionNet<dim,1>(dimsize,weight);
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);
  //weightNet.disp();
  //controlNet.disp();

  IGA::NURBSPatch<dim,dimworld> patch(knotSpans, controlNet, weightNet, order);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  std::size_t offset = 0;

  const auto validknotes = patch.validKnotSize();

  for (unsigned int j=0; j<validknotes[1]; ++j)
  {
    for (unsigned int i=0; i<validknotes[0]; i++)
    {
      auto geometry = patch.geometry({i,j});

      //Add vertex coordinates to the VTK file
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
          if (vtkFile.cellOffsets_.empty())
            vtkFile.cellOffsets_.push_back(4);
          else
            vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+4);

          // Element type: a 4-node quadrilateral
          vtkFile.cellTypes_.push_back(9);
        }
      }

      offset += ((1<<subSampling)+1) * ((1<<subSampling)+1);
    }
  }

  // Actually write the VTK file
  vtkFile.write("NURBSsurface");
 }

void testBSplineSurface()
{
  //parameters
  unsigned int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim      = 2;
  const int dimworld = 3;
  const std::array<int,dim> order = {2,2};

  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,1},{0,0,0,1,1,1}}};

  const std::vector<std::vector<FieldVector<double,dimworld> > > controlPoints = {{{{0,0,1},{1,0,1},{2,0,2}}, {{0,1,0},{1,1,0},{2,1,0}}, {{0,2,1},{1,2,2},{2,2,2}} }};

  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size()),static_cast<unsigned int>(controlPoints[0].size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);
  //controlNet.disp();

  IGA::BSplinePatch<dim,dimworld> patch(knotSpans, controlNet, order);


  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  IGA::VTKFile vtkFile;

  //  The number of vertices that have been inserted so far
  std::size_t offset = 0;

  const auto validknotes = patch.validKnotSize();

  for (unsigned int j=0; j<validknotes[1]; ++j)
  {
    for (unsigned int i=0; i<validknotes[0]; i++)
    {
      auto geometry = patch.geometry({i,j});

      //Add vertex coordinates to the VTK file
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
          if (vtkFile.cellOffsets_.empty())
            vtkFile.cellOffsets_.push_back(4);
          else
            vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+4);

          // Element type: a 4-node quadrilateral
          vtkFile.cellTypes_.push_back(9);
        }
      }

      offset += ((1<<subSampling)+1) * ((1<<subSampling)+1);
    }
  }
  // Actually write the VTK file
  vtkFile.write("bsplinesurface");
 }

void testNURBSCurve()
{
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int,dim> order = {2};
  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,2,3,4,4,5,5,5}}};

  const std::vector<FieldVector<double,dimworld> > controlPoints = {{1,3,4}, {2,2,2}, {3,4,5}, {5,1,7}, {4,7,2}, {8,6,2}, {2,9,9}, {1,4,3},{1,7,1}};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);

  //auto weightNet = MultiDimensionNet<dim,1>(dimsize,FieldVector<double,1>(1));
  const std::vector<FieldVector<double,1> > weight = {{2},{2},{1},{1},{4},{2},{1},{2},{4}};
  //auto weightNet = MultiDimensionNet<dim,1>(dimsize,FieldVector<double,1>(1));
  auto weightNet = MultiDimensionNet<dim,1>(dimsize,weight);
  //controlNet.disp();

  IGA::NURBSPatch<dim,dimworld> patch(knotSpans, controlNet, weightNet, order);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.                                      //
  //  The higher-order geometry is captured by subsampling.     //
  ////////////////////////////////////////////////////////////////

 IGA::VTKFile vtkFile;

 //  The number of vertices that have been inserted so far
 std::size_t offset = 0;

  const auto validknots = patch.validKnotSize();

  for (unsigned int i=0; i<validknots[0]; i++)
  {
    auto geometry = patch.geometry({i});

    FieldVector<double,dim> localPos;

    // Add vertex coordinates to the VTK file
    for (int ix=0; ix<=(1<<subSampling); ix++)
    {
      localPos[0] = ((double)ix)/(1<<subSampling);
      vtkFile.points_.push_back(geometry.global(localPos));
      std::cout << "Jacob transposed:" << geometry.jacobianTransposed(localPos);
    }

    // Add elements to the VTK file
    for (int l=0; l<(1<<subSampling); l++)
    {
      vtkFile.cellConnectivity_.push_back(offset + l);
      vtkFile.cellConnectivity_.push_back(offset + l+1);

      // 2 corners per element
      if (vtkFile.cellOffsets_.empty())
        vtkFile.cellOffsets_.push_back(2);
      else
        vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+2);

      // Element type: a line segment
      vtkFile.cellTypes_.push_back(3);
    }

    offset += ((1<<subSampling)+1);
  }

  // Actually write the VTK file
  vtkFile.write("bsplinecurve");
}

void testBSplineCurve()
{
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int,dim> order = {2};
  const std::array<std::vector<double>,dim> knotSpans = {{{0,0,0,1,1,2,3,4,4,5,5,5}}};

  const std::vector<FieldVector<double,dimworld> > controlPoints = {{1,3,4}, {2,2,2}, {3,4,5}, {5,1,7}, {4,7,2}, {8,6,2}, {2,9,9}, {1,4,3},{1,7,1}};
  std::array<unsigned int,dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet = MultiDimensionNet<dim,dimworld>(dimsize,controlPoints);
  //controlNet.disp();

  IGA::BSplinePatch<dim,dimworld> patch(knotSpans, controlNet, order);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.                                      //
  //  The higher-order geometry is captured by subsampling.     //
  ////////////////////////////////////////////////////////////////

 IGA::VTKFile vtkFile;

 //  The number of vertices that have been inserted so far
 std::size_t offset = 0;

  const auto validknots = patch.validKnotSize();

  for (unsigned int i=0; i<validknots[0]; i++)
  {
    auto geometry = patch.geometry({i});

    FieldVector<double,dim> localPos;

    // Add vertex coordinates to the VTK file
    for (int ix=0; ix<=(1<<subSampling); ix++)
    {
      localPos[0] = ((double)ix)/(1<<subSampling);
      vtkFile.points_.push_back(geometry.global(localPos));
      std::cout << "Jacob transposed:" << geometry.jacobianTransposed(localPos);
    }

    // Add elements to the VTK file
    for (int l=0; l<(1<<subSampling); l++)
    {
      vtkFile.cellConnectivity_.push_back(offset + l);
      vtkFile.cellConnectivity_.push_back(offset + l+1);

      // 2 corners per element
      if (vtkFile.cellOffsets_.empty())
        vtkFile.cellOffsets_.push_back(2);
      else
        vtkFile.cellOffsets_.push_back(vtkFile.cellOffsets_.back()+2);

      // Element type: a line segment
      vtkFile.cellTypes_.push_back(3);
    }

    offset += ((1<<subSampling)+1);
  }

  // Actually write the VTK file
  vtkFile.write("bsplinecurve");
}


int main(int argc, char** argv) try
{

  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  std::cout<<"test B-Spline grid Surface"<<std::endl;
  testBSplineGridSurface();
  //testBSplineCurve();
  std::cout<< "done with B-Spline grid Surface" << std::endl;

//  std::cout<<"test B-Spline grid Curve"<<std::endl;
  testBSplineGridCurve();
  //testBSplineCurve();
  std::cout<< "done with B-Spline grid Curve" << std::endl;

  std::cout<<"test NURBS grid surface" <<std::endl;
  testNURBSGridSurface();
  std::cout<< "done with NURBS grid surface" << std::endl;

  std::cout<<"test NURBS grid Curve"<<std::endl;
  testNURBSGridCurve();
  //testBSplineCurve();
  std::cout<< "done with NURBS grid Curve" << std::endl;

  testNURBSCurve();
  std::cout<< "done with NURBS curve" << std::endl;

  testBSplineSurface();
  std::cout<< "done with surface" << std::endl;

  testNURBSSurface();
  std::cout<< "done with NURBS surface" << std::endl;

  return 0;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << std::endl;
}
