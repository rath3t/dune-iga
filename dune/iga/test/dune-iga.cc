// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/iga/NURBSgrid.hh>
#include <dune/iga/NURBSpatch.hh>
#include <dune/iga/bsplinepatch.hh>
//#include <dune/iga/vtkfile.hh>

#include <dune/common/test/testsuite.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

using namespace Dune;
using namespace Dune::IGA;



#ifndef DUNE_IGA_VTKFILE_HH
#define DUNE_IGA_VTKFILE_HH

#include <vector>
#include <fstream>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/dataarraywriter.hh>

namespace Dune::IGA {

  /** \brief A class representing a VTK file, but independent from the Dune grid interface
     *
     * This file is supposed to represent an abstraction layer in between the pure XML used for VTK files,
     * and the VTKWriter from dune-grid, which knows about grids.  In theory, the dune-grid VTKWriter
     * could use this class to simplify its own code.  More importantly, the VTKFile class allows to
     * write files containing second-order geometries, which is currently not possible with the dune-grid
     * VTKWriter.
   */
  class VTKFile
  {

  public:

    /** \brief Write the file to disk */
    void write(const std::string& filename) const
    {
      std::ofstream outFile(filename + ".vtu");

      // Write header
      outFile << "<?xml version=\"1.0\"?>" << std::endl;
      outFile << R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">)" << std::endl;
      outFile << "  <UnstructuredGrid>" << std::endl;
      outFile << "    <Piece NumberOfCells=\"" << cellOffsets_.size() << "\" NumberOfPoints=\"" << points_.size() << "\">" << std::endl;

      // Write vertex coordinates
      outFile << "      <Points>" << std::endl;
      {  // extra parenthesis to control destruction of the pointsWriter object
        Dune::VTK::AsciiDataArrayWriter pointsWriter(outFile, "Coordinates", 3, Dune::Indent(4),VTK::Precision::float32);
        for (const auto & point : points_)
          for (int j=0; j<3; j++)
            pointsWriter.write(point[j]);
      }  // destructor of pointsWriter objects writes trailing </DataArray> to file
      outFile << "      </Points>" << std::endl;

      // Write elements
      outFile << "      <Cells>" << std::endl;
      {  // extra parenthesis to control destruction of the cellConnectivityWriter object
        Dune::VTK::AsciiDataArrayWriter cellConnectivityWriter(outFile, "connectivity", 1, Dune::Indent(4),VTK::Precision::int32);
        for (int i : cellConnectivity_)
          cellConnectivityWriter.write(i);
      }

      {  // extra parenthesis to control destruction of the writer object
        Dune::VTK::AsciiDataArrayWriter cellOffsetsWriter(outFile, "offsets", 1, Dune::Indent(4),VTK::Precision::int32);
        for (int cellOffset : cellOffsets_)
          cellOffsetsWriter.write(cellOffset);
      }

      {  // extra parenthesis to control destruction of the writer object
        Dune::VTK::AsciiDataArrayWriter cellTypesWriter(outFile, "types", 1, Dune::Indent(4),VTK::Precision::uint32);
        for (int cellType : cellTypes_)
          cellTypesWriter.write(cellType);
      }

      outFile << "      </Cells>" << std::endl;

      //////////////////////////////////////////////////
      //   Write footer
      //////////////////////////////////////////////////
      outFile << "    </Piece>" << std::endl;
      outFile << "  </UnstructuredGrid>" << std::endl;
      outFile << "</VTKFile>" << std::endl;


    }

    std::vector<Dune::FieldVector<double,3> > points_;

    std::vector<int> cellConnectivity_;

    std::vector<int> cellOffsets_;

    std::vector<int> cellTypes_;

  };

}

#endif








void testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};
//  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<ControlPoint > controlPoints
      = {{.p = {1, 3, 4}, .w = 1},
         {.p = {2, 2, 2}, .w = 3},
         {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 2},
         {.p = {4, 7, 2}, .w = 1},
         {.p = {8, 6, 2}, .w = 1},
         {.p = {2, 9, 9}, .w = 7},
         {.p = {1, 4, 3}, .w = 1},
         {.p = {1, 7, 1}, .w = 5}};

  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();


  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
//  vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-CurveNewFineResample_knotRefine");
}

void testNURBSGridSurface() {
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim                    = 2;
  const int dimworld               = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
//  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
//      = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
 using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0, 0, 1}, .w = 2},        {.p = {1, 0, 1}, .w = 2},     {.p = {2, 0, 2}, .w = 1}},
          {{.p = {0, 1, 0}, .w = 1},        {.p = {1, 1, 0}, .w = 4},     {.p = {2, 1, 0}, .w = 1}},
          {{.p = {0, 2, 1}, .w = 1},        {.p = {1, 2, 2}, .w = 2},     {.p = {2, 2, 2}, .w = 4}}};

  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");
}


void testNURBSSurface() {
  // parameters
  int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const int dim                    = 2;
  const int dimworld               = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0, 0, 1}, .w = 2},        {.p = {1, 0, 1}, .w = 2},     {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1},        {.p = {1, 1, 0}, .w = 4},     {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1},        {.p = {1, 2, 2}, .w = 2},     {.p = {2, 2, 2}, .w = 4}}};


  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
  //

//  auto weightNet  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 1);
  testSuite.check(patch.size(1) == 4);
  testSuite.check(patch.size(2) == 9);
}

void testNURBSCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
  const std::vector<ControlPoint > controlPoints
      = {{.p = {1, 3, 4}, .w = 2},
         {.p = {2, 2, 2}, .w = 2},
         {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 1},
         {.p = {4, 7, 2}, .w = 4},
         {.p = {8, 6, 2}, .w = 2},
         {.p = {2, 9, 9}, .w = 1},
         {.p = {1, 4, 3}, .w = 2},
         {.p = {1, 7, 1}, .w = 4}};


  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 5);
  testSuite.check(patch.size(1) == controlPoints.size());
}

void testBSplineCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::GlobalCoordinateType;
  const std::vector<ControlPoint > controlPoints
      = {{ {1, 3, 4}},
         { {2, 2, 2}},
         { {3, 4, 5}},
         { {5, 1, 7}},
         { {4, 7, 2}},
         { {8, 6, 2}},
         { {2, 9, 9}},
         { {1, 4, 3}},
         { {1, 7, 1}}};

  std::array<unsigned int, dim> dimsize = {static_cast<unsigned int>(controlPoints.size())};
  auto controlNet                       = MultiDimensionNet<dim,typename Dune::IGA::NURBSPatchData<dim,dimworld>::GlobalCoordinateType>(dimsize, controlPoints);

  IGA::BSplinePatch<dim, dimworld> patch(knotSpans, controlNet, order);
}

void testNurbsGridCylinder() {
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  const int dim                    = 2;
  const int dimworld               = 3;
  const std::array<int, dim> order = {2, 1};
  constexpr double invsqr2         = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l                                       = 10;
  const double rad                                     = 5;
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};


  using ControlPoint = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint >> controlPoints
      =  {{{.p = {0,   0, rad}, .w =       1},  {.p = {0,   l, rad}, .w = 1      }},
          {{.p = {rad, 0, rad}, .w = invsqr2},  {.p = {rad, l, rad}, .w = invsqr2}},
          {{.p = {rad, 0,   0}, .w =       1},  {.p = {rad, l,   0}, .w = 1      }}};

  std::array<unsigned int, dim> dimsize
      = {static_cast<unsigned int>(controlPoints.size()), static_cast<unsigned int>(controlPoints[0].size())};
  auto controlNet = Dune::IGA::NURBSPatchData<dim,dimworld>::ControlPointNetType(dimsize, controlPoints);

  IGA::NURBSGrid<dim, dimworld> grid(knotSpans, controlNet, order);
  grid.globalRefine(2);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, order);

  testSuite.check(nurbsPatch.size(0) == 1);
  testSuite.check(nurbsPatch.size(1) == 4);
  testSuite.check(nurbsPatch.size(2) == 6);

  ////////////////////////////////////////////////////////////////
  //  Write to a VTK file.
  //  The higher-order geometry is captured by subsampling.
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

//  vtkWriter.write("Zyl");
  vtkWriter.write("ZylRefine");


//  auto gridViewRefined  = grid.leafGridView();
//
//  SubsamplingVTKWriter<decltype(gridViewRefined)> vtkWriter2(gridViewRefined, refinementIntervals1);
//  vtkWriter2.write("ZylRefine");


//  IGA::VTKFile vtkFile;
//
//  //  The number of vertices that have been inserted so far
//  std::size_t offset = 0;
//
//  //Range-based for loop to get each element and its corresponding geometry
//  for (auto const &element: elements(gridView))
//  {
//    auto geometry = element.geometry();
//    std::cout<<"Element index: "<<indexSet.index(element)<<std::endl;
//    //Add vertex coordinates to the VTK file
//    for (int iy=0; iy<=(1<<subSampling); iy++)
//    {
//      FieldVector<double,dim> localPos;
//      localPos[1] = ((double)iy)/(1<<subSampling);
//
//
//      // Add vertex coordinates to the VTK file
//      for (int ix=0; ix<=(1<<subSampling); ix++)
//      {
//        localPos[0] = ((double)ix)/(1<<subSampling);
//        std::cout<<"localPos: "<<localPos<<std::endl;
//        vtkFile.points_.push_back(geometry.global(localPos));
//        std::cout<<"localPos: "<<geometry.global(localPos)<<std::endl;
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
//        if (vtkFile.cellOffsets_.empty())
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
//  vtkFile.write("ZylByHand");

  ////////////////////////////////////////////////////////////////
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

  std::cout << "test NURBS grid surface" << std::endl;
  testNURBSGridSurface();
  std::cout << "done with NURBS grid surface" << std::endl;

  std::cout << "test NURBS grid Curve" << std::endl;
  testNURBSGridCurve();
  std::cout << "done with NURBS grid Curve" << std::endl;

  std::cout << "test NURBS grid Curve" << std::endl;
  testNURBSCurve();
  std::cout << "done with NURBS curve" << std::endl;

  testNURBSSurface();
  std::cout << "done with NURBS surface" << std::endl;

  testNurbsGridCylinder();
  std::cout << "done with NURBS surface cylinder" << std::endl;
  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
