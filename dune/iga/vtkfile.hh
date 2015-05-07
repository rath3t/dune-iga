#ifndef DUNE_IGA_VTKFILE_HH
#define DUNE_IGA_VTKFILE_HH

#include <vector>
#include <fstream>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/dataarraywriter.hh>

namespace Dune {

  namespace IGA {

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
        outFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        outFile << "  <UnstructuredGrid>" << std::endl;
        outFile << "    <Piece NumberOfCells=\"" << cellOffsets_.size() << "\" NumberOfPoints=\"" << points_.size() << "\">" << std::endl;

        // Write vertex coordinates
        outFile << "      <Points>" << std::endl;
        {  // extra parenthesis to control destruction of the pointsWriter object
          Dune::VTK::AsciiDataArrayWriter<float> pointsWriter(outFile, "Coordinates", 3, Dune::Indent(4));
          for (size_t i=0; i<points_.size(); i++)
            for (int j=0; j<3; j++)
              pointsWriter.write(points_[i][j]);
        }  // destructor of pointsWriter objects writes trailing </DataArray> to file
        outFile << "      </Points>" << std::endl;

        // Write elements
        outFile << "      <Cells>" << std::endl;
        {  // extra parenthesis to control destruction of the cellConnectivityWriter object
          Dune::VTK::AsciiDataArrayWriter<int> cellConnectivityWriter(outFile, "connectivity", 1, Dune::Indent(4));
          for (size_t i=0; i<cellConnectivity_.size(); i++)
            cellConnectivityWriter.write(cellConnectivity_[i]);
        }

        {  // extra parenthesis to control destruction of the writer object
          Dune::VTK::AsciiDataArrayWriter<int> cellOffsetsWriter(outFile, "offsets", 1, Dune::Indent(4));
          for (size_t i=0; i<cellOffsets_.size(); i++)
            cellOffsetsWriter.write(cellOffsets_[i]);
        }

        {  // extra parenthesis to control destruction of the writer object
          Dune::VTK::AsciiDataArrayWriter<unsigned int> cellTypesWriter(outFile, "types", 1, Dune::Indent(4));
          for (size_t i=0; i<cellTypes_.size(); i++)
            cellTypesWriter.write(cellTypes_[i]);
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

}

#endif
