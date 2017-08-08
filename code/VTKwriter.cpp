#include "VTKwriter.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

// #include "Grid.hpp"
// #include "DataTypes.hpp"



void VTKwriter::write()
{
  std::stringstream ss;
	ss << basename_ << "_" << iteration_++ << ".vtk";

  std::ofstream file(ss.str());
	if(file.is_open()){
	  file.precision(6);

    //header
    file << "# vtk DataFile Version 4.0\n"
         << "Quadfree DG: " << basename_ << "\n"
         << "ASCII\n";

    // Points
    file << "\nDATASET UNSTRUCTURED_GRID\n"
         << "POINTS " << (mesh_.getRows()+1)*(mesh_.getColumns()+1)  << " float\n";
    for (unsigned int i=0; i<mesh_.getRows()+1; ++i)
      for (unsigned int j=0; j<mesh_.getColumns()+1; ++j)
          file << mesh_.getVertex(i,j) << " 0\n";

    // Topology
    unsigned int numberCells = mesh_.getRows() * mesh_.getColumns() * 2;
    file << "\nCELLS " << numberCells << " " << numberCells*4 << "\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          file << "3 " <<  i   *(mesh_.getColumns()+1) + j
               << " "  <<  i   *(mesh_.getColumns()+1) + j + 1
               << " "  << (i+1)*(mesh_.getColumns()+1) + j     << "\n"; // lower
          file << "3 " <<  i   *(mesh_.getColumns()+1) + j + 1
               << " "  << (i+1)*(mesh_.getColumns()+1) + j + 1
               << " "  << (i+1)*(mesh_.getColumns()+1) + j     << "\n"; // upper
        }

    file << "\nCELL_TYPES " << numberCells << "\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        file << "5 5\n";

    // Data
    file << "\nCELL_DATA " << numberCells << "\n";

    // C
    file << "SCALARS C float 1\nLOOKUP_TABLE default\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          file << (mesh_.getLower(i,j).C()[0]*pol::phi[0]).get(0,0) << "\n";
          file << (mesh_.getUpper(i,j).C()[0]*pol::phi[0]).get(0,0) << "\n";
        }

    // U
    file << "\nVECTORS U float\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          file << (mesh_.getLower(i,j).U1()[0]*pol::phi[0]).get(0,0) << " "
               << (mesh_.getLower(i,j).U2()[0]*pol::phi[0]).get(0,0) << " 0\n";
          file << (mesh_.getUpper(i,j).U1()[0]*pol::phi[0]).get(0,0) << " "
               << (mesh_.getUpper(i,j).U2()[0]*pol::phi[0]).get(0,0) << " 0\n";
        }

    // U
    file << "\nVECTORS F float\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          file << (mesh_.getLower(i,j).F1()[0]*pol::phi[0]).get(0,0) << " "
               << (mesh_.getLower(i,j).F2()[0]*pol::phi[0]).get(0,0) << " 0\n";
          file << (mesh_.getUpper(i,j).F1()[0]*pol::phi[0]).get(0,0) << " "
               << (mesh_.getUpper(i,j).F2()[0]*pol::phi[0]).get(0,0) << " 0\n";
        }


    // close file
		file.close();
  }
  else
    throw std::runtime_error("Couldn't open file " + ss.str());

}
