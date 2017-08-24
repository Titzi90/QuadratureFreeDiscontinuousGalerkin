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

  unsigned int numberCells = mesh_.getRows() * mesh_.getColumns() * 2;
  unsigned int localPoints=-1;
  int cellType=-1;
  switch (order_)
    {
    case 0:
    case 1:
      localPoints = 3;
      cellType = 5;
      break;
    case 2:
      localPoints = 6;
      cellType = 22;
      break;
    }
  unsigned int numberPoits = localPoints * numberCells;

  std::ofstream file(ss.str());
	if(file.is_open()){
	  file.precision(6);

    //header
    file << "# vtk DataFile Version 4.0\n"
         << "Quadfree DG: " << basename_ << "\n"
         << "ASCII\n";

    // Points
    file << "\nDATASET UNSTRUCTURED_GRID\n"
         << "POINTS " << numberPoits << " float\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          switch (order_)
            {
            case 0:
            case 1: // 0 and 1 -> 3 points
              //lower
              file << mesh_.getVertex(i,j) << " 0 "
                   << mesh_.getVertex(i+1,j) << " 0 "
                   << mesh_.getVertex(i, j+1) << " 0\n";
              //upper
              file << mesh_.getVertex(i+1,j) << " 0 "
                   << mesh_.getVertex(i+1,j+1) << " 0 "
                   << mesh_.getVertex(i, j+1) << " 0\n";
              break;
            case 2: // quadratic triangle -> 6 points
              //lower
              {
                Point const & A = mesh_.getVertex(i, j);
                Point const & B = mesh_.getVertex(i+1, j);
                Point const & C = mesh_.getVertex(i, j+1);
                file << A << " 0 " << B << " 0 " << C << " 0 "
                     << 0.5*(A+B) << " 0 " << 0.5*(B+C) << " 0 " << 0.5*(B+C) << " 0\n";
              }
              //upper
              {
                Point const & A = mesh_.getVertex(i+1, j);
                Point const & B = mesh_.getVertex(i+1, j+1);
                Point const & C = mesh_.getVertex(i, j+1);
                file << A << " 0 " << B << " 0 " << C << " 0 "
                     << 0.5*(A+B) << " 0 " << 0.5*(B+C) << " 0 " << 0.5*(B+C) << " 0\n";
              }
              break;
            }
        }

    // Topology
    file << "\nCELLS " << numberCells << " " << numberCells*(localPoints +1) << "\n";
    for (unsigned int c=0; c<numberCells; ++c)
      {
        file << localPoints;

        for (unsigned int p=c*localPoints; p<(c+1)*localPoints; ++p)
          file << " " << p;

        file << "\n";
      }
    file << "\nCELL_TYPES " << numberCells << "\n";
    for (unsigned int c=0; c<numberCells; ++c)
      file << cellType << "\n";

    // Data
    file << "\nPOINT_DATA " << numberPoits << "\n";
    // C
    file << "SCALARS C double 1\nLOOKUP_TABLE default\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          std::vector<double> data_l (lagrangeProjection(mesh_.getLower(i, j).C(), order_) );
          std::vector<double> data_u (lagrangeProjection(mesh_.getUpper(i, j).C(), order_) );

          switch (order_)
            {
            case 0:
              file << data_l[0] << " " << data_l[0] << " " << data_l[0] << "\n"
                   << data_u[0] << " " << data_u[0] << " " << data_u[0] << "\n";
              break;
            case 1:
              file << data_l[0] << " "
                   << data_l[1] << " "
                   << data_l[2] << "\n";
              file << data_u[0] << " "
                   << data_u[1] << " "
                   << data_u[2] << "\n";
              break;
            case 2:
              file << data_l[0] << " "
                   << data_l[1] << " "
                   << data_l[2] << " "
                   << data_l[3] << " "
                   << data_l[4] << " "
                   << data_l[5] << "\n";
              file << data_u[0] << " "
                   << data_u[1] << " "
                   << data_u[2] << " "
                   << data_u[3] << " "
                   << data_u[4] << " "
                   << data_u[5] << "\n";
              break;
            }
        }

    // U
    file << "\nVECTORS U double\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          std::vector<double> data1_l (lagrangeProjection(mesh_.getLower(i, j).U1(), order_) );
          std::vector<double> data2_l (lagrangeProjection(mesh_.getLower(i, j).U2(), order_) );
          std::vector<double> data1_u (lagrangeProjection(mesh_.getUpper(i, j).U1(), order_) );
          std::vector<double> data2_u (lagrangeProjection(mesh_.getUpper(i, j).U2(), order_) );

          switch (order_)
            {
            case 0:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[0] << " " << data2_l[0] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[0] << " " << data2_u[0] << " 0\n";
              break;
            case 1:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[1] << " " << data2_l[1] << " 0\n"
                   << data1_l[2] << " " << data2_l[2] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[1] << " " << data2_u[1] << " 0\n"
                   << data1_u[2] << " " << data2_u[2] << " 0\n";
              break;
            case 2:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[1] << " " << data2_l[1] << " 0\n"
                   << data1_l[2] << " " << data2_l[2] << " 0\n"
                   << data1_l[3] << " " << data2_l[3] << " 0\n"
                   << data1_l[4] << " " << data2_l[4] << " 0\n"
                   << data1_l[5] << " " << data2_l[5] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[1] << " " << data2_u[1] << " 0\n"
                   << data1_u[2] << " " << data2_u[2] << " 0\n"
                   << data1_u[3] << " " << data2_u[3] << " 0\n"
                   << data1_u[4] << " " << data2_u[4] << " 0\n"
                   << data1_u[5] << " " << data2_u[5] << " 0\n";
              break;
            }
        }

    // F
    file << "\nVECTORS F double\n";
    for (unsigned int i=0; i<mesh_.getRows(); ++i)
      for (unsigned int j=0; j<mesh_.getColumns(); ++j)
        {
          std::vector<double> data1_l (lagrangeProjection(mesh_.getLower(i, j).F1(), order_) );
          std::vector<double> data2_l (lagrangeProjection(mesh_.getLower(i, j).F2(), order_) );
          std::vector<double> data1_u (lagrangeProjection(mesh_.getUpper(i, j).F1(), order_) );
          std::vector<double> data2_u (lagrangeProjection(mesh_.getUpper(i, j).F2(), order_) );

          switch (order_)
            {
            case 0:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[0] << " " << data2_l[0] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[0] << " " << data2_u[0] << " 0\n";
              break;
            case 1:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[1] << " " << data2_l[1] << " 0\n"
                   << data1_l[2] << " " << data2_l[2] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[1] << " " << data2_u[1] << " 0\n"
                   << data1_u[2] << " " << data2_u[2] << " 0\n";
              break;
            case 2:
              file << data1_l[0] << " " << data2_l[0] << " 0\n"
                   << data1_l[1] << " " << data2_l[1] << " 0\n"
                   << data1_l[2] << " " << data2_l[2] << " 0\n"
                   << data1_l[3] << " " << data2_l[3] << " 0\n"
                   << data1_l[4] << " " << data2_l[4] << " 0\n"
                   << data1_l[5] << " " << data2_l[5] << " 0\n";
              file << data1_u[0] << " " << data2_u[0] << " 0\n"
                   << data1_u[1] << " " << data2_u[1] << " 0\n"
                   << data1_u[2] << " " << data2_u[2] << " 0\n"
                   << data1_u[3] << " " << data2_u[3] << " 0\n"
                   << data1_u[4] << " " << data2_u[4] << " 0\n"
                   << data1_u[5] << " " << data2_u[5] << " 0\n";
              break;
            }
        }


    // close file
		file.close();
  }
  else
    throw std::runtime_error("Couldn't open file " + ss.str());

}
