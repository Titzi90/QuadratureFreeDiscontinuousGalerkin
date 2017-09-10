
#include <cmath>
#include <cstdlib>

#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "stepper.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"



int main(int argc, char** argv)
{
  int order = 1;
  int orderF = 2*order;
  int refiment = 64;
  double tEnd = 1;
  unsigned int numSteps = 100000;

  if (argc > 1)
    refiment = std::atoi(argv[1]);

  auto u1 = [](double, double, double){return 1;};
  auto u2 = [](double, double, double){return 0.0;};
  auto f  = [](double, double, double){return 0;};
  auto cExact = [](double, double, double){return 1;};
  // auto f  = [](double, double, double t){return -std::exp(-t);};
  // auto cExact = [](double, double, double t){return std::exp(-t);};
  auto c0 = std::bind(cExact, _1, _2, 0);

  UniqueSquareGrid mesh (refiment /*,0.25*/);
  VTKwriter writer ("analyticalTest", mesh, order);
  auto bcHanderl = [order, orderF, &cExact](UniqueSquareGrid & mesh, double time)
    {
      // setBoundary_Periodic(mesh, Boundary::bottom);
      // setBoundary_Periodic(mesh, Boundary::top);
      // setBoundary_Periodic(mesh, Boundary::left);
      // setBoundary_Periodic(mesh, Boundary::right);

      setBoundary_Diriclet(mesh, Boundary::bottom, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::top, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::left, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::right, order, orderF, cExact, time);
    };

  Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl, tEnd, numSteps, writer, numSteps/100, true, true);

  std::cout << "refiment level: " << refiment
            << ", number of Triangles: " << mesh.getColumns()*mesh.getRows()*2
            << ", basic polynomial degree: " << order
            << ", number of DOFs per triangle: " << numberOf2DBasefunctions(order)
            << ", number of DOFs toal: " << mesh.getColumns()*mesh.getRows()*2*numberOf2DBasefunctions(order)
            << std::endl;

  stepper.go();

  std::cout << "L2 Error: " << stepper.l2error() << std::endl;

  return 0;
}
