
#include <cmath>
#include <functional>

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
  unsigned int numSteps = 5000;

  if (argc > 1)
    refiment = std::atoi(argv[1]);

  using std::placeholders::_1;
  using std::placeholders::_2;

  auto u1 = [](double x, double y, double){return std::exp((x+y)*0.5);};
  auto u2 = [](double x, double y, double){return std::exp((x-y)*0.5);};
  auto cExact = [](double x, double y, double){return std::cos(7.*x)*std::cos(7.*y);};
  auto c0 = std::bind(cExact, _1, _2, 0);
  auto f  = [&u1, &u2, &cExact](double x, double y, double t)
    {
      return 0.5*cExact(x,y,t)*u1(x,y,t) -7.*std::sin(7.*x)*std::cos(7.*y)*u1(x,y,t)
            -0.5*cExact(x,y,t)*u2(x,y,t) -7.*std::cos(7.*x)*std::sin(7.*y)*u2(x,y,t);
    };

  UniqueSquareGrid mesh (refiment);
  VTKwriter writer ("analyticalConvergenceTest", mesh, order);
  auto bcHanderl = [&cExact, order, orderF](UniqueSquareGrid & mesh, double t)
    {
      setBoundary_Diriclet(mesh, Boundary::bottom, order, orderF, cExact, t);
      setBoundary_Diriclet(mesh, Boundary::top, order, orderF, cExact, t);
      setBoundary_Diriclet(mesh, Boundary::left, order, orderF, cExact, t);
      setBoundary_Diriclet(mesh, Boundary::right, order, orderF, cExact, t);
    };

  Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl, tEnd, numSteps, writer, 1, true, true);

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