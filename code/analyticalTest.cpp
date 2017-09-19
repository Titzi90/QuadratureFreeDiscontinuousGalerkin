
#include <cmath>
#include <cstdlib>

#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "stepper.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"



int main(int argc, char** argv)
{
  int order = 3;
int refiment = 64;
  double tEnd = 1;
  unsigned int numSteps = tEnd * 1000;

  if (argc > 1)
    numSteps = std::atoi(argv[1]);
  if (argc > 2)
    refiment = std::atoi(argv[2]);
  if (argc > 3)
    order = std::atoi(argv[3]);

  int orderF = 2*order;

  auto u1 = [](double, double, double){return 1;};
  auto u2 = [](double, double, double){return 0.0;};
  // constant
  // auto f  = [](double, double, double){return 0;};
  // auto cExact = [](double, double, double){return 1;};

  // linear
  // auto f  = [](double, double, double){return 1;};
  // auto cExact = [](double x, double, double){return x;};

  // quadratic
  auto f  = [](double x, double, double){return 8.*x-4.;};
  auto cExact = [](double x, double, double){return (2.*x-1.)*(2.*x-1.)*100.;};

  // time dependent
  // auto f  = [](double, double, double t){return -std::exp(-t);};
  // auto cExact = [](double, double, double t){return std::exp(-t);};

  // balken
  // auto f  = [](double, double, double){return 0;};
  // auto cExact = [](double , double, double){return 1.;};

  auto c0 = std::bind(cExact, _1, _2, 0);

  UniqueSquareGrid mesh (refiment /*,0.25*/);
  VTKwriter writer ("analyticalTest_time", mesh, order);
  auto bcHanderl = [order, orderF, &cExact](UniqueSquareGrid & mesh, double time)
    {
      // setBoundary_Periodic(mesh, Boundary::bottom);
      // setBoundary_Periodic(mesh, Boundary::top);
      // setBoundary_Periodic(mesh, Boundary::left);
      // setBoundary_Periodic(mesh, Boundary::right);

      setBoundary_Diriclet(mesh, Boundary::bottom, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::top, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::right, order, orderF, cExact, time);
      setBoundary_Diriclet(mesh, Boundary::left, order, orderF, cExact, time);

      // setBoundary_Diriclet(mesh, Boundary::left, order, orderF, [](double , double y, double){
      //     if (y>0.2 && y<0.8) return 2.; else return 1.;}, time);
    };

  Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl,
                   tEnd, numSteps, writer, false, numSteps/100, true ); //, numSteps/100, true);

  // stepper.go();
  // for (int i=0; i<1; ++i)
    stepper.next();

  std::cout << "L2 Error: " << stepper.l2error() << std::endl;

  return 0;
}
