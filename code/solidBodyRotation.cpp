#include <cmath>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>

#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "stepper.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"



int main(int argc, char** argv)
{
  int order = 0;
  int fMultiplier = 2;
  int refiment = 64;
  double tEnd = 2*M_PI;
  unsigned int numSteps = tEnd * 1000;
  std::string outputName ("analyticalConvergenceTest");
  double disorder= 0.;

  if (argc > 1)
    numSteps = std::atoi(argv[1]) * tEnd;
  if (argc > 2)
    refiment = std::atoi(argv[2]);
  if (argc > 3)
    order = std::atoi(argv[3]);
  if (argc > 4)
    fMultiplier = std::atoi(argv[4]);
  if (argc > 5)
    outputName = argv[5];
  if (argc > 6)
    disorder = std::atof(argv[6]);

  int orderF = fMultiplier*order;

  auto u1 = [](double, double y, double){return 0.5 - y;};
  auto u2 = [](double x, double, double){return x - 0.5;};
  auto f  = [](double, double, double ){return 0;};
  double r = 0.0225;
  auto G = [](double x, double y, double x0, double y0)
    {return 1./0.15 * std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));};
  auto c0 = [r, &G](double x, double y)
    {
      // slotted cylinder
      if( (x-0.5)*(x-0.5) + (y-0.75)*(y-0.75) <= r
          && (x<=0.475 || x>=0.525 || y>=0.85) )
        return 1.;
      // sharp cone
      else if( (x-0.5)*(x-0.5) + (y-0.25)*(y-0.25) <= r )
        return 1. - G(x, y, 0.5, 0.25);
      // smooth hump
      else if( (x-0.25)*(x-0.25) + (y-0.5)*(y-0.5) <= r )
        return 0.25*(1 + std::cos(M_PI*G(x,y,0.25,0.5)));
      else
        return 0.;
    };
  auto cExact = [&c0](double x, double y, double){return c0(x,y);};

  UniqueSquareGrid mesh (refiment , disorder);
  VTKwriter writer (outputName, mesh, order);
  auto bcHanderl = [order, orderF](UniqueSquareGrid & mesh, double time)
    {
      setBoundary_Diriclet(mesh, Boundary::bottom, order, orderF, [](double,double,double){return 0.;}, time);
      setBoundary_Diriclet(mesh, Boundary::top, order, orderF, [](double,double,double){return 0.;}, time);
      setBoundary_Diriclet(mesh, Boundary::left, order, orderF, [](double,double,double){return 0.;}, time);
      setBoundary_Diriclet(mesh, Boundary::right, order, orderF, [](double,double,double){return 0.;}, time);
    };

  Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl,
                   tEnd, numSteps, writer, false, numSteps/100, true);

  stepper.go();

  std::cout << "L2 Error: " << stepper.l2error() << std::endl;

  return 0;
}
