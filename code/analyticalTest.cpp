
#include <cmath>

#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "stepper.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"


int main()
{
int order = 1;
int orderF = 2;
int refiment = 64;


auto u1 = [](double, double, double){return 0.;};
auto u2 = [](double, double, double){return 0.1;};
auto f  = [](double, double, double t){return -std::exp(-t);};
auto c0 = [](double, double){return 1.;};
auto cExact = [](double, double, double t){return std::exp(-t);};

 UniqueSquareGrid mesh (refiment /*,0.25*/);
 VTKwriter writer ("analyticalTest", mesh, order);
 auto bcHanderl = [](UniqueSquareGrid & mesh, double)
  {
    setBoundary_Periodic(mesh, Boundary::bottom);
    setBoundary_Periodic(mesh, Boundary::top);
    setBoundary_Periodic(mesh, Boundary::left);
    setBoundary_Periodic(mesh, Boundary::right);
  };

 Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl, writer, 1, true);


for (int i=0; i<1000; ++i)
  stepper.next();

 return 0;
}
