
#include <cmath>

#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "stepper.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"


int main()
{
int order = 1;
int orderF = 2;
int refiment = 16;


auto u1 = [](double x, double y, double t){return 0.;};
auto u2 = [](double x, double y, double t){return 0.;};
auto f  = [](double x, double y, double t){return -std::exp(-t);};
auto c0 = [](double x, double y){return 1.;};
auto cExact = [](double x, double y, double t){return std::exp(-t);};

 UniqueSquareGrid mesh (refiment, 0.25);
 VTKwriter writer ("analyticalTest", mesh, order);
auto bcHanderl = [](UniqueSquareGrid & mesh)
  {
    setBoundary_Periodic(mesh, Boundary::bottom);
    setBoundary_Periodic(mesh, Boundary::top);
    setBoundary_Periodic(mesh, Boundary::left);
    setBoundary_Periodic(mesh, Boundary::right);
  };

 Stepper stepper (mesh, order, orderF, u1, u2, f, c0, cExact, bcHanderl, writer, 1, true);


for (int i=0; i<100; ++i)
  stepper.next();

 return 0;
}
