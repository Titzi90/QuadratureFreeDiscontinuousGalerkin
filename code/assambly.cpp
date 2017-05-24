#include "assembly.hpp"
#include "Matrix.hpp"
#include "squareGrid.hpp"

#include <iostream>
#include <vector>


/**
 * Advection Problem
 * ∂ₜc + ∇*(uc) = f
 */
int main()
{
  /* TODO
  * weitere Matrizen aufatellen
  * Testes für verschiedene Matrizen
  */

  auto u1 = [](double x, double y)->double{ return 1.; };
  auto u2 = [](double x, double y)->double{ return 1.; };
  auto f  = [](double x, double y)->double{ return 1.; };

  int const order = 1;     // max polynomial degree
  int const refinment = 2; // order if mesh refinement

  GridOnSquer mesh(refinment);

  Coefficient U1 (mesh, order, u1);  // velocity x component
  Coefficient U2 (mesh, order, u2);  // velocity y component
  Coefficient F  (mesh, order, f);   // rhs vector

  Matrix hatM = assemblyHatM(order);
  Matrix M = assemblyM(hatM, mesh, order);

  auto hatG = assemblyHatG(order);
  Matrix A = assemblyA(hatG,mesh,{U1,U2}, order);

  Coefficient L = assemblyL(M,F);

  std::cout << "M hat:\n" << hatM
            << "\v\v\vG1 hat:\n" << hatG[0]
            << "\v\v\vG2 hat:\n" << hatG[1]
            << std::endl;

  std::cout << "\v\v\vL:\n" << L << std::endl;

  std::cout << "\v\v\vM:\n" << M
            << "\v\vG:\n" << A
            << std::endl;

}
