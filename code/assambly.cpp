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
  auto f  = [](double x, double y)->double{ return 0.; };  // rhs function
  auto c  = [](double x, double y)->double{ return 1.; };

  int const order = 1;     // max polynomial degree
  int const refinment = 2; // order if mesh refinement

  GridOnSquer mesh(refinment);

  Coefficient C  (mesh, order, c);   // c -> gesuchte größe
  Coefficient U1 (mesh, order, u1);  // velocity x component
  Coefficient U2 (mesh, order, u2);  // velocity y component
  Coefficient F  (mesh, order, f);   // rhs vector
  Coefficient F1 = C*U1;             // Flux x component
  Coefficient F2 = C*U2;             // Flux x component

  Matrix hatM = assemblyHatM(order);
  Matrix M = assemblyM(hatM, mesh, order);

  auto hatG = assemblyHatG(order);
  Matrix G = assemblyG(hatG, mesh,{U1,U2}, order);

  auto hatE = assemblyHatE(order);
  auto E = assemblyE(hatE, mesh, order);

  auto hatI = getHatI(order);
  auto Fr = assemblyRiemanFlux_UpWinding(mesh, {U1,U2}, {F1,F2}, hatI, {U1,U2}, order);

  Coefficient L = assemblyL(M,F);

  std::cout << "M hat:\n" << hatM
            << "\v\v\vG1 hat:\n" << hatG[0]
            << "\v\v\vG2 hat:\n" << hatG[1]
            << std::endl;

  std::cout << "\v\v\vL:\n" << L << std::endl;

  std::cout << "\v\v\vM:\n" << M
            << "\v\vG:\n" << G
            << std::endl;

  std::cout << "\v\v\vE1 hat:\n" << hatE[0]
            << "\v\v\vE2 hat:\n" << hatE[1]
            << "\v\v\vE3 hat:\n" << hatE[2]
            << "\v\v\vE1:\n" << E[0]
            << "\v\v\vE2:\n" << E[1]
            << "\v\v\vE3:\n" << E[2]
            << "\v\v\vRiemanFlux1:\n" << Fr[0]
            << "\v\v\vRiemanFlux2:\n" << Fr[1]
            << "\v\v\vRiemanFlux3:\n" << Fr[2]
            <<std::endl;

}
