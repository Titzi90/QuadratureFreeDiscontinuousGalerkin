#include "assembly.hpp"
#include "DataTypes.hpp"
#include "Grid.hpp"

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

  UniqueSquareGrid mesh(refinment);

  assamblyC(mesh, order, c);          // gesuchte Größe c
  assamblyU(mesh, order, u1, u2);     // velocity field U
  assamblyF(mesh, order, 2*order, assamblyLocalLinearF);  // Flux field f=cu

  auto hatM (assemblyHatM(order));
  auto hatG (assemblyHatG(order));
  auto hatE (assemblyHatE(order));
  auto hatT (getLinearTrasformationToRefEdge(order));
  auto hatI (getHatI(order));

  assemblyM(mesh, hatM);
  assemblyG(mesh, hatG);
  assemblyE(mesh, hatE);
  assemblyFr(mesh, order, riemanSolver_UpWinding, hatT, hatI);
  assamblyL(mesh, order, f);          // RHS vector



  printC(mesh, std::cout);
  std::cout << "\v\v\v";
  printU(mesh, std::cout);
  std::cout << "\v\v\v";
  printF(mesh, std::cout);
  std::cout << std::endl;

  std::cout << "hatM:\n" << hatM
            << "\v\v\vhatG1:\n" << hatG[0]
            << "\v\v\vhatG2:\n" << hatG[1]
            << "\v\v\vhatEa:\n" << hatE[0]
            << "\v\v\vhatEb:\n" << hatE[1]
            << "\v\v\vhatEc:\n" << hatE[2]
            << "\v\v\vhatI:\n"  << hatI
            << std::endl;

  std::cout << "\v\v\v";
  printM(mesh, std::cout);
  std::cout << "\v\v\v";
  printG(mesh, std::cout);
  std::cout << "\v\v\v";
  printE(mesh, std::cout);
  std::cout << "\v\v\v";
  printFr(mesh, std::cout);
  std::cout << "\v\v\v";
  printL(mesh, std::cout);
  std::cout << std::endl;

}
