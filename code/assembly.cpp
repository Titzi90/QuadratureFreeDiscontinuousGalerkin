#include "assembly.hpp"
#include "DataTypes.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"

#include <iostream>
#include <vector>


/**
 * Advection Problem
 * ∂ₜc + ∇*(uc) = f
 */
int main()
{
  int const order = 1;     // max polynomial degree
  int const refinment = 2; // order if mesh refinement

  UniqueSquareGrid mesh(refinment);

  auto u1 = [](double x, double y)->double{ return 1.; };
  auto u2 = [](double x, double y)->double{ return 1.; };
  auto f  = [](double x, double y)->double{ return 0.; };
  auto c  = [](double x, double y)->double{ return 1.; };

  auto hatM (assemblyHatM(order));
  auto hatG (assemblyHatG(order));
  auto hatE (assemblyHatE(order, 2*order));
  auto hatI (getHatI(2*order));

  assamblyC(mesh, order, c);          // gesuchte Größe c
  assamblyU(mesh, order, u1, u2);     // velocity field U
  assamblyF(mesh, order, 2*order, assamblyLocalLinearF);  // Flux field f=cu
  setBoundary_Periodic(mesh, Boundary::bottom);
  setBoundary_Periodic(mesh, Boundary::top);
  setBoundary_Periodic(mesh, Boundary::left);
  setBoundary_Periodic(mesh, Boundary::right);

  // assemblyM(mesh, hatM);
  // assemblyG(mesh, hatG);
  assemblyMGaus(mesh, order);
  assemblyGgaus(mesh, order);
  assemblyE(mesh, hatE);
  assemblyFr(mesh, order, order*2, riemanSolver_UpWinding, hatI);
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
  printFr(mesh, std::cout); //TODO wiso anderes als bei test?!
  std::cout << "\v\v\v";
  printL(mesh, std::cout);
  std::cout << std::endl;



  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        auto & T_l = mesh.getLower(row, col);
        auto & T_u = mesh.getUpper(row, col);

        //TODO durch M teilen
        auto dc_l = T_l.M()*T_l.L() + T_l.G()*T_l.C()
          - T_l.E_a()*T_l.F_a() - T_l.E_b()*T_l.F_b() - T_l.E_c()*T_l.F_c();
        auto dc_u = T_u.M()*T_u.L() + T_u.G()*T_u.C()
          - T_u.E_a()*T_u.F_a() - T_u.E_b()*T_u.F_b() - T_u.E_c()*T_u.F_c();

        std::cout << dc_l << std::endl;
        std::cout << dc_u << std::endl;
      }

  VTKwriter writer ("test", mesh, order);
  writer.write();

}
