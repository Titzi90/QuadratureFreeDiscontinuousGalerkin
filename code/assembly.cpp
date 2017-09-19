#include "assembly.hpp"
#include "DataTypes.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"
#include "monmomials_and_basefunctions.hpp"
#include <iostream>
#include <vector>


/**
 * Advection Problem
 * ∂ₜc + ∇*(uc) = f
 */
int main()
{
  UniqueSquareGrid mesh(3);

  auto t=mesh.getUpper(2, 2);
  auto const & B  ( t.getJakobian() );
  auto const A = t.getA();

  auto const f  = [](double x, double y)->double{ return (2.*x-1.)*(2.*x-1.)*100.; };

  std::cout << l2Projection(3,f,B,A);






  std::vector<double> w,x1,x2;
    x1 = {
      0.054830900955589,   0.175654279195255,   0.343651813106453,   0.533230731173959,
      0.715527432866568,   0.862793031223432,   0.952646581185227,   0.048991501878362,
      0.156947392786903,   0.307053470832875,   0.476442551784230,   0.639324960202548,
      0.770907019092335,   0.851191316541618,   0.039548223967455,   0.126695251279609,
      0.247867874404688,   0.384606636317686,   0.516092908865112,   0.622312080263295,
      0.687121307473297,   0.028131280268461,   0.090120345868446,   0.176312358556585,
      0.273576813165278,   0.367105088607705,   0.442660473419548,   0.488760306780644,
      0.016714336569468,   0.053545440457283,   0.104756842708482,   0.162546990012870,
      0.218117268350298,   0.263008866575801,   0.290399306087990,   0.007271058658560,
      0.023293298949990,   0.045571246280295,   0.070711074546325,   0.094885217012863,
      0.114413927746761,   0.126329297019669,   0.001431659581333,   0.004586412541638,
      0.008972904006717,   0.013922895156596,   0.018682744348843,   0.022527915615664,
      0.024874032376061};
    x2 = {
      0.001431659581333,   0.004586412541638,   0.008972904006717,   0.013922895156596,   0.018682744348843,
      0.022527915615664,   0.024874032376061,   0.007271058658560,   0.023293298949990,   0.045571246280295,
      0.070711074546325,   0.094885217012863,   0.114413927746761,   0.126329297019669,   0.016714336569468,
      0.053545440457283,   0.104756842708482,   0.162546990012870,   0.218117268350298,   0.263008866575801,
      0.290399306087990,   0.028131280268461,   0.090120345868446,   0.176312358556585,   0.273576813165278,
      0.367105088607705,   0.442660473419548,   0.488760306780644,   0.039548223967455,   0.126695251279609,
      0.247867874404688,   0.384606636317686,   0.516092908865112,   0.622312080263295,   0.687121307473297,
      0.048991501878362,   0.156947392786903,   0.307053470832875,   0.476442551784230,   0.639324960202548,
      0.770907019092335,   0.851191316541618,   0.054830900955589,   0.175654279195255,   0.343651813106453,
      0.533230731173959,   0.715527432866568,   0.862793031223432,   0.952646581185227    };

  // for (int i=0; i<10; ++i)
  //   {
  //     std::cout << i << ":\n";
    // for (int j=0; j<x1.size(); ++j)
    //   {
    //     double x1_ = B[0][0]*x1[j] + B[0][1]*x2[j] + A.x;
    //     double x2_ = B[1][0]*x1[j] + B[1][1]*x2[j] + A.y;
    //     std::cout << f(x1_, x2_) << std::endl;
    //   }
        // std::cout << pol::phi[i](x1[j],x2[j]) << "\n";
        // std::cout << 1[1][0]*x1[j] + B[1][1]*x2[j] +A.y << "\n";
      // std::cout << "\n\n\n\n" << std::endl;
    // }



















  /*
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
  */
}
