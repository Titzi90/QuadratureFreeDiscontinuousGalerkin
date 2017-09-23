#include "assembly.hpp"
#include "DataTypes.hpp"
#include "Grid.hpp"
#include "VTKwriter.hpp"
#include "monmomials_and_basefunctions.hpp"
#include <iostream>
#include <vector>

std::vector<double> l2Error(UniqueSquareGrid const & mesh,
                            unsigned int polynomialDegree,
                            unsigned int integragradeDegree,
                            std::function<double(double,double)> c_ex,
                            std::function<double(double,double)> u1_ex,
                            std::function<double(double,double)> u2_ex
                            )
{
  double err_c = 0.;
  double err_u1 = 0.;
  double err_u2 = 0.;

// #pragma omp parallel for reduction(+: err)
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle t = mesh.getLower(row,col);
        Jakobian const & B_k (t.getJakobian());
        Point A_k = t.getA();
        Polynomial2D c_aprox = reconstructFunction2D(polynomialDegree, t.C());
        Polynomial2D u1_aprox = reconstructFunction2D(polynomialDegree, t.U1());
        Polynomial2D u2_aprox = reconstructFunction2D(polynomialDegree, t.U2());
        double area_k = t.getArea();

        err_c += 2*area_k * integradeOverRefTriangle_gaus([&c_ex, &c_aprox, &B_k, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c_aprox(x1_hat, x2_hat);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err_u1 += 2*area_k * integradeOverRefTriangle_gaus([&u1_ex, &u1_aprox, &B_k, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c_ex = u1_ex(x1,x2);
                                                          double c_aprox = u1_aprox(x1_hat, x2_hat);

                                                          return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                        },
                                                        integragradeDegree);

        err_u2 += 2*area_k * integradeOverRefTriangle_gaus([&u2_ex, &u2_aprox, &B_k, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c_ex = u2_ex(x1,x2);
                                                          double c_aprox = u2_aprox(x1_hat, x2_hat);

                                                          return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                        },
                                                        integragradeDegree);


        t = mesh.getLower(row,col);
        Jakobian const & B_u = t.getJakobian();
        A_k = t.getA();
        c_aprox = reconstructFunction2D(polynomialDegree, t.C());
        u1_aprox = reconstructFunction2D(polynomialDegree, t.U1());
        u2_aprox = reconstructFunction2D(polynomialDegree, t.U2());
        area_k = t.getArea();

        err_c += 2*area_k * integradeOverRefTriangle_gaus([&c_ex, &c_aprox, &B_u, &A_k](double x1_hat, double x2_hat)
                                                          {
                                                            double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                            double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                            double c = c_ex(x1,x2);
                                                            double c_ = c_aprox(x1_hat, x2_hat);

                                                            return (c-c_)*(c-c_);
                                                          },
                                                          integragradeDegree);
        err_u1 += 2*area_k * integradeOverRefTriangle_gaus([&u1_ex, &u1_aprox, &B_u, &A_k](double x1_hat, double x2_hat)
                                                          {
                                                            double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                            double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                            double c_ex = u1_ex(x1,x2);
                                                            double c_aprox = u1_aprox(x1_hat, x2_hat);

                                                            return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                          },
                                                          integragradeDegree);
        err_u2 += 2*area_k * integradeOverRefTriangle_gaus([&u2_ex, &u2_aprox, &B_u, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                          double c_ex = u2_ex(x1,x2);
                                                          double c_aprox = u2_aprox(x1_hat, x2_hat);

                                                          return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                        },
                                                        integragradeDegree);
      }

  return {std::sqrt(err_c), std::sqrt(err_u1), std::sqrt(err_u2)};
}



std::vector<double> l2ErrorEdge(UniqueSquareGrid const & mesh,
                                unsigned int polynomialDegree,
                                unsigned int integragradeDegree,
                                std::function<double(double,double)> c_ex,
                                std::function<double(double,double)> u1_ex,
                                std::function<double(double,double)> u2_ex
                                )
{
  double err1_c = 0.;
  double err2_c = 0.;
  double err3_C = 0.;
  double err1_f = 0.;
  double err2_f = 0.;
  double err3_f = 0.;

  std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);

// #pragma omp parallel for reduction(+: err)
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle t = mesh.getLower(row,col);
        Jakobian const & B_k (t.getJakobian());
        Point A_k = t.getA();
        auto norm1 = t.getNormalA();
        auto norm2 = t.getNormalB();
        auto norm3 = t.getNormalC();

        auto f_ex1 = [&c_ex, &u1_ex, &u2_ex, &norm1](double x, double y)
          {
            return (norm1.x*u1_ex(x,y) + norm1.y*u2_ex(x,y)) * c_ex(x,y);
          };
        auto f_ex2 = [&c_ex, &u1_ex, &u2_ex, &norm2](double x, double y)
          {
            return (norm2.x*u1_ex(x,y) + norm2.y*u2_ex(x,y)) * c_ex(x,y);
          };
        auto f_ex3 = [&c_ex, &u1_ex, &u2_ex, &norm3](double x, double y)
          {
            return (norm3.x*u1_ex(x,y) + norm3.y*u2_ex(x,y)) * c_ex(x,y);
          };



        Polynomial1D c1_aprox = reconstructFunction1D(polynomialDegree, T[0]*t.C());
        Polynomial1D c2_aprox = reconstructFunction1D(polynomialDegree, T[1]*t.C());
        Polynomial1D c3_aprox = reconstructFunction1D(polynomialDegree, T[2]*t.C());
        /*
        Polynomial1D f1_aprox = reconstructFunction1D(polynomialDegree,
                                                      l2Projection_edge(polynomialDegree,f_ex,B_k,A_k,0));
        Polynomial1D f2_aprox = reconstructFunction1D(polynomialDegree,
                                                      l2Projection_edge(polynomialDegree,f_ex,B_k,A_k,1));
        Polynomial1D f3_aprox = reconstructFunction1D(polynomialDegree,
                                                      l2Projection_edge(polynomialDegree,f_ex,B_k,A_k,2));
        */
        Polynomial1D f1_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_a());
        Polynomial1D f2_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_b());
        Polynomial1D f3_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_c());

        double l1 = t.getLengthA();
        double l2 = t.getLengthB();
        double l3 = t.getLengthC();

        err1_c += l1 * integradeOverRefEdge_gaus([&c_ex, &c1_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 1-x_bar;
                                                          double x2_hat = x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c1_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err2_c += l2 * integradeOverRefEdge_gaus([&c_ex, &c2_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 0;
                                                          double x2_hat = 1-x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c2_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err3_C += l3 * integradeOverRefEdge_gaus([&c_ex, &c3_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = x_bar;
                                                          double x2_hat = 0;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c3_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err1_f += l1 * integradeOverRefEdge_gaus([&f_ex1, &f1_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 1-x_bar;
                                                          double x2_hat = x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex1(x1,x2);
                                                          double c_ = f1_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);

        err2_f += l2 * integradeOverRefEdge_gaus([&f_ex2, &f2_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 0;
                                                          double x2_hat = 1-x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex2(x1,x2);
                                                          double c_ = f2_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);

        err3_f += l3 * integradeOverRefEdge_gaus([&f_ex3, &f3_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = x_bar;
                                                          double x2_hat = 0;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex3(x1,x2);
                                                          double c_ = f3_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);


        t = mesh.getLower(row,col);
        Jakobian const & B_u = t.getJakobian();
        A_k = t.getA();
        c1_aprox = reconstructFunction1D(polynomialDegree, T[0]*t.C());
        c2_aprox = reconstructFunction1D(polynomialDegree, T[1]*t.C());
        c3_aprox = reconstructFunction1D(polynomialDegree, T[2]*t.C());
        f1_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_a());
        f2_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_b());
        f3_aprox = reconstructFunction1D(polynomialDegree*2, t.Fn_c());
        l1 = t.getLengthA();
        l2 = t.getLengthB();
        l3 = t.getLengthC();

        err1_c += l1 * integradeOverRefEdge_gaus([&c_ex, &c1_aprox, &B_u, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 1-x_bar;
                                                          double x2_hat = x_bar;
                                                          double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c1_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err2_c += l2 * integradeOverRefEdge_gaus([&c_ex, &c2_aprox, &B_u, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 0;
                                                          double x2_hat = 1-x_bar;
                                                          double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c2_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err3_C += l3 * integradeOverRefEdge_gaus([&c_ex, &c3_aprox, &B_u, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = x_bar;
                                                          double x2_hat = 0;
                                                          double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                          double c = c_ex(x1,x2);
                                                          double c_ = c3_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree);

        err1_f += l1 * integradeOverRefEdge_gaus([&f_ex1, &f1_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 1-x_bar;
                                                          double x2_hat = x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex1(x1,x2);
                                                          double c_ = f1_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);

        err2_f += l2 * integradeOverRefEdge_gaus([&f_ex2, &f2_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = 0;
                                                          double x2_hat = 1-x_bar;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex2(x1,x2);
                                                          double c_ = f2_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);

        err3_f += l3 * integradeOverRefEdge_gaus([&f_ex3, &f3_aprox, &B_k, &A_k](double x_bar)
                                                        {
                                                          double x1_hat = x_bar;
                                                          double x2_hat = 0;
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c = f_ex3(x1,x2);
                                                          double c_ = f3_aprox(x_bar);

                                                          return (c-c_)*(c-c_);
                                                        },
                                                        integragradeDegree*2);

      }


  std::cout << "c error on edges A: " << std::sqrt(err1_c)
            <<" B: " << std::sqrt(err2_c)
            <<" C: " << std::sqrt(err3_C)
            << std::endl;
  std::cout << "f error on edges A: " << std::sqrt(err1_f)
            <<" B: " << std::sqrt(err2_f)
            <<" C: " << std::sqrt(err3_f)
            << std::endl;

  return {std::sqrt(err1_c+err2_c+err3_C), std::sqrt(err1_f+err2_f+err3_f)};
}

/**
 * Advection Problem
 * ∂ₜc + ∇*(uc) = f
 */
int main(int argc, char** argv)
{
  int order = 0;
  int refiment = 64;

  if (argc > 2)
    refiment = std::atoi(argv[2]);
  if (argc > 3)
    order = std::atoi(argv[3]);

  // constant
  // auto f = [](double, double){return 1.;};

  // linear
  // auto f = [](double x, double){return x;};

  // quadratic
  // auto f = [](double x, double){return (2.*x-1.)*(2.*x-1.);};

  // more complex
  auto f = [](double x, double y) {return std::cos(7.*x)*std::cos(7.*y);};
  auto u1 = [](double x, double y){return std::exp((x+y)*0.5);};
  auto u2 = [](double x, double y){return std::exp((x-y)*0.5);};


  UniqueSquareGrid mesh(refiment);
  assamblyC(mesh, order, f);
  assamblyU(mesh, order, u1, u2);
  assamblyF(mesh, order, order*2, assamblyLocalLinearF);

  std::vector<double> err_h= l2Error(mesh, order,order+1, f, u1, u2);
  std::vector<double> err_e = l2ErrorEdge(mesh, order, order+1, f,u1,u2);





  std::cout << "Error analysis of mappings:\n"
            << "c in element: " << err_h[0] << "\n"
            << "u1 in element: " << err_h[1] << "\n"
            << "u2 in element: " << err_h[2] << "\n"
            << "c on edge   :  " << err_e[0] << "\n"
            << "f:             " << err_e[1]
            << std::endl;


  /*

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


    */
















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
