#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "Grid.hpp"
#include "DataTypes.hpp"
#include "assembly.hpp"
#include "VTKwriter.hpp"

#include <functional>
#include <vector>
#include <iostream>

using namespace std::placeholders;


/*** forward declaration ***********/
double l2Error(UniqueSquareGrid const &, unsigned int, unsigned int, std::function<double(double,double)> f_ex);
class Stepper;


/*** class implementation ***********/

class Stepper
{
public:
  Stepper(UniqueSquareGrid & mesh, unsigned int order, unsigned int orderF,
          std::function<double(double,double, double)> u1, std::function<double(double,double, double)> u2,
          std::function<double(double,double, double)> f,
          std::function<double(double,double)> c0, std::function<double(double,double, double)> cExact,
          std::function<void(UniqueSquareGrid &, double)> boundaryHandler,
          VTKwriter & writer, unsigned int writeInterval,
          bool writeInitialData)
    :mesh_(mesh), order_(order), orderF_(orderF), u1_(u1), u2_(u2), f_(f), c0_(c0), cExact_(cExact),
     hatM_(assemblyHatM(order_)), hatG_(assemblyHatG(order_)), hatE_(assemblyHatE(order_, orderF_)),
     hatI_(getHatI(orderF_)), boundaryHandler_(boundaryHandler),
     writer_(writer), writeInterval_(writeInterval),
     t_(0.), step_(0)
  {
    // set initial data
    assamblyC(mesh_, order_, c0_);

    // assembly constant data: mass matrix and rhs
    assemblyM(mesh_, hatM_);            // mass matrix:W

    // set initial U and F
    assamblyU(mesh_, order_, std::bind(u1_, _1, _2, t_), std::bind(u2_, _1, _2, t_));
    assamblyF(mesh_, order_, orderF_, assamblyLocalLinearF);

    // write initial data
    if (writeInitialData)
      writer_.write();

    //error
    std::cout << "initial L2 error: " << this->l2error() << std::endl;
  }

  /**
   * explicit Euler step
   */
  void next()
  {
    std::cout << "STEP " << step_ << ": time=" << t_ << std::endl;

    // update boundary
    boundaryHandler_(mesh_, t_);

    // Assembly matrices and vectors for computation
    assamblyL(mesh_, order_, std::bind(f_, _1, _2, t_));          // RHS vector
    assemblyG(mesh_, hatG_);
    assemblyE(mesh_, hatE_);
    assemblyFr(mesh_, order_, orderF_, riemanSolver_UpWinding, hatI_);

    //TODO set deltaT;
    deltaT_ = 0.001;

    //increase time
    ++step_;
    t_ += deltaT_;

    // update c
    for (unsigned int row=0; row<mesh_.getRows(); ++row)
      for (unsigned int col=0; col<mesh_.getColumns(); ++col)
        {
          auto & T_l = mesh_.getLower(row, col);
          auto & T_u = mesh_.getUpper(row, col);

          std::vector<double> tmp_l = T_l.G()*T_l.C()
            - T_l.E_a()*T_l.F_a() - T_l.E_b()*T_l.F_b() - T_l.E_c()*T_l.F_c();
          std::vector<double> tmp_u = T_u.G()*T_u.C()
            - T_u.E_a()*T_u.F_a() - T_u.E_b()*T_u.F_b() - T_u.E_c()*T_u.F_c();


          // std::cout << row  << "," << col << " lower:\n"
                    // << "C_old:\n" << T_l.C()
                    // << "\vG:\n" << T_l.G()
                    // << "\vU:\n" << T_l.U1() << "\n" << T_l.U2()
                    // << "\vE:\n" << T_l.E_a() << "\n" << T_l.E_b() << "\n" << T_l.E_c()
                    // << "\vFr:\n" << T_l.F_a() << "\n" << T_l.F_b() << "\n" << T_l.F_c()
                    // << "\vF:\n" << T_l.F1() << "\n" << T_l.F2()
                    // << "\vL:\n" << T_l.L()
                    // << "\vdc:\n" << tmp_l
                    // << "\vM^-1:\n" << invertM(T_l.M())
                    // << std::endl;

          // std::cout << row  << "," << col << " upper:\n"
                    // << "C_old:\n" << T_u.C()
                    // << "\vG:\n" << T_u.G()
                    // << "\vU:\n" << T_u.U1() << "\n" << T_u.U2()
                    // << "\vE:\n" << T_u.E_a() << "\n" << T_u.E_b() << "\n" << T_u.E_c()
                    // << "\vFr:\n" << T_u.F_a() << "\n" << T_u.F_b() << "\n" << T_u.F_c()
                    // << "\vF:\n" << T_u.F1() << "\n" << T_u.F2()
                    // << "\vL:\n" << T_u.L()
                    // << "\vdc:\n" << tmp_u
                    // << "\vM^-1:\n" << invertM(T_u.M())
                    // << std::endl;





          //TODO umstellen um M-1*M bei L zu sparen
          T_l.C() += deltaT_ * ( T_l.L() + invertM(T_l.M()) * tmp_l);
          T_u.C() += deltaT_ * ( T_u.L() + invertM(T_u.M()) * tmp_u);
          // T_l.C() += deltaT_ * ( T_l.L() /*+ invertM(T_l.M()) * tmp_l*/);
          // T_u.C() += deltaT_ * ( T_u.L() /*+ invertM(T_u.M()) * tmp_u*/);

          // std::cout << "C_neu:\n" << T_l.C() << "\n" << T_u.C() << std::endl;
        }

    // update U and F
    assamblyU(mesh_, order_, std::bind(u1_, _1, _2, t_), std::bind(u2_, _1, _2, t_));
    assamblyF(mesh_, order_, orderF_, assamblyLocalLinearF);

    // write data
    if (step_ % writeInterval_ == 0)
      writer_.write();

    // error
    std::cout << "L2 error: " << this->l2error() << std::endl;;

  }

  /**
   * get L2 error of current time step
   */
  double l2error()
  {
    return l2Error(mesh_, order_, order_*2, std::bind(cExact_, _1, _2, t_));
  }


  double getTime() const {return t_;}
  unsigned int getSteps() const { return step_;}


private:
  UniqueSquareGrid & mesh_;
  unsigned int const order_, orderF_;
  std::function<double(double, double, double)> u1_, u2_, f_;
  std::function<double(double,double)> c0_;
  std::function<double(double,double,double)> cExact_;
  BlockMatrix hatM_;
  std::vector<Tensor> hatG_;
  std::vector<BlockMatrix> hatE_;
  std::vector<double> hatI_;
  std::function<void(UniqueSquareGrid &, double)> boundaryHandler_;
  VTKwriter & writer_;
  unsigned int const writeInterval_;

  double deltaT_; //TODO
  double t_;
  unsigned int step_;
};


/*** free function implementation ***********/

double l2Error(UniqueSquareGrid const & mesh,
               unsigned int polynomialDegree,
               unsigned int integragradeDegree,
               std::function<double(double,double)> f_ex)
{
  double err = 0;

  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle t = mesh.getLower(row,col);
        Jakobian const & B_k (t.getJakobian());
        Point A_k = t.getA();
        Polynomial2D f_aprox = reconstructFunction2D(polynomialDegree, t.C());
        double area_k = t.getArea();

        err += 2*area_k * integradeOverRefTriangle_gaus([&f_ex, &f_aprox, &B_k, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                          double c_ex = f_ex(x1,x2);
                                                          double c_aprox = f_aprox(x1_hat, x2_hat);

                                                          return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                        },
                                                        integragradeDegree);

        t = mesh.getLower(row,col);
        Jakobian const & B_u = t.getJakobian();
        A_k = t.getA();
        f_aprox = reconstructFunction2D(polynomialDegree, t.C());
        area_k = t.getArea();

        err += 2*area_k * integradeOverRefTriangle_gaus([&f_ex, &f_aprox, &B_u, &A_k](double x1_hat, double x2_hat)
                                                        {
                                                          double x1 = B_u[0][0]*x1_hat + B_u[0][1]*x2_hat + A_k.x;
                                                          double x2 = B_u[1][0]*x2_hat + B_u[1][1]*x2_hat + A_k.y;
                                                          double c_ex = f_ex(x1,x2);
                                                          double c_aprox = f_aprox(x1_hat, x2_hat);

                                                          return (c_ex-c_aprox)*(c_ex-c_aprox);
                                                        },
                                                        integragradeDegree);

      }

  return std::sqrt(err);
}




#endif
