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

class Stepper
{
public:
  Stepper(UniqueSquareGrid & mesh, unsigned int order, unsigned int orderF,
          std::function<double(double,double, double)> u1, std::function<double(double,double, double)> u2,
          std::function<double(double,double, double)> f,
          std::function<double(double,double)> c0, std::function<double(double,double, double)> cExact,
          std::function<void(UniqueSquareGrid &)> boundaryHandler,
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
  }

  /**
   * explicit Euler step
   */
  void next()
  {
    std::cout << "STEP " << step_ << ": time=" << t_ << std::endl;

    // update boundary
    boundaryHandler_(mesh_);

    // Assembly matrices and vectors for computation
    assamblyL(mesh_, order_, std::bind(f_, _1, _2, t_));          // RHS vector
    assemblyG(mesh_, hatG_);
    assemblyE(mesh_, hatE_);
    assemblyFr(mesh_, order_, orderF_, riemanSolver_UpWinding, hatI_);

    //TODO set deltaT;
    deltaT_ = 0.01;

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

  }


  double getTime() const {return t_;}
  unsigned int getSteps() const { return step_;}

  double l2Error()
  {
    //TODO
    return -1.;
  }



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
  std::function<void(UniqueSquareGrid &)> boundaryHandler_;
  VTKwriter & writer_;
  unsigned int const writeInterval_;

  double deltaT_; //TODO
  double t_;
  unsigned int step_;
};






#endif
