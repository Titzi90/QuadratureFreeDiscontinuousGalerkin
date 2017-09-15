#ifndef STEPPER_HPP
#define STEPPER_HPP

#include "Grid.hpp"
#include "DataTypes.hpp"
#include "assembly.hpp"
#include "VTKwriter.hpp"

#include <functional>
#include <vector>
#include <iostream>
#include <limits>
#include <omp.h>

// include likwid
extern "C"
{
#ifdef USE_LIKWID
#include <likwid.h>
  // This block enables compilation of the code with and without LIKWID in place
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif
}

#ifndef _GITVERSION
#define _GITVERSION "000000"
#endif

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
          double tEnd, unsigned int numSteps,
          VTKwriter & writer, bool writeError=false,
          unsigned int writeInterval=std::numeric_limits<unsigned int>::max(),
          bool writeInitialData=false)
    :mesh_(mesh), order_(order), orderF_(orderF), u1_(u1), u2_(u2), f_(f), c0_(c0), cExact_(cExact),
     hatM_(assemblyHatM(order_)), hatG_(assemblyHatG(order_)), hatE_(assemblyHatE(order_, orderF_)),
     hatI_(getHatI(orderF_)), boundaryHandler_(boundaryHandler),
     numSteps_(numSteps), deltaT_(tEnd/numSteps_),
     writer_(writer), writeInterval_(writeInterval), writeError_(writeError),
     t_(0.), step_(0)
  {

    // initialize likwid
    LIKWID_MARKER_INIT;

#pragma omp parallel
    {
      LIKWID_MARKER_THREADINIT;
    } // Implicit barrier after thread init

#pragma omp parallel
    {
      // set initial data
      assamblyC(mesh_, order_, c0_);

      // assembly constant data: mass matrix and rhs
      assemblyM(mesh_, hatM_);            // mass matrix:W

      // set initial U and F
      assamblyU(mesh_, order_, std::bind(u1_, _1, _2, t_), std::bind(u2_, _1, _2, t_));
      assamblyF(mesh_, order_, orderF_, assamblyLocalLinearF);

      // write initial data
#pragma omp master
      {
        if (writeInitialData)
          writer_.write();

        // get some OMP information
        int chunkSize;
        omp_sched_t schedType;
        std::string schedName;
        omp_get_schedule(&schedType, &chunkSize);
        switch (schedType)
          {
          case omp_sched_static  : schedName = "static"; break;
          case omp_sched_dynamic : schedName = "dynamic"; break;
          case omp_sched_guided  : schedName = "guided"; break;
          case omp_sched_auto    : schedName = "auto"; break;
          }

        std::cout << "Computing with basic polynomial order " << order_
                  << " (" << numberOf2DBasefunctions(order_) <<" local DOFs) "
                  << "and numerical flux aproximation of order "  << orderF_ << " on "
                  << mesh_.getRows()*mesh_.getColumns()*2 << " triangles ("
               // << mesh.getColumns()*mesh.getRows()*2*numberOf2DBasefunctions(order) << " DOFs total)\n"
                  << "refiment level: " << mesh_.getColumns() <<")\n"
                  << "Data written to " << writer_.getName() << ".vtk\n"
                  << "Starting time integration from 0 to " << tEnd
                  << " using time step size " << deltaT_
                  << " (" << numSteps_ << ")\n"
                  << omp_get_num_threads() << " OpenMP treads are used (chunk size " << chunkSize << ") "
                  << "with " << schedName << " scheduler.\n"
                  << "This is git version " << _GITVERSION
                  << std::endl;

        //error
        std::cout << "initial L2 error: " << this->l2error() << std::endl;
      }
    }
  }

  // close likwid
  ~Stepper()
  {
    std::cerr << "END" << std::endl;
    LIKWID_MARKER_CLOSE;
  }

  /**
   * explicit Euler step
   */
  void next()
  {
    //increase time
    ++step_;
    t_ += deltaT_;

#pragma omp parallel
    {
      // update boundary
      boundaryHandler_(mesh_, t_);

      // Assembly matrices and vectors for computation
      assamblyL(mesh_, order_, std::bind(f_, _1, _2, t_));          // RHS vector
      assemblyG(mesh_, hatG_);
      // assemblyGquadFree(mesh_, order_);
      // assemblyGgaus(mesh_, order_);
      assemblyE(mesh_, hatE_);
      // assemblyEquadFree(mesh_, order_, orderF_);
      assemblyFr(mesh_, order_, orderF_, riemanSolver_UpWinding, hatI_);

      // update c
#pragma omp for
      for (unsigned int row=0; row<mesh_.getRows(); ++row)
        for (unsigned int col=0; col<mesh_.getColumns(); ++col)
          {
            auto & T_l = mesh_.getLower(row, col);
            auto & T_u = mesh_.getUpper(row, col);

            std::vector<double> tmp_l = T_l.G()*T_l.C()
              - T_l.E_a()*T_l.Fr_a() - T_l.E_b()*T_l.Fr_b() - T_l.E_c()*T_l.Fr_c();
            std::vector<double> tmp_u = T_u.G()*T_u.C()
              - T_u.E_a()*T_u.Fr_a() - T_u.E_b()*T_u.Fr_b() - T_u.E_c()*T_u.Fr_c();


            // std::cout << row  << "," << col << " lower:\n"
            //           << "C_old:\n" << T_l.C()
            //           << "\vG:\n" << T_l.G()
            //           << "\vU:\n" << T_l.U1() << "\n" << T_l.U2()
            //           << "\vE:\n" << T_l.E_a() << "\n" << T_l.E_b() << "\n" << T_l.E_c()
            //           // << "\vF:\na: " << T_l.Fn_a() << "\nb: " << T_l.Fn_b() << "\nc: " << T_l.Fn_c()
            //           << "\vFr:\na: " << T_l.Fr_a() << "\nb: " << T_l.Fr_b() << "\nc: " << T_l.Fr_c()
            //           << "\vL:\n" << T_l.L()
            //           << "\vdc:\n" << tmp_l
            //           << "\vM^-1:\n" << invertM(T_l.M())
                      // << std::endl;

            // std::cout << row  << "," << col << " upper:\n"
            //           << "C_old:\n" << T_u.C()
            //           << "\vG:\n" << T_u.G()
            //           << "\vU:\n" << T_u.U1() << "\n" << T_u.U2()
            //           << "\vE:\n" << T_u.E_a() << "\n" << T_u.E_b() << "\n" << T_u.E_c()
            //           // << "\vF:\na: " << T_u.Fn_a() << "\nb: " << T_u.Fn_b() << "\nc: " << T_u.Fn_c()
            //           << "\vFr:\na: " << T_u.Fr_a() << "\nb: " << T_u.Fr_b() << "\nc: " << T_u.Fr_c()
            //           << "\vL:\n" << T_u.L()
            //           << "\vdc:\n" << tmp_u
            //           << "\vM^-1:\n" << invertM(T_u.M())
                      // << std::endl;

            T_l.C() += deltaT_ * ( T_l.L() + invertM(T_l.M()) * tmp_l);
            T_u.C() += deltaT_ * ( T_u.L() + invertM(T_u.M()) * tmp_u);

            // std::cout << "C_neu:\n" << T_l.C() << "\n" << T_u.C() << std::endl;
          }

      // update U and F
      assamblyU(mesh_, order_, std::bind(u1_, _1, _2, t_), std::bind(u2_, _1, _2, t_));
      assamblyF(mesh_, order_, orderF_, assamblyLocalLinearF);
    }

    // write data
    if (step_ % writeInterval_ == 0)
      writer_.write();

    // error
    if (writeError_)
      std::cout << t_ << " " << this->l2error() << std::endl;;

  }

  /**
   * get L2 error of current time step
   */
  double l2error()
  {
    return l2Error(mesh_, order_, order_+1, std::bind(cExact_, _1, _2, t_));
  }

  void go()
  {
    for ( ;step_<numSteps_;)
      this->next();
  }

  double getTime() const {return t_;}
  unsigned int getSteps() const { return step_;}


private:
  UniqueSquareGrid & mesh_;
  unsigned int const order_, orderF_;
  std::function<double(double, double, double)> u1_, u2_, f_;
  std::function<double(double,double)> c0_;
  std::function<double(double,double,double)> cExact_;
  BlockMatrix const hatM_;
  std::vector<Tensor> const  hatG_;
  std::vector<BlockMatrix> const hatE_;
  std::vector<double> const hatI_;
  std::function<void(UniqueSquareGrid &, double)> boundaryHandler_;
  unsigned int const numSteps_;
  double const deltaT_;
  VTKwriter & writer_;
  unsigned int const writeInterval_;
  bool const writeError_;

  double t_;
  unsigned int step_;
};


/*** free function implementation ***********/

double l2Error(UniqueSquareGrid const & mesh,
               unsigned int polynomialDegree,
               unsigned int integragradeDegree,
               std::function<double(double,double)> f_ex)
{
  double err = 0.;

#pragma omp parallel for reduction(+: err)
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
