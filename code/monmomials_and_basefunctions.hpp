#ifndef MOMOMIALS_AND_BASISFUNCRIONS_HPP
#define MOMOMIALS_AND_BASISFUNCRIONS_HPP

#include <cmath>
#include <functional>
#include <algorithm>
#include "polynomial.hpp"
#include "Grid.hpp"
#include "DataTypes.hpp"

/**
 * gives the number of base function for a given order in 1D and 2D
 */
inline unsigned int numberOf2DBasefunctions(unsigned int polynomialDegree)
{
  assert(polynomialDegree <= 3); // not implemented
  return (polynomialDegree+1)*(polynomialDegree+2)/2;
}
inline unsigned int numberOf1DBasefunctions(int polynomialDegree)
{
  assert(polynomialDegree <= 2); // not implemented
  return polynomialDegree+1;
}

namespace pol
{
  using std::sqrt;

  //monomials 1D
  Polynomial1D const s0 = Polynomial1D(0,1.);
  Polynomial1D const s1 = Polynomial1D(1,1.);
  Polynomial1D const s2 = Polynomial1D(2,1.);
  Polynomial1D const s3 = Polynomial1D(3,1.);

  //monomials 2D
  // Order 0
  Polynomial2D const c = Polynomial2D(0,0,1.);
  // Order 1
  Polynomial2D const x = Polynomial2D(1,0,1.);
  Polynomial2D const y = Polynomial2D(0,1,1.);
  // Order 2
  Polynomial2D const x2 = Polynomial2D(2,0,1.);
  Polynomial2D const xy = Polynomial2D(1,1,1.);
  Polynomial2D const y2 = Polynomial2D(0,2,1.);
  // Order 3
  Polynomial2D const x3  = Polynomial2D(3,0,1.);
  Polynomial2D const x2y = Polynomial2D(2,1,1.);
  Polynomial2D const xy2 = Polynomial2D(1,2,1.);
  Polynomial2D const y3  = Polynomial2D(0,3,1.);

  // integral over reference triangle
  double const monomialIntegralsRefTriangle[4][4] = {{ 1./2,  1./6 , 1./12, 1./20},
                                                     { 1./6,  1./24, 1./60, 0    },
                                                     { 1./12, 1./60, 0    , 0    },
                                                     { 1./20, 0    , 0    , 0    }};

  // integral over reference edge
  double const monomialIntegralsRefEdge[4] = {1., 1./2., 1./3., 1./4.};

  // 2D base Functions
  Polynomial2D const phi [] = {
    // Order 0
    sqrt(2)*c,
    // Order 1
    -6*x + 2,
    (-2*sqrt(3))*(x + 2*y - 1),
    // Order 2
    sqrt(6)*(x * ( 10*x - 8) + 1),
    sqrt(3)*(x * (5*x - 4) + y*(-15*y +12) -1),
    (3*sqrt(5))*(x*(3*x + 8*y - 4) + y*(3*y - 4) +1),
    // Order 3
    (2*sqrt(2))*(x*(x*(35*x - 45) + 15) - 1),
    (2*sqrt(6))*(x*(x*(21*x -33) + 13) + y*(x*(42*x -24) + 2) -1),
    (2*sqrt(10))*(x*(x*(7*x - 15) + 9) + y*(x*(42*x - 48) + y*(42*x - 6) +6) - 1),
    (2*sqrt(14))*(x*(x*(x - 3) + 3) + y*(x*(12*x - 24) +y*(30*x + 20*y -30) + 12) -1)
  };

  // X derivatives of 2D base functions
  Polynomial2D const dXphi [] = {
    // Order 0
    derive( phi[0], Variable::X ),
    // Order 1
    derive( phi[1], Variable::X),
    derive( phi[2], Variable::X),
    // Order 2
    derive( phi[3], Variable::X ),
    derive( phi[4], Variable::X ),
    derive( phi[5], Variable::X ),
    // Order 3
    derive( phi[6], Variable::X ),
    derive( phi[7], Variable::X ),
    derive( phi[8], Variable::X ),
    derive( phi[9], Variable::X ),
  };

  // Y derivatives of 2D base functions
  Polynomial2D const dYphi [] = {
    // Order 0
    derive( phi[0], Variable::Y ),
    // Order 1
    derive( phi[1], Variable::Y ),
    derive( phi[2], Variable::Y ),
    // Order 2
    derive( phi[3], Variable::Y ),
    derive( phi[4], Variable::Y ),
    derive( phi[5], Variable::Y ),
    // Order 3
    derive( phi[6], Variable::Y ),
    derive( phi[7], Variable::Y ),
    derive( phi[8], Variable::Y ),
    derive( phi[9], Variable::Y ),
  };

  // 1D base functions
  Polynomial1D const phi1D [] = {
    // Order 0
    s0,
    // Order 1
    sqrt(3) * (1-2*s1),
    // Order 2
    sqrt(5) * ( (6 * s1 - 6) * s1 + 1 )
  };
}  // end namespace pol


/**
 * calculates transformation matrices from reference triangle to reference edge
 * returns transposed matrix!
 * linear case!
 */
inline std::vector<BlockMatrix> getLinearTrasformationToRefEdge(int const polynomialDegree)
{
  unsigned int dof1D = numberOf1DBasefunctions(polynomialDegree);
  unsigned int dof2D = numberOf2DBasefunctions(polynomialDegree);

  BlockMatrix t1(dof1D,dof2D);
  BlockMatrix t2(dof1D,dof2D);
  BlockMatrix t3(dof1D,dof2D);

  using std::sqrt;
  switch(polynomialDegree)
    {
    case 2:
      t1(0,3) = sqrt(6)/3    ; t1(1,3) = sqrt(2)       ; t1(2,3) = sqrt(30)/3;
      t1(0,4) = 19*sqrt(3)/3 ; t1(1,4) =-8             ; t1(2,4) = sqrt(15)/3;
      t1(0,5) = sqrt(5)      ; t1(1,5) = 0             ; t1(2,5) = -1;

      t2(0,3) = sqrt(6);
      t2(0,4) =-2*sqrt(3)    ; t2(1,4) =-2/3           ; t2(2,4) = -sqrt(15)/2;
      t2(0,5) =-3*sqrt(5)/5  ; t2(1,5) =-3*sqrt(15)/10 ; t2(2,5) = 3/2;

      t2(0,3) =-2*sqrt(6)/3  ; t3(1,3) =-sqrt(2)       ; t2(2,3) = sqrt(30)/3;
      t2(0,4) =-4*sqrt(3)/3  ; t2(1,4) =-1/2           ; t2(2,4) = sqrt(15)/6;
      t2(0,5) = 0            ; t2(1,5) = sqrt(15)/2    ; t2(2,5) = 3/2;
    case 1:
      t1(0,0) =  sqrt(2)     ; t1(1,0) = 0.;
      t1(0,1) = -1.          ; t1(1,1) = -sqrt(3);
      t1(0,2) = -sqrt(3)     ; t1(1,2) = 1.;

      t2(0,0) =  sqrt(2)     ; t2(1,0) = 0.;
      t2(0,1) =  2.          ; t2(1,1) = 0.;
      t2(0,2) =  0.          ; t2(1,2) = -2.;

      t3(0,0) =  sqrt(2)     ; t3(1,0) = 0.;
      t3(0,1) = -1.          ; t3(1,1) = sqrt(3);
      t3(0,2) =  sqrt(3)     ; t3(1,2) = 1.;

      break;
    default:
      assert(-1);
    }



  return {t1,t2,t3};
}

/**
 * Transformation Matrix from 2D Polynomial to 2D Monomial base space
 */
inline BlockMatrix getPolynomialMapping (unsigned int polynomialDegree)
{
  using std::sqrt;
  BlockMatrix T;
  switch(polynomialDegree)
    {
    case 0:
      T = BlockMatrix(1);
      T(0,0) = 1./sqrt(2);
      break;
    case 1:
      T = BlockMatrix(3);
      T(0,0) = 1./sqrt(2);
      T(0,1) = 1./(3.*sqrt(2));
      T(0,2) = 1./(3.*sqrt(2));
      T(1,1) = -1./6.;
      T(1,2) = 1./12.;
      T(2,2) = -1./(4.*sqrt(3));
      break;
    case 2:
      T = BlockMatrix(6);
      T(0,0) = 1./sqrt(2);
      T(0,1) = 1./(3.*sqrt(2));
      T(0,2) = 1./(3.*sqrt(2));
      T(0,3) = sqrt(2)/12.;
      T(0,5) = sqrt(2)/12.;
      T(0,4) = sqrt(2)/24.;
      T(1,1) = -1./6.;
      T(1,2) = 1./12.;
      T(1,3) = -2./15;
      T(1,5) = 1./15;
      T(1,4) = -1./60;
      T(2,2) = -1./(4.*sqrt(3));
      T(2,5) = -sqrt(3)/15.;
      T(2,4) = -sqrt(3)/60.;
      T(3,3) = sqrt(6)/60.;
      T(3,5) = sqrt(6)/180.;
      T(3,4) = sqrt(6)/120.;
      T(4,5) = -sqrt(3)/45.;
      T(4,4) = sqrt(3)/120.;
      T(5,4) = sqrt(5)/120.;
    default:
      assert(-1);
    }

  return T;
}




inline double integradeOverRefTriangle(Polynomial2D const & pol)
{
  assert(pol.getOrder() <= 3); // not implemented!

  double integral = 0.;

  for (unsigned int x=0; x<=pol.getOrder(); ++x)
    for (unsigned int y=0; y<=pol.getOrder(); ++y)
      integral += pol.get(x,y)*pol::monomialIntegralsRefTriangle[x][y];
  return integral;
}

inline double integradeOverRefEdge(Polynomial1D const & pol)
{
  assert(pol.getOrder() <= 3); // not implemented!
  double integral = 0.;
  for (unsigned int s=0; s<=pol.getOrder(); ++s)
    integral += pol.get(s) * pol::monomialIntegralsRefEdge[s];
  return integral;
}


/**
 * weighted sum of a function evaluated at q1,q2
 * returns value + sum(weight * f(q1, q2))
 * helper function for Gaus integration
 */
template<class InputIt1, class InputIt2, class InputIt3, class T1, class T2, class T3>
T1 weighted_sum(InputIt1 weight_first, InputIt1 weight_last,
                InputIt2 q1, InputIt3 q2,
                std::function<T1(T2,T3)> f,
                T1 value)
{
  while (weight_first != weight_last) {
    value += *weight_first * f(*q1, *q2);
    ++weight_first;
    ++q1;
    ++q2;
  }
  return value;
}

inline Polynomial1D reconstructFunction1D (unsigned int polynomialDegree, std::vector<double> coeficents)
{
  return std::inner_product(std::begin(coeficents), std::end(coeficents),
                            std::begin(pol::phi1D), Polynomial1D(polynomialDegree));
}
inline Polynomial2D reconstructFunction2D (unsigned int polynomialDegree, std::vector<double> coeficents)
{
  return std::inner_product(std::begin(coeficents), std::end(coeficents),
                            std::begin(pol::phi), Polynomial2D(polynomialDegree));
}

/**
 * unsing gausion quadrature rule to integrate a function f over the reference triangle
 * is exact for polynomials f of order 'order'
 */
inline double integradeOverRefTriangle_gaus(std::function<double(double,double)> f, int order)
{
  using std::begin;
  using std::end;

  std::vector<double> w,x1,x2;

  switch(order){
  case 0:
  case 1:
    return f(1./3,1./3) * 0.5;
    break;
  case 2:
    w  = {1./6, 1./6, 1./6};
    x1 = {1./6, 2./3, 1./6};
    x2 = {1./6, 1./6, 2./3};
    break;
  default:
    std::cerr << "gaus integration for given order is not implemented!" << std::endl;
    return -1;
  }

  return weighted_sum(begin(w), end(w), begin(x1), begin(x2), f, 0.);
}

/**
 * use L2-Projection to represent a continues function f to the discreet polynomial space on element T_k
 */
inline std::vector<double> l2Projection (unsigned int polynomialDegree,
                                         std::function<double(double,double)> const & f,
                                         Jakobian const & B_k,
                                         Point const & A_k)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);

  std::vector<double> F (dof);

  // hatM(i,:) F_k(i) = int_hatT {phi(i)*f(x)} d hatx
  // x = T_k(x^) = B_k * x^ + a_k
  // here hatM == I
  // =>  F_k(i) = int_hatT { phi(i)*f(t_k(x^)) } = int_hatT { phi(i) * f(B_k*x^+a_k) }
  for (unsigned int i=0; i<dof; ++i)
    {
      F[i] = integradeOverRefTriangle_gaus([&f, i, &B_k, &A_k](double x1_hat, double x2_hat)
                                           {
                                             double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                             double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                             return pol::phi[i](x1_hat,x2_hat) * f(x1,x2);
                                           },
                                           polynomialDegree);
    }

  return F;
}

/**
 * f = F^T*B = sum{F_j * B_j}
 * int{b_i * f} = sum{int{F_j * B_j * B_i} }= sum {F_j * int{δ_{ij}} } = F_i
 */
inline std::vector<double> l2Projection (unsigned int polynomialDegree,
                                         Polynomial2D const & f)
{
  std::cout << "l2Projetion: " << polynomialDegree << std::endl;
  unsigned int dof (numberOf2DBasefunctions(polynomialDegree));

  std::vector<double> F;
  F.reserve(dof);

  for (unsigned int i=0; i<dof; ++i)
    F.push_back(integradeOverRefTriangle(pol::phi[i] * f));
  //TODO poly hier ist natürlich nochmal ein gard höher

  return F;
}

inline std::vector<double> lagrangeProjection(std::vector<double> const & c, unsigned int polynomialDegree)
{
  Polynomial2D pol = reconstructFunction2D(polynomialDegree, c);
  std::vector<double> l;

  switch (polynomialDegree)
    {
    case 0: // locally constant -> evaluate in ref. triangle center
      l.push_back(pol(1./3., 1./3.));
      break;
    case 1: // locally linear -> evaluate at corners
      l.push_back(pol(0,0));
      l.push_back(pol(1,0));
      l.push_back(pol(0,1));
      break;
    default: // locally quadratic -> evaluate at corners and edge mid points
      if (polynomialDegree > 2)
        std::cerr << "ploting is limited to 2 order" << std::endl;

      l.push_back(pol(0,0));
      l.push_back(pol(1,0));
      l.push_back(pol(0,1));
      l.push_back(pol(0.5,0));
      l.push_back(pol(0.5,0.5));
      l.push_back(pol(0,0.5));
    }
  return l;
}

#endif
