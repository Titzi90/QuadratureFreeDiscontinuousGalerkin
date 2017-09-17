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
  assert(polynomialDegree <= 4); // not implemented
  return (polynomialDegree+1)*(polynomialDegree+2)/2;
}
inline unsigned int numberOf1DBasefunctions(int polynomialDegree)
{
  assert(polynomialDegree <= 8); // not implemented
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
  Polynomial1D const s4 = Polynomial1D(4,1.);
  Polynomial1D const s5 = Polynomial1D(5,1.);
  Polynomial1D const s6 = Polynomial1D(6,1.);
  Polynomial1D const s7 = Polynomial1D(7,1.);
  Polynomial1D const s8 = Polynomial1D(8,1.);

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
  // Order 4
  Polynomial2D const x4  = Polynomial2D(4,0,1.);
  Polynomial2D const x3y1= Polynomial2D(3,1,1.);
  Polynomial2D const x2y2= Polynomial2D(2,2,1.);
  Polynomial2D const x1y3= Polynomial2D(1,3,1.);
  Polynomial2D const   y4= Polynomial2D(0,4,1.);


  // integral over reference triangle
  /*
  double const monomialIntegralsRefTriangle[5][5] = {{ 1./2,  1./6  , 1./12 , 1./20 , 1./30},
                                                     { 1./6,  1./24 , 1./60 , 1./120, 0    },
                                                     { 1./12, 1./60 , 1./180, 0     , 0    },
                                                     { 1./20, 1./120, 0     , 0     , 0    },
                                                     { 1./30, 0     , 0     , 0     , 0    }};
  */
  double const monomialIntegralsRefTriangle[12][12] =
    {{ 1./2, 1./6, 1./12, 1./20, 1./30, 1./42, 1./56, 1./72, 1./90, 1./110, 1./132, 1./156, },
     { 1./6, 1./24, 1./60, 1./120, 1./210, 1./336, 1./504, 1./720, 1./990, 1./1320, 1./1716, 0, },
     { 1./12, 1./60, 1./180, 1./420, 1./840, 1./1512, 1./2520, 1./3960, 1./5940, 1./8580, 0, 0, },
     { 1./20, 1./120, 1./420, 1./1120, 1./2520, 1./5040, 1./9240, 1./15840, 1./25740, 0, 0, 0, },
     { 1./30, 1./210, 1./840, 1./2520, 1./6300, 1./13860, 1./27720, 1./51480, 0, 0, 0, 0, },
     { 1./42, 1./336, 1./1512, 1./5040, 1./13860, 1./33264, 1./72072, 0, 0, 0, 0, 0, },
     { 1./56, 1./504, 1./2520, 1./9240, 1./27720, 1./72072, 0, 0, 0, 0, 0, 0, },
     { 1./72, 1./720, 1./3960, 1./15840, 1./51480, 0, 0, 0, 0, 0, 0, 0, },
     { 1./90, 1./990, 1./5940, 1./25740, 0, 0, 0, 0, 0, 0, 0, 0, },
     { 1./110, 1./1320, 1./8580, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
     { 1./132, 1./1716, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
     { 1./156, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }
    };



  // integral over reference edge
  double const monomialIntegralsRefEdge[13] = {1., 1./2., 1./3., 1./4., 1./5, 1./6, 1./7, 1./8, 1./9,
                                               1./10, 1./11, 1./12, 1./13};

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
    (2*sqrt(14))*(x*(x*(x - 3) + 3) + y*(x*(12*x - 24) +y*(30*x + 20*y -30) + 12) -1),
    // Order 4
    sqrt(10)*(1+(-24+(126+(-224+126*x)*x)*x)*x),
    sqrt(30)*(1+(-22+(105+(-168+84*x)*x)*x)*x+(-2+(42+(-168+168*x)*x)*x)*y),
    5*sqrt(2)*(1+(-18+(69+(-88+36*x)*x)*x)*x+(-6+(102+(-312+216*x)*x)*x
                                              +(6+(-96+216*x)*x)*y)*y),
    sqrt(70)*(1+(-12+(30+(-28+9*x)*x)*x)*x+(-12+(132+(-228+108*x)*x)*x
                                            +(30+(-300+270*x)*x+(-20+180*x)*y)*y)*y),
    3*sqrt(10)*(1+(-4+(6+(-4+x)*x)*x)*x+(-20+(60+(-60+20*x)*x)*x
                                         +(90+(-180+90*x)*x+(-140+140*x+70*y)*y)*y)*y)
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
    // Order 4
    derive( phi[10], Variable::X ),
    derive( phi[11], Variable::X ),
    derive( phi[12], Variable::X ),
    derive( phi[13], Variable::X ),
    derive( phi[14], Variable::X )
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
    // Order 4
    derive( phi[10], Variable::Y ),
    derive( phi[11], Variable::Y ),
    derive( phi[12], Variable::Y ),
    derive( phi[13], Variable::Y ),
    derive( phi[14], Variable::Y )
  };

  // 1D base functions
  Polynomial1D const phi1D [] = {
    s0,
    sqrt(3) * (1-2*s1),
    sqrt(5) * ( (6 * s1 - 6) * s1 + 1 ),
    20*sqrt(7)*s3 - 30*sqrt(7)*s2 + 12*sqrt(7)*s1 - sqrt(7),
    210*s4 - 420*s3 + 270*s2 - 60*s1 + 3,
    252*sqrt(11)*s5 - 630*sqrt(11)*s4 + 560*sqrt(11)*s3 - 210*sqrt(11)*s2 + 30*sqrt(11)*s1 - sqrt(11),
    924*sqrt(13)*s6 - 2772*sqrt(13)*s5 + 3150*sqrt(13)*s4 - 1680*sqrt(13)*s3 + 420*sqrt(13)*s2 - 42*sqrt(13)*s1 + sqrt(13),
    3432*sqrt(15)*s7 - 12012*sqrt(15)*s6 + 16632*sqrt(15)*s5 - 11550*sqrt(15)*s4 + 4200*sqrt(15)*s3 - 756*sqrt(15)*s2 + 56*sqrt(15)*s1 - sqrt(15),
    12870*sqrt(17)*s8 - 51480*sqrt(17)*s7 + 84084*sqrt(17)*s6 - 72072*sqrt(17)*s5 + 34650*sqrt(17)*s4 - 9240*sqrt(17)*s3 + 1260*sqrt(17)*s2 - 72*sqrt(17)*s1 + sqrt(17)
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
    case 4:
      t1(0,10) = sqrt(10)/5;  t1(1,10) = sqrt(30)/5;  t1(2,10) = sqrt(2);  t1(3,10) = -sqrt(70)/5;  t1(4,10) = 3.*sqrt(10)/5;
      t1(0,11) = sqrt(30)/5;  t1(1,11) = sqrt(10)/2;  t1(2,11) = sqrt(6)/2;  t1(4,11) = -2.*sqrt(30)/5;
      t1(0,12) = sqrt(2);  t1(1,12) = sqrt(6)/2;  t1(2,12) = -3.*sqrt(10)/14;  t1(3,12) = 4.*sqrt(14)/7;  t1(4,12) = 6.*sqrt(2)/7;
      t1(0,13) = sqrt(70)/5;  t1(2,13) = -4.*sqrt(14)/7;  t1(3,13) = -sqrt(10)/2;  t1(4,13) = -3.*sqrt(70)/70;
      t1(0,14) = 3.*sqrt(10)/5;  t1(1,14) = -2.*sqrt(30)/5;  t1(2,14) = 6.*sqrt(2)/7;  t1(3,14) = 3.*sqrt(70)/70;  t1(4,14) = sqrt(10)/70;

      t2(0,10) = sqrt(10);
      t2(1,11) = -sqrt(10);
      t2(2,12) = sqrt(10);
      t2(3,13) = sqrt(10);
      t2(4,14) = sqrt(10);

      t3(0,10) = sqrt(10)/5;  t3(1,10) = -sqrt(30)/5;  t3(2,10) = sqrt(2);  t3(3,10) = sqrt(70)/5;  t3(4,10) = 3.*sqrt(10)/5;
      t3(0,11) = -sqrt(30)/5;  t3(1,11) = sqrt(10)/2;  t3(2,11) = -sqrt(6)/2;  t3(4,11) = 2.*sqrt(30)/5;
      t3(0,12) = sqrt(2);  t3(1,12) = -sqrt(6)/2;  t3(2,12) = -3.*sqrt(10)/14;  t3(3,12) = -4.*sqrt(14)/7;  t3(4,12) = 6.*sqrt(2)/7;
      t3(0,13) = -sqrt(70)/5;  t3(2,13) = 4.*sqrt(14)/7;  t3(3,13) = -sqrt(10)/2;  t3(4,13) = 3.*sqrt(70)/70;
      t3(0,14) = 3.*sqrt(10)/5;  t3(1,14) = 2.*sqrt(30)/5;  t3(2,14) = 6.*sqrt(2)/7;  t3(3,14) = -3.*sqrt(70)/70;  t3(4,14) = sqrt(10)/70;

    case 3:
      t1(0,6) = sqrt(2)/2;  t1(1,6) = sqrt(6)/2;  t1(2,6) = sqrt(10)/2;  t1(3,6) = -sqrt(14)/2;
      t1(0,7) = sqrt(6)/2;  t1(1,7) = 11.*sqrt(2)/10;  t1(2,7) = sqrt(30)/10;  t1(3,7) = 3.*sqrt(42)/10;
      t1(0,8) = sqrt(10)/2;  t1(1,8) = sqrt(30)/10;  t1(2,8) = -3.*sqrt(2)/2;  t1(3,8) = -sqrt(70)/10;
      t1(0,9) = sqrt(14)/2;  t1(1,9) = -3.*sqrt(42)/10;  t1(2,9) = sqrt(70)/10;  t1(3,9) = sqrt(2)/10;

      t2(0,6) = -2.*sqrt(2);
      t2(1,7) = 2.*sqrt(2);
      t2(2,8) = -2.*sqrt(2);
      t2(3,9) = -2.*sqrt(2);

      t3(0,6) = sqrt(2)/2;  t3(1,6) = -sqrt(6)/2;  t3(2,6) = sqrt(10)/2;  t3(3,6) = sqrt(14)/2;
      t3(0,7) = -sqrt(6)/2;  t3(1,7) = 11.*sqrt(2)/10;  t3(2,7) = -sqrt(30)/10;  t3(3,7) = 3.*sqrt(42)/10;
      t3(0,8) = sqrt(10)/2;  t3(1,8) = -sqrt(30)/10;  t3(2,8) = -3.*sqrt(2)/2;  t3(3,8) = sqrt(70)/10;
      t3(0,9) = -sqrt(14)/2;  t3(1,9) = -3.*sqrt(42)/10;  t3(2,9) = -sqrt(70)/10;  t3(3,9) = sqrt(2)/10;

    case 2:
      t1(0,3) = sqrt(6)/3    ; t1(1,3) = sqrt(2)       ; t1(2,3) = sqrt(30)/3;
      t1(0,4) = -sqrt(3)/3   ; t1(1,4) = 2             ; t1(2,4) = -sqrt(15)/3;
      t1(0,5) = sqrt(5)      ; t1(1,5) = 0             ; t1(2,5) = -1;

      t2(0,3) = sqrt(6);
      t2(1,4) = -3./2;  t2(2,4) = -sqrt(15)/2;
      t2(1,5) = -sqrt(15)/2;  t2(2,5) = 3./2;

      t3(0,3) = sqrt(6)/3;  t3(1,3) = -sqrt(2);  t3(2,3) = sqrt(30)/3;
      t3(0,4) = -4.*sqrt(3)/3;  t3(1,4) = -1./2;  t3(2,4) = sqrt(15)/6;
      t3(1,5) = sqrt(15)/2;  t3(2,5) = 3./2;

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
 * Transformation Matrix (transposed) from 2D Polynomial to 2D Monomial base space
 */
inline BlockMatrix get1DPolynomialMapping (unsigned int polynomialDegree)
{
  using std::sqrt;
  if(polynomialDegree>8)
    {
      std::cerr << "maping from monomial into hirachical space for given order not implemented" << std::endl;
      assert(false);
    }

  BlockMatrix T(polynomialDegree+1);
  switch(polynomialDegree)
    {
    case 8:
      T(0,8) = 1./9.;
      T(1,8) = -4.*sqrt(3)/45.;
      T(2,8) = 28.*sqrt(5)/495.;
      T(3,8) = 14.*sqrt(7)/495.;
      T(4,8) = 14./429.;
      T(5,8) = 4.*sqrt(11)/1287.;
      T(6,8) = 4.*sqrt(13)/6435.;
      T(7,8) = sqrt(15)/12870.;
      T(8,8) = sqrt(17)/218790.;
    case 7:
      T(0,7) = 1./8.;
      T(1,7) = -7.*sqrt(3)/72.;
      T(2,7) = 7.*sqrt(5)/120.;
      T(3,7) = 7.*sqrt(7)/264.;
      T(4,7) = 7./264.;
      T(5,7) = 7.*sqrt(11)/3432.;
      T(6,7) = sqrt(13)/3432.;
      T(7,7) = sqrt(15)/51480.;
    case 6:
      T(0,6) = 1./7.;
      T(1,6) = -3.*sqrt(3)/28.;
      T(2,6) = 5.*sqrt(5)/84.;
      T(3,6) = sqrt(7)/42.;
      T(4,6) = 3./154.;
      T(5,6) = sqrt(11)/924.;
      T(6,6) = sqrt(13)/12012.;
    case 5:
      T(0,5) = 1./6.;
      T(1,5) = -5.*sqrt(3)/42.;
      T(2,5) = 5.*sqrt(5)/84.;
      T(3,5) = 5.*sqrt(7)/252.;
      T(4,5) = 1./84.;
      T(5,5) = sqrt(11)/2772.;
    case 4:
      T(0,4) = 1./5.;
      T(1,4) = -2.*sqrt(3)/15.;
      T(2,4) = 2.*sqrt(5)/35.;
      T(3,4) = sqrt(7)/70.;
      T(4,4) = 1./210.;
    case 3:
      T(0,3) = 1./4.;
      T(1,3) = -3.*sqrt(3)/20.;
      T(2,3) = sqrt(5)/20.;
      T(3,3) = sqrt(7)/140.;
    case 2:
      T(0,2) = 1./3.;
      T(1,2) = -sqrt(3)/6.;
      T(2,2) = sqrt(5)/30.;
    case 1:
      T(0,1) = 1./2.;
      T(1,1) = -sqrt(3)/6.;
    case 0:
      T(0,0) = 1.;
    }

  return T;
}



inline double integradeOverRefTriangle(Polynomial2D const & pol)
{
  assert(pol.getOrder() <= 11); // not implemented!

  double integral = 0.;

  for (unsigned int x=0; x<=pol.getOrder(); ++x)
    for (unsigned int y=0; y<=pol.getOrder(); ++y)
      integral += pol.get(x,y)*pol::monomialIntegralsRefTriangle[x][y];
  return integral;
}

inline double integradeOverRefEdge(Polynomial1D const & pol)
{
  assert(pol.getOrder() <= 12); // not implemented!
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
template<class InputIt1, class InputIt2, class InputIt3, class T1, typename Function>
T1 weighted_sum(InputIt1 weight_first, InputIt1 weight_last,
                InputIt2 q1, InputIt3 q2,
                Function f,
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
template<class InputIt1, class InputIt2, class T1, typename Function>
T1 weighted_sum(InputIt1 weight_first, InputIt1 weight_last,
                InputIt2 q,
                Function f,
                T1 value)
{
  while (weight_first != weight_last) {
    value += *weight_first * f(*q);
    ++weight_first;
    ++q;
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
template<typename Function>
double integradeOverRefTriangle_gaus(Function f, int order)
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
  case 3:
    w  = {0.159020691, 0.159020691, 0.090979309, 0.090979309};
    x1 = {0.666390246, 0.178558728, 0.280019915, 0.075031109};
    x2 = {0.178558728, 0.666390246, 0.075031109, 0.280019915};
    break;
  case 4:
    w  = {0.111690794839005, 0.111690794839005, 0.111690794839005,
          0.054975871827661, 0.054975871827661, 0.054975871827661};
    x1 = {0.445948490915965, 0.108103018168070, 0.445948490915965,
          0.091576213509771, 0.816847572980458, 0.091576213509771};
    x2 = {0.108103018168070, 0.445948490915965, 0.445948490915965,
          0.816847572980458, 0.091576213509771, 0.091576213509771};
    break;
  case 5:
    {
    double a1 = (6.-std::sqrt(15))/21;     double a2 = (6.+std::sqrt(15))/21;
    double w1 = (155.-std::sqrt(15))/2400; double w2 = (155.+std::sqrt(15))/2400;

    w  = {9./80,       w1,       w1,     w1,       w2,       w2,     w2};
    x1 = {1./3 ,       a1, 1.-2.*a1,     a1,       a2, 1.-2.*a2,     a2};
    x2 = {1./3 , 1.-2.*a1,       a1,     a1, 1.-2.*a2,       a2,     a2};
    }
    break;
  case 6:
    w  = {0.025422453185103, 0.025422453185103, 0.025422453185103,
          0.058393137863189, 0.058393137863189, 0.058393137863189,
          0.041425537809187, 0.041425537809187, 0.041425537809187,
          0.041425537809187, 0.041425537809187, 0.041425537809187};
    x1 = {0.063089014491502, 0.873821971016996, 0.063089014491502,
          0.249286745170910, 0.501426509658179, 0.249286745170910,
          0.310352451033785, 0.053145049844816, 0.636502499121399,
          0.053145049844816, 0.636502499121399, 0.310352451033785};
    x2 = {0.063089014491502, 0.063089014491502, 0.873821971016996,
          0.249286745170910, 0.249286745170910, 0.501426509658179,
          0.053145049844816, 0.310352451033785, 0.053145049844816,
          0.636502499121399, 0.310352451033785, 0.636502499121399};
    break;
    // komisches werte aus matlab
  case 7:
    x1 = {
      0.054830900955589,0.175654279195255,0.343651813106453,0.533230731173959,0.715527432866568,0.862793031223432,0.952646581185227,0.048991501878362,
      0.156947392786903,0.307053470832875,0.476442551784230,0.639324960202548,0.770907019092335,0.851191316541618,0.039548223967455,0.126695251279609,
      0.247867874404688,0.384606636317686,0.516092908865112,0.622312080263295,0.687121307473297,0.028131280268461,0.090120345868446,0.176312358556585,
      0.273576813165278,0.367105088607705,0.442660473419548,0.488760306780644,0.016714336569468,0.053545440457283,0.104756842708482,0.162546990012870,
      0.218117268350298,0.263008866575801,0.290399306087990,0.007271058658560,0.023293298949990,0.045571246280295,0.070711074546325,0.094885217012863,
      0.114413927746761,0.126329297019669,0.001431659581333,0.004586412541638,0.008972904006717,0.013922895156596,0.018682744348843,0.022527915615664,
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
    w = {
      0.000337590756711,   0.001774485071438,   0.004297910087982,   0.006935542753734,   0.008247603013529,
      0.007154643779096,   0.003623466079726,   0.000729242610652,   0.003833132573485,   0.009284078756889,
      0.014981729219389,   0.017815960400676,   0.015455017662734,   0.007827186648495,   0.000995500091625,
      0.005232667115688,   0.012673836002093,   0.020451784622510,   0.024320836374897,   0.021097877818152,
      0.010685010601315,   0.001089695284832,   0.005727787200653,   0.013873046771564,   0.022386952504607,
      0.026622097721383,   0.023094179670909,   0.011696036764419,   0.000995500091625,   0.005232667115688,
      0.012673836002093,   0.020451784622510,   0.024320836374897,   0.021097877818152,   0.010685010601315,
      0.000729242610652,   0.003833132573485,   0.009284078756889,   0.014981729219389,   0.017815960400676,
      0.015455017662734,   0.007827186648495,   0.000337590756711,   0.001774485071438,   0.004297910087982,
      0.006935542753734,   0.008247603013529,   0.007154643779096,   0.003623466079726    };
    break;
  case 8:
    x1 = {
      0.043747744905146,   0.141499854650117,   0.281129831011298,   0.445782964189930,   0.615597503479888,
      0.770091559088435,   0.890634557136406,   0.962718034592387,   0.040096165611934,   0.129689007248651,
      0.257664213026852,   0.408573918447520,   0.564214212717487,   0.705812808329004,   0.816294206251800,
      0.882360949948550,   0.034045272688803,   0.110117702008052,   0.218780231495104,   0.346916226396932,
      0.479068919277106,   0.599298939439204,   0.693107643137351,   0.749204286556732,   0.026410684460876,
      0.085424014895553,   0.169719176965060,   0.269120916535959,   0.371638617134706,
      0.464907281899148,   0.537679560614672,   0.581196637484813,   0.018223270829094,   0.058942242146592,
      0.117105580179370,   0.185692398660614,   0.256429218282022,   0.320784238705222,
      0.370996831485534,   0.401023447367823,   0.010588682601167,   0.034248555034093,   0.068044525649326,
      0.107897088799642,   0.148998916139621,   0.186392581165165,   0.215568748962855,
      0.233015798295905,   0.004537789678036,   0.014677249793495,   0.029160544117579,   0.046239396749053,
      0.063853622699241,   0.079878712275365,   0.092382185848406,   0.099859134904086,
      0.000886210384824,   0.002866402392029,   0.005694926133132,   0.009030351006644,   0.012470331936840,
      0.015599961515934,   0.018041834963800,   0.019502050260250    };
    x2 = {
      0.000886210384824,   0.002866402392029,   0.005694926133132,   0.009030351006644,   0.012470331936840,
      0.015599961515934,   0.018041834963800,   0.019502050260250,   0.004537789678036,
      0.014677249793495,   0.029160544117579,   0.046239396749053,   0.063853622699241,   0.079878712275365,
      0.092382185848406,   0.099859134904086,   0.010588682601167,   0.034248555034093,
      0.068044525649326,   0.107897088799642,   0.148998916139621,   0.186392581165165,   0.215568748962855,
      0.233015798295905,   0.018223270829094,   0.058942242146592,   0.117105580179370,
      0.185692398660614,   0.256429218282022,   0.320784238705222,   0.370996831485534,   0.401023447367823,
      0.026410684460876,   0.085424014895553,   0.169719176965060,   0.269120916535959,
      0.371638617134706,   0.464907281899148,   0.537679560614672,   0.581196637484814,   0.034045272688803,
      0.110117702008052,   0.218780231495104,   0.346916226396932,   0.479068919277106,
      0.599298939439204,   0.693107643137351,   0.749204286556732,   0.040096165611934,   0.129689007248651,
      0.257664213026852,   0.408573918447521,   0.564214212717487,   0.705812808329004,
      0.816294206251800,   0.882360949948550,   0.043747744905146,   0.141499854650117,   0.281129831011298,
      0.445782964189930,   0.615597503479888,   0.770091559088435,   0.890634557136406,
      0.962718034592387    };
    w = {
      0.000166783703248,   0.000903105459519,   0.002299877901746,   0.004008629765696,   0.005367509486579,
      0.005694398702308,   0.004611922695459,   0.002254906358040,   0.000366394040825,
      0.001983961575145,   0.005052421438156,   0.008806244431697,   0.011791460746205,   0.012509578034170,
      0.010131571367319,   0.004953626979826,   0.000516861727437,   0.002798718572468,
      0.007127308256396,   0.012422720355804,   0.016633880716426,   0.017646908496911,   0.014292321640317,
      0.006987941703713,   0.000597556249615,   0.003235665720862,   0.008240052156051,
      0.014362205192963,   0.019230828768754,   0.020402014502054,   0.016523696115091,   0.008078927139199,
      0.000597556249615,   0.003235665720862,   0.008240052156051,   0.014362205192963,
      0.019230828768754,   0.020402014502054,   0.016523696115091,   0.008078927139199,   0.000516861727437,
      0.002798718572468,   0.007127308256396,   0.012422720355804,   0.016633880716426,
      0.017646908496911,   0.014292321640317,   0.006987941703713,   0.000366394040825,   0.001983961575145,
      0.005052421438156,   0.008806244431697,   0.011791460746205,   0.012509578034170,
      0.010131571367319,   0.004953626979826,   0.000166783703248,   0.000903105459519,   0.002299877901746,
      0.004008629765696,   0.005367509486579,   0.005694398702308,   0.004611922695459,
      0.002254906358040    };
    break;
  default:
    std::cerr << "2d gaus integration for given order(" << order << ") is not implemented!" << std::endl;
    return -1;
  }

  return weighted_sum(begin(w), end(w), begin(x1), begin(x2), f, 0.);
}

template<typename Function>
double integradeOverRefEdge_gaus(Function f, int order)
{
  using std::begin;
  using std::end;
  using std::sqrt;

  std::vector<double> w, x;

  switch(order){
  case 0:
  case 1:
    x = {0};
    w = {2};
    break;
  case 2:
    x = {-1./sqrt(3), 1./sqrt(3)};
    w = {1. , 1.};
    break;
  case 3:
    x = {-sqrt(3./5.), 0., sqrt(3./5.)};
    w = {5./9., 8./9., 5./9.};
    break;
  case 4:
    x = {-sqrt(3./7.+2./7.*sqrt(6./5.)), -sqrt(3./7.-2./7.*sqrt(6./5.)),
          sqrt(3./7.-2./7.*sqrt(6./5.)),  sqrt(3./7.+2./7.*sqrt(6./5.)) };
    w = {(18.-sqrt(30.))/36., (18.+sqrt(30.))/36.,
         (18.+sqrt(30.))/36., (18.-sqrt(30.))/36. };
    break;
  case 5:
    x = {-1./3*sqrt(5+2*sqrt(10./7)), -1./3*sqrt(5-2*sqrt(10./7)), 0,
         1./3*sqrt(5-2*sqrt(10./7)),  1./3*sqrt(5+2*sqrt(10./7)) };
    w = {(322.-13.*sqrt(70))/900, (322.+13.*sqrt(70))/900, 128./225,
         (322.+13.*sqrt(70))/900, (322.-13.*sqrt(70))/900};
    break;
  case 6:
    x = { 0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
          0.2386191860831969, -0.9324695142031521,  0.9324695142031521};
    w = { 0.3607615730481386,  0.3607615730481386,  0.4679139345726910,
          0.4679139345726910,  0.1713244923791704,  0.171324492379170};
    break;
  case 7:
    x = { 0.0000000000000000,  0.4058451513773972, -0.4058451513773972,
          -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,
          0.9491079123427585};
    w = { 0.4179591836734694,  0.3818300505051189,  0.3818300505051189,
          0.2797053914892766,  0.2797053914892766,  0.1294849661688697,
          0.1294849661688697};
    break;
  case 8:
    x = {-0.1834346424956498,  0.1834346424956498, -0.5255324099163290,
         0.5255324099163290, -0.7966664774136267,  0.7966664774136267,
         -0.9602898564975363,  0.9602898564975363};
    w = { 0.3626837833783620,  0.3626837833783620,  0.3137066458778873,
          0.3137066458778873,  0.2223810344533745,  0.2223810344533745,
          0.1012285362903763,  0.1012285362903763};
    break;
  case 9:
    x = { 0.0000000000000000, -0.8360311073266358,  0.8360311073266358,
          -0.9681602395076261,  0.9681602395076261, -0.3242534234038089,
          0.3242534234038089, -0.6133714327005904,  0.6133714327005904};
    w = { 0.3302393550012598,  0.1806481606948574,  0.1806481606948574,
          0.0812743883615744,  0.0812743883615744,  0.3123470770400029,
          0.3123470770400029,  0.2606106964029354,  0.2606106964029354};
    break;
  default:
    std::cerr << "1d gaus integration for given order(" << order << ") is not implemented!" << std::endl;
    assert(-1);
    return -1;
    }
  // map gaus quad role to interfal (0,1) (ref edge)
  x = (x+1)*0.5;
  w = w*0.5;

  return weighted_sum(begin(w), end(w), begin(x), f, 0.);
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

  std::vector<double> F;
  F.reserve(dof);

  // hatM(i,:) F_k(i) = int_hatT {phi(i)*f(x)} d hatx
  // x = T_k(x^) = B_k * x^ + a_k
  // here hatM == I
  // =>  F_k(i) = int_hatT { phi(i)*f(t_k(x^)) } = int_hatT { phi(i) * f(B_k*x^+a_k) }
  for (unsigned int i=0; i<dof; ++i)
    F.push_back(integradeOverRefTriangle_gaus([&f, i, &B_k, &A_k](double x1_hat, double x2_hat)
                                              {
                                                double x1 = B_k[0][0]*x1_hat + B_k[0][1]*x2_hat + A_k.x;
                                                double x2 = B_k[1][0]*x2_hat + B_k[1][1]*x2_hat + A_k.y;
                                                return pol::phi[i](x1_hat,x2_hat) * f(x1,x2);
                                              },
                                              2*polynomialDegree)
                                              // polynomialDegree) //TODO Balti fragen
                );

  return F;
}
/**
 * l2Projection onto edge
 */
inline std::vector<double> l2Projection_edge (unsigned int polynomialDegree,
                                              std::function<double(double,double)> const & f,
                                              Jakobian const & B_k,
                                              Point const & A_k,
                                              int edgeID)
{
  unsigned int dof = numberOf1DBasefunctions(polynomialDegree);
  std::vector<double> F;
  F.reserve(dof);

  for (unsigned int i=0; i<dof; ++i)
    F.push_back(integradeOverRefEdge_gaus([&f, i, &B_k, &A_k, edgeID](double x_bar)
                                          {
                                            double x_hat, y_hat;
                                            // gama_e: refEdge(1D) -> refElement(2D)
                                            switch(edgeID)
                                              {
                                              case 0: // edge a
                                                x_hat = 1-x_bar;
                                                y_hat = x_bar;
                                                break;
                                              case 1: // edge b
                                                x_hat = 0;
                                                y_hat = 1-x_bar;
                                                break;
                                              case 2: // edge c
                                                x_hat = x_bar;
                                                y_hat = 0;
                                                break;
                                              default: // no more warnings
                                                x_hat =0;
                                                y_hat =0;
                                              }

                                            // F_k: refElemnt -> phyikal Elemnt
                                            double x1 = B_k[0][0]*x_hat + B_k[0][1]*y_hat + A_k.x;
                                            double x2 = B_k[1][0]*x_hat + B_k[1][1]*y_hat + A_k.y;

                                            // evaluate function
                                            return pol::phi1D[i](x_bar) * f(x1,x2);
                                          },
                                          2*polynomialDegree)

                );



  return F;
}

/**
 * trim polynomial to lower order
 * make use of herachical base function to project polynomial in lower sapce
 + f_h apoximate f
 * int{f_h \phi_i} = int{f \phi_i}  -> both f and f_h are polynomials in basis expansion (differnet order)
 * F_h int{B_h \phi_i} = F int{B \phi_i}  with i \in {1,..N_h} -> phi_i aus ziel raum
 * F_h M(=I)   =  F \hat{M}        -> M ist I (ortogonale basisvectoren)
 * and \hat{M} ist I_h mit 0 spalten hinten dran (= unten abscgnittenes "göseres" I)
 * -> F unten abschneiden
 */
inline std::vector<double> trim1D (std::vector<double> const & pol,
                                   unsigned int polynomialDegree)
{
  unsigned int dof (numberOf1DBasefunctions(polynomialDegree));

  assert(pol.size() >= dof);

  std::vector<double> F;
  F.reserve(dof);

  for (unsigned int i=0; i<dof; ++i)
    F.push_back(pol[i]);

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
