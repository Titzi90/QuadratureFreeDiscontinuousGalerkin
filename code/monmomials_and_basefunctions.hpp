#ifndef MOMOMIALS_AND_BASISFUNCRIONS_HPP
#define MOMOMIALS_AND_BASISFUNCRIONS_HPP

#include <cmath>
#include "polynomial.hpp"

namespace pol
{
  using std::sqrt;

  //monomials 1D
  Polynomial1D const c1D = Polynomial1D(0,1.);
  Polynomial1D const s   = Polynomial1D(1,1.);
  Polynomial1D const s2  = Polynomial1D(2,1.);
  Polynomial1D const s3  = Polynomial1D(3,1.);

  //monomials 2D
  // Order 0
  Polynomial2D const c2D = Polynomial2D(0,0,1.);
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

  // int const polynomialGad [] = {1,3,6,10};
  inline int numberOf2DBasefunctions(int polynomialDegree)
  {
    assert(polynomialDegree <= 3); // not implemented
    return (polynomialDegree+1)*(polynomialDegree+2)/2;
  }
  inline int numberOf1DBasefunctions(int polynomialDegree)
  {
    assert(polynomialDegree <= 2); // not implemented
    return polynomialDegree+1;
  }

  // 2D base Functions
  Polynomial2D const phi [] = {
    // Order 0
    sqrt(2)*c2D,
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
    c1D,
    // Order 1
    sqrt(3) * (1-2*s),
    // Order 2
    sqrt(5) * ( (6 * s - 6) * s + 1 )
  };
}  // end namespace pol

#endif
