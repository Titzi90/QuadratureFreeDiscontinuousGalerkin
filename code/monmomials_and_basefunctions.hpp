#ifndef MOMOMIALS_AND_BASISFUNCRIONS_HPP
#define MOMOMIALS_AND_BASISFUNCRIONS_HPP

#include <cmath>
#include "polynomial.hpp"

namespace pol
{
  using std::sqrt;

  //monomials
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


  int const polynomialGad [] = {1,3,6,10};

  // Base Functions
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

  // X derivatives of base functions
  Polynomial2D const dXphi [] = {
    // Order 0
    derive( sqrt(2)*c, Variable::X ),
    // Order 1
    derive( -6*x + 2, Variable::X),
    derive( (-2*sqrt(3))*(x + 2*y - 1), Variable::X),
    // Order 2
    derive( sqrt(6)*(x * ( 10*x - 8) + 1), Variable::X ),
    derive( sqrt(3)*(x * (5*x - 4) + y*(-15*y +12) -1), Variable::X ),
    derive( (3*sqrt(5))*(x*(3*x + 8*y - 4) + y*(3*y - 4) +1), Variable::X ),
    // Order 3
    derive( (2*sqrt(2))*(x*(x*(35*x - 45) + 15) - 1), Variable::X ),
    derive( (2*sqrt(6))*(x*(x*(21*x -33) + 13) + y*(x*(42*x -24) + 2) -1), Variable::X ),
    derive( (2*sqrt(10))*(x*(x*(7*x - 15) + 9) + y*(x*(42*x - 48) + y*(42*x - 6) +6) - 1), Variable::X ),
    derive( (2*sqrt(14))*(x*(x*(x - 3) + 3) + y*(x*(12*x - 24) +y*(30*x + 20*y -30) + 12) -1), Variable::X ),
  };

  // Y derivatives of base functions
  Polynomial2D const dYphi [] = {
    // Order 0
    derive( sqrt(2)*c, Variable::Y ),
    // Order 1
    derive( -6*x + 2, Variable::Y),
    derive( (-2*sqrt(3))*(x + 2*y - 1), Variable::Y),
    // Order 2
    derive( sqrt(6)*(x * ( 10*x - 8) + 1), Variable::Y ),
    derive( sqrt(3)*(x * (5*x - 4) + y*(-15*y +12) -1), Variable::Y ),
    derive( (3*sqrt(5))*(x*(3*x + 8*y - 4) + y*(3*y - 4) +1), Variable::Y ),
    // Order 3
    derive( (2*sqrt(2))*(x*(x*(35*x - 45) + 15) - 1), Variable::Y ),
    derive( (2*sqrt(6))*(x*(x*(21*x -33) + 13) + y*(x*(42*x -24) + 2) -1), Variable::X ),
    derive( (2*sqrt(10))*(x*(x*(7*x - 15) + 9) + y*(x*(42*x - 48) + y*(42*x - 6) +6) - 1), Variable::Y ),
    derive( (2*sqrt(14))*(x*(x*(x - 3) + 3) + y*(x*(12*x - 24) +y*(30*x + 20*y -30) + 12) -1), Variable::Y )
  };
}  // end namespace pol

#endif
