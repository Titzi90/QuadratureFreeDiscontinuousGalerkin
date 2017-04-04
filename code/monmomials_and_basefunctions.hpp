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

  //TODO testen
  // Base Functions
  // Order 0
  Polynomial2D const phi1  = sqrt(2)*c;
  // Order 1
  Polynomial2D const phi2  = -6*x + 2;
  Polynomial2D const phi3  = (-2*sqrt(3))*(x + 2*y - 1);
  // Order 2
  Polynomial2D const phi4  = sqrt(6)*(x * ( 10*x - 8) + 1);
  Polynomial2D const phi5  = sqrt(3)*(x * (5*x - 4) + y*(-15*y +12) -1);
  Polynomial2D const phi6  = (3*sqrt(5))*(x*(3*x + 8*y - 4) + y*(3*y - 4) +1);
  // Order 3
  Polynomial2D const phi7  = (2*sqrt(2))*(x*(x*(35*x - 45) + 15) - 1);
  Polynomial2D const phi8  = (2*sqrt(6))*(x*(x*(21*x -33) + 13) + y*(x*(42*x -24) + 2) -1);
  Polynomial2D const phi9  = (2*sqrt(10))*(x*(x*(7*x - 15) + 9) + y*(x*(42*x - 48) + y*(42*x - 6) +6) - 1);
  Polynomial2D const phi10 = (2*sqrt(14))*(x*(x*(x - 3) + 3) + y*(x*(12*x - 24) +y*(30*x + 20*y -30) + 12) -1);

}  // end namespace pol

#endif
