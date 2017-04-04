#include "polynomial.hpp"
#include "monmomials_and_basefunctions.hpp"

#include <iostream>

template<typename T>
int test(T a, T b, std::string msg)
{
  if ( a != b){
    std::cerr << "Error!\n" << msg
              << "\n\tExpackted: " << a
              << "\n\tgot: " << b
              << std::endl;
    return -1;
  }
  return 0;
}

int test(Polynomial2D pol, std::string pol_str, std::string msg)
{
  if ( pol_str != to_string(pol)){
    std::cerr << "Error!\n" << msg
              << "\n\tExpackted: " << pol_str
              << "\n\tgot: " << to_string(pol)
              << std::endl;
    return -1;
  }
  return 0;
}
int test(Polynomial2D pol1, Polynomial2D pol2, std::string msg){ return test(pol1, to_string(pol2), msg); }


int main()
{
  int status = 0;

  // Test create 1
  Polynomial2D pol1 = Polynomial2D(2);
  pol1(2,0) = 3.2;
  pol1(1,1) = 2.2;
  pol1(1,0) = 1.2;
  pol1(0,2) = -2.2;
  pol1(0,1) = -1.2;
  pol1(0,0) = 5.2;
  std::string const pol1_str ("3.2x^2y^0 + 2.2x^1y^1 + 1.2x^1y^0 + -2.2x^0y^2 + -1.2x^0y^1 + 5.2x^0y^0");
  status = status | test(pol1, pol1_str, "creating complex 2D polynomial");

  // Test create 2
  Polynomial2D pol2 = Polynomial2D(2);
  pol2(2,0) = -2;
  pol2(1,1) = 3;
  std::string const pol2_str ("-2x^2y^0 + 3x^1y^1");
  status = status | test(pol2, pol2_str, "creating complex 2D polynomial");

  // Test create 2
  Polynomial2D pol3 = Polynomial2D(2,1,-2.5);
  std::string const pol3_str ("-2.5x^2y^1");
  status = status | test(pol3, pol3_str, "creating simple 2D polynomial");

  // Test plus
  Polynomial2D pol4 = pol1 + pol2;
  std::string const pol4_str ("1.2x^2y^0 + 5.2x^1y^1 + 1.2x^1y^0 + -2.2x^0y^2 + -1.2x^0y^1 + 5.2x^0y^0");
  status = status | test(pol4, pol4_str, "add two polynomials");

  // Test plus of different orders
  Polynomial2D pol5 = pol2 + pol3;
  std::string const pol5_str ("-2.5x^2y^1 + -2x^2y^0 + 3x^1y^1");
  status = status | test(pol5, pol5_str, "add two polynomials of different order");

  // Test plus real
  Polynomial2D pol6 = pol5 + 3.7;
  std::string const pol6_str ("-2.5x^2y^1 + -2x^2y^0 + 3x^1y^1 + 3.7x^0y^0");
  status = status | test(pol6, pol6_str, "add real number to polynomial");

  // Test ...
  // Polynomial2D pol7 =

  // Test Multiplication
  Polynomial2D pol8 = pol1 * pol2;
  std::string const pol8_str ("-6.4x^4y^0 + 5.2x^3y^1 + -2.4x^3y^0 + 11x^2y^2 + 6x^2y^1 + -10.4x^2y^0 + -6.6x^1y^3 + -3.6x^1y^2 + 15.6x^1y^1");
  status = status | test(pol8, pol8_str, "multiply two polynomials");

  // Test scale
  Polynomial2D pol9 = pol1 * (-1);
  std::string const pol9_str ("-3.2x^2y^0 + -2.2x^1y^1 + -1.2x^1y^0 + 2.2x^0y^2 + 1.2x^0y^1 + -5.2x^0y^0");
  status = status | test(pol9, pol9_str, "multiply real with polynomial (scale)");

  // Test minus
  Polynomial2D pol10 = pol2 - pol5;
  status = status | test(pol10, -1*pol3, "subtract two polynomials");

  // Test integrate
  // double integral1 = integradeOverRefTriangle(pol::c);
  // status = status | test(integral1, std::sqrt(2)/2, "integral 1");

  /****************************************************************************/

  status = status | test(pol::c  , "1x^0y^0", "constant monomial");
  status = status | test(pol::x  , "1x^1y^0", "x monomial");
  status = status | test(pol::y  , "1x^0y^1", "y monomial");
  status = status | test(pol::x2 , "1x^2y^0", "x^2 monomial");
  status = status | test(pol::xy , "1x^1y^1", "xy monomial");
  status = status | test(pol::y2 , "1x^0y^2", "y^2 monomial");
  status = status | test(pol::x3 , "1x^3y^0", "x^3 monomial");
  status = status | test(pol::x2y, "1x^2y^1", "x^2y monomial");
  status = status | test(pol::xy2, "1x^1y^2", "xy^2 monomial");
  status = status | test(pol::y3 , "1x^0y^3", "y^3 monomial");

  return status;
}
