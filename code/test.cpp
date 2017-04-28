#include "polynomial.hpp"
#include "monmomials_and_basefunctions.hpp"
#include "squareGrid.hpp"

#include <iostream>
#include <cmath>

template<typename T>
int test(T const & a, T const & b, std::string msg)
{
  double const EPSILON = 1e-9;

  if ( EPSILON < std::abs(a-b)){
    std::cerr << "Error!\n" << msg
              << "\n\tExpackted: " << a
              << "\n\tgot: " << b
              << std::endl;
    return -1;
  }
  return 0;
}

int test(Polynomial2D const & pol, std::string const & pol_str, std::string const & msg)
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
int test(Polynomial2D const & pol1, Polynomial2D const & pol2, std::string const & msg)
{
  return test(pol1, to_string(pol2), msg);
}

int test(Coordinate2D const & cor, double x, double y, std::string const & msg)
{
  return test<double>(cor.x, x, "x coordinate of "+msg) | test<double>(cor.y, y, "y coordinate of "+msg);
}

int main()
{
  int status = 0;

  // Test create 1
  Polynomial2D pol1 (2);
  pol1(2,0) = 3.2;
  pol1(1,1) = 2.2;
  pol1(1,0) = 1.2;
  pol1(0,2) = -2.2;
  pol1(0,1) = -1.2;
  pol1(0,0) = 5.2;
  std::string const pol1_str ("3.2x^2y^0 + 2.2x^1y^1 + 1.2x^1y^0 + -2.2x^0y^2 + -1.2x^0y^1 + 5.2x^0y^0");
  status = status | test(pol1, pol1_str, "creating complex 2D polynomial");

  // Test create 2
  Polynomial2D pol2 (2);
  pol2(2,0) = -2;
  pol2(1,1) = 3;
  std::string const pol2_str ("-2x^2y^0 + 3x^1y^1");
  status = status | test(pol2, pol2_str, "creating complex 2D polynomial");

  // Test create 2
  Polynomial2D pol3 (2,1,-2.5);
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

  // Test integrate 1
  double integral1 = integradeOverRefTriangle(pol::phi1);
  status = status | test(integral1, std::sqrt(2)/2, "integral 1");

  // Test integrate 2
  double integral2 = integradeOverRefTriangle(pol2);
  status = status | test(integral2, -1./24, "integral 2");

  // Test derive constant
  Polynomial2D pol11 = derive(Polynomial2D(0,0,5), Variable::X);
  status = status | test(pol11, "", "derive constant polynomial"); //TODO das geht nicht !!!!

  // // Test derive X
  Polynomial2D pol12 = derive(pol2, Variable::X);
  std::string const pol12_str ("-4x^1y^0 + 3x^0y^1");
  status = status | test(pol12, pol12_str, "simple derivative (dPol/dX)");

  // // Test derive Y (more complex then the ones before)
  Polynomial2D pol13 = derive(pol1, Variable::Y);
  std::string const pol13_str ("2.2x^1y^0 + -4.4x^0y^1 + -1.2x^0y^0");
  status = status | test(pol13, pol13_str, "more complex derivative (dPol/dY)");



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


  /****************************************************************************/

  Vector vec1 (1,0);
  Vector vec2 (0,1);
  Vector vec3 (3,5);
  status = status | test(vec3,3,5, "Vector constructor failed");

  // test Vector - Vector
  Vector vec4 = vec1-vec2;
  status = status | test(vec4, 1., -1., "Vector - operator failed");

  // test scale operator
  Vector vec5 = 3*vec4;
  status = status | test(vec5, 3.,-3., "Vector scaling failed");

  // test length function
  double l = length(vec5);
  status |= test(4.242640687119285,l, "vector length failed");

  // test nomalizing
  Vector vec6 = normalize(vec5);
  status |= test(1.,length(vec6), "nomalizing of vectr failed");
  status |= test(vec6,0.7071067811865476,-0.7071067811865476,"nomalizing of vectr failed");

  // test distance
  double d = distance(vec5, vec3);
  status |= test (8., d, "calculating distance failed");

  // test dot product
  double dot1 = dot(vec1, vec2);
  status |= test(0., dot1, "dot product of unit vectors failed");

  double dot2 = dot(vec3, vec5);
  status |= test(-6.,dot2, "more complex dot product failed");

  // test getNormal vector
  Vector vec7 = getNormal(vec1, vec2);
  status |= test(1., length(vec7), "get normale vector failed (not normalized)");
  status |= test(0., dot(vec7,vec1-vec2), "get normale vector failed  (not perpendicular)");
  // status |= test(vec7, 0.7071067811865475,0.7071067811865475, "gte normalized vector failed"); // vec7*(-1) is also ok

  // Test construct a triangle
  Triangle t1 (vec1, vec2, Vector(0,0));
  status |= test(t1.getPointA(), 1, 0, "triangle inizailisation Point A failed");
  status |= test(t1.getPointB(), 0, 1, "triangle inizailisation Point B failed");
  status |= test(t1.getPointC(), 0, 0, "triangle initialization Point C failed");

  status |= test(1., t1.getLengthA(), "triangle initialization Length a failed");
  status |= test(1., t1.getLengthB(), "triangle initialization Length b failed");
  status |= test(1.4142135623730951, t1.getLengthC(), "triangle initialization Length c failed");

  status |= test(t1.getNormalA(), -1,0, "triangle initialization nomal a failed");
  status |= test(t1.getNormalB(), 0,-1, "triangle initialization nomal b failed");
  status |= test(t1.getNormalC(), 0.7071067811865476,0.7071067811865476, "triangle initialization nomal c failed");

  status |= test(0.5, t1.getArea(), "triangle initialization of area failed");




  return status;
}
