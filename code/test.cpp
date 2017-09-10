#include "polynomial.hpp"
#include "monmomials_and_basefunctions.hpp"
#include "Grid.hpp"
#include "assembly.hpp"
#include "DataTypes.hpp"

#include <iostream>
#include <cmath>
#include <vector>

template<typename T>
int test (T const & a, T const & b, std::string msg)
{
  double const EPSILON = 1e-3;

  if ( EPSILON < std::abs(a-b))
    {
      std::cerr << "Error!\n" << msg
                << "\n\tExpackted: " << b
                << "\n\tgot: " << a
                << std::endl;
      return -1;
    }
  return 0;
}

int test (std::string const & str1, std::string const & str2, std::string const & msg)
{
  if ( str1 != str2)
    {
      std::cerr << "Error!\n" << msg
                << "\n\tExpackted: " << str2
                << "\n\tgot: " << str1
                << std::endl;
      return -1;
    }
  return 0;
}

int test (Coordinates2D const & cor, double x, double y, std::string const & msg)
{
  return test<double>(cor.x, x, "x coordinate of "+msg) | test<double>(cor.y, y, "y coordinate of "+msg);
}







int polynomialTest()
{
  int status = 0;

  // Test creating polynomials
  Polynomial1D pol1Dc (0, 3);
  Polynomial1D pol1Ds (1, 5);
  Polynomial1D pol1D0 (1);

  std::string const pol1Dc_str ("3x^0");
  std::string const pol1Ds_str ("5x^1");
  std::string const pol1D0_str ("");
  status |= test(to_string(pol1Dc), pol1Dc_str, "construct simple 1D polynomial");
  status |= test(to_string(pol1Ds), pol1Ds_str, "construct simple 1D polynomial");
  status |= test(to_string(pol1D0), pol1D0_str, "construct null 1D polynomial");

  Polynomial2D pol1 (2);
  pol1.get(2,0) = 3.2;
  pol1.get(1,1) = 2.2;
  pol1.get(1,0) = 1.2;
  pol1.get(0,2) = -2.2;
  pol1.get(0,1) = -1.2;
  pol1.get(0,0) = 5.2;
  std::string const pol1_str ("3.2x^2y^0 + 2.2x^1y^1 + 1.2x^1y^0 + -2.2x^0y^2 + -1.2x^0y^1 + 5.2x^0y^0");
  status |= test(to_string(pol1), pol1_str, "creating complex 2D polynomial");

  Polynomial2D pol2 (2);
  pol2.get(2,0) = -2;
  pol2.get(1,1) = 3;
  std::string const pol2_str ("-2x^2y^0 + 3x^1y^1");
  status |= test(to_string(pol2), pol2_str, "creating complex 2D polynomial");

  Polynomial2D pol3 (2,1,-2.5);
  std::string const pol3_str ("-2.5x^2y^1");
  status |= test(to_string(pol3), pol3_str, "creating simple 2D polynomial");

  // Test evaluation
  double eval1 = pol3(2,3);
  status |= test(eval1, -30., "evaluation of simple polynomial");

  double eval2 = pol2(3,-1.5);
  status |= test(eval2, -31.5, "evaluation of more complex polynomial");

  // Test plus operator
  Polynomial2D pol4 = pol1 + pol2;
  std::string const pol4_str ("1.2x^2y^0 + 5.2x^1y^1 + 1.2x^1y^0 + -2.2x^0y^2 + -1.2x^0y^1 + 5.2x^0y^0");
  status |= test(to_string(pol4), pol4_str, "add two polynomials");

  Polynomial2D pol5 = pol2 + pol3;
  std::string const pol5_str ("-2.5x^2y^1 + -2x^2y^0 + 3x^1y^1");
  status |= test(to_string(pol5), pol5_str, "add two polynomials of different order");

  Polynomial2D pol6 = pol5 + 3.7;
  std::string const pol6_str ("-2.5x^2y^1 + -2x^2y^0 + 3x^1y^1 + 3.7x^0y^0");
  status |= test(to_string(pol6), pol6_str, "add real number to polynomial");

  Polynomial1D pol1D_1 = pol1Dc + pol1Ds;
  std::string const pol1D_1_str ("5x^1 + 3x^0");
  status |= test(to_string(pol1D_1), pol1D_1_str, "add two simple 1D polynomials");

  // Test assignment-add
  pol1D_1 += pol1Dc;
  std::string pol1D_2_str ("5x^1 + 6x^0");
  status |= test(to_string(pol1D_1), pol1D_2_str, "assigment-add two simple 1D polynomials");

  pol1D0 += pol1D_1;
  status |= test(to_string(pol1D0), pol1D_2_str, "assigment-add two simple 1D polynomials to zero pol");

  // Test Multiplication
  Polynomial2D pol8 = pol1 * pol2;
  std::string const pol8_str ("-6.4x^4y^0 + 5.2x^3y^1 + -2.4x^3y^0 + 11x^2y^2 + 6x^2y^1 + -10.4x^2y^0 + -6.6x^1y^3 + -3.6x^1y^2 + 15.6x^1y^1");
  status |= test(to_string(pol8), pol8_str, "multiply two polynomials");

  // Test scale
  Polynomial2D pol9 = pol1 * (-1);
  std::string const pol9_str ("-3.2x^2y^0 + -2.2x^1y^1 + -1.2x^1y^0 + 2.2x^0y^2 + 1.2x^0y^1 + -5.2x^0y^0");
  status |= test(to_string(pol9), pol9_str, "multiply real with polynomial (scale)");

  // Test minus
  Polynomial2D pol10 = pol2 - pol5;
  status |= test(to_string(pol10), to_string(-1*pol3), "subtract two polynomials");

  // Test integrate
  double integral1 = integradeOverRefTriangle(pol::phi[0]);
  status |= test(integral1, std::sqrt(2)/2, "integral 1");

  double integral2 = integradeOverRefTriangle(pol2);
  status |= test(integral2, -1./24, "integral 2");

  // Test derive constant
  Polynomial2D pol11 = derive(Polynomial2D(0,0,5), Variable::X);
  status |= test(to_string(pol11), "", "derive constant polynomial");

  Polynomial2D pol12 = derive(pol2, Variable::X);
  std::string const pol12_str ("-4x^1y^0 + 3x^0y^1");
  status |= test(to_string(pol12), pol12_str, "simple derivative (dPol/dX)");

  Polynomial2D pol13 = derive(pol1, Variable::Y);
  std::string const pol13_str ("2.2x^1y^0 + -4.4x^0y^1 + -1.2x^0y^0");
  status |= test(to_string(pol13), pol13_str, "more complex derivative (dPol/dY)");

  // TEST monomials
  status |= test(to_string(pol::c  ), "1x^0y^0", "constant monomial");
  status |= test(to_string(pol::x  ), "1x^1y^0", "x monomial");
  status |= test(to_string(pol::y  ), "1x^0y^1", "y monomial");
  status |= test(to_string(pol::x2 ), "1x^2y^0", "x^2 monomial");
  status |= test(to_string(pol::xy ), "1x^1y^1", "xy monomial");
  status |= test(to_string(pol::y2 ), "1x^0y^2", "y^2 monomial");
  status |= test(to_string(pol::x3 ), "1x^3y^0", "x^3 monomial");
  status |= test(to_string(pol::x2y), "1x^2y^1", "x^2y monomial");
  status |= test(to_string(pol::xy2), "1x^1y^2", "xy^2 monomial");
  status |= test(to_string(pol::y3 ), "1x^0y^3", "y^3 monomial");

  return status;
}

int vectorTest()
{
  int status = 0;

  Vector vec1 (1,0);
  Vector vec2 (0,1);
  Vector vec3 (3,5);
  status |= test(vec3,3,5, "Vector constructor failed");

  // test Vector - Vector
  Vector vec4 = vec1-vec2;
  status |= test(vec4, 1., -1., "Vector - operator failed");

  // test scale operator
  Vector vec5 = 3*vec4;
  status |= test(vec5, 3.,-3., "Vector scaling failed");

  // test length function
  double l = length(vec5);
  status |= test(l, 4.242640687119285, "vector length failed");

  // test nomalizing
  Vector vec6 = normalize(vec5);
  status |= test(length(vec6), 1., "nomalizing of vectr failed");
  status |= test(vec6,0.7071067811865476,-0.7071067811865476,"nomalizing of vectr failed");

  // test distance
  double d = distance(vec5, vec3);
  status |= test (d, 8., "calculating distance failed");

  // test dot product
  double dot1 = dot(vec1, vec2);
  status |= test(dot1, 0., "dot product of unit vectors failed");

  double dot2 = dot(vec3, vec5);
  status |= test(dot2, -6., "more complex dot product failed");

  // test getNormal vector
  Vector vec7 = getNormal(vec1, vec2);
  status |= test(length(vec7), 1., "get normale vector failed (not normalized)");
  status |= test(dot(vec7,vec1-vec2), 0., "get normale vector failed  (not perpendicular)");
  // status |= test(vec7, 0.7071067811865475,0.7071067811865475, "gte normalized vector failed"); // vec7*(-1) is also ok

  return status;
}

int triangleMeshTest()
{
  int status = 0;

  // Test construct a triangle
  Triangle t1 (Vector(1,0), Vector(0,1), Vector(0,0));
  status |= test(t1.getA(), 1, 0, "triangle inizailisation Point A failed");

  status |= test(t1.getLengthA(), 1., "triangle initialization Length a failed");
  status |= test(t1.getLengthB(), 1., "triangle initialization Length b failed");
  status |= test(t1.getLengthC(), 1.4142135623730951, "triangle initialization Length c failed");

  status |= test(t1.getArea(), 0.5, "triangle initialization of area failed");

  status |= test(t1.getNormalA(), -1,0, "triangle initialization nomal a failed");
  status |= test(t1.getNormalB(), 0,-1, "triangle initialization nomal b failed");
  status |= test(t1.getNormalC(), 0.7071067811865476,0.7071067811865476, "triangle initialization nomal c failed");


  // Test mesh class
  UniqueSquareGrid grid1(4);
  Triangle const & t2 = grid1.getLower(0,0);
  Triangle const & t3 = grid1.getUpper(0,0);
  Triangle const & t4 = grid1.getLower(1,2);
  Triangle const & t5 = grid1.getUpper(3,3);

  status |= test(t2.getA(),0.,0., "first Triangle in mesh is incorrect (Point A)");
  status |= test(t2.getArea(), 0.03125, "Area of first triangle is wrong");

  status |= test(t3.getA(),0.25, 0., "second Triangle in mesh is incorrect (Point A)");

  status |= test(t4.getA(),0.5,0.25, "Triangle in mesh is incorrect (Point A)");

  status |= test(t5.getA(),1.,0.75, "last Triangle in mesh is incorrect (Point A)");


  // test Jakobian
  Jakobian const & b1 ( t2.getJakobian() );
  status |= test(b1[0][0], 0.25, "Jakobia[1][1] of first (lower) triangle is wrong");
  status |= test(b1[1][0], 0.  , "Jakobia[2][1] of first (lower) triangle is wrong");
  status |= test(b1[0][1], 0.  , "Jakobia[1][2] of first (lower) triangle is wrong");
  status |= test(b1[1][1], 0.25, "Jakobia[2][2] of first (lower) triangle is wrong");

  Jakobian const & b2 ( t5.getJakobian() );
  status |= test(b2[0][0],  0.  , "Jakobia[1][1] of last (upper) triangle is wrong");
  status |= test(b2[1][0],  0.25, "Jakobia[2][1] of last (upper) triangle is wrong");
  status |= test(b2[0][1], -0.25, "Jakobia[1][2] of last (upper) triangle is wrong");
  status |= test(b2[1][1],  0.25, "Jakobia[2][2] of last (upper) triangle is wrong");

  return status;
}

int assamblyTest()
{
  int status=0;

  int const order = 1;     // max polynomial degree
  int const refinment = 2; // order if mesh refinement

  UniqueSquareGrid mesh(refinment);

  auto u1 = [](double , double )->double{ return 1.; };
  auto u2 = [](double , double )->double{ return 1.; };
  auto f  = [](double , double )->double{ return 0.; };
  auto c  = [](double , double )->double{ return 1.; };

  assamblyC(mesh, order, c);          // gesuchte Größe c
  assamblyU(mesh, order, u1, u2);     // velocity field U
  assamblyF(mesh, order, 2*order, assamblyLocalLinearF);  // Flux field f=cu
  setBoundary_Periodic(mesh, Boundary::bottom);
  setBoundary_Periodic(mesh, Boundary::top);
  setBoundary_Periodic(mesh, Boundary::left);
  setBoundary_Periodic(mesh, Boundary::right);

  auto hatM (assemblyHatM(order));
  auto hatG (assemblyHatG(order));
  auto hatE (assemblyHatE(order, order*2));
  auto hatI (getHatI(2*order));

  /*
  assemblyM(mesh, hatM);
  assemblyG(mesh, hatG);
  assemblyE(mesh, hatE);
  assemblyFr(mesh, order, order*2, riemanSolver_UpWinding, hatI);
  assamblyL(mesh, order, f);          // RHS vector
  */
  /*
  assemblyMquadFree(mesh, order);
  assemblyGquadFree(mesh, order);
  assemblyEquadFree(mesh, order, order*2);
  assemblyFr(mesh, order, order*2, riemanSolver_UpWinding, hatI);
  assamblyL(mesh, order, f);          // RHS vector
  */
  assemblyMGaus(mesh, order);
  assemblyGgaus(mesh, order);
  assemblyE(mesh, hatE);
  assemblyFr(mesh, order, order*2, riemanSolver_UpWinding, hatI);
  assamblyL(mesh, order, f);          // RHS vector

  auto t1 (mesh.getLower(0, 0));
  auto t2 (mesh.getUpper(1, 1));
  auto t3 (mesh.getLower(1, 1));
  auto t4 (mesh.getUpper(0, 0));

  //TODO hatE

  // test hatI
  auto hatI_test (getHatI(2));
  status |= test(hatI_test[0], 1., "hatI[0]");
  status |= test(hatI_test[1],-1., "hatI[0]");
  status |= test(hatI_test[2], 1., "hatI[0]");
  status |= test(hatI_test.size()*1., 3., "size of hatI (order 2)");

  //test Mhat
  status |= test(1., hatM(0,0), "hatM(0,0)");
  status |= test(1., hatM(1,1), "hatM(1,1)");
  status |= test(1., hatM(2,2), "hatM(2,2)");
  status |= test(0., hatM(0,1), "hatM(0,1)");
  status |= test(0., hatM(0,2), "hatM(0,2)");
  status |= test(0., hatM(1,0), "hatM(1,0)");
  status |= test(0., hatM(1,2), "hatM(1,2)");
  status |= test(0., hatM(2,0), "hatM(2,0)");
  status |= test(0., hatM(2,1), "hatM(2,1)");

  // test hatG
  status |= test(0., hatG[0](0,0,0), "hatG1(0,0,0)");
  status |= test(0., hatG[0](0,1,0), "hatG1(0,1,0)");
  status |= test(0., hatG[0](0,2,0), "hatG1(0,2,0)");
  status |= test(0., hatG[0](0,0,1), "hatG1(0,0,1)");
  status |= test(0., hatG[0](0,1,1), "hatG1(0,1,1)");
  status |= test(0., hatG[0](0,2,1), "hatG1(0,2,1)");
  status |= test(0., hatG[0](0,0,2), "hatG1(0,0,2)");
  status |= test(0., hatG[0](0,1,2), "hatG1(0,1,2)");
  status |= test(0., hatG[0](0,2,2), "hatG1(0,2,2)");

  status |= test(-6., hatG[0](1,0,0), "hatG1(1,0,0)");
  status |= test(0., hatG[0](1,1,0), "hatG1(1,1,0)");
  status |= test(0., hatG[0](1,2,0), "hatG1(1,2,0)");
  status |= test(0., hatG[0](1,0,1), "hatG1(1,0,1)");
  status |= test(-6., hatG[0](1,1,1), "hatG1(1,1,1)");
  status |= test(0., hatG[0](1,2,1), "hatG1(1,2,1)");
  status |= test(0., hatG[0](1,0,2), "hatG1(1,0,2)");
  status |= test(0., hatG[0](1,1,2), "hatG1(1,1,2)");
  status |= test(-6., hatG[0](1,2,2), "hatG1(1,2,2)");

  status |= test(-3.4641, hatG[0](2,0,0), "hatG1(2,0,0)");
  status |= test(0., hatG[0](2,1,0), "hatG1(2,1,0)");
  status |= test(0., hatG[0](2,2,0), "hatG1(2,2,0)");
  status |= test(0., hatG[0](2,0,1), "hatG1(2,0,1)");
  status |= test(-3.4641, hatG[0](2,1,1), "hatG1(2,1,1)");
  status |= test(0., hatG[0](2,2,1), "hatG1(2,2,1)");
  status |= test(0., hatG[0](2,0,2), "hatG1(2,0,2)");
  status |= test(0., hatG[0](2,1,2), "hatG1(2,1,2)");
  status |= test(-3.4641, hatG[0](2,2,2), "hatG1(2,2,2)");

  status |= test(0., hatG[1](0,0,0), "hatG2(0,0,0)");
  status |= test(0., hatG[1](0,1,0), "hatG2(0,1,0)");
  status |= test(0., hatG[1](0,2,0), "hatG2(0,2,0)");
  status |= test(0., hatG[1](0,0,1), "hatG2(0,0,1)");
  status |= test(0., hatG[1](0,1,1), "hatG2(0,1,1)");
  status |= test(0., hatG[1](0,2,1), "hatG2(0,2,1)");
  status |= test(0., hatG[1](0,0,2), "hatG2(0,0,2)");
  status |= test(0., hatG[1](0,1,2), "hatG2(0,1,2)");
  status |= test(0., hatG[1](0,2,2), "hatG2(0,2,2)");

  status |= test(0., hatG[1](1,0,0), "hatG2(1,0,0)");
  status |= test(0., hatG[1](1,1,0), "hatG2(1,1,0)");
  status |= test(0., hatG[1](1,2,0), "hatG2(1,2,0)");
  status |= test(0., hatG[1](1,0,1), "hatG2(1,0,1)");
  status |= test(0., hatG[1](1,1,1), "hatG2(1,1,1)");
  status |= test(0., hatG[1](1,2,1), "hatG2(1,2,1)");
  status |= test(0., hatG[1](1,0,2), "hatG2(1,0,2)");
  status |= test(0., hatG[1](1,1,2), "hatG2(1,1,2)");
  status |= test(0., hatG[1](1,2,2), "hatG2(1,2,2)");

  status |= test(-6.9282, hatG[1](2,0,0), "hatG2(2,0,0)");
  status |= test(0., hatG[1](2,1,0), "hatG2(2,1,0)");
  status |= test(0., hatG[1](2,2,0), "hatG2(2,2,0)");
  status |= test(0., hatG[1](2,0,1), "hatG2(2,0,1)");
  status |= test(-6.9282, hatG[1](2,1,1), "hatG2(2,1,1)");
  status |= test(0., hatG[1](2,2,1), "hatG2(2,2,1)");
  status |= test(0., hatG[1](2,0,2), "hatG2(2,0,2)");
  status |= test(0., hatG[1](2,1,2), "hatG2(2,1,2)");
  status |= test(-6.9282, hatG[1](2,2,2), "hatG2(2,2,2)");

  // test global M
  status |= test(t1.M()(1,1), 0.25, "M(1,1) @ lower(0,0)");
  status |= test(t2.M()(0,0), 0.25, "M(0,0) @ upper(1,1)");
  status |= test(t3.M()(2,2), 0.25, "M(2,2) @ lower(1,1)");
  status |= test(t1.M()(0,1), 0.  , "M(0,1) @ lower(0,0)");
  status |= test(t2.M()(1,0), 0.  , "M(1,0) @ upper(1,1)");
  status |= test(t3.M()(1,2), 0.  , "M(1,2) @ lower(1,1)");


  // test global G
  status |= test(t1.G()(0,0), 0., "G(0,0) @ lower(0,0)");
  status |= test(t1.G()(1,0), -2.1213, "G(1,0) @ lower(0,0)");
  status |= test(t1.G()(2,0), -3.6742, "G(2,0) @ lower(0,0)");
  status |= test(t1.G()(0,1), 0., "G(0,1) @ lower(0,0)");
  status |= test(t1.G()(1,1), 0., "G(1,1) @ lower(0,0)");
  status |= test(t1.G()(2,1), 0., "G(2,1) @ lower(0,0)");
  status |= test(t1.G()(0,0), 0., "G(0,0) @ lower(0,0)");
  status |= test(t1.G()(1,2), 0., "G(1,2) @ lower(0,0)");
  status |= test(t1.G()(2,2), 0., "G(2,2) @ lower(0,0)");

  status |= test(t3.G()(0,0), 0., "G(0,0) @ lower(1,1)");
  status |= test(t3.G()(1,0), -2.1213, "G(1,0) @ lower(1,1)");
  status |= test(t3.G()(2,0), -3.6742, "G(2,0) @ lower(1,1)");
  status |= test(t3.G()(0,1), 0., "G(0,1) @ lower(1,1)");
  status |= test(t3.G()(1,1), 0., "G(1,1) @ lower(1,1)");
  status |= test(t3.G()(2,1), 0., "G(2,1) @ lower(1,1)");
  status |= test(t3.G()(0,0), 0., "G(0,0) @ lower(1,1)");
  status |= test(t3.G()(1,2), 0., "G(1,2) @ lower(1,1)");
  status |= test(t3.G()(2,2), 0., "G(2,2) @ lower(1,1)");

  status |= test(t2.G()(0,0), 0., "G(0,0) @ upper(1,1)");
  status |= test(t2.G()(1,0), -4.2426, "G(1,0) @ upper(1,1)");
  status |= test(t2.G()(2,0), 0., "G(2,0) @ upper(1,1)");
  status |= test(t2.G()(0,1), 0., "G(0,1) @ upper(1,1)");
  status |= test(t2.G()(1,1), 0., "G(1,1) @ upper(1,1)");
  status |= test(t2.G()(2,1), 0., "G(2,1) @ upper(1,1)");
  status |= test(t2.G()(0,0), 0., "G(0,0) @ upper(1,1)");
  status |= test(t2.G()(1,2), 0., "G(1,2) @ upper(1,1)");
  status |= test(t2.G()(2,2), 0., "G(2,2) @ upper(1,1)");

  //test C (=1)
  status |= test(t1.C()[0], 0.7071, "C[0] @ lower(0,0)");
  status |= test(t1.C()[1], 0., "C[1] @ lower(0,0)");
  status |= test(t1.C()[2], 0., "C[2] @ lower(0,0)");

  status |= test(t2.C()[0], 0.7071, "C[0] @ upper(1,1)");
  status |= test(t2.C()[1], 0., "C[1] @ upper(1,1)");
  status |= test(t2.C()[2], 0., "C[2] @ upper(1,1)");

  //test U[1,2] (=1)
  status |= test(t1.U1()[0], 0.7071, "C[0] @ lower(0,0)");
  status |= test(t1.U1()[1], 0., "C[1] @ lower(0,0)");
  status |= test(t1.U1()[2], 0., "C[2] @ lower(0,0)");

  status |= test(t2.U2()[0], 0.7071, "C[0] @ upper(1,1)");
  status |= test(t2.U2()[1], 0., "C[1] @ upper(1,1)");
  status |= test(t2.U2()[2], 0., "C[2] @ upper(1,1)");

  //test F[1,2] (=1)
  status |= test(t1.F1()[0], 0.7071, "C[0] @ lower(0,0)");
  status |= test(t1.F1()[1], 0., "C[1] @ lower(0,0)");
  status |= test(t1.F1()[2], 0., "C[2] @ lower(0,0)");
  status |= test(t1.F1()[3], 0., "C[3] @ lower(0,0)");
  status |= test(t1.F1()[4], 0., "C[4] @ lower(0,0)");
  status |= test(t1.F1()[5], 0., "C[5] @ lower(0,0)");

  status |= test(t2.F2()[0], 0.7071, "C[0] @ upper(1,1)");
  status |= test(t2.F2()[1], 0., "C[1] @ upper(1,1)");
  status |= test(t2.F2()[2], 0., "C[2] @ upper(1,1)");
  status |= test(t2.F2()[3], 0., "C[3] @ upper(1,1)");
  status |= test(t2.F2()[4], 0., "C[4] @ upper(1,1)");
  status |= test(t2.F2()[5], 0., "C[5] @ upper(1,1)");

  // test L(=rhs= 0)
  status |= test(t1.L()[0], 0., "C[0] @ lower(0,0)");
  status |= test(t1.L()[1], 0., "C[1] @ lower(0,0)");
  status |= test(t1.L()[2], 0., "C[2] @ lower(0,0)");

  status |= test(t2.L()[0], 0., "C[0] @ upper(1,1)");
  status |= test(t2.L()[1], 0., "C[1] @ upper(1,1)");
  status |= test(t2.L()[2], 0., "C[2] @ upper(1,1)");

  // test F at the edges (=Fr)
  status |= test(t1.F_a()[0], 1.4142, "F_a[0] @ lower(0,0)");
  status |= test(t1.F_a()[1], 0., "F_a[1] @ lower(0,0)");

  status |= test(t1.F_b()[0], -1., "F_b[0] @ lower(0,0)");
  status |= test(t1.F_b()[1], 0., "F_b[1] @ lower(0,0)");

  status |= test(t1.F_c()[0], -1., "F_c[0] @ lower(0,0)");
  status |= test(t1.F_c()[1], 0., "F_c[1] @ lower(0,0)");

  status |= test(t2.F_a()[0], 1., "F_a[0] @ upper(1,1)");
  status |= test(t2.F_a()[1], 0., "F_a[1] @ upper(1,1)");

  status |= test(t2.F_b()[0], -1.4142, "F_b[0] @ upper(1,1)");
  status |= test(t2.F_b()[1], 0., "F_b[1] @ upper(1,1)");

  status |= test(t2.F_c()[0], 1., "F_c[0] @ upper(1,1)");
  status |= test(t2.F_c()[1], 0., "F_c[1] @ upper(1,1)");

  //TODO E


  // test dC = M^{-1}*(L + GC - sum_e{E_e*F_e})
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        auto & T_l = mesh.getLower(row, col);
        auto & T_u = mesh.getUpper(row, col);

        //TODO durch M teilen
        auto dc_l = T_l.M()*T_l.L() + T_l.G()*T_l.C()
          - T_l.E_a()*T_l.F_a() - T_l.E_b()*T_l.F_b() - T_l.E_c()*T_l.F_c();
        auto dc_u = T_u.M()*T_u.L() + T_u.G()*T_u.C()
          - T_u.E_a()*T_u.F_a() - T_u.E_b()*T_u.F_b() - T_u.E_c()*T_u.F_c();

        status |= test(dc_l[0], 0., "lower dc");
        status |= test(dc_l[1], 0., "lower dc");
        status |= test(dc_l[2], 0., "lower dc");

        status |= test(dc_u[0], 0., "upper dc");
        status |= test(dc_u[1], 0., "upper dc");
        status |= test(dc_u[2], 0., "upper dc");
      }



  /*
  std::cout << "\v\vtriangle position:\n";
  for ( unsigned int row=0; row<mesh.getRows(); ++row)
    for ( unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        auto & l (mesh.getLower(row, col));
        std::cout << row << "," << col <<" lower:\n\tA: "
                  << "(" << l.getA().x <<","<< l.getA().y << ")" << "\n\tnormale a: "
                  << "(" << l.getNormalA().x << "," << l.getNormalA().y << ")\n\tnormal b: "
                  << "(" << l.getNormalB().x << "," << l.getNormalB().y << ")\n\tnormal c: "
                  << "(" << l.getNormalC().x << "," << l.getNormalC().y << ")\n\tlength a: "
                  << l.getLengthA() << "\n\tlength b: "
                  << l.getLengthB() << "\n\tlength c: "
                  << l.getLengthC() << "\n";
      }

  for ( unsigned int row=0; row<mesh.getRows(); ++row)
    for ( unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        auto & u (mesh.getUpper(row, col));
        std::cout << row << "," << col <<" upper:\n\tA: "
                  << "(" << u.getA().x <<","<< u.getA().y << ")" << "\n\tnormale a: "
                  << "(" << u.getNormalA().x << "," << u.getNormalA().y << ")\n\tnormal b: "
                  << "(" << u.getNormalB().x << "," << u.getNormalB().y << ")\n\tnormal c: "
                  << "(" << u.getNormalC().x << "," << u.getNormalC().y << ")\n\tlength a: "
                  << u.getLengthA() << "\n\tlength b: "
                  << u.getLengthB() << "\n\tlength c: "
                  << u.getLengthC() << "\n";
      }
  std::cout << std::endl;
  */
  return status;
}











int main()
{
  int status = 0;

  std::cout << "Polynomial test";
  if (0==polynomialTest())
    std::cout << " successful" << std::endl;
  else
    {
      std::cout << " failed" << std::endl;
      --status;
    }

  std::cout << "Vector test ";
  if (0==vectorTest() )
    std::cout << "successful" << std::endl;
  else
    {
      std::cout << "failed" << std::endl;
      --status;
    }

  std::cout << "Mesh test ";
  if (0==triangleMeshTest() )
    std::cout << "successful" << std::endl;
  else
    {
      std::cout << "failed" << std::endl;
      --status;
    }

  std::cout << "Assembly test ";
  if (0==assamblyTest() )
    std::cout << "successful" << std::endl;
  else
  {
  std::cout << "failed" << std::endl;
  --status;
  }

  return status;
}
