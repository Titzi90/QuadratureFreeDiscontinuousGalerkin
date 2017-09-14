#ifndef SQUARGRID_HPP
#define SQUARGRID_HPP

#include <cmath>
#include <cassert>
#include <vector>

#include "monmomials_and_basefunctions.hpp"
#include "DataTypes.hpp"


/*****************************************************************************/
// Forward declaration
class UniqueSquareGrid;
class Triangle;
class Coordinates2D;
typedef Coordinates2D Point;
typedef Coordinates2D Vector;
typedef double (Jakobian)[2][2];



/*****************************************************************************/
// Class definition

/**
 * Grid Class for a square domain with unique distributed triangles
 * physical cells have indices [0 ... rows_/columns_[
 * ghost cells have indices [-1, rows_+1/columns_+1[
 *    x_2
 *    ^
 *    |
 *    +--->x_1
 */
class UniqueSquareGrid
{
public:
  /**
   * Construct an uniform triangle mesh on the unique square with 2*n^2 triangles
   * n is the refiment level in x and y direction
   * variance the max shift of the points in % from the uniform mesh; variance=[0,1]
   */
  UniqueSquareGrid(unsigned int n, double variance=0.);

  /**
   * get the lower triangle from cell defined by 'row' and 'col'
   */
  Triangle const & getLower(unsigned int row, unsigned int col) const
  {
    assert (row < rows_+1); assert (col < columns_+1);
    return lower_[row*(columns_+1) + col];
  }
  Triangle & getLower(unsigned int row, unsigned int col)
  {
    assert (row < rows_+1); assert (col < columns_+1);
    return lower_[row*(columns_+1) + col];
  }

  /**
   * get the upper triangle from cell defined by 'row' and 'col'
   */
  Triangle const & getUpper(unsigned int row, unsigned int col) const
  {
    ++row; ++col; // ghost layer in front of physical domain
    assert (row < rows_+1); assert (col < columns_+1);
    return upper_[row*(columns_+1) + col];
  }
  Triangle & getUpper(unsigned int row, unsigned int col)
  {
    ++row; ++col; // ghostlayer in front of physical domain
    assert (row < rows_+1); assert (col < columns_+1);
    return upper_[row*(columns_+1) + col];
  }

  unsigned int getRows() const { return rows_; }
  unsigned int getColumns() const { return columns_; }

  Point const & getVertex(unsigned int row, unsigned int col) const
  {
    assert (row < rows_+1); assert (col < columns_+1);
    return vertices_[row*(columns_+1) + col];
  }

private:
  std::vector<Triangle> lower_, upper_;
  std::vector<Point> vertices_;
  unsigned int rows_, columns_; // without ghost layer!!!!
};


/**
 * Struct holding 2D coordinates
 */
class Coordinates2D
{
public:
  Coordinates2D(){}
  Coordinates2D(double x, double y)
    :x(x), y(y)
  {}

  double x;
  double y;
};


/**
   Class holding all necessary information of an physical triangle for the simulation
*/
class Triangle
{
public:
  Triangle (Point const & a, Point const & b, Point const & c);

  Point const & getA() const { return A_; }

  double getLengthA() const { return length1_; }
  double getLengthB() const { return length2_; }
  double getLengthC() const { return length3_; }

  double getArea() const { return area_; }

  Vector const & getNormalA() const { return normale1_; }
  Vector const & getNormalB() const { return normale2_; }
  Vector const & getNormalC() const { return normale3_; }

  Jakobian const & getJakobian() const { return B_; }

  BlockMatrix const & M() const { return M_; }
  BlockMatrix const & G() const { return G_; }
  BlockMatrix const & E_a() const { return E_a_; }
  BlockMatrix const & E_b() const { return E_b_; }
  BlockMatrix const & E_c() const { return E_c_; }
  std::vector<double> const & C() const { return C_; }
  std::vector<double> const & U1() const { return U2_; }
  std::vector<double> const & U2() const { return U1_; }
  std::vector<double> const & L() const { return L_; }
  std::vector<double> const & Fn_a() const { return Fn_a_; }
  std::vector<double> const & Fn_b() const { return Fn_b_; }
  std::vector<double> const & Fn_c() const { return Fn_c_; }
  std::vector<double> const & Fr_a() const { return Fr_a_; }
  std::vector<double> const & Fr_b() const { return Fr_b_; }
  std::vector<double> const & Fr_c() const { return Fr_c_; }

  BlockMatrix & M() { return M_; }
  BlockMatrix & G() { return G_; }
  BlockMatrix & E_a() { return E_a_; }
  BlockMatrix & E_b() { return E_b_; }
  BlockMatrix & E_c() { return E_c_; }
  std::vector<double> & C() { return C_; }
  std::vector<double> & U1() { return U2_; }
  std::vector<double> & U2() { return U1_; }
  std::vector<double> & L() { return L_; }
  std::vector<double> & Fn_a() { return Fn_a_; }
  std::vector<double> & Fn_b() { return Fn_b_; }
  std::vector<double> & Fn_c() { return Fn_c_; }
  std::vector<double> & Fr_a() { return Fr_a_; }
  std::vector<double> & Fr_b() { return Fr_b_; }
  std::vector<double> & Fr_c() { return Fr_c_; }

private:
  double length1_, length2_, length3_, area_ ;       // TODO hier könnte man 1.5 double sparen
  Point A_;
  Vector normale1_, normale2_, normale3_;          // TODO auch hier muss nur 1.5 normalen gespeichert werden
  Jakobian B_;
  BlockMatrix M_, G_, E_a_, E_b_, E_c_;
  std::vector<double> C_, U1_, U2_, L_, Fn_a_, Fn_b_, Fn_c_, Fr_a_, Fr_b_, Fr_c_;
};



/*****************************************************************************/
// Free function definition

/**
 * subtract a two Vectors from each other
 */
inline Vector operator-(Vector const & lhs, Vector const & rhs)
{
  return Point(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline Point operator+(Point lhs, Point const & rhs)
{
  lhs.x+=rhs.x;
  lhs.y+=lhs.y;
  return lhs;
}

/**
 * scaling a vector
 */
inline Vector operator*(double lhs, Vector const & rhs)
{
  return Vector(lhs*rhs.x, lhs*rhs.y);
}

/**
 * computes the length of the vector
 */
inline double length(Vector const & vec)
{
  return std::sqrt(vec.x*vec.x + vec.y*vec.y);
}

/**
 * return a nomelized Vector of 'vec'
 */
inline Vector normalize(Vector const & vec)
{
  return 1./length(vec) * vec;
}
inline Vector normalize(Vector && vec)
{
  double l = length(vec);
  vec.x *= ( 1./l );
  vec.y *= ( 1./l );

  return vec;
}

/**
 * calculates the distance of two Points
 */
inline double distance (Point const & a, Point const &b)
{
  return length(b-a);
}

/**
 * calculates one (normalized) perpendicular vector to the edge given by p1 and p2
 */
inline Vector getNormal(Point const & p1, Point const & p2)
{

  return normalize( Vector(p1.y-p2.y, p2.x-p1.x) );
}

/**
 * calculate the inner product of two vectors
 */
inline double dot(Vector const & a, Vector const & b)
{
  return a.x*b.x + a.y*b.y;
}








std::ostream& printC(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printU(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printF(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printM(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printG(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printL(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printE(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& printFr(UniqueSquareGrid const & mesh, std::ostream& os);
std::ostream& operator<<(std::ostream & os, Point const & p);


#endif
