#ifndef SQUARGRID_HPP
#define SQUARGRID_HPP

#include <cmath>
#include <cassert>
#include <vector>

#include "monmomials_and_basefunctions.hpp"
#include "DataTypes.hpp"


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
typedef Coordinates2D Point;
typedef Coordinates2D Vector;

/**
 * subtract a two Vectors from each other
 */
inline Vector operator-(Vector const & lhs, Vector const & rhs)
{
  return Point(lhs.x-rhs.x, lhs.y-rhs.y);
}

/**
 * scaling a vector
 */
inline Vector operator*(double lhs, Vector const & rhs)
{
  return Vector(lhs*rhs.x, lhs*rhs.y);
}
inline Vector operator*(Vector const & lhs, double rhs)
{
  return rhs*lhs;
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

// typedef double (*Jakobian)[2]; // 2D Jacobian Matrix -> 2x2
typedef double (Jakobian)[2][2];


/**
   Class holding all necessary information of an physical triangle for the simulation
*/
class Triangle
{
public:
  Triangle (Point const & a, Point const & b, Point const & c)
    :length1_(distance(b, c)), length2_(distance(c, a)), length3_(distance(a,b)),
     A_(a)
  {
    // calculate area of the triangle using Heron's formula
    double s = 0.5*(length1_ + length2_ + length3_);
    area_ = std::sqrt(s * (s-length1_) * (s-length2_) * (s-length3_));

    // calculate the outer normal of each edge
    Vector na = getNormal(c,b);
    Vector nb = getNormal(a,c);
    Vector nc = getNormal(b,a);
    normale1_ = dot(na,a-b)>0 ? -1*na : na;
    normale2_ = dot(nb,b-c)>0 ? -1*nb : nb;
    normale3_ = dot(nc,c-a)>0 ? -1*nc : nc;

    // calculate Jakobian matrix [B-A | C-A]
    B_[0][0] = b.x - a.x;
    B_[1][0] = b.y - a.y;
    B_[0][1] = c.x - a.x;
    B_[1][1] = c.y - a.y;
  }

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
  std::vector<double> const & F1() const { return F1_; }
  std::vector<double> const & F2() const { return F2_; }
  std::vector<double> const & F_a() const { return Fr_a_; }
  std::vector<double> const & F_b() const { return Fr_b_; }
  std::vector<double> const & F_c() const { return Fr_c_; }

  BlockMatrix & M() { return M_; }
  BlockMatrix & G() { return G_; }
  BlockMatrix & E_a() { return E_a_; }
  BlockMatrix & E_b() { return E_b_; }
  BlockMatrix & E_c() { return E_c_; }
  std::vector<double> & C() { return C_; }
  std::vector<double> & U1() { return U2_; }
  std::vector<double> & U2() { return U1_; }
  std::vector<double> & L() { return L_; }
  std::vector<double> & F1() { return F1_; }
  std::vector<double> & F2() { return F2_; }
  std::vector<double> & F_a() { return Fr_a_; }
  std::vector<double> & F_b() { return Fr_b_; }
  std::vector<double> & F_c() { return Fr_c_; }



private:
  double length1_, length2_, length3_, area_ ;       // TODO hier könnte man 1.5 double sparen
  Point A_;
  Vector normale1_, normale2_, normale3_;          // TODO auch hier muss nur 1.5 normalen gespeichert werden
  Jakobian B_;
  BlockMatrix M_, G_, E_a_, E_b_, E_c_;
  std::vector<double> C_, U1_, U2_, L_, F1_, F2_, Fr_a_, Fr_b_, Fr_c_;
};




/**
 * Grid Class for a square domain with unique distributed triangles
 * physical cells have indices [0 ... rows_/columns_[
 * ghost cells have indices [-1, rows_+1/columns_+1[
 *
 *    x_2
 *    ^
 *    |
 *    |
 *    |
 *    |
 *    |
 *    +----------->x_1
 */
class UniqueSquareGrid
{
public:
  /**
   * Construct an uniform triangle mesh on the unique square with 2*n^2 triangles
   */
  UniqueSquareGrid(unsigned int n)
    :rows_(n), columns_(n)
  {
    lower_.reserve((n+1)*(n+1));
    upper_.reserve((n+1)*(n+1));
    vertices_.reserve((n+1)*(n+1));

    double const h = 1./n;

    // POINTS
    for (unsigned int row=0; row<rows_+1; ++row)
      for (unsigned int col=0; col<columns_+1; ++col)
        vertices_.push_back(Point(row*h, col*h));





    //TRIANGLES
    //TODO neu aus punkt menge generieren

    // ghost layer bottom row
    upper_.push_back( Triangle({0.,-h}, {0.,0.}, {-h,0.}) );
    for (unsigned int col=0; col<columns_; ++col)
      {
        Point lr ((col+1)*h, -h);
        Point ur ((col+1)*h, 0.);
        Point ul ( col   *h, 0.);
        upper_.push_back( Triangle(lr,ur,ul) );
      }
    for (unsigned int row=0; row<rows_; ++row)
      {
        // add ghost triangle (left site)
        {
          Point lr (0.,  row   *h); // lower right corner
          Point ur (0., (row+1)*h); // upper right corner
          Point ul (-h, (row+1)*h); // upper left corner
          upper_.push_back( Triangle(lr,ur,ul) );
        }
        for (unsigned int col=0; col<columns_; ++col)
          {
            // add physical triangles
            Point ll ( col   *h,  row   *h); // lower left corner
            Point lr ((col+1)*h,  row   *h); // lower right corner
            Point ur ((col+1)*h, (row+1)*h); // upper right corner
            Point ul ( col   *h, (row+1)*h); // upper left corner

            lower_.push_back( Triangle(ll,lr,ul) );
            upper_.push_back( Triangle(lr,ur,ul) );
          }
        // add ghost triangle (right site)
        {
          Point ll ( 1.  ,  row   *h); // lower left corner
          Point lr ( 1.+h,  row   *h); // lower right corner
          Point ul ( 1.  , (row+1)*h); // upper left corner
          lower_.push_back( Triangle(ll,lr,ul) );
        }
      }
    // ghost layer top row
    for (unsigned int col=0; col<columns_; ++col)
      {
        Point ll ( col   *h, 1.  ); // lower left corner
        Point lr ((col+1)*h, 1.  ); // lower right corner
        Point ul ( col   *h, 1.+h); // upper left corner
        lower_.push_back( Triangle(ll,lr,ul) );
      }
  }

  /**
   * get the lower triangle from cell defined by 'row' and 'col'
   */
  Triangle const & getLower(unsigned int row, unsigned int col) const
  {
    assert (row < rows_+1);
    assert (col < columns_+1);
    return lower_[row*(columns_+1) + col];
  }
  Triangle & getLower(unsigned int row, unsigned int col)
  {
    assert (row < rows_+1);
    assert (col < columns_+1);
    return lower_[row*(columns_+1) + col];
  }

  /**
   * get the upper triangle from cell defined by 'row' and 'col'
   */
  Triangle const & getUpper(unsigned int row, unsigned int col) const
  {
    ++row; ++col; // ghost layer in front of physical domain

    assert (row < rows_+1);
    assert (col < columns_+1);
    return upper_[row*(columns_+1) + col];
  }
  Triangle & getUpper(unsigned int row, unsigned int col)
  {
    ++row; ++col; // ghostlayer in front of physical domain

    assert (row < rows_+1);
    assert (col < columns_+1);
    return upper_[row*(columns_+1) + col];
  }

  unsigned int getRows() const { return rows_; }
  unsigned int getColumns() const { return columns_; }

  Point const & getVertex(unsigned int row, unsigned int col) const
  {
    assert (row < rows_+1);
    assert (col < columns_+1);
    return vertices_[row*(columns_+1) + col];
  }

private:
  std::vector<Triangle> lower_, upper_;
  std::vector<Point> vertices_;

  unsigned int rows_, columns_; // without ghost layer!!!!
};

inline std::ostream& printC(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "C:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).C();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).C();

  return os;
}

inline std::ostream& printU(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "U1:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).U1();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).U1();
  os << "\vU2:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).U2();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).U2();

  return os;
}

inline std::ostream& printF(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "F1:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).F1();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).F1();
  os << "\vF2:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).F2();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).F2();

  return os;
}

inline std::ostream& printM(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "M:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).M();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).M();

  return os;
}

inline std::ostream& printG(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "G:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).G();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).G();

  return os;
}

inline std::ostream& printL(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "L:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).L();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).L();

  return os;
}

inline std::ostream& printE(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "Ea:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).E_a();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).E_a();

  os << "\vEb:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).E_b();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).E_b();

  os << "\vEc:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).E_c();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).E_c();

  return os;
}

inline std::ostream& printFr(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "Fra:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).F_a();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).F_a();

  os << "\vFrb:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).F_b();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).F_b();

  os << "\vFrc:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).F_c();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).F_c();

  return os;
}

inline std::ostream& operator<<(std::ostream & os, Point const & p)
{
  os << p.x << " " << p.y;
  return os;
}


#endif
