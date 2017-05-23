#ifndef SQUARGRID_HPP
#define SQUARGRID_HPP

#include <cmath>
#include <cassert>
#include <vector>



class Coordinate2D
{
public:
  Coordinate2D(){}
  Coordinate2D(double x, double y)
    :x(x), y(y)
  {}

  double x;
  double y;
};
typedef Coordinate2D Point;
typedef Coordinate2D Vector;

/**
 * subtract a two Vectors from each other
 */
inline Vector operator-(Vector const & lhs, Vector const & rhs)
{
  return Point(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline Vector operator*(double lhs, Vector const & rhs)
{
  return Vector(lhs*rhs.x, lhs*rhs.y);
}
inline Vector operator*(Vector const & lhs, double rhs)
{
  return rhs*lhs;
}

/**
 * computes the length of the vector 'vec'
 */
inline double length(Vector const & vec)
{
  return std::sqrt(vec.x*vec.x + vec.y*vec.y);
}

/**
 * return a nomelized Vector auf 'vec'
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
 * calculates one (normalize) perpendicular vector to the edge given by p1 and p2
 */
inline Vector getNormal(Point const & p1, Point const & p2)
{

  return normalize( Vector(p1.y-p2.y, p2.x-p1.x) );
}

/**
 * calculate the dot product of two vectors
 */
inline double dot(Vector const & a, Vector const & b)
{
  return a.x*b.x + a.y*b.y;
}

typedef double const (*Jakobian)[2];

class Triangle
{
public:
  Triangle(Point const & a, Point const & b, Point const & c)
  {
    coordinates_[0]=a;
    coordinates_[1]=b;
    coordinates_[2]=c;

    edgeLength_[0] = distance(b,c);
    edgeLength_[1] = distance(c,a);
    edgeLength_[2] = distance(a,b);

    // calculate the outer normal of each edge
    Vector na = getNormal(c,b);
    Vector nb = getNormal(a,c);
    Vector nc = getNormal(b,a);
    edgeNormale_[0] = dot(na,a-b)>0 ? -1*na : na;
    edgeNormale_[1] = dot(nb,b-c)>0 ? -1*nb : nb;
    edgeNormale_[2] = dot(nc,c-a)>0 ? -1*nc : nc;

    // calculate area of the triangle using Heron's formula
    double s = 0.5*(edgeLength_[0]+edgeLength_[1]+edgeLength_[2]);
    area_ = std::sqrt(s * (s-edgeLength_[0]) * (s-edgeLength_[1]) * (s-edgeLength_[2]));

    // calculate Jakobian matrix [B-A | C-A]
    jakobian_[0][0] = getPointB().x - getPointA().x;
    jakobian_[1][0] = getPointB().y - getPointA().y;
    jakobian_[0][1] = getPointC().x - getPointA().x;
    jakobian_[1][1] = getPointC().y - getPointA().y;
  }

  Point const & getPointA() const { return coordinates_[0]; }
  Point const & getPointB() const { return coordinates_[1]; }
  Point const & getPointC() const { return coordinates_[2]; }
  // Point       & getPointA()       { return coordinates_[0]; }
  // Point       & getPointB()       { return coordinates_[1]; }
  // Point       & getPointC()       { return coordinates_[2]; }

  double   getLengthA() const { return edgeLength_[0]; }
  double   getLengthB() const { return edgeLength_[1]; }
  double   getLengthC() const { return edgeLength_[2]; }
  // double & getLengthA()       { return edgeLength_[0]; }
  // double & getLengthB()       { return edgeLength_[1]; }
  // double & getLengthC()       { return edgeLength_[2]; }

  Vector const & getNormalA() const { return edgeNormale_[0]; }
  Vector const & getNormalB() const { return edgeNormale_[1]; }
  Vector const & getNormalC() const { return edgeNormale_[2]; }
  // Vector       & getNormalA()       { return edgeNormale_[0]; }
  // Vector       & getNormalB()       { return edgeNormale_[1]; }
  // Vector       & getNormalC()       { return edgeNormale_[2]; }

  double   getArea() const { return area_; }
  // double & getArea()       { return area_; }


  Jakobian getJakobian() const { return jakobian_; }


private:
  Point coordinates_[3];
  double edgeLength_[3];
  Vector edgeNormale_[3];
  double area_;
  double jakobian_[2][2];
};




/**
 *    x_2
 *    ^
 *    |
 *    |
 *    |
 *    |
 *    |
 *    +----------->x_1
 */
class GridOnSquer
{
public:
  /**
   * Construct an uniform triangle mesh on the unique square with 2*n^2 triangles
   */
  GridOnSquer(unsigned int n)
    : numberOfRows_(n), numberOfColums_(n)
  {
    double const h = 1./n;

    for (unsigned int row=0; row<n; ++row)
      for (unsigned int col=0; col<n; ++col)
      {
        Point ll ( col   *h,  row   *h); // lower left corner
        Point lr ((col+1)*h,  row   *h); // lower right corner
        Point ul ( col   *h, (row+1)*h); // upper left corner
        Point ur ((col+1)*h, (row+1)*h); // upper right corner

        gridData_.push_back(Triangle(ll,lr,ul)); // lower triangle
        gridData_.push_back(Triangle(lr,ur,ul)); // upper triangle
      }
  }

  /**
   * returns the triangle defined by 'row', 'col' and the 'position'
   * 'position' is 0 for lower triangle and 1 for upper triangle
   */
  Triangle & get(unsigned int row, unsigned int col, unsigned int position)
  {
    assert (row < numberOfRows_);
    assert (col < numberOfColums_);
    assert (position < 2);
    return gridData_[row*numberOfRows_*2 + col*2 + position];
  }
  // Triangle const & get(unsigned int row, unsigned int col, unsigned int position) const
  // {
  //   assert (row < numberOfRows_);
  //   assert (col < numberOfColums_);
  //   assert (position < 2);
  //   return gridData_[row*numberOfRows_*2 + col*2 + position];
  // }

  unsigned int getRows() const { return numberOfRows_; }
  unsigned int getColums() const { return numberOfColums_; }
  unsigned int getSize() const { return gridData_.size(); }

private:
  std::vector<Triangle> gridData_;
  unsigned int numberOfRows_;
  unsigned int numberOfColums_;

  // friends declarations
  friend auto begin(GridOnSquer &) ->decltype(gridData_.begin());
  friend auto begin(GridOnSquer const &) ->decltype(const_cast<const std::vector<Triangle>&>(gridData_).begin());

  friend auto end(GridOnSquer &) ->decltype(gridData_.end());
  friend auto end(GridOnSquer const &) ->decltype(const_cast<const std::vector<Triangle>&>(gridData_).end());
};



// iterator access function for grid elements
inline auto begin(GridOnSquer & mesh)
  ->decltype(mesh.gridData_.begin()){ return mesh.gridData_.begin(); }
inline auto begin(GridOnSquer const & mesh)
  ->decltype(const_cast<const std::vector<Triangle>&>(mesh.gridData_).begin()){ return mesh.gridData_.begin(); }

inline auto end(GridOnSquer & mesh) ->decltype(mesh.gridData_.end()){ return mesh.gridData_.end(); }
inline auto end(GridOnSquer const & mesh)
  ->decltype(const_cast<const std::vector<Triangle>&>(mesh.gridData_).end()){ return mesh.gridData_.end(); }


#endif
