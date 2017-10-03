#include "monmomials_and_basefunctions.hpp"
#include "Grid.hpp"
#include "DataTypes.hpp"

#include  <cstdlib>
#include <ctime>

UniqueSquareGrid::UniqueSquareGrid(unsigned int n, double maxVariance)
  :rows_(n), columns_(n)
{
  lower_.reserve((n+1)*(n+1));
  upper_.reserve((n+1)*(n+1));
  vertices_.reserve((n+1)*(n+1));

  double const h = 1./n;
  // std::srand(std::time(0));


  // POINTS
  for (unsigned int row=0; row<rows_+1; ++row)
    for (unsigned int col=0; col<columns_+1; ++col)
      {
        // randome shift
        double varX = 2.*std::rand()/RAND_MAX  - 1;
        double varY = 2.*std::rand()/RAND_MAX  - 1;
        double dhx = 0., dhy = 0.;
        if (0 != col && columns_ != col)
          dhx = varX*maxVariance * h;
        if (0 != row && rows_ != row)
          dhy = varY*maxVariance * h;
        vertices_.push_back(Point(col*h +dhx, row*h +dhy));
      }

  //TRIANGLES
  // lower left ghost corner (just for indexing)
  upper_.push_back( Triangle({0.,-h}, getVertex(0, 0), {-h,0.}) );

  // ghost layer bottom row
  for (unsigned int col=0; col<columns_; ++col)
    {
      Point lr ((col+1)*h, -h);
      Point const & ur = getVertex(0, col+1);
      Point const & ul = getVertex(0, col  );

      upper_.push_back( Triangle(lr,ur,ul) );
    }

  for (unsigned int row=0; row<rows_; ++row)
    {
      { // add ghost triangle (left site)
        Point const & lr = getVertex(row  , 0); // lower right corner
        Point const & ur = getVertex(row+1, 0); // upper right corner
        Point ul (-h, (row+1)*h);               // upper left corner

          upper_.push_back( Triangle(lr,ur,ul) );
      }

      for (unsigned int col=0; col<columns_; ++col)
        { // add physical triangles
          Point const & ll = getVertex(row  , col  ); // lower left corner
            Point const & lr = getVertex(row  , col+1); // lower right corner
            Point const & ul = getVertex(row+1, col  ); // upper left corner
            Point const & ur = getVertex(row+1, col+1); // upper right corner

            lower_.push_back( Triangle(ll,lr,ul) );
            upper_.push_back( Triangle(lr,ur,ul) );
        }

      { // add ghost triangle (right site)
        Point const & ll = getVertex(row  , columns_); // lower left corner
        Point lr ( 1.+h,  row*h);                        // lower right corner
        Point const & ul = getVertex(row+1, columns_); // upper left corner

        lower_.push_back( Triangle(ll,lr,ul) );
      }
    }

    // ghost layer top row
  for (unsigned int col=0; col<columns_; ++col)
    {
      Point const & ll = getVertex(rows_, col  ); // lower left corner
      Point const & lr = getVertex(rows_, col+1); // lower right corner
      Point ul ( col*h, 1.+h);                      // upper left corner

      lower_.push_back( Triangle(ll,lr,ul) );
    }

  // top right ghost corner (just for indexing)
    lower_.push_back( Triangle(getVertex(rows_, columns_), {1.+h, 1.}, {1., 1.+h}) );
}



Triangle::Triangle (Point const & a, Point const & b, Point const & c)
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


std::ostream& printC(UniqueSquareGrid const & mesh, std::ostream& os)
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

std::ostream& printU(UniqueSquareGrid const & mesh, std::ostream& os)
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

 std::ostream& printF(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "Fn_a:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fn_a();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fn_a();
  os << "\vFn_b:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fn_b();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fn_b();
  os << "\vFn_c:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fn_c();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fn_c();

  return os;
}

 std::ostream& printM(UniqueSquareGrid const & mesh, std::ostream& os)
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

 std::ostream& printG(UniqueSquareGrid const & mesh, std::ostream& os)
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

 std::ostream& printL(UniqueSquareGrid const & mesh, std::ostream& os)
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

 std::ostream& printE(UniqueSquareGrid const & mesh, std::ostream& os)
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

 std::ostream& printFr(UniqueSquareGrid const & mesh, std::ostream& os)
{
  os << "Fra:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fr_a();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fr_a();

  os << "\vFrb:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fr_b();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fr_b();

  os << "\vFrc:\n";
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getLower(row, col).Fr_c();
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      os << mesh.getUpper(row, col).Fr_c();

  return os;
}

 std::ostream& operator<<(std::ostream & os, Point const & p)
{
  os << p.x << " " << p.y;
  return os;
}
