#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "polynomial.hpp"
#include "monmomials_and_basefunctions.hpp"

#include <vector>
#include <cassert>
#include <ostream>
#include <iomanip>
#include <iterator>
#include <functional>
#include <algorithm>



/**
 * Matrix class to hold dense data for a local Block-Matrix
 */
class BlockMatrix
{
public:
  BlockMatrix(){}
  BlockMatrix(unsigned int N)
    :data_(N*N, 0.), N_(N),M_(N) {}
  BlockMatrix(unsigned int N, unsigned int M)
    :data_(N*M, 0.), N_(N),M_(M) {}

  double operator()(unsigned int i, unsigned int j) const
  {
    assert(i<N_);
    assert(j<M_);

    return data_[i*M_+j];
  }
  double & operator()(unsigned int i, unsigned int j)
  {
    assert(i<N_);
    assert(j<M_);

    return data_[i*M_+j];
  }

  unsigned int getN() const { return N_; }
  unsigned int getM() const { return M_; }

private:
  std::vector<double> data_;
  unsigned int N_, M_;
};

/**
 * dense cubic 3D-Tensor class
 */
class Tensor
{
public:
  Tensor(unsigned int N)
    :data_(N*N*N, 0.), N_(N) {}

  double operator()(unsigned int i, unsigned int j, unsigned int z) const
  {
    assert(i<N_);
    assert(j<N_);
    assert(z<N_);

    return data_[i*N_*N_+j*N_+z];
  }
  double & operator()(unsigned int i, unsigned int j, unsigned int z)
  {
    assert(i<N_);
    assert(j<N_);
    assert(z<N_);

    return data_[i*N_*N_+j*N_+z];
  }

  int getN() const { return N_; }

private:
  std::vector<double> data_;
  unsigned int N_;
};


// c = A*b
inline std::vector<double> matVecMul(BlockMatrix const & A, std::vector<double> const & b)
{
  assert(b.size() == A.getM());

  std::vector<double> c;
  c.reserve(A.getN());

  for (unsigned int i=0; i<A.getN(); ++i)
    {
      double tmp = 0.;
      for (unsigned int j=0; j<A.getM(); ++j)
        tmp += A(i,j) * b[j];
      c.push_back(tmp);
    }

  return c;
}
inline std::vector<double> operator* (BlockMatrix const & A, std::vector<double> const & b)
{
  return matVecMul(A,b);
}

inline std::vector<double> operator*(BlockMatrix const & M, Polynomial2D const & pol)
{
  return matVecMul(M,serialize(pol));
}

inline std::vector<double> operator* (double c, std::vector<double> const & v)
{
  std::vector<double> res;
  res.reserve(v.size());

  std::transform(v.begin(), v.end(), std::back_inserter(res), [c](double v){return c*v;});

  return res;
}

inline std::vector<double> operator+ (std::vector<double> const & lhs, std::vector<double> const & rhs)
{
  assert(lhs.size() == rhs.size());
  std::vector<double> res;
  res.reserve(lhs.size());

  std::transform(lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(res), [](double l, double r)
                 {return l+r;});

  return res;
}


/**
 * element wise multiplication of two vectors
 */
inline std::vector<double> elementWiseMul (std::vector<double> const & lhs,
                                           std::vector<double> const & rhs)
{
  assert (lhs.size() == rhs.size());

  std::vector<double> c;
  c.reserve(lhs.size());

  std::transform( lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(c),
                  std::multiplies<double>() );

  return c;
}



inline std::ostream& operator<< (std::ostream& os, std::vector<double> const & v)
{
  os <<  std::fixed << std::setprecision(3);
  for (auto i : v)
    os << std::setw(6) << i << "\n";

  return os;
}

inline std::ostream& operator<< (std::ostream& os, BlockMatrix const & m)
{
  os <<  std::fixed << std::setprecision(3);
  for (unsigned int i=0; i<m.getN(); ++i)
    {
      for (unsigned int j=0; j<m.getM();++j)
        os << std::setw(6) << m(i,j) << " ";

      os << "\n" ;
    }

  return os;
}

inline std::ostream& operator<< (std::ostream& os, Tensor const & t)
{
  os <<  std::fixed << std::setprecision(3);
  for (int i=0; i<t.getN(); ++i)
  {
    os << "i = " << i <<":\n" ;
    for (int j=0; j<t.getN();++j)
    {
      for (int z=0; z<t.getN(); ++z)
        os << std::setw(6) << t(i,j,z) << " ";
      os << "\n" ;
    }
    os << "\v" ;
  }

  return os;
}



#endif
