#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <cassert>
#include <ostream>
#include <iomanip>

/**
 * dense square Matrix class
 * TODO testen
 */
class Matrix
{
public:
  Matrix(unsigned int N)
    :data_(N*N, 0.), N_(N) {}

  double operator()(unsigned int i, unsigned int j) const
  {
    assert(i<N_);
    assert(j<N_);

    return data_[i*N_+j];
  }
  double & operator()(unsigned int i, unsigned int j)
  {
    assert(i<N_);
    assert(j<N_);

    return data_[i*N_+j];
  }

  int getSize() const { return N_; }

private:
  std::vector<double> data_;
  unsigned int N_;
};

/**
 * dense cubic 3D-Tensor class
 * TODO testen
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

  int getSize() const { return N_; }

private:
  std::vector<double> data_;
  unsigned int N_;
};

std::ostream& operator<< (std::ostream& os, Matrix const & m)
{
  os <<  std::fixed << std::setprecision(3);
  for (int i=0; i<m.getSize(); ++i)
  {
    for (int j=0; j<m.getSize();++j)
    os << std::setw(6) << m(i,j) << " ";
    os << "\n" ;
  }

  return os;
}

std::ostream& operator<< (std::ostream& os, Tensor const & t)
{
  os <<  std::fixed << std::setprecision(3);
  for (int i=0; i<t.getSize(); ++i)
  {
    os << "i = " << i <<":\n" ;
    for (int j=0; j<t.getSize();++j)
    {
      for (int z=0; z<t.getSize(); ++z)
        os << std::setw(6) << t(i,j,z) << " ";
      os << "\n" ;
    }
    os << "\v" ;
  }

  return os;
}


#endif
