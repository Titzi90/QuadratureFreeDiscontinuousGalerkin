#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "monmomials_and_basefunctions.hpp"
#include "squareGrid.hpp"

#include <vector>
#include <cassert>
#include <ostream>
#include <iomanip>
#include <functional>
#include <iostream>
#include <iterator>

// weighted sum of a function evaluated at q1,q2
// val = value + sum(weight * f(q1, q2))
template<class InputIt1, class InputIt2, class InputIt3, class T1, class T2, class T3>
T1 weighted_sum(InputIt1 weight_first, InputIt1 weight_last,
                InputIt2 q1, InputIt3 q2,
                std::function<T1(T2,T3)> f, T1 value)
{
  while (weight_first != weight_last) {
    value += *weight_first * f(*q1, *q2);
    ++weight_first;
    ++q1;
    ++q2;
  }
  return value;
}

// unsing gausion quatrature role to intefgare a function f over the reference triangle
// is exact for polynomials of order 'Order'
double integradeOverRefTriangle_gaus(std::function<double(double,double)> f, int Order)
{
  using std::begin;
  using std::end;

  std::vector<double> w,x1,x2;

  switch(Order){
  case 1:
    return f(1./3,1./3) * 0.5;
    break;
  case 2:
    w  = {1./6, 1./6, 1./6};
    x1 = {1./6, 2./3, 1./6};
    x2 = {1./6, 1./6, 2./3};
    break;
  default:
    std::cerr << "gaus integration not implemented!" << std::endl;
    return -1;
  }

  return weighted_sum(begin(w), end(w), begin(x1), begin(x2), f, 0.);
    }

/**
 * vector class to hold coefficients functions in local base representation
 */
class Coefficient
{
public:
  // Construct an empty coefficent vector
  Coefficient(int size, int numberOfBasefunctions)
    :data_(size), numberOfBasefunctions_(numberOfBasefunctions)
  {}

  // create a discreet Coefficient d_h using L2 projection in to the polynomial space
  // from a continues coeficentFunction d
  //TODO testen
  Coefficient(GridOnSquer mesh, int polynomialDegree, std::function<double(double,double)> coeficentFunction)
    :data_(mesh.getSize()*pol::polynomialGad[polynomialDegree]),
     numberOfBasefunctions_(pol::polynomialGad[polynomialDegree])
  {
    int t_id = 0;
    for (Triangle const & t : mesh)
    {
      int offset = t_id * numberOfBasefunctions_;

      // hatM(i,:) D_k(i) = int_hatT {phi(i)*d(F_k(hatx))} d hatx
      // here hatM == I     =>  D_k(i) = int...
      // int_hatT { phi*d(f_k(hatx)) }= int_hatT { phi * d(B_k*hatx+a_k) }
      // TODO algemein -> solver for SoE
      for (int i=0; i<numberOfBasefunctions_; ++i)
      {
        data_[offset+i] = integradeOverRefTriangle_gaus([&coeficentFunction, i,&t](double x1_hat, double x2_hat)
                                                        {
                                                          Jakobian B = t.getJakobian();
                                                          double x1 = B[0][0]*x1_hat + B[0][1]*x2_hat + t.getPointA().x;
                                                          double x2 = B[1][0]*x2_hat + B[1][1]*x2_hat + t.getPointA().y;
                                                          return pol::phi[i](x1_hat,x2_hat) * coeficentFunction(x1,x2);
                                                        },
                                                        polynomialDegree);
      }

      ++t_id;
    }
  }
  //TODO iterator

  double & get(int triangle_id, int base_id)
  {
    return data_[numberOfBasefunctions_*triangle_id + base_id];
  }
  double get(int triangle_id, int base_id) const
  {
    return data_[numberOfBasefunctions_*triangle_id + base_id];
  }
  double   operator[](int i) const { return data_[i]; }
  double & operator[](int i)       { return data_[i]; }

  int size() const { return data_.size(); }

  int getNumberOfbasisFunctions() const { return numberOfBasefunctions_; }

private:
  std::vector<double> data_;
  int const numberOfBasefunctions_;
};



/**
 * dense square Matrix class
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

std::ostream& operator<< (std::ostream& os, Coefficient const & c)
{
  os << std::fixed << std::setprecision(3);
  for (int i=0; i<c.size(); ++i)
    os << c[i] << "\n";

  return os;
}

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


// c = A*b
inline Coefficient matVecMul (Matrix const & A, Coefficient const & b)
{
  Coefficient c (b.size(), b.getNumberOfbasisFunctions());

  for (int i=0; i<b.size(); ++i)
  {
    double tmp = 0;
    for (int j=0; j<b.size(); ++j)
        tmp += A(i,j) * b[j];
    c[i] = tmp;
  }

  return c;
}

inline Coefficient operator* (Matrix const & A, Coefficient const & b) { return matVecMul(A,b); }

#endif
