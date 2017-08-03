#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <ostream>
#include <limits>
#include <string>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <numeric>
#include <iostream>

// #include "monmomials_and_basefunctions.hpp"

enum class Variable {X, Y};

/**
 * Class holding an 1 dimensional polynomial of order 'order'
 */
class Polynomial1D
{
public:
  // standard Constructor creates 0 as a polynomial
  Polynomial1D()
    :order_(0), coeficents_(order_+1,0.)
  {}
  // Constructor getting the order (-> size) of the polynomial
  Polynomial1D(unsigned int order)
    :order_(order), coeficents_(order_+1, 0.)
  {}
  // Constructor for simple polynomials with one coefficient non equal to 0
  Polynomial1D(unsigned int xEx, double coeficent)
    :order_(xEx), coeficents_(order_+1,0.)
  {
    (*this).get(xEx) = coeficent;
  }

  // geter for the order of the polynomial
  unsigned int getOrder() const {return order_;}

  // getter / setter for coefficients
  double get(unsigned int xExponent) const { return coeficents_[xExponent]; }
  double & get(unsigned int xExponent) { return coeficents_[xExponent]; }

  // evaluate polynomial
  double operator()(double x) const
  {
    double val = 0.;

    for (unsigned int xExponent=0; xExponent<=order_; ++xExponent)
        val += (*this).get(xExponent) * std::pow(x,xExponent);

    return val;
  }


private:
  unsigned int order_;
  std::vector<double> coeficents_;

  // friends declarations
  friend auto begin(Polynomial1D & pol)      ->decltype(coeficents_.begin());
  friend auto begin(Polynomial1D const & pol)->decltype(const_cast<const std::vector<double>&>(coeficents_).begin());
  friend auto end(Polynomial1D & pol)        ->decltype(coeficents_.end());
  friend auto end(Polynomial1D const & pol)  ->decltype(const_cast<const std::vector<double>&>(coeficents_).end());
};

/**
 * Class holding an 2 dimensional polynomial of order 'order'
 */
class Polynomial2D
{
public:
  // standard Constructor creates 0 as a polynomial
  Polynomial2D()
    :order_(0), coeficents_(( order_+1 )*( order_+1 ), 0.)
  {}
  // Constructor getting the order (-> size) of the polynomial
  Polynomial2D(unsigned int order)
    :order_(order), coeficents_(( order_+1 )*( order_+1 ), 0.)
  {}
  // Constructor for simple polynomials with one coefficient non equal to 0
  Polynomial2D(unsigned int xEx, unsigned int yEx, double coeficent)
    :order_(xEx+yEx), coeficents_((order_+1)*(order_+1),0.)
  {
    (*this).get(xEx,yEx) = coeficent;
  }

  // geter for the order of the polynomial
  unsigned int getOrder() const { return order_; }

  // getter / setter for coefficients
  double get(unsigned int xExponent, unsigned int yExponent) const
  {
    return coeficents_[xExponent * (order_+1) + yExponent];
  }
  double & get(unsigned int xExponent, unsigned int yExponent)
  {
    return coeficents_[xExponent * (order_+1) + yExponent];
  }

  // evaluate polynomial
  double operator()(double x, double y) const
  {
    double val = 0.;

    for (unsigned int xExponent=0; xExponent<=order_; ++xExponent)
      for (unsigned int yExponent=0; yExponent<=order_; ++yExponent)
        val += (*this).get(xExponent,yExponent) * std::pow(x,xExponent) * std::pow(y,yExponent);

    return val;
  }


private:
  unsigned int order_;
  std::vector<double> coeficents_;

  // friends declarations
  friend auto begin(Polynomial2D & pol)      ->decltype(coeficents_.begin());
  friend auto begin(Polynomial2D const & pol)->decltype(const_cast<const std::vector<double>&>(coeficents_).begin());
  friend auto end(Polynomial2D & pol)        ->decltype(coeficents_.end());
  friend auto end(Polynomial2D const & pol)  ->decltype(const_cast<const std::vector<double>&>(coeficents_).end());
};

inline std::vector<double> serialize(Polynomial2D const & pol)
{
  //TODO testen
  std::vector<double> c;
  // c.reserve(numberOf2DBasefunctions(pol.getOrder()));

  for (unsigned int i=0; i<=pol.getOrder(); ++i)
    for (unsigned int x=i; x!=std::numeric_limits<unsigned int>::max(); --x)
      {
        unsigned int y = i-x;
        c.push_back(pol.get(x,y));
      }

  return c;
}

// plus assignment operator for polynomials
inline Polynomial1D& operator+=(Polynomial1D & lhs, Polynomial1D const & rhs)
{
  assert(lhs.getOrder()>=rhs.getOrder());
  for (unsigned int s=0; s<=rhs.getOrder(); ++s)
    lhs.get(s) += rhs.get(s);
  return lhs;
}
inline Polynomial2D& operator+=(Polynomial2D & lhs, Polynomial2D const & rhs)
{
  assert(lhs.getOrder()>=rhs.getOrder());
  for (unsigned int x=0; x<=rhs.getOrder(); ++x)
    for (unsigned int y=0; y<=rhs.getOrder(); ++y)
      lhs.get(x,y) += rhs.get(x,y);
  return lhs;
}

// plus operator for polynomials
template <typename Polynomial> Polynomial operator+(Polynomial const & lhs, Polynomial const & rhs);
inline Polynomial1D operator+(Polynomial1D const & lhs, double rhs){ return lhs + Polynomial1D(0,rhs); }
inline Polynomial2D operator+(Polynomial2D const & lhs, double rhs){ return lhs + Polynomial2D(0,0,rhs); }
template<typename Polynomial> inline Polynomial operator+(double lhs, Polynomial const & rhs){ return rhs+lhs; }

// multiplication operator for polynomials
Polynomial1D operator*(Polynomial1D const & lhs, Polynomial1D const & rhs);
Polynomial2D operator*(Polynomial2D const & lhs, Polynomial2D const & rhs);
template<typename Polynomial> Polynomial operator*(double lhs, Polynomial const & rhs);
template<typename Polynomial> inline Polynomial operator*(Polynomial const & lhs, double rhs){ return rhs*lhs; }

// minus operator for polynomials
template<typename Polynomial> Polynomial operator-(Polynomial const & lhs, Polynomial const & rhs)
{ return lhs + (-1*rhs); }
template<typename Polynomial> Polynomial operator-(Polynomial const & lhs, double rhs)
{ return lhs + (-1*rhs); }
template<typename Polynomial> Polynomial operator-(double lhs, Polynomial const & rhs)
{ return lhs + (-1*rhs); }


//TODO division


Polynomial2D derive(Polynomial2D const & pol, Variable const var);

// iterator access functions
inline auto begin(Polynomial1D & pol)      ->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto begin(Polynomial1D const & pol)->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto end(Polynomial1D & pol)        ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }
inline auto end(Polynomial1D const & pol)  ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }

inline auto begin(Polynomial2D & pol)      ->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto begin(Polynomial2D const & pol)->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto end(Polynomial2D & pol)        ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }
inline auto end(Polynomial2D const & pol)  ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }

// converting a polynomial to a string
std::string to_string(Polynomial1D const & pol);
std::string to_string(Polynomial2D const & pol);

// output operator for polynomials
inline std::ostream& operator<< (std::ostream& os, Polynomial1D const & pol)
{
  os << to_string(pol);
  return os;
}
inline std::ostream& operator<< (std::ostream& os, Polynomial2D const & pol)
{
  os << to_string(pol);
  return os;
}














#endif
