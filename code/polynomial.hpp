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


enum class Variable {X, Y};

/**
 * Class holding an 2 dimensional polynomial of order 'order'
 */
class Polynomial2D
{
public:
  // Constructor getting the order (-> size) of the polynomial
  Polynomial2D(unsigned int order);
  // Constructor for simple polynomials with one coefficient non equal to 0
  Polynomial2D(unsigned int xEx, unsigned int yEx, double coeficent);

  // geter for the order of the polynomial
  unsigned int getOrder() const {return order_;}

  // getter / setter for coefficients
  double get(unsigned int xExponent, unsigned int yExponent) const
  {
    return coeficents_[xExponent * (order_+1) + yExponent];
  }

  double& get(unsigned int xExponent, unsigned int yExponent)
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
  unsigned int const order_;
  std::vector<double> coeficents_;

  // friends declarations
  friend auto begin(Polynomial2D & pol)      ->decltype(coeficents_.begin());
  friend auto begin(Polynomial2D const & pol)->decltype(const_cast<const std::vector<double>&>(coeficents_).begin());
  friend auto end(Polynomial2D & pol)        ->decltype(coeficents_.end());
  friend auto end(Polynomial2D const & pol)  ->decltype(const_cast<const std::vector<double>&>(coeficents_).end());
};

// plus operator for polynomials
Polynomial2D operator+(Polynomial2D const & lhs, Polynomial2D const & rhs);
inline Polynomial2D operator+(Polynomial2D const & lhs, double rhs){ return lhs + Polynomial2D(0,0,rhs); }
inline Polynomial2D operator+(double lhs, Polynomial2D const & rhs){ return rhs+lhs; }

// multiplication operator for polynomials
Polynomial2D operator*(Polynomial2D const & lhs, Polynomial2D const & rhs);
Polynomial2D operator*(double const lhs, Polynomial2D const & rhs);
inline Polynomial2D operator*(Polynomial2D const & lhs, double rhs){ return rhs*lhs; }

// minus operator for polynomials
inline Polynomial2D operator-(Polynomial2D const & lhs, Polynomial2D const & rhs){ return lhs + (-1*rhs); }
inline Polynomial2D operator-(Polynomial2D const & lhs, double rhs){ return lhs-Polynomial2D(0,0,rhs); }
inline Polynomial2D operator-(double lhs, Polynomial2D const & rhs){ return Polynomial2D(0,0,lhs)-rhs; }

//TODO division

double integradeOverRefTriangle(Polynomial2D const & pol);

//TODO integrieren über kante

Polynomial2D derive(Polynomial2D const & pol, Variable const var);

// iterator access functions
inline auto begin(Polynomial2D & pol)      ->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto begin(Polynomial2D const & pol)->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
inline auto end(Polynomial2D & pol)        ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }
inline auto end(Polynomial2D const & pol)  ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }

// converting a polynomial to a string
std::string to_string(Polynomial2D const & pol);

// output operator for polynomials
std::ostream& operator<< (std::ostream& os, Polynomial2D const & pol);

#endif
