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

/**
 * Class holding an 2 dimensional polynomial of order 'order'
 */
class Polynomial2D
{
public:
  // Constructor getting the order (-> size) of the polynomial
  Polynomial2D(unsigned int order)
    :order_(order), coeficents_(( order_+1 )*( order_+1 ), 0.)
  {}
  // Constructor for simple polynomials with one coefficient non equal to 0
  Polynomial2D(unsigned int xEx, unsigned int yEx, double coeficent)
    :order_(xEx+yEx), coeficents_((order_+1)*(order_+1),0.)
  {
    (*this)(xEx,yEx) = coeficent;
  }

  // geter for the order of the polynomial
  unsigned int getOrder() const {return order_;}

  // getter / setter for coefficients
  double const & operator()(unsigned int xExponent, unsigned int yExponent) const
  {
    // assert (xExponent + yExponent <= order_);
    return coeficents_[xExponent * (order_+1) + yExponent];
  }
  double& operator() (unsigned int xExponent, unsigned int yExponent)
  {
    // assert (xExponent + yExponent <= order_);
    return coeficents_[xExponent * (order_+1) + yExponent];
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

template <class BinaryOperation>
Polynomial2D operatorPlusMinusHelper(Polynomial2D const & lhs, Polynomial2D const & rhs, BinaryOperation op)
{
  Polynomial2D result (lhs);
  for (unsigned int x=0; x<=rhs.getOrder(); ++x)
    for (unsigned int y=0; y<=rhs.getOrder(); ++y)
      result(x,y) = op(result(x,y), rhs(x,y));
  return result;
}

// plus operator for polynomials
Polynomial2D operator+(Polynomial2D const & lhs, Polynomial2D const & rhs)
{
  return lhs.getOrder()>rhs.getOrder() ? operatorPlusMinusHelper(lhs, rhs, std::plus<double>())
                                       : operatorPlusMinusHelper(rhs, lhs, std::plus<double>());
}
Polynomial2D operator+(Polynomial2D const & lhs, double rhs){ return lhs + Polynomial2D(0,0,rhs); }
Polynomial2D operator+(double lhs, Polynomial2D const & rhs){ return rhs+lhs; }

// multiplication operator for polynomials
Polynomial2D operator*(Polynomial2D const & lhs, Polynomial2D const & rhs)
{
  Polynomial2D result = Polynomial2D(lhs.getOrder()+rhs.getOrder());

  //loop over lhs
  for (unsigned int lhs_x=0; lhs_x<=lhs.getOrder(); ++lhs_x)
    for (unsigned int lhs_y=0; lhs_y<=lhs.getOrder(); ++lhs_y)
    {
      double lhs_coefficent = lhs(lhs_x, lhs_y);
      if (0 == lhs_coefficent)
        continue;
      //loop over rhs
      for (unsigned int rhs_x=0; rhs_x<=rhs.getOrder(); ++rhs_x)
        for (unsigned int rhs_y=0; rhs_y<=rhs.getOrder(); ++rhs_y)
        {
          double rhs_coefficent = rhs(rhs_x, rhs_y);
          if (0 == rhs_coefficent)
            continue;
          result(lhs_x+rhs_x, lhs_y+rhs_y) += lhs_coefficent*rhs_coefficent;
        }
    }
  return result;
}
Polynomial2D operator*(double const lhs, Polynomial2D const & rhs){
  Polynomial2D result = Polynomial2D(rhs.getOrder());
  std::transform(begin(rhs), end(rhs), begin(result),
                 [lhs](decltype(*begin(rhs)) x){ return lhs * x; });
  return result;
}
Polynomial2D operator*(Polynomial2D const & lhs, double rhs){ return rhs*lhs; }

// minus operator for polynomials
Polynomial2D operator-(Polynomial2D const & lhs, Polynomial2D const & rhs){ return lhs + (-1*rhs); }
Polynomial2D operator-(Polynomial2D const & lhs, double rhs){ return lhs-Polynomial2D(0,0,rhs); }
Polynomial2D operator-(double lhs, Polynomial2D const & rhs){ return Polynomial2D(0,0,lhs)-rhs; }

//TODO division

// iterator access functions
auto begin(Polynomial2D & pol)      ->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
auto begin(Polynomial2D const & pol)->decltype(pol.coeficents_.begin()){ return pol.coeficents_.begin(); }
auto end(Polynomial2D & pol)        ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }
auto end(Polynomial2D const & pol)  ->decltype(pol.coeficents_.end())  { return pol.coeficents_.end(); }

// converting a polynomial to a string
std::string to_string(Polynomial2D const & pol)
{
  std::stringstream ss;
  unsigned int order = pol.getOrder();

  for(unsigned int x=order; x!=std::numeric_limits<unsigned int>::max(); --x)
    {
      for(unsigned int y=order-x; y!=std::numeric_limits<unsigned int>::max(); --y)
        {
          if(0 != pol(x,y))
            {
              //TODO plus nur wenn ned minus
              ss << pol(x,y) << "x^" << x << "y^" << y << " + ";
            }
        }
    }
  std::string str = ss.str();
  return str.erase(str.size()-3);
}

// output operator for polynomials
std::ostream& operator<< (std::ostream& os, Polynomial2D const & pol)
{
  os << to_string(pol);
  return os;
}
#endif
