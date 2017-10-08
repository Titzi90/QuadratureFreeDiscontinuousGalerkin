
#include "polynomial.hpp"
#include "monmomials_and_basefunctions.hpp"

#include <vector>
#include <ostream>
#include <limits>
#include <string>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <numeric>

#include <iostream> // debugung



// template <class BinaryOperation>
// Polynomial1D operatorPlusMinusHelper(Polynomial1D  lhs, Polynomial1D const & rhs, BinaryOperation op)
// {

//   for (unsigned int x=0; x<=rhs.getOrder(); ++x)
//       lhs.get(x) = op(lhs.get(x), rhs.get(x));

//   return lhs;
// }
// template <class BinaryOperation>
// Polynomial2D operatorPlusMinusHelper(Polynomial2D lhs, Polynomial2D const & rhs, BinaryOperation op)
// {

//   for (unsigned int x=0; x<=rhs.getOrder(); ++x)
//     for (unsigned int y=0; y<=rhs.getOrder(); ++y)
//       lhs.get(x,y) = op(lhs.get(x,y), rhs.get(x,y));

//   return lhs;
// }
template <typename Polynomial>
Polynomial operatorPlusHelper(Polynomial lhs, Polynomial const & rhs)
{
  return lhs += rhs;
}
// plus operator for polynomials
template <typename Polynomial>
Polynomial operator+(Polynomial const & lhs, Polynomial const & rhs)
{
  return lhs.getOrder()>rhs.getOrder() ? operatorPlusHelper(lhs, rhs)
                                       : operatorPlusHelper(rhs, lhs);
}
template Polynomial1D operator+(Polynomial1D const & lhs, Polynomial1D const & rhs);
template Polynomial2D operator+(Polynomial2D const & lhs, Polynomial2D const & rhs);

// multiplication operator for polynomials
Polynomial1D operator*(Polynomial1D const & lhs, Polynomial1D const & rhs)
{
  Polynomial1D result = Polynomial1D(lhs.getOrder()+rhs.getOrder());

  //loop over lhs
  for (unsigned int lhs_x=0; lhs_x<=lhs.getOrder(); ++lhs_x)
  {
    double lhs_coefficent = lhs.get(lhs_x);
    if (0 == lhs_coefficent)
      continue;
    //loop over rhs
    for (unsigned int rhs_x=0; rhs_x<=rhs.getOrder(); ++rhs_x)
    {
      double rhs_coefficent = rhs.get(rhs_x);
      result.get(lhs_x+rhs_x) += lhs_coefficent*rhs_coefficent;
    }
  }
  return result;
}
Polynomial2D operator*(Polynomial2D const & lhs, Polynomial2D const & rhs)
{
  Polynomial2D result = Polynomial2D(lhs.getOrder()+rhs.getOrder());

  //loop over lhs
  for (unsigned int lhs_x=0; lhs_x<=lhs.getOrder(); ++lhs_x)
    for (unsigned int lhs_y=0; lhs_y<=lhs.getOrder(); ++lhs_y)
      {
        double lhs_coefficent = lhs.get(lhs_x, lhs_y);
        if (0 == lhs_coefficent)
          continue;
        //loop over rhs
        for (unsigned int rhs_x=0; rhs_x<=rhs.getOrder(); ++rhs_x)
          for (unsigned int rhs_y=0; rhs_y<=rhs.getOrder(); ++rhs_y)
            {
              double rhs_coefficent = rhs.get(rhs_x, rhs_y);
              result.get(lhs_x+rhs_x, lhs_y+rhs_y) += lhs_coefficent*rhs_coefficent;
            }
      }
  return result;
}
template <typename Polynomial>
Polynomial operator*(double const lhs, Polynomial const & rhs)
{
  Polynomial result = Polynomial(rhs.getOrder());
  std::transform(begin(rhs), end(rhs), begin(result),
                 [lhs](decltype(*begin(rhs)) x){ return lhs * x; });
  return result;
}
template Polynomial1D operator*(double lhs, Polynomial1D const & rhs);
template Polynomial2D operator*(double lhs, Polynomial2D const & rhs);



Polynomial2D derive(Polynomial2D const & pol, Variable const var)
{
  if (0 == pol.getOrder() )
  {
    return Polynomial2D(0);
  }

  Polynomial2D result (pol.getOrder()-1);

  for (unsigned int x=0; x<=result.getOrder(); ++x)
    for (unsigned int y=0; y<=result.getOrder(); ++y)
    {
      double exponentMul;
      unsigned int exponentX, exponentY;
      switch(var) {
      case Variable::X:
        exponentX   = x+1;
        exponentY   = y;
        exponentMul = exponentX;
        break;
      case Variable::Y:
        exponentX   = x;
        exponentY   = y+1;
        exponentMul = exponentY;
        break;
      default:
        assert(-1); // should never happen! But removed compiler warnings
        exponentY=-1;
        exponentX=-1;
        exponentMul=0;
      }
      result.get(x,y) = exponentMul*pol.get(exponentX,exponentY);
    }

  return result;
}

// converting a polynomial to a string
std::string to_string(Polynomial1D const & pol)
{
  std::stringstream ss;
  for(unsigned int x=pol.getOrder(); x!=std::numeric_limits<unsigned int>::max(); --x)
    if (0!=pol.get(x))
      //TODO plus nur wenn ned minus
      ss << pol.get(x) << "x^" << x << " + ";
  std::string str = ss.str();
  return str.size()>3 ? str.erase(str.size()-3) : str;
}
std::string to_string(Polynomial2D const & pol)
{
  std::stringstream ss;
  unsigned int order = pol.getOrder();

  for(unsigned int x=order; x!=std::numeric_limits<unsigned int>::max(); --x)
    {
      for(unsigned int y=order-x; y!=std::numeric_limits<unsigned int>::max(); --y)
        {
          if(0 != pol.get(x,y))
            {
              //TODO plus nur wenn ned minus
              ss << pol.get(x,y) << "x^" << x << "*y^" << y << " + ";
            }
        }
    }
  std::string str = ss.str();
  return str.size()>3 ? str.erase(str.size()-3) : str;
}
