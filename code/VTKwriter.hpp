#ifndef VTKwriter_HPP
#define VTKwriter_HPP


#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "Grid.hpp"

#include <string>
#include <algorithm>


class VTKwriter
{
public:
  VTKwriter(std::string const & basename, UniqueSquareGrid const & mesh, unsigned int polynomialDegree)
    :basename_(basename), mesh_(mesh), order_(std::min(polynomialDegree, 2u)), iteration_(0)
  {}

  void write();

  std::string const & getName() const {return basename_;}

private:
  std::string const basename_;
  UniqueSquareGrid const & mesh_;
  unsigned int const order_;
  unsigned int iteration_;
};

#endif
