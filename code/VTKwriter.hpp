#ifndef VTKWRITER_HPP
#define VTKwriter_HPP


#include "monmomials_and_basefunctions.hpp"  //muss hier als ertets includet werden wegen seltsammen abhägigkeitssscheis
#include "Grid.hpp"

#include <string>


class VTKwriter
{
public:
  VTKwriter(std::string const & basename, UniqueSquareGrid const & mesh)
    :basename_(basename), mesh_(mesh), iteration_(0)
  {}

  void write();

private:
  std::string const basename_;
  UniqueSquareGrid const & mesh_;
  unsigned int iteration_;
};

#endif
