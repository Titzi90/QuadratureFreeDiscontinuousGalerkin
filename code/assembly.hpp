#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include "monmomials_and_basefunctions.hpp"
#include "squareGrid.hpp"
#include "Matrix.hpp"

#include <vector>






/**
 * Generate local Mass-Matrix for ref triangle hatM
 */
inline Matrix assemblyHatM (int polynomialDegree)
{
  Matrix hatM (pol::polynomialGad[polynomialDegree]);
  for (int i=0; i<pol::polynomialGad[polynomialDegree]; ++i)
    for (int j=0; j<pol::polynomialGad[polynomialDegree]; ++j)
      hatM(i,j) = integradeOverRefTriangle(pol::phi[i]*pol::phi[j]);

  return hatM;
}

/**
 * Assembly global Mass-Matrix M
 */
inline Matrix assemblyM (Matrix const & hatM, GridOnSquer const & mesh, int polynomialDegree)
{
  Matrix M (pol::polynomialGad[polynomialDegree] * mesh.getSize());
  int t_id = 0;
  for (Triangle const & t : mesh)
  {
    int offset = t_id * pol::polynomialGad[polynomialDegree];

    for (int i=0; i<pol::polynomialGad[polynomialDegree]; ++i)
      for (int j=0; j<pol::polynomialGad[polynomialDegree]; ++j)
        M(offset+i,offset+j) = t.getArea() * hatM(i,j);

    ++t_id;
  }

  return M;
}

inline Matrix assemblyM (GridOnSquer const & mesh, int polynomialDegree)
{
  return assemblyM(assemblyHatM(polynomialDegree), mesh, polynomialDegree);
}


/**
 * Generate local tensor representing hatG_l on ref triangle
 */
inline std::vector<Tensor> assemblyHatG (int const polynomialDegree)
{
  Tensor hatG1 (pol::polynomialGad[polynomialDegree]);
  Tensor hatG2 (pol::polynomialGad[polynomialDegree]);
  for (int i=0; i<pol::polynomialGad[polynomialDegree]; ++i)
    for (int j=0; j<pol::polynomialGad[polynomialDegree]; ++j)
      for (int z=0; z<pol::polynomialGad[polynomialDegree]; ++z)
        {
          hatG1(i,j,z) = integradeOverRefTriangle(pol::dXphi[i]*pol::phi[j]*pol::phi[z]);
          hatG2(i,j,z) = integradeOverRefTriangle(pol::dYphi[i]*pol::phi[j]*pol::phi[z]);
        }
  return {hatG1, hatG2};
}

/**
 * Assembly global matrix G
 */
inline Matrix assemblyG (std::vector<Tensor> const & hatG,
                         GridOnSquer const & mesh,
                         std::vector< std::vector<double> > const & u,
                         int polynomialDegree)
{
  Matrix G (pol::polynomialGad[polynomialDegree] * mesh.getSize());
  int t_id = 0;
  for (Triangle const & t : mesh)
    {
      int offset = t_id*pol::polynomialGad[polynomialDegree];
      Jakobian B = t.getJakobian();

      for (int i=0; i<pol::polynomialGad[polynomialDegree]; ++i)
        for (int j=0; j<pol::polynomialGad[polynomialDegree]; ++j)
          {
            // G = -G1*u1 -G2*u2
            // with G1 = sum_l {  B_{2,2}*hatG1 - B_{2,1}*hatG2 }
            // and  G2 = sum_l { -B_{1,2}*hatG1 + B_{1,1}*harG2 }
            double g1 = 0.;
            double g2 = 0.;
            for (int l=0; l<pol::polynomialGad[polynomialDegree]; ++l)
              {
                g1 +=  B[1][1]*hatG[0](i,j,l) - B[1][0]*hatG[1](i,j,l);
                g2 += -B[0][1]*hatG[0](i,j,l) + B[0][0]*hatG[1](i,j,l);
              }
            G(offset+i,offset+j) =  -g1*u[0][t_id] - g2*u[1][t_id];
          }
      ++t_id;
    }

  return G;
}

inline Matrix assemblyG (GridOnSquer const & mesh,
                         std::vector< std::vector<double> > const & u,
                         int const polynomialDegree)
{
  return assemblyG(assemblyHatG(polynomialDegree), mesh, u, polynomialDegree);
}






#endif
