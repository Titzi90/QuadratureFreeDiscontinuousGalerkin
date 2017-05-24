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
        M(offset+i,offset+j) = 2. * t.getArea() * hatM(i,j);

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
 * Assembly global matrix A
 * A = -G1 -G2 +R
 */
inline Matrix assemblyA (std::vector<Tensor> const & hatG,
                         GridOnSquer const & mesh,
                         std::vector<Coefficient> const & u,
                         int polynomialDegree)
{
  Matrix A (pol::polynomialGad[polynomialDegree] * mesh.getSize());
  int t_id = 0;
  for (Triangle const & t : mesh)
    {
      int offset = t_id*pol::polynomialGad[polynomialDegree];
      Jakobian B = t.getJakobian();

      for (int i=0; i<pol::polynomialGad[polynomialDegree]; ++i)
        for (int j=0; j<pol::polynomialGad[polynomialDegree]; ++j)
          {
            // G = -G1*u1 -G2*u2
            // with G1 = sum_l {  B_{2,2}*hatG1_l - B_{2,1}*hatG2_l }
            // and  G2 = sum_l { -B_{1,2}*hatG1_l + B_{1,1}*hatG2_l }
            double g1 = 0.;
            double g2 = 0.;
            double r =0.; //TODO
            for (int l=0; l<pol::polynomialGad[polynomialDegree]; ++l)
              {
                g1 += u[0].get(t_id, l) * (  B[1][1]*hatG[0](i,j,l) - B[1][0]*hatG[1](i,j,l) );
                g2 += u[1].get(t_id, l) * ( -B[0][1]*hatG[0](i,j,l) + B[0][0]*hatG[1](i,j,l) );
              }
            A(offset+i,offset+j) =  -g1 - g2 +r;
          }
      ++t_id;
    }

  return A;
}

inline Matrix assemblyA (GridOnSquer const & mesh,
                         std::vector<Coefficient> const & u,
                         int const polynomialDegree)
{
  return assemblyA(assemblyHatG(polynomialDegree), mesh, u, polynomialDegree);
}

inline Coefficient assemblyL (Matrix const & M, Coefficient const & F) { return M*F; }







#endif
