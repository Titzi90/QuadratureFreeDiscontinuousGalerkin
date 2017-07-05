#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include "monmomials_and_basefunctions.hpp"
#include "squareGrid.hpp"
#include "Matrix.hpp"

#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>

#include <iostream> //TODO debuging





/**
 * Generate local Mass-Matrix for ref triangle hatM
 */
inline Matrix assemblyHatM (int polynomialDegree)
{
  Matrix hatM (pol::numberOf2DBasefunctions(polynomialDegree));
  for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
    for (int j=0; j<pol::numberOf2DBasefunctions(polynomialDegree); ++j)
      hatM(i,j) = integradeOverRefTriangle(pol::phi[i]*pol::phi[j]);

  return hatM;
}

/**
 * Assembly global Mass-Matrix M
 */
inline Matrix assemblyM (Matrix const & hatM, GridOnSquer const & mesh, int polynomialDegree)
{
  Matrix M (pol::numberOf2DBasefunctions(polynomialDegree) * mesh.getSize());
  int t_id = 0;
  for (Triangle const & t : mesh)
  {
    int offset = t_id * pol::numberOf2DBasefunctions(polynomialDegree);

    for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
      for (int j=0; j<pol::numberOf2DBasefunctions(polynomialDegree); ++j)
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
  Tensor hatG1 (pol::numberOf2DBasefunctions(polynomialDegree));
  Tensor hatG2 (pol::numberOf2DBasefunctions(polynomialDegree));
  for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
    for (int j=0; j<pol::numberOf2DBasefunctions(polynomialDegree); ++j)
      for (int z=0; z<pol::numberOf2DBasefunctions(polynomialDegree); ++z)
        {
          hatG1(i,j,z) = integradeOverRefTriangle(pol::dXphi[i]*pol::phi[j]*pol::phi[z]);
          hatG2(i,j,z) = integradeOverRefTriangle(pol::dYphi[i]*pol::phi[j]*pol::phi[z]);
        }
  return {hatG1, hatG2};
}

/**
 * Assembly global matrix G
 * G = G1 + G2
 * including transformation
 */
inline Matrix assemblyG (std::vector<Tensor> const & hatG,
                         GridOnSquer const & mesh,
                         std::vector<Coefficient> const & u,
                         int polynomialDegree)
{
  Matrix G (pol::numberOf2DBasefunctions(polynomialDegree) * mesh.getSize());
  int t_id = 0;
  for (Triangle const & t : mesh)
    {
      int offset = t_id*pol::numberOf2DBasefunctions(polynomialDegree); //TODO iterator und offset müssen zussammen passen!!!
      Jakobian B = t.getJakobian();

      for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
        for (int j=0; j<pol::numberOf2DBasefunctions(polynomialDegree); ++j)
          {
            // G = G1*u1 +G2*u2
            // with G1 = sum_l {  B_{2,2}*hatG1_l - B_{2,1}*hatG2_l }
            // and  G2 = sum_l { -B_{1,2}*hatG1_l + B_{1,1}*hatG2_l }
            double g1 = 0.;
            double g2 = 0.;
            for (int l=0; l<pol::numberOf2DBasefunctions(polynomialDegree); ++l)
              {
                g1 += u[0].get(t_id, l) * (  B[1][1]*hatG[0](i,j,l) - B[1][0]*hatG[1](i,j,l) );
                g2 += u[1].get(t_id, l) * ( -B[0][1]*hatG[0](i,j,l) + B[0][0]*hatG[1](i,j,l) );
              }
            G(offset+i,offset+j) =  g1 + g2;
          }
      ++t_id;
    }

  return G;
}

inline Matrix assemblyG (GridOnSquer const & mesh,
                         std::vector<Coefficient> const & u,
                         int const polynomialDegree)
{
  return assemblyG(assemblyHatG(polynomialDegree), mesh, u, polynomialDegree);
}


/**
 * gernate matrix hat E
 * hatE_i = int_refEdge (T_i*barB*barB^t)
 */
inline std::vector<Matrix> assemblyHatE (int const polynomialDegree)
{
  std::vector<T_linear> T = getLinearTrasformationToRefEdge(polynomialDegree);

  Matrix E1(pol::numberOf2DBasefunctions(polynomialDegree), pol::numberOf1DBasefunctions(polynomialDegree));
  Matrix E2(pol::numberOf2DBasefunctions(polynomialDegree), pol::numberOf1DBasefunctions(polynomialDegree));
  Matrix E3(pol::numberOf2DBasefunctions(polynomialDegree), pol::numberOf1DBasefunctions(polynomialDegree));

  for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
    {
      Polynomial1D tmp1 (polynomialDegree);
      Polynomial1D tmp2 (polynomialDegree);
      Polynomial1D tmp3 (polynomialDegree);
      // MatVecMul tmp=T*B
      for (int ii=0; ii<pol::numberOf1DBasefunctions(polynomialDegree); ++ii)
        {
          tmp1 += T[0][i][ii] * pol::phi1D[ii];
          tmp2 += T[1][i][ii] * pol::phi1D[ii];
          tmp3 += T[2][i][ii] * pol::phi1D[ii];
        }

      // MatMatMul tmp*B^t
      for (int j=0; j<pol::numberOf1DBasefunctions(polynomialDegree); ++j)
        {
          E1(i,j) = integradeOverRefEdge(tmp1 * pol::phi1D[j]);
          E2(i,j) = integradeOverRefEdge(tmp2 * pol::phi1D[j]);
          E3(i,j) = integradeOverRefEdge(tmp3 * pol::phi1D[j]);
        }
    }
  return {E1, E2, E3};
}

/**
 * gernate matrix E
 * E = transformation_to_refTriangle * length_egdge * hatE
 */
inline std::vector<Matrix> assemblyE (std::vector<Matrix> const & hatE,
                                      GridOnSquer const & mesh,
                                      int polynomialDegree)
{
  std::vector<Matrix> E (3, Matrix(pol::numberOf2DBasefunctions(polynomialDegree)* mesh.getSize(),
                                   pol::numberOf1DBasefunctions(polynomialDegree)* mesh.getSize()));

  int t_id=0;
  for (Triangle const & t : mesh)
    {
      int offset_1d = t_id * pol::numberOf1DBasefunctions(polynomialDegree);
      int offset_2d = t_id * pol::numberOf2DBasefunctions(polynomialDegree);
      double length[3] = { t.getLengthA(), t.getLengthB(), t.getLengthC()};

      for (int i=0; i<pol::numberOf2DBasefunctions(polynomialDegree); ++i)
        for (int j=0; j<pol::numberOf1DBasefunctions(polynomialDegree); ++j)
          for (int e=0; e<3; ++e)
            E[e](offset_2d+i, offset_1d+j) = length[e] * hatE[e](i,j);

      ++t_id;
    }

  return E;
}


/**
 * convert coefficient from ref. Element to ref. edge
 * V_ = T^T * V
 * TODO als scalere operation nicht benötigt -> nur Vector verion
 */
//TODO für alle vectoren alle drei dreicke erzeugen
//oder individuell nur für die Teile wo ich es brauche
inline std::vector<double> convertToEdge (T_linear const & T,
                                          Coefficient const & v,
                                          int elementID,
                                          int polynomialDegree)
{
  int bases1D = pol::numberOf1DBasefunctions(polynomialDegree);
  int bases2D = pol::numberOf2DBasefunctions(polynomialDegree);

  std::vector<double> v_ (bases1D, 0.);

  for (int i=0; i<bases1D; ++i)
    for (int j=0; j<bases2D; ++j)
      {
        v_[i] += T[j][i] * v[elementID*bases2D + j];
      }

  return v_;
}
inline std::vector< std::vector<double> > convertToEdge (T_linear const & T,
                                                         std::vector<Coefficient> const & v,
                                                         int elementID,
                                                         int polynomialDegree)
{
  //TODO mit normalen multiplizieren
  return {convertToEdge(T, v[0], elementID, polynomialDegree),
          convertToEdge(T, v[1], elementID, polynomialDegree)};
}

/**
 * Generate Matrix hatI reprsentig the differences in left and right hand site coordinate system
 * as hatI is a diagonal Matrix in OUR case(!) this is just a vector with the diagonal elements
 */
std::vector<double> getHatI (int const polynomialDegree)
{
  int sig = -1;
  std::vector<double> hatI(pol::numberOf1DBasefunctions(polynomialDegree));
  std::generate(hatI.begin(), hatI.end(), [&sig]{ sig*=-1; return sig;});
  return hatI;
}

/**
 * Return the solution of the RiemanProblem for the Flux on the edges multiplied with the edge normal
 * in edge coordinate system
 * using the UpWinding method
 * F_r_local = upwinding(F_l,F_r) * n
 */
inline std::vector< std::vector<double> > getRiemenFlux_upWinding (GridOnSquer const & mesh,
                                                                   unsigned int triangleRow,
                                                                   unsigned int triangleCol,
                                                                   unsigned int trianglePos,
                                                                   std::vector<Coefficient> const & F, // Fluss in Refelemnt Basen expansion
                                                                   std::vector<Coefficient> const & U,
                                                                   std::vector<double> const & hatI,
                                                                   std::vector<Coefficient> const & boundaryFlux,
                                                                   int const polynomialDegree)
{
  unsigned int triangleID = mesh.getID(triangleRow, triangleCol, trianglePos);

  //Nachbar IDs und kanten Typ:
  //TODO duch kleveres bennenen der kanten if sparen (a<->a, ...)
  unsigned int triangleID_neigbour [3];
  unsigned int edgeType_neigbour [3];
  if (0 == trianglePos) // lower triangle
    {
      triangleID_neigbour[0] = mesh.getID(triangleRow  , triangleCol  , 1);
      triangleID_neigbour[1] = mesh.getID(triangleRow  , triangleCol-1, 1);
      triangleID_neigbour[2] = mesh.getID(triangleRow+1, triangleCol  , 1);

      edgeType_neigbour[0] = 2; //=b
      edgeType_neigbour[1] = 3; //=c
      edgeType_neigbour[2] = 1; //=a
    }
  else
    {
      triangleID_neigbour[0]  = mesh.getID(triangleRow-1, triangleCol  , 0);
      triangleID_neigbour[1]  = mesh.getID(triangleRow  , triangleCol  , 0);
      triangleID_neigbour[2]  = mesh.getID(triangleRow  , triangleCol+1, 0);

      edgeType_neigbour[0] = 3; //=c
      edgeType_neigbour[1] = 1; //=a
      edgeType_neigbour[2] = 2; //=b
    }

  Triangle const triangle = mesh.get(triangleRow, triangleCol, trianglePos);
  std::vector<T_linear> const T = getLinearTrasformationToRefEdge(polynomialDegree);
  Vector const normale[3] = {triangle.getNormalA(), triangle.getNormalB(), triangle.getNormalC()};

  std::vector< std::vector<double> > f_r(3, std::vector<double>(pol::numberOf1DBasefunctions(polynomialDegree)));

  for (unsigned int e=0; e<3; ++e)
    {
      // u at edge
      //TODO -> eigene funktion
      std::vector< std::vector<double> > U_bar = convertToEdge(T[e], U, triangleID, polynomialDegree);
      Polynomial1D u1_bar = std::inner_product(U_bar[0].begin(), U_bar[0].end(),
                                               std::begin(pol::phi1D),
                                               Polynomial1D(polynomialDegree));
      Polynomial1D u2_bar = std::inner_product(U_bar[1].begin(), U_bar[1].end(),
                                               std::begin(pol::phi1D),
                                               Polynomial1D(polynomialDegree));
      double u1_int = integradeOverRefEdge(u1_bar);
      double u2_int = integradeOverRefEdge(u2_bar);

      Vector u (u1_int, u2_int); // TODO hier eigentlich mittelwerrt von beiden seiten

      if ( dot(u, normale[e]) >= 0 )
        {
          std::vector< std::vector<double> > vecF = convertToEdge(T[e], F, triangleID, polynomialDegree);
          std::transform(vecF.begin(), vecF.end(), f_r[e].begin(),
                         [&normale, e](std::vector<double> & f){ return f[0]*normale[e].x + f[1]*normale[e].y; });
        }
      else
        {
          // RAND
          if ( (e==0 && trianglePos==1 && triangleRow==0) ||  // Top
               (e==1 && trianglePos==0 && triangleCol==0) ||  // Left
               (e==2 && trianglePos==0 && triangleRow==mesh.getRows()-1)  ||  // Bottom
               (e==2 && trianglePos==1 && triangleCol==mesh.getColums()-1))   // Right
            {
              std::vector < std::vector<double> > vecF = convertToEdge(T[e],boundaryFlux, triangle, e);
              std::transform(vecF.begin(), vecF.end(), f_r[e].begin(),
                             [&normale, e](std::vector<double> & f)
                             { return f[0]*normale[e].x + f[1]*normale[e].y; });
            }
          else  // nachbarzelle
            {
              std::vector< std::vector<double> > vecF = convertToEdge(T[edgeType_neigbour[e]],
                                                                      F,
                                                                      triangleID_neigbour[e],
                                                                      polynomialDegree);
              std::transform(vecF.begin(), vecF.end(), hatI.begin(), f_r[e].begin(),
                             [&normale,e](std::vector<double> &f, double i)
                             {return i * (f[0]*normale[e].x + f[1]*normale[e].y); });

            }
        }
    }

  return f_r;
}


/**
 * Riemanflux auf Kanten 1,2,3 von Referenzdreieck in Edges kordinaten
 * TODO diese funktion algemein machen, das sie eine funktion für riemensolver übergeeben bekommt
 */
inline std::vector<Coefficient> assemblyRiemanFlux_UpWinding (GridOnSquer const & mesh,
                                                              std::vector<Coefficient> const & U,
                                                              std::vector<Coefficient> const & F, // Fluss in Refelemnt Basen expantion
                                                              std::vector<double> const & hatI,
                                                              std::vector<Coefficient> const & boundaryU,
                                                              int const polynomialDegree)

{

  unsigned int numberBaseFunctions = pol::numberOf1DBasefunctions(polynomialDegree);
  std::vector<Coefficient> F_r (3, Coefficient(numberBaseFunctions*mesh.getSize(), numberBaseFunctions));


  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColums(); ++col)
      for (unsigned int pos=0; pos<2; ++pos)
        {
          unsigned int t_id = mesh.getID(row, col, pos);
          unsigned int offset = t_id * numberBaseFunctions;
          std::vector< std::vector<double> > F_r_local = getRiemenFlux_upWinding(mesh,
                                                                                 row, col, pos,
                                                                                 F, U, hatI, boundaryU,
                                                                                 polynomialDegree);

          for (unsigned int i=0; i<numberBaseFunctions; ++i)
            for (unsigned int e=0; e<3; ++e)
              F_r[e][offset + i] = F_r_local[e][i];
        }

  return F_r;
}



// Assembly rhs Vector L = M*F
inline Coefficient assemblyL (Matrix const & M, Coefficient const & F) { return M*F; }







#endif
