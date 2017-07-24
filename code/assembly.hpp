#ifndef ASSEMBLY_HPP
#define ASSEMBLY_HPP

#include "monmomials_and_basefunctions.hpp"
#include "DataTypes.hpp"
#include "Grid.hpp"

#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <iostream> //TODO debuging




/******** Mass Matrix M *******************************************************/

/**
 * Generate Mass-Matrix hatM for reference element t^
 */
inline BlockMatrix assemblyHatM (unsigned int polynomialDegree)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);

  BlockMatrix hatM (dof);

  for (unsigned int i=0; i<dof; ++i)
    for (unsigned int j=0; j<dof; ++j)
      hatM(i,j) = integradeOverRefTriangle(pol::phi[i] * pol::phi[j]);

  return hatM;
}

/**
 * Calculate local Mass-Matrix M_k for element t_k
 */
inline BlockMatrix assemblyLocalM_k (BlockMatrix hatM,
                                     double area_k)
{
  for (unsigned int i=0; i<hatM.getN(); ++i)
    for (unsigned int j=0; j<hatM.getM(); ++j)
      hatM(i,j) *= 2. * area_k;

  return hatM;
}

inline void assemblyM (UniqueSquareGrid & mesh, BlockMatrix const & hatM)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
          Triangle & l = mesh.getLower(row, col);
          l.M() = assemblyLocalM_k(hatM, l.getArea());
          Triangle & u = mesh.getUpper(row, col);
          u.M() = assemblyLocalM_k(hatM, u.getArea());
      }
}

/******** Volume integral Matrix part G ***************************************/

/**
 * Generate tensors [hatG_1, hatG_2] for the reference element t^
 */
inline std::vector<Tensor> assemblyHatG (unsigned int polynomialDegree)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);

  Tensor hatG1 (dof);
  Tensor hatG2 (dof);
  for (unsigned int i=0; i<dof; ++i)
    for (unsigned int j=0; j<dof; ++j)
      for (unsigned int z=0; z<dof; ++z)
        {
          hatG1(i,j,z) = integradeOverRefTriangle(pol::dXphi[i] * pol::phi[j] * pol::phi[z]);
          hatG2(i,j,z) = integradeOverRefTriangle(pol::dYphi[i] * pol::phi[j] * pol::phi[z]);
        }

  return {hatG1, hatG2};
}

/**
 * Calculate local Matrix G_k for element t_k
 */
inline BlockMatrix assemblyLocalG_k (std::vector<Tensor> const & hatG,
                                     std::vector<double> const & u1_k,
                                     std::vector<double> const & u2_k,
                                     Jakobian const & B_k)
{
  unsigned int dof = hatG[1].getN();

  BlockMatrix G_k (dof);

  for (unsigned int i=0; i<dof; ++i)
    for (unsigned int j=0; j<dof; ++j)
      {
        // G = G1*u1 + G2*u2
        // with G1 = sum_l {  B_{2,2}*hatG1_l - B_{2,1}*hatG2_l }
        // and  G2 = sum_l { -B_{1,2}*hatG1_l + B_{1,1}*hatG2_l }
        double g1 = 0.;
        double g2 = 0.;
        for (unsigned int l=0; l<dof; ++l)
          {
            g1 += u1_k[l] * (  B_k[1][1]*hatG[0](i,j,l) - B_k[1][0]*hatG[1](i,j,l) );
            g2 += u2_k[l] * ( -B_k[0][1]*hatG[0](i,j,l) + B_k[0][0]*hatG[1](i,j,l) );
          }
        G_k(i,j) = g1 + g2;
      }

  return G_k;
}

inline void assemblyG (UniqueSquareGrid & mesh,
                       std::vector<Tensor> const & hatG)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);
        l.G() = assemblyLocalG_k(hatG, l.U1(), l.U2(), l.getJakobian());
        u.G() = assemblyLocalG_k(hatG, u.U1(), u.U2(), u.getJakobian());
      }
}


/******** Edge integral Matrix part E *****************************************/

/**
 * Generate matrices[hatE_1, hatE_2, hatE_3] for the reference element t^
 */
inline std::vector<BlockMatrix> assemblyHatE (unsigned int polynomialDegree)
{
  unsigned int dof1D      = numberOf1DBasefunctions(polynomialDegree);
  unsigned int dof2D      = numberOf2DBasefunctions(polynomialDegree);
  std::vector<T_linear> T = getLinearTrasformationToRefEdge(polynomialDegree);

  BlockMatrix E1(dof2D, dof1D);
  BlockMatrix E2(dof2D, dof1D);
  BlockMatrix E3(dof2D, dof1D);

  // E = (T*B)*b^t
  for (unsigned int i=0; i<dof2D; ++i)
    {
      Polynomial1D tmp1 (polynomialDegree);
      Polynomial1D tmp2 (polynomialDegree);
      Polynomial1D tmp3 (polynomialDegree);

      // MatRow-Vec-Mul tmp=T_i*B^
      for (unsigned int ii=0; ii<dof1D; ++ii)
        {
          tmp1 += T[0][i][ii] * pol::phi1D[ii];
          tmp2 += T[1][i][ii] * pol::phi1D[ii];
          tmp3 += T[2][i][ii] * pol::phi1D[ii];
        }

      // Vec-VecTrans-Mul E=tmp*B^t
      for (unsigned int j=0; j<dof1D; ++j)
        {
          E1(i,j) = integradeOverRefEdge(tmp1 * pol::phi1D[j]);
          E2(i,j) = integradeOverRefEdge(tmp2 * pol::phi1D[j]);
          E3(i,j) = integradeOverRefEdge(tmp3 * pol::phi1D[j]);
        }
    }

  return {E1, E2, E3};
}

/**
 * Calculate local matrices [E1_k, E2_k, E3_k] for element t_k
 */
inline std::vector<BlockMatrix> assemblyLocalE_k (std::vector<BlockMatrix> hatE,
                                                  std::vector<double> const & length_k)
{
  unsigned int dof1D = hatE[0].getM();
  unsigned int dof2D = hatE[0].getN();

  for (unsigned int i=0; i<dof2D; ++i)
    for (unsigned int j=0; j<dof1D; ++j)
      for (int e=0; e<3; ++e)
        hatE[e](i, j) *= length_k[e];

  return hatE;
}

inline void assemblyE (UniqueSquareGrid & mesh,
                       std::vector<BlockMatrix> const & hatE)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
    {
      Triangle & l = mesh.getLower(row, col);
      Triangle & u = mesh.getUpper(row, col);

      auto E_l = assemblyLocalE_k(hatE, {l.getLengthA(), l.getLengthB(), l.getLengthC()});
      l.E_a() = E_l[0];
      l.E_b() = E_l[1];
      l.E_c() = E_l[2];

      auto E_u = assemblyLocalE_k(hatE, {u.getLengthA(), u.getLengthB(), u.getLengthC()});
      u.E_a() = E_u[0];
      u.E_b() = E_u[1];
      u.E_c() = E_u[2];
    }
}

/******** Vector representing the Flux F **************************************/

typedef std::function<std::vector<double>(unsigned int,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          Vector const &,
                                          T_linear const &,
                                          T_linear const &,
                                          std::vector<double> const &)> RiemanSolver;

/**
 * Calculate local Riemanflues [F_k1, F_k2, F_k3] on edges of element t_k in edge coordinate system
 */
inline void assemblyFr(UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       RiemanSolver const & riemanSolver,
                       std::vector<T_linear> const & hatT,
                       std::vector<double> const & hatI)
{
  for ( unsigned int row=0; row<mesh.getRows(); ++row)
    for ( unsigned int col=0; col<mesh.getColumns(); ++col)
    {
      { // Lower Triangle
        Triangle & k   = mesh.getLower(row, col);
        Triangle & n_a = mesh.getUpper(row, col);
        Triangle & n_b = mesh.getUpper(row, col-1);
        Triangle & n_c = mesh.getUpper(row-1, col);

        k.F_a() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_a.F1(), n_a.F2(),
                               k.U1(), k.U2(), n_a.U1(), n_a.U2(), k.getNormalA(),
                               hatT[0],  hatT[1], hatI);
        k.F_b() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_b.F1(), n_b.F2(),
                               k.U1(), k.U2(), n_b.U1(), n_b.U2(), k.getNormalB(),
                               hatT[1],  hatT[2], hatI);
        k.F_c() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_c.F1(), n_c.F2(),
                               k.U1(), k.U2(), n_c.U1(), n_c.U2(), k.getNormalC(),
                               hatT[2],  hatT[0], hatI);
      }

      { // Upper Triangle
        Triangle & k   = mesh.getUpper(row, col);
        Triangle & n_a = mesh.getLower(row+1, col);
        Triangle & n_b = mesh.getLower(row, col);
        Triangle & n_c = mesh.getLower(row, col+1);

        k.F_a() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_a.F1(), n_a.F2(),
                               k.U1(), k.U2(), n_a.U1(), n_a.U2(), k.getNormalA(),
                               hatT[0],  hatT[2], hatI);
        k.F_b() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_b.F1(), n_b.F2(),
                               k.U1(), k.U2(), n_b.U1(), n_b.U2(), k.getNormalB(),
                               hatT[1],  hatT[0], hatI);
        k.F_c() = riemanSolver(polynomialDegree, k.F1(), k.F2(), n_c.F1(), n_c.F2(),
                               k.U1(), k.U2(), n_c.U1(), n_c.U2(), k.getNormalC(),
                               hatT[2],  hatT[1], hatI);
      }
    }
}

/**
 * Riemansolver: UP_WINDING
 *
 * Return the solution of the RiemanProblem for the Flux F_ke on the elemnt t_e and the edge e
 * multiplied with the edge normal n in edge-coordinate-system
 *
 * Fr_ke = upwinding(F_k,F_n) * n_ke
 */
inline std::vector<double> riemanSolver_UpWinding (unsigned int polynomialDegree,
                                                   std::vector<double> const & F1_k,
                                                   std::vector<double> const & F2_k,
                                                   std::vector<double> const & F1_n,
                                                   std::vector<double> const & F2_n,
                                                   std::vector<double> const & U1_k,
                                                   std::vector<double> const & U2_k,
                                                   std::vector<double> const & U1_n,
                                                   std::vector<double> const & U2_n,
                                                   Vector const & n_ke,
                                                   T_linear const & T_ke,
                                                   T_linear const & T_ne,
                                                   std::vector<double> const & hatI)
{
  std::vector<double> Frn_e;


  // normal u at edge from both sides
  std::vector<double> ubar1_k = convertToEdge(T_ke, U1_k, polynomialDegree);
  std::vector<double> ubar2_k = convertToEdge(T_ke, U2_k, polynomialDegree);
  std::vector<double> ubar1_n = convertToEdge(T_ne, U1_n, polynomialDegree);
  std::vector<double> ubar2_n = convertToEdge(T_ne, U2_n, polynomialDegree);

  Polynomial1D uNormal_k = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_k + n_ke.y*ubar2_k);
  Polynomial1D uNormal_n = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_n + n_ke.y*ubar2_n);

  double u_k = integradeOverRefEdge(uNormal_k);
  double u_n = integradeOverRefEdge(uNormal_n);

  // up_wind direction
  if ( 0 >= u_k+u_n)
    { // from elemnt k to n
      auto Fr1_e = convertToEdge(T_ke, F1_k, polynomialDegree);
      auto Fr2_e = convertToEdge(T_ke, F2_k, polynomialDegree);

      Frn_e = Fr1_e*n_ke.x + Fr2_e*n_ke.y;
    }
  else
    { // from elemnt n to k
      auto Fr1_e = convertToEdge(T_ne, F1_n, polynomialDegree);
      auto Fr2_e = convertToEdge(T_ne, F2_n, polynomialDegree);

      Frn_e = elementWiseMul(hatI, Fr1_e*n_ke.x + Fr2_e*n_ke.y);
    }

  return Frn_e;
}


/******** Vector representing the Matrix hatI **********************************/

/**
 * Generate Matrix hatI reprsentig the differences in left and right hand site coordinate system
 * as hatI is a diagonal Matrix in OUR case(!) this is just a vector with the diagonal elements
 */
std::vector<double> getHatI (int const polynomialDegree)
{
  unsigned int dof = numberOf1DBasefunctions(polynomialDegree);

  std::vector<double> hatI;
  int sig = -1;

  std::generate_n(std::back_inserter(hatI), dof, [&sig]{sig*=-1; return sig;});

  return hatI;
}


/******** RHS Vector L ********************************************************/

/**
 * Assembly rhs Vector L_k = M_k*f_k for element t_k
 */
inline std::vector<double> assemblyLocalL_k (unsigned int polynomialDegree,
                                             std::function<double(double,double)> const & f,
                                             BlockMatrix const & M_k,
                                             Jakobian const & B_k,
                                             Point const & A_k)
{
  return M_k * l2Projection(polynomialDegree, f, B_k, A_k);
}

inline void assamblyL (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       std::function<double(double,double)> const & f)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);

        l.L() = assemblyLocalL_k(polynomialDegree, f, l.M(), l.getJakobian(), l.getA());
        u.L() = assemblyLocalL_k(polynomialDegree, f, u.M(), u.getJakobian(), u.getA());
      }
}


/******** Vectors representing the velocity field U_1,2 *************************/

/**
 *
 */
inline void assamblyU (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       std::function<double(double,double)> const & u1,
                       std::function<double(double,double)> const & u2)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        l.U1() = l2Projection(polynomialDegree, u1, l.getJakobian(), l.getA());
        l.U2() = l2Projection(polynomialDegree, u2, l.getJakobian(), l.getA());

        Triangle & u = mesh.getUpper(row, col);
        u.U1() = l2Projection(polynomialDegree, u1, u.getJakobian(), u.getA());
        u.U2() = l2Projection(polynomialDegree, u2, u.getJakobian(), u.getA());
    }
}


/******** Vectors representing the Conscentration C *****************************/

/**
 *
 */
inline void assamblyC (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       std::function<double(double,double)> const & c)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);

        l.C() = l2Projection(polynomialDegree, c, l.getJakobian(), l.getA());
        u.C() = l2Projection(polynomialDegree, c, u.getJakobian(), u.getA());
      }
}

/******** Vectors representing the Flux F ***************************************/

/**
 *
 */
inline void assamblyF (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree, unsigned int targetPolynomialDegree,
                       std::function<std::vector< std::vector<double> >
                         (Triangle const &, unsigned int, unsigned int)> f)
{
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        auto F_l = f(l, targetPolynomialDegree, polynomialDegree);
        l.F1() = F_l[0];
        l.F2() = F_l[1];

        Triangle & u = mesh.getUpper(row, col);
        auto F_u = f(u, targetPolynomialDegree, polynomialDegree);
        u.F1() = F_u[0];
        u.F2() = F_u[1];
      }
}

inline std::vector< std::vector<double> > assamblyLocalLinearF (Triangle const & t_k,
                                                                unsigned int targetPolynomialDegree,
                                                                unsigned int polynomialDegree)
{
  //TODO hier mit etwas rumspeilen

  auto c_pol ( reconstructFunction2D(polynomialDegree, t_k.C()) );
  auto u1_pol ( reconstructFunction2D(polynomialDegree, t_k.U1()) );
  auto u2_pol ( reconstructFunction2D(polynomialDegree, t_k.U2()) );

  auto F1_pol = c_pol * u1_pol;
  auto F2_pol = c_pol * u2_pol;

  auto F1 (l2Projection(targetPolynomialDegree, F1_pol));
  auto F2 (l2Projection(targetPolynomialDegree, F2_pol));

  return {F1, F2};
}

//TODO Rand machen für U[1,2] und F[1,2]
















//TODO print methode die globale maxtix prinded


#endif
