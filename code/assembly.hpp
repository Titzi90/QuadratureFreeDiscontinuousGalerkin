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

// include likwid
extern "C"
{
#ifdef USE_LIKWID
#include <likwid.h>
  // This block enables compilation of the code with and without LIKWID in place
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif
}


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
  LIKWID_MARKER_START("ASSEMBLY_M_PRECOMPUTED");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
          Triangle & l = mesh.getLower(row, col);
          l.M() = assemblyLocalM_k(hatM, l.getArea());
          Triangle & u = mesh.getUpper(row, col);
          u.M() = assemblyLocalM_k(hatM, u.getArea());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_M_PRECOMPUTED");
}





inline BlockMatrix assemblyLocalM_kQuadFree(unsigned int polynomialDegree, double area_k)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);
  BlockMatrix M_k (dof);

  for (unsigned int i=0; i<dof; ++i)
    for (unsigned int j=0; j<dof; ++j)
      M_k(i,j) = 2. * area_k * integradeOverRefTriangle(pol::phi[i]*pol::phi[j]);

  return M_k;
}


inline void assemblyMquadFree (UniqueSquareGrid & mesh, unsigned int polynomialDegree)
{
  LIKWID_MARKER_START("ASSEMBLY_M_QUADFREE");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        l.M() = assemblyLocalM_kQuadFree(polynomialDegree, l.getArea());
        Triangle & u = mesh.getUpper(row, col);
        u.M() = assemblyLocalM_kQuadFree(polynomialDegree, u.getArea());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_M_QUADFREE");
}


inline BlockMatrix assemblyLocalM_kGaus(unsigned int polynomialDegree, double area_k)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);
  BlockMatrix M_k (dof);

  for (unsigned int i=0; i<dof; ++i)
    for (unsigned int j=0; j<dof; ++j)
      M_k(i,j) = 2. * area_k * integradeOverRefTriangle_gaus(pol::phi[i]*pol::phi[j], 2*polynomialDegree);

  return M_k;
}


inline void assemblyMGaus (UniqueSquareGrid & mesh, unsigned int polynomialDegree)
{
  LIKWID_MARKER_START("ASSEMBLY_M_GAUS");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        l.M() = assemblyLocalM_kGaus(polynomialDegree, l.getArea());
        Triangle & u = mesh.getUpper(row, col);
        u.M() = assemblyLocalM_kGaus(polynomialDegree, u.getArea());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_M_GAUS");
}






/**
 * invert mass matrix M
 * M have to be a quadratic diagonal matrix!!!
 */
BlockMatrix invertM(BlockMatrix M)
{
  for (unsigned int i=0; i<M.getM(); ++i)
    M(i,i) = 1./M(i,i);

  return M;
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
  LIKWID_MARKER_START("ASSEMBLY_G_PRECOMPUTED");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);
        l.G() = assemblyLocalG_k(hatG, l.U1(), l.U2(), l.getJakobian());
        u.G() = assemblyLocalG_k(hatG, u.U1(), u.U2(), u.getJakobian());
      }

    LIKWID_MARKER_STOP("ASSEMBLY_G_PRECOMPUTED");
}




inline BlockMatrix assemblyLocalG_kquadFree (unsigned int polynomialDegree,
                                             std::vector<double> const & u1_k,
                                             std::vector<double> const & u2_k,
                                             Jakobian const & B_k)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);

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
            double hatG1 = integradeOverRefTriangle(pol::dXphi[i] * pol::phi[j] * pol::phi[l]);
            double hatG2 = integradeOverRefTriangle(pol::dYphi[i] * pol::phi[j] * pol::phi[l]);
            g1 += u1_k[l] * ( B_k[1][1]*hatG1 - B_k[1][0]*hatG2);
            g2 += u2_k[l] * (-B_k[0][1]*hatG1 + B_k[0][0]*hatG2);
          }
        G_k(i,j) = g1 + g2;
      }

  return G_k;
}

inline void assemblyGquadFree (UniqueSquareGrid & mesh,
                               unsigned int polynomialDegree)
{
  LIKWID_MARKER_START("ASSEMBLY_G_QUADFREE");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);
        l.G() = assemblyLocalG_kquadFree(polynomialDegree, l.U1(), l.U2(), l.getJakobian());
        u.G() = assemblyLocalG_kquadFree(polynomialDegree, u.U1(), u.U2(), u.getJakobian());
      }

  LIKWID_MARKER_START("ASSEMBLY_G_QUADFREE");
}


inline BlockMatrix assemblyLocalG_kgaus (unsigned int polynomialDegree,
                                             std::vector<double> const & u1_k,
                                             std::vector<double> const & u2_k,
                                             Jakobian const & B_k)
{
  unsigned int dof = numberOf2DBasefunctions(polynomialDegree);

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
            //TODO order ist sehr ungenau
            double hatG1 = integradeOverRefTriangle_gaus(pol::dXphi[i] * pol::phi[j] * pol::phi[l],2*polynomialDegree);
            double hatG2 = integradeOverRefTriangle_gaus(pol::dYphi[i] * pol::phi[j] * pol::phi[l],2*polynomialDegree);
            g1 += u1_k[l] * ( B_k[1][1]*hatG1 - B_k[1][0]*hatG2);
            g2 += u2_k[l] * (-B_k[0][1]*hatG1 + B_k[0][0]*hatG2);
          }
        G_k(i,j) = g1 + g2;
      }

  return G_k;
}

inline void assemblyGgaus (UniqueSquareGrid & mesh,
                               unsigned int polynomialDegree)
{
  LIKWID_MARKER_START("ASSEMBLY_G_GAUS");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);
        l.G() = assemblyLocalG_kgaus(polynomialDegree, l.U1(), l.U2(), l.getJakobian());
        u.G() = assemblyLocalG_kgaus(polynomialDegree, u.U1(), u.U2(), u.getJakobian());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_G_GAUS");
}


/******** Edge integral Matrix part E *****************************************/

/**
 * Generate matrices[hatE_1, hatE_2, hatE_3] for the reference element t^
 */
inline std::vector<BlockMatrix> assemblyHatE (unsigned int polynomialDegree,
                                              unsigned int polynomialDegreeF)
{
  unsigned int dof1D      = numberOf1DBasefunctions(polynomialDegree);
  unsigned int dof2D      = numberOf2DBasefunctions(polynomialDegree);
  unsigned int dofF       = numberOf1DBasefunctions(polynomialDegreeF);
  std::vector<BlockMatrix> T = getLinearTrasformationToRefEdge(polynomialDegree);

  BlockMatrix E1(dof2D, dofF);
  BlockMatrix E2(dof2D, dofF);
  BlockMatrix E3(dof2D, dofF);

  // E = (T*B)*b^t
  for (unsigned int i=0; i<dof2D; ++i)
    {
      Polynomial1D tmp1 (polynomialDegree);
      Polynomial1D tmp2 (polynomialDegree);
      Polynomial1D tmp3 (polynomialDegree);

      // MatRow-Vec-Mul tmp=T_i*B^
      for (unsigned int ii=0; ii<dof1D; ++ii)
        {
          tmp1 += T[0](ii, i) * pol::phi1D[ii];
          tmp2 += T[1](ii, i) * pol::phi1D[ii];
          tmp3 += T[2](ii, i) * pol::phi1D[ii];
        }

      // Vec-VecTrans-Mul E=tmp*B^t
      for (unsigned int j=0; j<dofF; ++j)
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
  unsigned int N = hatE[0].getN();
  unsigned int M = hatE[0].getM();

  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<M; ++j)
      for (int e=0; e<3; ++e)
        hatE[e](i, j) *= length_k[e];

  return hatE;
}

inline void assemblyE (UniqueSquareGrid & mesh,
                       std::vector<BlockMatrix> const & hatE)
{
  LIKWID_MARKER_START("ASSEMBLY_E_PRECOMPUTED");

#pragma omp for
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

  LIKWID_MARKER_STOP("ASSEMBLY_E_PRECOMPUTED");
}





inline std::vector<BlockMatrix> assemblyLocalE_kquadFree (unsigned int polynomialDegree,
                                                          unsigned int polynomialDegreeF,
                                                          std::vector<double> const & length_k)
{
  unsigned int dof1D      = numberOf1DBasefunctions(polynomialDegree);
  unsigned int dof2D      = numberOf2DBasefunctions(polynomialDegree);
  unsigned int dofF       = numberOf1DBasefunctions(polynomialDegreeF);
  std::vector<BlockMatrix> T = getLinearTrasformationToRefEdge(polynomialDegree);

  std::vector<BlockMatrix> E(3, BlockMatrix(dof2D, dofF));

  for (unsigned int i=0; i<dof2D; ++i)
    {
      Polynomial1D tmp[3] = {Polynomial1D(polynomialDegree),
                             Polynomial1D(polynomialDegree),
                             Polynomial1D(polynomialDegree)};
      // MatRow-Vec-Mul tmp=T_i*B^
      for (unsigned int ii=0; ii<dof1D; ++ii)
        for (int e=0;e<3;++e)
            tmp[e] += T[e](ii, i) * pol::phi1D[ii];

      for (unsigned int j=0; j<dofF; ++j)
        for (int e=0; e<3; ++e)
          E[e](i, j) = length_k[e] * integradeOverRefEdge(tmp[e] * pol::phi1D[j]);
    }

  return E;
}

inline void assemblyEquadFree (UniqueSquareGrid & mesh,
                               unsigned int polynomialDegree,
                               unsigned int polynomialDegreeF)
{
  LIKWID_MARKER_START("ASSEMBLY_E_QUADFREE");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
    {
      Triangle & l = mesh.getLower(row, col);
      Triangle & u = mesh.getUpper(row, col);

      auto E_l = assemblyLocalE_kquadFree(polynomialDegree, polynomialDegreeF,
                                          {l.getLengthA(), l.getLengthB(), l.getLengthC()});
      l.E_a() = E_l[0];
      l.E_b() = E_l[1];
      l.E_c() = E_l[2];

      auto E_u = assemblyLocalE_kquadFree(polynomialDegree, polynomialDegreeF,
                                          {u.getLengthA(), u.getLengthB(), u.getLengthC()});
      u.E_a() = E_u[0];
      u.E_b() = E_u[1];
      u.E_c() = E_u[2];
    }

  LIKWID_MARKER_STOP("ASSEMBLY_E_QUADFREE");
}




/******** Vector representing the Flux F **************************************/

typedef std::function<std::vector<double>(unsigned int, unsigned int,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
                                          Vector const &,
                                          unsigned int const,
                                          unsigned int const,
                                          std::vector<double> const &)> RiemanSolver;

/**
 * Calculate local Riemanflues [F_k1, F_k2, F_k3] on edges of element t_k in edge coordinate system
 */
inline void assemblyFr(UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       unsigned int polynomialDegreeF,
                       RiemanSolver const & riemanSolver,
                       std::vector<double> const & hatI)
{
  LIKWID_MARKER_START("ASSEMBLY_Fr");

#pragma omp for
  for ( unsigned int row=0; row<mesh.getRows(); ++row)
    for ( unsigned int col=0; col<mesh.getColumns(); ++col)
    {
      { // Lower Triangle
        Triangle & k   = mesh.getLower(row, col);
        Triangle & n_a = mesh.getUpper(row, col);
        Triangle & n_b = mesh.getUpper(row, col-1);
        Triangle & n_c = mesh.getUpper(row-1, col);

        k.Fr_a() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_a(), n_a.Fn_b(),
                                k.U1(), k.U2(), n_a.U1(), n_a.U2(),
                                k.getNormalA(), 0, 1, hatI);
        k.Fr_b() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_b(), n_b.Fn_c(),
                                k.U1(), k.U2(), n_b.U1(), n_b.U2(),
                                k.getNormalB(), 1,  2, hatI);
        k.Fr_c() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_c(), n_c.Fn_a(),
                                k.U1(), k.U2(), n_c.U1(), n_c.U2(),
                                k.getNormalC(), 2,  0, hatI);
      }

      { // Upper Triangle
        Triangle & k   = mesh.getUpper(row, col);
        Triangle & n_a = mesh.getLower(row+1, col);
        Triangle & n_b = mesh.getLower(row, col);
        Triangle & n_c = mesh.getLower(row, col+1);

        k.Fr_a() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_a(), n_a.Fn_c(),
                                k.U1(), k.U2(), n_a.U1(), n_a.U2(),
                                k.getNormalA(), 0,  2, hatI);
        k.Fr_b() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_b(), n_b.Fn_a(),
                                k.U1(), k.U2(), n_b.U1(), n_b.U2(),
                                k.getNormalB(), 1, 0, hatI);
        k.Fr_c() = riemanSolver(polynomialDegree, polynomialDegreeF,
                                k.Fn_c(), n_c.Fn_b(),
                                k.U1(), k.U2(), n_c.U1(), n_c.U2(),
                                k.getNormalC(), 2, 1, hatI);
      }
    }

  LIKWID_MARKER_STOP("ASSEMBLY_Fr");
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
                                                   unsigned int polynomialDegreeF,
                                                   std::vector<double> const & Fn_k,
                                                   std::vector<double> const & Fn_n,
                                                   std::vector<double> const & U1_k,
                                                   std::vector<double> const & U2_k,
                                                   std::vector<double> const & U1_n,
                                                   std::vector<double> const & U2_n,
                                                   Vector const & n_ke,
                                                   unsigned int const edgeID_k,
                                                   unsigned int const edgeID_n,
                                                   std::vector<double> const & hatI)
{

  std::vector<double> Frn_e;

  // u on the edge from both sides
  std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
  std::vector<double> ubar1_k = T[edgeID_k] * U1_k;
  std::vector<double> ubar2_k = T[edgeID_k] * U2_k;
  std::vector<double> ubar1_n = T[edgeID_n] * U1_n;
  std::vector<double> ubar2_n = T[edgeID_n] * U2_n;

  // u normal to edge as polynomial
  Polynomial1D uNormal_k = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_k + n_ke.y*ubar2_k);
  Polynomial1D uNormal_n = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_n + n_ke.y*ubar2_n);

  double u_k = integradeOverRefEdge(uNormal_k);
  double u_n = integradeOverRefEdge(uNormal_n);

  // up_wind direction
  if ( 0 <= u_k+u_n)
    { // from elemnt k to n
      Frn_e = Fn_k;
    }
  else
    { // from elemnt n to k
      Frn_e = -1 * elementWiseMul(hatI, Fn_n);
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
                                             // BlockMatrix const & M_k,
                                             Jakobian const & B_k,
                                             Point const & A_k)
{
  // return M_k * l2Projection(polynomialDegree, f, B_k, A_k);
  return l2Projection(polynomialDegree, f, B_k, A_k);
}

inline void assamblyL (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       std::function<double(double,double)> const & f)
{
  LIKWID_MARKER_START("ASSEMBLY_L");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);

        // l.L() = assemblyLocalL_k(polynomialDegree, f, l.M(), l.getJakobian(), l.getA());
        // u.L() = assemblyLocalL_k(polynomialDegree, f, u.M(), u.getJakobian(), u.getA());
        l.L() = assemblyLocalL_k(polynomialDegree, f, l.getJakobian(), l.getA());
        u.L() = assemblyLocalL_k(polynomialDegree, f, u.getJakobian(), u.getA());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_L");
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
  LIKWID_MARKER_START("ASSEMBLY_U");

#pragma omp for
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

  LIKWID_MARKER_STOP("ASSEMBLY_U");
}


/******** Vectors representing the Conscentration C *****************************/

/**
 *
 */
inline void assamblyC (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree,
                       std::function<double(double,double)> const & c)
{
  LIKWID_MARKER_START("ASSEMBLY_C");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        Triangle & u = mesh.getUpper(row, col);

        l.C() = l2Projection(polynomialDegree, c, l.getJakobian(), l.getA());
        u.C() = l2Projection(polynomialDegree, c, u.getJakobian(), u.getA());
      }

  LIKWID_MARKER_STOP("ASSEMBLY_C");
}

/******** Vectors representing the Flux F ***************************************/

/**
 * calculate expansion vector for \bar{f}=uc on t
 */
inline void assamblyF (UniqueSquareGrid & mesh,
                       unsigned int polynomialDegree, unsigned int targetPolynomialDegree,
                       std::function<std::vector< std::vector<double> >
                         (Triangle const &, unsigned int, unsigned int)> f)
{
  LIKWID_MARKER_START("ASSEMBLY_F");

#pragma omp for
  for (unsigned int row=0; row<mesh.getRows(); ++row)
    for (unsigned int col=0; col<mesh.getColumns(); ++col)
      {
        Triangle & l = mesh.getLower(row, col);
        auto F_l = f(l, targetPolynomialDegree, polynomialDegree);
        l.Fn_a() = F_l[0];
        l.Fn_b() = F_l[1];
        l.Fn_c() = F_l[2];

        Triangle & u = mesh.getUpper(row, col);
        auto F_u = f(u, targetPolynomialDegree, polynomialDegree);
        u.Fn_a() = F_u[0];
        u.Fn_b() = F_u[1];
        u.Fn_c() = F_u[2];
      }

  LIKWID_MARKER_STOP("ASSEMBLY_F");
}

inline std::vector< std::vector<double> > assamblyLocalLinearF (Triangle const & t_k,
                                                                unsigned int targetPolynomialDegree,
                                                                unsigned int polynomialDegree)
{
  Vector const n[3] = {t_k.getNormalA(),
                       t_k.getNormalB(),
                       t_k.getNormalC()};
  std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
  BlockMatrix Te = get1DPolynomialMapping(polynomialDegree*2);
  std::vector< std::vector<double> > Fn;
  for (int e=0; e<3; ++e)
    {
      // normal*u and c on edge
      auto Un = n[e].x * t_k.U1() + n[e].y * t_k.U2();
      auto Un_e = T[e] * Un;
      auto C_e  = T[e] * t_k.C();

      // reconstruct polynomials in monomial space to multiply them
      Polynomial1D fn_e = reconstructFunction1D(polynomialDegree, Un_e) *
                          reconstructFunction1D(polynomialDegree, C_e);

      // tranform f into basis expantion
      auto Fn_e = Te * fn_e;
      // trim F to order polynomialDegreeF
      Fn.push_back ( trim1D(Fn_e, targetPolynomialDegree) );
    }

  return Fn;
}




enum class Boundary {left, right, top, bottom};
/**
 * Set boundary
 */
inline void setBoundary_Periodic (UniqueSquareGrid & mesh, Boundary b)
{
  LIKWID_MARKER_START("BOUNDARY_PERIODIC");
  switch (b)
    {
    case Boundary::left :
#pragma omp for
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          mesh.getUpper(i, -1).U1() = mesh.getUpper(i,mesh.getColumns()-1).U1();
          mesh.getUpper(i, -1).U2() = mesh.getUpper(i,mesh.getColumns()-1).U2();
          mesh.getUpper(i, -1).Fn_c() = mesh.getUpper(i,mesh.getColumns()-1).Fn_c();
        }
      break;
    case Boundary::right :
#pragma omp for
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          mesh.getLower(i, mesh.getColumns()).U1() = mesh.getLower(i, 0).U1();
          mesh.getLower(i, mesh.getColumns()).U2() = mesh.getLower(i, 0).U2();
          mesh.getLower(i, mesh.getColumns()).Fn_b() = mesh.getLower(i, 0).Fn_b();
        }
      break;
    case Boundary::bottom :
#pragma omp for
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          mesh.getUpper(-1, i).U1() = mesh.getUpper(mesh.getRows()-1, i).U1();
          mesh.getUpper(-1, i).U2() = mesh.getUpper(mesh.getRows()-1, i).U2();
          mesh.getUpper(-1, i).Fn_a() = mesh.getUpper(mesh.getRows()-1, i).Fn_a();
        }
      break;
    case Boundary::top :
#pragma omp for
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          mesh.getLower(mesh.getRows(), i).U1() = mesh.getLower(0, i).U1();
          mesh.getLower(mesh.getRows(), i).U2() = mesh.getLower(0, i).U2();
          mesh.getLower(mesh.getRows(), i).Fn_c() = mesh.getLower(0, i).Fn_c();
        }
      break;
    }

  LIKWID_MARKER_STOP("BOUNDARY_PERIODIC");
}

inline void setBoundary_Diriclet (UniqueSquareGrid & mesh, Boundary b,
                                  unsigned int polynomialDegree,
                                  unsigned int polynomialDegreeF,
                                  std::function<double(double,double,double)> const & c,
                                  double time)
{
  LIKWID_MARKER_START("BOUNDARY_DIRICLET");

  using namespace std::placeholders;
  switch (b)
    {
    case Boundary::left :
#pragma omp for
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          unsigned int dof  = numberOf2DBasefunctions(polynomialDegree);
          Triangle & t_g = mesh.getUpper(i, -1);  // element in the ghost layer
          Triangle & t_i = mesh.getLower(i,  0);  // element on the insite of the BC

          // set u in ghost layer 0 -> no influnce on upwind direction
          t_g.U1() = std::vector<double> (dof, 0.);
          t_g.U2() = std::vector<double> (dof, 0.);

          // set F=CU*n on the edge
          //TODO nihct nur linar F
          // -> use U from the inner site; multiply it with normal vector and  map it onto the edge
          Vector const & n = t_i.getNormalB();
          std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
          auto Un = n.x*t_i.U1() + n.y*t_i.U2();
          auto Un_e = T[1] * Un;
          // -> project c onto the edge
          auto C_e = l2Projection_edge(polynomialDegree, std::bind(c,_1,_2,time),
                                       t_i.getJakobian(), t_i.getA(), 1);
          // reconstruct polynomials and multiply them to f
          auto fn_e = reconstructFunction1D(polynomialDegree, Un_e) *
                      reconstructFunction1D(polynomialDegree, C_e);
          // tranform f into basis expantion
          BlockMatrix Te = get1DPolynomialMapping(polynomialDegree*2);
          auto Fn_e = Te * fn_e;
          // trim F to order polynomialDegreeF
          t_g.Fn_c() =-1 * trim1D(Fn_e, polynomialDegreeF);
        }
      break;
    case Boundary::right :
#pragma omp for
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          unsigned int dof  = numberOf2DBasefunctions(polynomialDegree);
          Triangle & t_g = mesh.getLower(i, mesh.getColumns());  // element in the ghost layer
          Triangle & t_i = mesh.getUpper(i, mesh.getColumns()-1);  // element on the insite of the BC

          // set u in ghost layer 0 -> no influnce on upwind direction
          t_g.U1() = std::vector<double> (dof, 0.);
          t_g.U2() = std::vector<double> (dof, 0.);

          // set F=CU*n on the edge
          //TODO nihct nur linar F
          // -> use U from the inner site; multiply it with normal vector and  map it onto the edge
          Vector const & n = t_i.getNormalC();
          std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
          auto Un = n.x*t_i.U1() + n.y*t_i.U2();
          auto Un_e = T[2] * Un;
          // -> project c onto the edge
          auto C_e = l2Projection_edge(polynomialDegree, std::bind(c,_1,_2,time),
                                       t_i.getJakobian(), t_i.getA(), 2);
          // reconstruct polynomials and multiply them to f
          auto fn_e = reconstructFunction1D(polynomialDegree, Un_e) *
                      reconstructFunction1D(polynomialDegree, C_e);
          // tranform f into basis expantion
          BlockMatrix Te = get1DPolynomialMapping(polynomialDegree*2);
          auto Fn_e = Te * fn_e;
          // trim F to order polynomialDegreeF
          t_g.Fn_b() = -1 * trim1D(Fn_e, polynomialDegreeF);
        }
      break;
    case Boundary::bottom :
#pragma omp for
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          unsigned int dof  = numberOf2DBasefunctions(polynomialDegree);
          Triangle & t_g = mesh.getUpper(-1, i);  // element in the ghost layer
          Triangle & t_i = mesh.getLower(0,  i);  // element on the insite of the BC

          // set u in ghost layer 0 -> no influnce on upwind direction
          t_g.U1() = std::vector<double> (dof, 0.);
          t_g.U2() = std::vector<double> (dof, 0.);

          // set F=CU*n on the edge
          //TODO nihct nur linar F
          // -> use U from the inner site; multiply it with normal vector and  map it onto the edge
          Vector const & n = t_i.getNormalC();
          std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
          auto Un = n.x*t_i.U1() + n.y*t_i.U2();
          auto Un_e = T[2] * Un;
          // -> project c onto the edge
          auto C_e = l2Projection_edge(polynomialDegree, std::bind(c,_1,_2,time),
                                       t_i.getJakobian(), t_i.getA(), 2);
          // reconstruct polynomials and multiply them to f
          auto fn_e = reconstructFunction1D(polynomialDegree, Un_e) *
                      reconstructFunction1D(polynomialDegree, C_e);
          // tranform f into basis expantion
          BlockMatrix Te = get1DPolynomialMapping(polynomialDegree*2);
          auto Fn_e = Te * fn_e;
          // trim F to order polynomialDegreeF
          t_g.Fn_a() = -1 * trim1D(Fn_e, polynomialDegreeF);
        }
      break;
    case Boundary::top :
#pragma omp for
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          unsigned int dof  = numberOf2DBasefunctions(polynomialDegree);
          Triangle & t_g = mesh.getLower(mesh.getRows(), i);  // element in the ghost layer
          Triangle & t_i = mesh.getUpper(mesh.getRows()-1, i);  // element on the insite of the BC

          // set u in ghost layer 0 -> no influnce on upwind direction
          t_g.U1() = std::vector<double> (dof, 0.);
          t_g.U2() = std::vector<double> (dof, 0.);

          // set F=CU*n on the edge
          //TODO nihct nur linar F
          // -> use U from the inner site; multiply it with normal vector and  map it onto the edge
          Vector const & n = t_i.getNormalA();
          std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
          auto Un = n.x*t_i.U1() + n.y*t_i.U2();
          auto Un_e = T[0] * Un;
          // -> project c onto the edge
          auto C_e = l2Projection_edge(polynomialDegree, std::bind(c,_1,_2,time),
                                       t_i.getJakobian(), t_i.getA(), 0);
          // reconstruct polynomials and multiply them to f
          auto fn_e = reconstructFunction1D(polynomialDegree, Un_e) *
                      reconstructFunction1D(polynomialDegree, C_e);
          // tranform f into basis expantion
          BlockMatrix Te = get1DPolynomialMapping(polynomialDegree*2);
          auto Fn_e = Te * fn_e;
          // trim F to order polynomialDegreeF
          t_g.Fn_c() = -1 * trim1D(Fn_e, polynomialDegreeF);
        }
      break;
    }

  LIKWID_MARKER_STOP("BOUNDARY_DIRICLET");
}
















#endif
