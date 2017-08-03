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
inline std::vector<BlockMatrix> assemblyHatE (unsigned int polynomialDegree,
                                              unsigned int polynomialDegreeF)
{
  unsigned int dof1D      = numberOf1DBasefunctions(polynomialDegree);
  unsigned int dof2D      = numberOf2DBasefunctions(polynomialDegree);
  unsigned int dofF       = numberOf1DBasefunctions(polynomialDegreeF);
  std::vector<BlockMatrix> T = getLinearTrasformationToRefEdge(polynomialDegree);
  //TODO T muss hier ganz anderes sein

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
//TODO größe von E ist noch falsch
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

typedef std::function<std::vector<double>(unsigned int, unsigned int,
                                          std::vector<double>const &,
                                          std::vector<double>const &,
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
  for ( unsigned int row=0; row<mesh.getRows(); ++row)
    for ( unsigned int col=0; col<mesh.getColumns(); ++col)
    {
      { // Lower Triangle
        //TODO periodic boundary
        Triangle & k   = mesh.getLower(row, col);
        Triangle & n_a = mesh.getUpper(row, col);
        Triangle & n_b = mesh.getUpper(row, col-1);
        Triangle & n_c = mesh.getUpper(row-1, col);

        k.F_a() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_a.F1(), n_a.F2(),
                               k.U1(), k.U2(), n_a.U1(), n_a.U2(), k.getNormalA(),
                               0, 1, hatI);
        k.F_b() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_b.F1(), n_b.F2(),
                               k.U1(), k.U2(), n_b.U1(), n_b.U2(), k.getNormalB(),
                               1,  2, hatI);
        k.F_c() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_c.F1(), n_c.F2(),
                               k.U1(), k.U2(), n_c.U1(), n_c.U2(), k.getNormalC(),
                               2,  0, hatI);
      }

      { // Upper Triangle
        Triangle & k   = mesh.getUpper(row, col);
        Triangle & n_a = mesh.getLower(row+1, col);
        Triangle & n_b = mesh.getLower(row, col);
        Triangle & n_c = mesh.getLower(row, col+1);

        k.F_a() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_a.F1(), n_a.F2(),
                               k.U1(), k.U2(), n_a.U1(), n_a.U2(), k.getNormalA(),
                               0,  2, hatI);
        k.F_b() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_b.F1(), n_b.F2(),
                               k.U1(), k.U2(), n_b.U1(), n_b.U2(), k.getNormalB(),
                               1, 0, hatI);
        k.F_c() = riemanSolver(polynomialDegree, polynomialDegreeF, k.F1(), k.F2(),
                               n_c.F1(), n_c.F2(),
                               k.U1(), k.U2(), n_c.U1(), n_c.U2(), k.getNormalC(),
                               2, 1, hatI);
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
                                                   unsigned int polynomialDegreeF,
                                                   std::vector<double> const & F1_k,
                                                   std::vector<double> const & F2_k,
                                                   std::vector<double> const & F1_n,
                                                   std::vector<double> const & F2_n,
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

  std::vector<BlockMatrix> T   = getLinearTrasformationToRefEdge(polynomialDegree);
  std::vector<BlockMatrix> T_f = getLinearTrasformationToRefEdge(polynomialDegreeF);

  // normal u at edge from both sides
  std::vector<double> ubar1_k = convertToEdge(T[edgeID_k], U1_k, polynomialDegree);
  std::vector<double> ubar2_k = convertToEdge(T[edgeID_k], U2_k, polynomialDegree);
  std::vector<double> ubar1_n = convertToEdge(T[edgeID_n], U1_n, polynomialDegree);
  std::vector<double> ubar2_n = convertToEdge(T[edgeID_n], U2_n, polynomialDegree);

  Polynomial1D uNormal_k = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_k + n_ke.y*ubar2_k);
  Polynomial1D uNormal_n = reconstructFunction1D(polynomialDegree, n_ke.x*ubar1_n + n_ke.y*ubar2_n);

  double u_k = integradeOverRefEdge(uNormal_k);
  double u_n = integradeOverRefEdge(uNormal_n);

  // up_wind direction
  if ( 0 >= u_k+u_n)
    { // from elemnt k to n
      auto Fr1_e = convertToEdge(T_f[edgeID_k], F1_k, polynomialDegreeF);
      auto Fr2_e = convertToEdge(T_f[edgeID_k], F2_k, polynomialDegreeF);

      Frn_e = Fr1_e*n_ke.x + Fr2_e*n_ke.y;
    }
  else
    { // from elemnt n to k
      auto Fr1_e = convertToEdge(T_f[edgeID_n], F1_n, polynomialDegreeF);
      auto Fr2_e = convertToEdge(T_f[edgeID_n], F2_n, polynomialDegreeF);

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
  //TODO hier mit targetDegree etwas rumspeilen

  auto c_pol ( reconstructFunction2D(polynomialDegree, t_k.C()) );
  auto u1_pol ( reconstructFunction2D(polynomialDegree, t_k.U1()) );
  auto u2_pol ( reconstructFunction2D(polynomialDegree, t_k.U2()) );

  auto F1_pol = c_pol * u1_pol;
  auto F2_pol = c_pol * u2_pol;


  BlockMatrix T = getPolynomialMapping(targetPolynomialDegree);

  auto F1 (T * F1_pol);
  auto F2 (T * F2_pol);

  return {F1, F2};
}

enum class Boundary {left, right, top, bottom};
/**
 * Set boundary
 */
inline void setBoundary_Periodic (UniqueSquareGrid & mesh, Boundary b)
{
  switch (b)
    {
    case Boundary::left :
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          mesh.getUpper(i, -1).U1() = mesh.getUpper(i,mesh.getColumns()-1).U1();
          mesh.getUpper(i, -1).U2() = mesh.getUpper(i,mesh.getColumns()-1).U2();
          mesh.getUpper(i, -1).F1() = mesh.getUpper(i,mesh.getColumns()-1).F1();
          mesh.getUpper(i, -1).F2() = mesh.getUpper(i,mesh.getColumns()-1).F2();
        }
      break;
    case Boundary::right :
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          mesh.getLower(i, mesh.getColumns()).U1() = mesh.getLower(i, 0).U1();
          mesh.getLower(i, mesh.getColumns()).U2() = mesh.getLower(i, 0).U2();
          mesh.getLower(i, mesh.getColumns()).F1() = mesh.getLower(i, 0).F1();
          mesh.getLower(i, mesh.getColumns()).F2() = mesh.getLower(i, 0).F2();
        }
      break;
    case Boundary::bottom :
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          mesh.getUpper(-1, i).U1() = mesh.getUpper(mesh.getRows()-1, i).U1();
          mesh.getUpper(-1, i).U2() = mesh.getUpper(mesh.getRows()-1, i).U2();
          mesh.getUpper(-1, i).F1() = mesh.getUpper(mesh.getRows()-1, i).F1();
          mesh.getUpper(-1, i).F2() = mesh.getUpper(mesh.getRows()-1, i).F2();
        }
      break;
    case Boundary::top :
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          mesh.getLower(mesh.getRows(), i).U1() = mesh.getLower(0, i).U1();
          mesh.getLower(mesh.getRows(), i).U2() = mesh.getLower(0, i).U2();
          mesh.getLower(mesh.getRows(), i).F1() = mesh.getLower(0, i).F1();
          mesh.getLower(mesh.getRows(), i).F2() = mesh.getLower(0, i).F2();
        }
      break;
    }
}

inline void setBoundary_Diriclet (UniqueSquareGrid & mesh, Boundary b,
                                  unsigned int polynomialDegree,
                                  std::function<double(double,double)> const & f1,
                                  std::function<double(double,double)> const & f2)
{
  switch (b)
    {
    case Boundary::left :
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          Triangle & t = mesh.getLower(i, -1);
          t.U1() = mesh.getUpper(i,0).U1();
          t.U2() = mesh.getUpper(i,0).U2();
          t.F1() = l2Projection(polynomialDegree, f1, t.getJakobian(), t.getA());
          t.F2() = l2Projection(polynomialDegree, f2, t.getJakobian(), t.getA());
        }
      break;
    case Boundary::right :
      for (unsigned int i=0; i<mesh.getRows(); ++i)
        {
          Triangle & t = mesh.getLower(i, mesh.getColumns());

          t.U1() = mesh.getLower(i, mesh.getColumns()-1).U1();
          t.U2() = mesh.getLower(i, mesh.getColumns()-1).U2();
          t.F1() = l2Projection(polynomialDegree, f1, t.getJakobian(), t.getA());
          t.F2() = l2Projection(polynomialDegree, f2, t.getJakobian(), t.getA());
        }
      break;
    case Boundary::bottom :
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          Triangle & t = mesh.getUpper(-1, i);
          t.U1() = mesh.getUpper(0, i).U1();
          t.U2() = mesh.getUpper(0, i).U2();
          t.F1() = l2Projection(polynomialDegree, f1, t.getJakobian(), t.getA());
          t.F2() = l2Projection(polynomialDegree, f2, t.getJakobian(), t.getA());
        }
      break;
    case Boundary::top :
      for (unsigned int i=0; i<mesh.getColumns(); ++i)
        {
          Triangle & t = mesh.getLower(mesh.getRows(), i);
          t.U1() = mesh.getLower(mesh.getRows()-1, i).U1();
          t.U2() = mesh.getLower(mesh.getRows()-1, i).U2();
          t.F1() = l2Projection(polynomialDegree, f1, t.getJakobian(), t.getA());
          t.F2() = l2Projection(polynomialDegree, f2, t.getJakobian(), t.getA());
        }
      break;
    }
}
















#endif
