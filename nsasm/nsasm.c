# include <stdlib.h>
# include <math.h>
# include <string.h>

# include "mex.h"
# include "matrix.h"

/*
  Functions defined in this file:
*/
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
void assemble ( double *p, int *t, double *u0, int nt, int np,
  double *Pr, mwIndex *Ir, mwIndex *Jc, double *L, double nu, int ndof );
void assemble_constraints ( double *e, double *u0, int np, int ne, int n0,
  double *Pr, mwIndex *Ir, mwIndex *Jc, double *L, int ndof );
void i4vec_heap_d ( int n, int a[] );
void i4vec_sort_heap_a ( int n, int a[] );
void init_shape ( int ngp, double gp[][3] );
void localKL ( int ngp, double w[], double *p, int *tt, double u0[], int np,
  double nu, double lK[15][15], double lL[15] );
void quad_rule ( int ngp, double w[], double gp[][3] );
mxArray* sparse_create ( int *t, double *e, int nt, int np, int np0, int ne,
  int ndof );
void sparse_set ( double *Pr, mwIndex *Ir, mwIndex *Jc, int ri, int rj, 
  int ndof, double val );
/*
  Global data:
*/

/*
  Order of Gauss quadrature points.
*/
# define NGP 7
/*
  P1 shape, P1 r and S derivatives, evaluated at Gauss points.
*/
static double sh1[NGP][3];
static double sh1r[NGP][3];
static double sh1s[NGP][3];
/*
  P2 shape, P2 r and S derivatives, evaluated at Gauss points.
*/
static double sh2[NGP][6];
static double sh2r[NGP][6];
static double sh2s[NGP][6];

/******************************************************************************/

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )

/******************************************************************************/
/*
  Purpose:

    mexFunction is the interface between MATLAB and the C code.

  Discussion:

    This file should be called "nsasm.c".  It should be placed in the
    MATLAB user's path.  It can either be compiled externally, with
    a command like

      mex nsasm.c

    creating a compiled MEX file, or, inside of MATLAB, the command

      mex nsasm.c

    accomplishes the same task.  Once the file has been compiled,
    the MATLAB user can invoke the function by typing:

      [ K, L ] = nsasm ( p, t, np0, e, u0, nu )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2012

  Author:

    This C version by John Burkardt.
    Original C version by Per-Olof Persson.

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, int NLHS, the number of items on the left hand side of the
    MATLAB function call.

    Output, mxArray *PLHS[], pointers to the items on the left hand side
    of the MATLAB function call.

    Input, int NRHS, the number of items on the right hand side of the
    MATLAB function call.

    Input, const mxArray *PRHS[], pointers to the items on the right
    hand side of the MATLAB function call.

  Hidden Parameters:

    Input, double P[3][NP], the coordinates of the mesh nodes.

    Input, int T[6][NT], the indices of the nodes in each element.
    Note that these indices are 1-based, not 0-based!

    Input, int NP0, the number of pressure nodes.

    Input, double E[3][NE], the node, variable type (0,1,2) and
    value of each constraint.

    Input, double U0[2*NP+NP0+NE], the current solution vector.

    Input, double NU, the kinematic viscosity.

    Output, double sparse K[2*NP+NP0+NE,2*NP+NP0+NE], the stiffness matrix,
    in compressed column format.

    Output, double L[2*NP+NP0+NE], the residual vector.
*/
{
  int debug = 1;
  double *e;
  int i;
  mwIndex *Ir;
  mwIndex *Jc;
  double *L;
  int ndof;
  int ne;
  int np;
  int np0;
  int nt;
  double nu;
  int nv;
  double *p;
  double *Pr;
  double *t;
  int *t_int;
  double *u0;
/*
  Get input data ( P, T, NP0, E, U, NU ) from Matlab;
*/
  p = mxGetPr ( prhs[0] );
  np = mxGetN ( prhs[0] );
/*
  T is a pointer to a double.
  We need to make T_INT, an integer copy of this data.
*/
  t = mxGetPr ( prhs[1] );
  nt = mxGetN ( prhs[1] );

  t_int = ( int * ) malloc ( 6 * nt * sizeof ( int ) );
  for ( i = 0; i < 6 * nt; i++ )
  {
    t_int[i] = ( int ) t[i];
  }

  np0 = ( int ) mxGetScalar ( prhs[2] );

  e = mxGetPr ( prhs[3] );
  ne = mxGetN ( prhs[3] );

  ndof = 2 * np + np0 + ne;

  u0 = mxGetPr ( prhs[4] );

  nu = mxGetScalar ( prhs[5] );
/*
  Create the sparse matrix K, and retrieve pointers to the
  values, row indices, and compressed column vector.
*/
  plhs[0] = sparse_create ( t_int, e, nt, np, np0, ne, ndof );

  Pr = mxGetPr ( plhs[0] );
  Ir = mxGetIr ( plhs[0] );
  Jc = mxGetJc ( plhs[0] );
/*
  Create the residual vector L.
*/
  plhs[1] = mxCreateDoubleMatrix ( ndof, 1, mxREAL );

  L = mxGetPr ( plhs[1] );
/*
  Assemble the PDE.
*/
  assemble ( p, t_int, u0, nt, np, Pr, Ir, Jc, L, nu, ndof );
/*
  Assemble constraints.
*/
  nv = 2 * np + np0;
  assemble_constraints ( e, u0, np, ne, nv, Pr, Ir, Jc, L, ndof );

  free ( t_int );

  return;
}
/******************************************************************************/

void assemble ( double *p, int *t, double *u0, int nt, int np, double *Pr,
  mwIndex *Ir, mwIndex *Jc, double *L, double nu, int ndof )

/******************************************************************************/
/*
  Purpose:

    ASSEMBLE assembles the local stiffness and residual into global arrays.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2014

  Author:

    This C version by John Burkardt.
    Original C version by Per-Olof Persson.

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, double P[2][NP], the coordinates of the mesh nodes.

    Input, int T[6][NT], the indices of the nodes in each element.

    Input, double U0[2*NP+NP0+NE], the current solution vector.

    Input, int NT, the number of elements.

    Input, int NP, the number of nodes.

    Input/output, double Pr[NNZ], the values of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Ir[NNZ], the row indices of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Jc[N+1], points to the first element of each column.

    Input/output, double L[2*NP+NP0+NE], the residual vector.

    Input, double NU, the kinematic viscosity.

    Input, int NDOF, the number of degrees of freedom.
*/
{
  double gp[7][3];
  int i;
  int i2;
  int it;
  int j;
  int j2;
  double lK[15][15];
  double lL[15];
  const int ngp = 7;
  int ri;
  int rj;
  double w[7];
/*
  Get the quadrature rule.
*/
  quad_rule ( ngp, w, gp );
/*
  Initialize the shape functions.
*/
  init_shape ( ngp, gp );
/*
  For each element, determine the local matrix and right hand side.
*/
  for ( it = 0; it < nt; it++ )
  {
    localKL ( NGP, w, p, t + it * 6, u0, np, nu, lK, lL );
/*
  Add the local right hand side and matrix to the global data.
*/
    for ( i = 0; i < 15; i++ )
    {
      i2 = i % 6;
      ri = t[i2+it*6] + np * ( i / 6 ) - 1;

      L[ri] = L[ri] + lL[i];

      for ( j = 0; j < 15; j++ )
      {
        j2 = j % 6;
        rj = t[j2+it*6] + np * ( j / 6 ) - 1;
        sparse_set ( Pr, Ir, Jc, ri, rj, ndof, lK[i][j] );
      }
    }
  }
  return;
}
/******************************************************************************/

void assemble_constraints ( double *e, double *u0, int np, int ne, int nv,
  double *Pr, mwIndex *Ir, mwIndex *Jc, double *L, int ndof )

/******************************************************************************/
/*
  Purpose:

    ASSEMBLE_CONSTRAINTS assembles the constraints.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2014

  Author:

    This C version by John Burkardt.
    Original C version by Per-Olof Persson.

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, double E[3][NE], the node, variable type (0,1,2) and
    value of each constraint.

    Input, double U0[2*NP+NP0+NE], the current solution vector.

    Input, int NP, the number of nodes.

    Input, int NE, the number of constraints.

    Input, int NV, the number of variables.

    Input/output, double Pr[NNZ], the values of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Ir[NNZ], the row indices of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Jc[N+1], points to the first element of each column.

    Input/output, double L[2*NP+NP0+NE], the residual vector.  On output,
    the residual vector has been updated to include the constraints.

    Input, int NDOF, the number of degrees of freedom.
*/
{
  int ie;
  const double one = 1.0;
  int ri;
  int rj;

  for ( ie = 0; ie < ne; ie++ )
  {
    ri = nv + ie;
    rj = ( int ) e[1+ie*3] * np + ( int ) e[0+ie*3] - 1;
    sparse_set ( Pr, Ir, Jc, ri, rj, ndof, one );
    sparse_set ( Pr, Ir, Jc, rj, ri, ndof, one );
    L[rj] = L[rj] + u0[ri];
    L[ri] = u0[rj] - e[2+ie*3];
  }

  return;
}
/******************************************************************************/

void i4vec_heap_d ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_HEAP_D reorders an I4VEC into a descending heap.

  Discussion:

    A heap is an array A with the property that, for every index J,
    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
    2*J+1 and 2*J+2 are legal).

  Diagram:

                  A(0)
                /      \
            A(1)         A(2)
          /     \        /  \
      A(3)       A(4)  A(5) A(6)
      /  \       /   \
    A(7) A(8)  A(9) A(10)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 April 1999

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the size of the input array.

    Input/output, int A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
  int i;
  int ifree;
  int key;
  int m;
/*
  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
*/
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
/*
  Copy the value out of the parent node.
  Position IFREE is now "open".
*/
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
/*
  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
  IFREE.  (One or both may not exist because they equal or exceed N.)
*/
      m = 2 * ifree + 1;
/*
  Does the first position exist?
*/
      if ( n <= m )
      {
        break;
      }
      else
      {
/*
  Does the second position exist?
*/
        if ( m + 1 < n )
        {
/*
  If both positions exist, take the larger of the two values,
  and update M if necessary.
*/
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
/*
  If the large descendant is larger than KEY, move it up,
  and update IFREE, the location of the free position, and
  consider the descendants of THIS position.
*/
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
/*
  When you have stopped shifting items up, return the item you
  pulled out back to the heap.
*/
    a[ifree] = key;
  }

  return;
}
/********************************************************************/

void i4vec_sort_heap_a ( int n, int a[] )

/********************************************************************/
/*
  Purpose:

    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 April 1999

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN: 0-12-519260-6.

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, int A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
/*
  1: Put A into descending heap form.
*/
  i4vec_heap_d ( n, a );
/*
  2: Sort A.

  The largest object in the heap is in A[0].
  Move it to position A[N-1].
*/
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
/*
  Consider the diminished heap of size N1.
*/
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
/*
  Restore the heap structure of the initial N1 entries of A.
*/
    i4vec_heap_d ( n1, a );
/*
  Take the largest object from A[0] and move it to A[N1-1].
*/
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
/******************************************************************************/

void init_shape ( int ngp, double gp[][3] )

/******************************************************************************/
/*
  Purpose:

    INIT_SHAPE evaluates the shape functions at the quadrature points.

  Discussion:

    The shape functions and their spatial derivatives are evaluated at the
    Gauss points in the reference triangle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2006

  Author:

    Per-Olof Persson

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, int NGP, the number of quadrature points.

    Input, double GP[NGP][3], the quadrature points.

    Global, static double SH1[NGP][3], SH1R[NGP][3], SH1S[NGP][3], the
    P1 (pressure) shape functions, and R and S derivatives, evaluated at the
    Gauss points.

    Global, static double SH2[NGP][6], SH2R[NGP][6], SH2S[NGP][6], the
    P2 (velocity) shape functions, and R and S derivatives, evaluated at the
    Gauss points.
*/
{
  int i;

  for ( i = 0; i < ngp; i++ )
  {
    sh1[i][0] = gp[i][2];
    sh1[i][1] = gp[i][0];
    sh1[i][2] = gp[i][1];

    sh1r[i][0] = -1.0;
    sh1r[i][1] =  1.0;
    sh1r[i][2] =  0.0;

    sh1s[i][0] = -1.0;
    sh1s[i][1] =  0.0;
    sh1s[i][2] =  1.0;

    sh2[i][0] = 1.0 - 3.0 * gp[i][0] - 3.0 * gp[i][1]
      + 2.0 * gp[i][0] * gp[i][0] + 4.0 * gp[i][0] * gp[i][1]
      + 2.0 * gp[i][1] * gp[i][1];
    sh2[i][1] = - gp[i][0] + 2.0 * gp[i][0] * gp[i][0];
    sh2[i][2] = - gp[i][1] + 2.0 * gp[i][1] * gp[i][1];
    sh2[i][3] = 4.0 * gp[i][0] * gp[i][1];
    sh2[i][4] = 4.0 * gp[i][1] - 4.0 * gp[i][0] * gp[i][1]
      - 4.0 * gp[i][1] * gp[i][1];
    sh2[i][5] = 4.0 * gp[i][0] - 4.0 * gp[i][0] * gp[i][1]
      - 4.0 * gp[i][0] * gp[i][0];

    sh2r[i][0] = - 3.0 + 4.0 * gp[i][0] + 4.0 * gp[i][1];
    sh2r[i][1] = - 1.0 + 4.0 * gp[i][0];
    sh2r[i][2] = 0.0;
    sh2r[i][3] = 4.0 * gp[i][1];
    sh2r[i][4] = - 4.0 * gp[i][1];
    sh2r[i][5] = 4.0 - 8.0 * gp[i][0] - 4.0 * gp[i][1];

    sh2s[i][0] = - 3.0 + 4.0 * gp[i][0] + 4.0 * gp[i][1];
    sh2s[i][1] =   0.0;
    sh2s[i][2] = - 1.0 + 4.0 * gp[i][1];
    sh2s[i][3] =   4.0 * gp[i][0];
    sh2s[i][4] =   4.0 - 8.0 * gp[i][1] - 4.0 * gp[i][0];
    sh2s[i][5] = - 4.0 * gp[i][0];
  }
  return;
}
/******************************************************************************/

void localKL ( int ngp, double w[], double *p, int *tt, double u0[], int np, 
  double nu, double lK[15][15], double lL[15] )

/******************************************************************************/
/*
  Purpose:

    LOCALKL assembles the local stiffness matrix and residual.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2006

  Author:

    Per-Olof Persson

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, int NGP, the number of quadrature points.

    Input, double W[NGP], the quadrature weights.

    Input, double P[2][NP], the coordinates of the nodes.

    Input, int *TT, a pointer to the list of nodes for the current element.

    Input, double U0[2*NP+NP0+NE], called "U" elsewhere in the program,
    but renamed "U0" here because we want to use "U" for the horizontal
    velocity.  U0 contains the current finite element solution.

    Input, int NP, the number of nodes.

    Input, double NU, the kinematic viscosity.

    Output, double lK[15][15], the local stiffness matrix.

    Output, double lL[15], the local residual vector.
*/
{
  int i;
  int igp;
  int j;
  double Jdet;
  double Jinv[2][2];
  double mul;
  double px;
  double py;
  double sh1x[3];
  double sh1y[3];
  double sh2x[6];
  double sh2y[6];
  double u;
  double ux;
  double uy;
  double v;
  double vx;
  double vy;
  double xr;
  double xs;
  double yr;
  double ys;
/*
  Zero out lK and lL.
*/
  memset ( lK, 0, 15 * 15 * sizeof ( double ) );

  memset ( lL, 0, 15 * sizeof ( double ) );

  for ( igp = 0; igp < ngp; igp++ )
  {
/*
  Jacobian
*/
    xr = 0.0;
    xs = 0.0;
    yr = 0.0;
    ys = 0.0;

    for ( i = 0; i < 6; i++ )
    {
      xr = xr + sh2r[igp][i] * p[(tt[i]-1)*2+0];
      xs = xs + sh2s[igp][i] * p[(tt[i]-1)*2+0];
      yr = yr + sh2r[igp][i] * p[(tt[i]-1)*2+1];
      ys = ys + sh2s[igp][i] * p[(tt[i]-1)*2+1];
    }

    Jdet = xr * ys - xs * yr;

    Jinv[0][0] =  ys / Jdet;
    Jinv[0][1] = -xs / Jdet;
    Jinv[1][0] = -yr / Jdet;
    Jinv[1][1] =  xr / Jdet;
/*
  Set the X and Y derivatives of the shape functions.
*/
    for ( i = 0; i < 3; i++ )
    {
      sh1x[i] = sh1r[igp][i] * Jinv[0][0] + sh1s[igp][i] * Jinv[1][0];
      sh1y[i] = sh1r[igp][i] * Jinv[0][1] + sh1s[igp][i] * Jinv[1][1];
    }
    for ( i = 0; i < 6; i++ )
    {
      sh2x[i] = sh2r[igp][i] * Jinv[0][0] + sh2s[igp][i] * Jinv[1][0];
      sh2y[i] = sh2r[igp][i] * Jinv[0][1] + sh2s[igp][i] * Jinv[1][1];
    }
/*
  Solution and derivatives.
*/
    u  = 0.0;
    ux = 0.0;
    uy = 0.0;

    v  = 0.0;
    vx = 0.0;
    vy = 0.0;
    for ( i = 0; i < 6; i++ )
    {
      u  = u  + sh2[igp][i] * u0[tt[i]-1];
      ux = ux + sh2x[i]     * u0[tt[i]-1];
      uy = uy + sh2y[i]     * u0[tt[i]-1];

      v  = v  + sh2[igp][i] * u0[np+tt[i]-1];
      vx = vx + sh2x[i]     * u0[np+tt[i]-1];
      vy = vy + sh2y[i]     * u0[np+tt[i]-1];
    }

    px = 0.0;
    py = 0.0;
    for ( i = 0; i < 3; i++ )
    {
      px = px + sh1x[i] * u0[2*np+tt[i]-1];
      py = py + sh1y[i] * u0[2*np+tt[i]-1];
    }
/*
  Local K.
*/
    mul = w[igp] * Jdet / 2.0;

    for ( i = 0; i < 6; i++ )
    {
      lL[i]   = lL[i]   + mul * (
        ( u * ux + v * uy + px ) * sh2[igp][i]
        + nu * ( ux * sh2x[i] + uy * sh2y[i] ) );

      lL[6+i] = lL[6+i] + mul * (
        ( u * vx + v * vy + py ) * sh2[igp][i]
        + nu * ( vx * sh2x[i] + vy * sh2y[i] ) );

      for ( j = 0; j < 6; j++ )
      {
        lK[i][j] +=     nu * ( sh2x[i] * sh2x[j] + sh2y[i] * sh2y[j] ) * mul;
        lK[6+i][6+j] += nu * ( sh2x[i] * sh2x[j] + sh2y[i] * sh2y[j] ) * mul;

        lK[i][j] +=     ( u  * sh2[igp][i] * sh2x[j]
                        + v  * sh2[igp][i] * sh2y[j] ) * mul;
        lK[6+i][6+j] += ( u  * sh2[igp][i] * sh2x[j]
                        + v  * sh2[igp][i] * sh2y[j] ) * mul;

        lK[i][j] +=     ( ux * sh2[igp][i] * sh2[igp][j] ) * mul;
        lK[i][6+j] +=   ( uy * sh2[igp][i] * sh2[igp][j] ) * mul;
        lK[6+i][j] +=   ( vx * sh2[igp][i] * sh2[igp][j] ) * mul;
        lK[6+i][6+j] += ( vy * sh2[igp][i] * sh2[igp][j] ) * mul;
      }
      for ( j = 0; j < 3; j++ )
      {
        lK[i][12+j] +=   ( sh2[igp][i] * sh1x[j] ) * mul;
        lK[6+i][12+j] += ( sh2[igp][i] * sh1y[j] ) * mul;
      }
    }
/*
  Local L.
*/
    for ( i = 0; i < 3; i++ )
    {
      lL[12+i] = lL[12+i] + ( ux + vy ) * sh1[igp][i] * mul;
      for ( j = 0; j < 6; j++ )
      {
        lK[12+i][j] =   lK[12+i][j]   + ( sh1[igp][i] * sh2x[j] ) * mul;
        lK[12+i][6+j] = lK[12+i][6+j] + ( sh1[igp][i] * sh2y[j] ) * mul;
      }
    }

  }
  return;
}
/******************************************************************************/

void quad_rule ( int ngp, double w[], double gp[][3] )

/******************************************************************************/
/*
  Purpose:

    QUAD_RULE returns the points and weights of a quadrature rule.

  Discussion:

    At the moment, only a 7-point rule is available.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2011

  Author:

    John Burkardt.

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, int NGP, the number of quadrature points.

    Output, double W[NGP], the quadrature weights.

    Output, double GP[NGP][3], the quadrature points.
*/
{
  static double gp_save[7][3] = {
    0.333333333333333, 0.333333333333333, 0.333333333333334,
    0.059715871789770, 0.470142064105115, 0.470142064105115,
    0.470142064105115, 0.059715871789770, 0.470142064105115,
    0.470142064105115, 0.470142064105115, 0.059715871789770,
    0.797426985353087, 0.101286507323457, 0.101286507323457,
    0.101286507323457, 0.797426985353087, 0.101286507323457,
    0.101286507323457, 0.101286507323457, 0.797426985353087 };
  int i;
  int j;
  static double w_save[7] = {
    0.225000000000000,
    0.132394152788506,
    0.132394152788506,
    0.132394152788506,
    0.125939180544827,
    0.125939180544827,
    0.125939180544827 };

  for ( i = 0; i < ngp; i++ )
  {
    w[i] = w_save[i];
  }

  for ( i = 0; i < ngp; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      gp[i][j] = gp_save[i][j];
    }
  }

  return;
}
/******************************************************************************/

mxArray* sparse_create ( int *t, double *e, int nt, int np, int np0, int ne,
  int ndof )

/******************************************************************************/
/*
  Purpose:

    SPARSE_CREATE creates a sparse matrix.

  Discussion:

    The invocation of the routine SORT was replaced by a call to
    I4VEC_SORT_HEAP_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    Per-Olof Persson

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input, int T[6][NT], the indices of the nodes in each element.

    Input, double E[3][NE], the node, variable type (0,1,2) and
    value of each constraint.

    Input, int NT, the number of elements.

    Input, int NP, the number of nodes.

    Input, int NP0, the number of pressure nodes.

    Input, int NE, the number of constraints.

    Input, int NDOF, the number of degrees of freedom, = 2 * NP + NP0 + NE.

    Output, mxArray *K, the structure that stores the sparse Jacobian.

  Local Parameters:

    Local, int COL[NDOF], the number of entries in each column.

    Local, int COLSTART[NDOF], the start index of each column.

    Local, int IDXI[NIDX], the row indices.

    Local, int IDXJ[NIDX], the column indices.

    Local, int NIDX, the upper limit on the number of entries.
*/
{
  int *col;
  int *colstart;
  int i;
  int *idxi;
  int *idxj;
  int ie;
  int ii;
  mwIndex *Ir;
  int it;
  int j;
  mwIndex *Jc;
  mxArray *K;
  int n;
  int *n1;
  double *n2;
  int nidx;
  int nnz;
  int *p;
  int ri;
  int rj;

  nidx = nt * 15 * 15 + 2 * ne;

  idxi =     ( int * ) mxCalloc ( nidx, sizeof ( int ) );
  idxj =     ( int * ) mxCalloc ( nidx, sizeof ( int ) );
  col =      ( int * ) mxCalloc ( ndof, sizeof ( int ) );
  colstart = ( int * ) mxCalloc ( ndof, sizeof ( int ) );
/*
  Count the elements in each column, including duplicates.
*/
  for ( it = 0, n1 = t; it < nt; it++, n1+=6 )
  {
    for ( j = 0; j < 6; j++ )
    {
      col[n1[j]-1]    = col[n1[j]-1]    + 15;
      col[n1[j]-1+np] = col[n1[j]-1+np] + 15;
    }

    for ( j = 0; j < 3; j++ )
    {
      col[n1[j]-1+2*np] = col[n1[j]-1+2*np] + 15;
    }
  }
/*
  Constraints.
*/
  for ( ie = 0, n2 = e; ie < ne; ie++, n2+=3 )
  {
    col[2*np+np0+ie]++;
    col[ ( int ) n2[1]*np+ ( int ) n2[0]-1]++;
  }
/*
  Form cumulative sum (start of each column)
*/
  for ( i = 0; i < ndof-1; i++ )
  {
    colstart[i+1] = colstart[i] + col[i];
  }
/*
  Insert row/column indices, sorted columns.
*/
  for ( it = 0, n1 = t; it < nt; it++, n1+=6 )
  {
    for ( i = 0; i < 15; i++ )
    {
      ri = n1[i%6] + np * (i/6) - 1;
      for ( j = 0; j < 15; j++ )
      {
        rj = n1[j%6] + np * (j/6) - 1;
        idxi[colstart[rj]] = ri;
        idxj[colstart[rj]] = rj;
        colstart[rj]++;
      }
    }
  }
/*
  Constraints.
*/
  for ( ie = 0, n2 = e; ie < ne; ie++, n2+=3 )
  {
    ri = ( int ) n2[1] * np + (int) n2[0] - 1;
    rj = 2 * np + np0 + ie;
    idxi[colstart[rj]] = ri;
    idxj[colstart[rj]] = rj;
    colstart[rj]++;

    ri = 2 * np + np0 + ie;
    rj = ( int ) n2[1] * np + (int) n2[0] - 1;
    idxi[colstart[rj]] = ri;
    idxj[colstart[rj]] = rj;
    colstart[rj]++;
  }

  mxFree ( colstart );
/*
  Sort the row indices in each column.
  NOTE: The original code called SORT, a C library function.
  However, the MATLAB Mex compiler did not seem to provide an interface
  to this library routine, so we had to supply the source code for
  a sort routine doing our best to make a proper patch.
*/
  p = idxi;
  for ( i = 0; i < ndof; i++ )
  {
    n = col[i];
/*
    sort ( p, p+col[i] );
*/
    i4vec_sort_heap_a ( n, p );

    p = p + col[i];
  }
/*
  Zero out COL.
*/
  memset ( col, 0, ndof * sizeof ( int ) );
/*
  Count the unique entries.
*/
  nnz = 1;
  col[idxj[0]]++;

  for ( ii = 1; ii < nidx; ii++ )
  {
    if ( idxi[ii] != idxi[ii-1] || idxj[ii] != idxj[ii-1] )
    {
      nnz = nnz + 1;
      col[idxj[ii]]++;
    }
  }
/*
  Create the sparse matrix.
*/
  K = mxCreateSparse ( ndof, ndof, nnz, mxREAL );

  if ( !K )
  {
    printf ( "\n" );
    printf ( "SPARSE_CREATE: Fatal error!\n" );
    printf ( "  mxCreateSparse failed while creating the matrix K.\n" );
    exit ( 1 );
  }

  Ir = mxGetIr ( K );
  Jc = mxGetJc ( K );
/*
  Form the Jc vector.
*/
  Jc[0] = 0;
  for ( i = 1; i <= ndof; i++ )
  {
    Jc[i] = Jc[i-1] + col[i-1];
  }
  mxFree ( col );
/*
  Form the Ir vector.
*/
  *(Ir++) = idxi[0];
  for ( ii = 1; ii < nidx; ii++ )
  {
    if ( idxi[ii] != idxi[ii-1] || idxj[ii] != idxj[ii-1] )
    {
      *(Ir++) = idxi[ii];
    }
  }
  mxFree ( idxj );
  mxFree ( idxi );

  return K;
}
/******************************************************************************/

void sparse_set ( double *Pr, mwIndex *Ir, mwIndex *Jc, int ri, int rj, 
  int ndof, double val )

/******************************************************************************/
/*
  Purpose:

    SPARSE_SET increments an entry of the sparse matrix.

  Discussion:

    We know RJ, the column in which the entry occurs.  We know that entries
    in column RJ occur in positions K1 = Jc[RJ] through K2 = Jc[RJ+1]-1.
    We may assume that the corresponding entries Ir[K1] through Ir[K2]
    are ascending sorted and that one of these entries is equal to RI.

    We now simply use binary search on Ir to locate the index K for which
    Ir[K] = RI, and then increment Pr[K].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2012

  Author:

    This C version by John Burkardt.
    Original C version by Per-Olof Persson.

  Reference:

    Per-Olof Persson,
    Implementation of Finite Element-Based Navier-Stokes Solver,
    April 2002.

  Parameters:

    Input/output, double Pr[NNZ], the values of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Ir[NNZ], the row indices of the nonzero entries
    of the sparse matrix.

    Input, mwIndex Jc[N+1], points to the first element of each column.

    Input, int RI, RJ, the row and column index of the matrix entry
    that is to be incremented.

    Input, double VAL, the amount that is to be added to the matrix entry.
*/
{
  mwIndex cr;
  mwIndex k;
  mwIndex k1;
  mwIndex k2;

  cr = ri;
/*
  Set the range K1 <= K2 of entries in Ir to be searched.
*/
  k1 = Jc[rj];
  k2 = Jc[rj+1] - 1;

  if ( ri < 0 || ndof <= ri )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPARSE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of index i = %d!\n", ri );
    fprintf ( stderr, "  Nominal increment is A(%d,%d) = %g\n", ri, rj, val );
    exit ( 1 );
  }

  if ( rj < 0 || ndof <= rj )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPARSE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of index j = %d!\n", rj );
    fprintf ( stderr, "  Nominal increment is A(%d,%d) = %g\n", ri, rj, val );
    exit ( 1 );
  }
/*
  Find K so that Ir[K] = ri;
*/
  do
  {
    k = ( k1 + k2 ) / 2;

    if ( cr < Ir[k] )
    {
      k2 = k - 1;
    }
    else
    {
      k1 = k + 1;
    }

    if ( cr == Ir[k] )
    {
      break;
    }

  } while ( k1 <= k2 );
/*
  Increment the matrix entry.
*/
  Pr[k] = Pr[k] + val;

  return;
}
