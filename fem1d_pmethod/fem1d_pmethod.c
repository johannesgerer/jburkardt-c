# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( void );
void alpbet ( double a[], double alpha[], double beta[], int np, int problem,
  int quad_num, double quad_w[], double quad_x[] );
void exact ( double alpha[], double beta[], double f[], int np, int nprint,
  int problem, int quad_num, double quad_w[], double quad_x[] );
double ff ( double x, int problem );
void ortho ( double a[], double alpha[], double beta[], int np, int problem,
  int quad_num, double quad_w[], double quad_x[] );
void out ( double alpha[], double beta[], double f[], int np, int nprint );
void phi ( double alpha[], double beta[], int i, int np, double *phii,
  double *phiix, double x );
double pp ( double x, int problem );
void quad ( int quad_num, double quad_w[], double quad_x[] );
double qq ( double x, int problem );
void sol ( double a[], double alpha[], double beta[], double f[], int np,
  int problem, int quad_num, double quad_w[], double quad_x[] );
void timestamp ( void );
double uex ( double x, int problem );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM1D_PMETHOD.

  Discussion:

    FEM1D_PMETHOD implements the P-version of the finite element method.

    Program to solve the one dimensional problem:

      - d/dX (P dU/dX) + Q U  =  F

    by the finite-element method using a sequence of polynomials
    which satisfy the boundary conditions and are orthogonal
    with respect to the inner product:

      (U,V)  =  Integral (-1 to 1) P U' V' + Q U V dx

    Here U is an unknown scalar function of X defined on the
    interval [-1,1], and P, Q and F are given functions of X.

    The boundary values are U(-1) = U(1)=0.

    Sample problem #1:

      U=1-x^4,        P=1, Q=1, F=1.0+12.0*x^2-x^4

    Sample problem #2:

      U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x)

    The program should be able to get the exact solution for
    the first problem, using NP = 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Local Parameters:

    Local, double A[NP+1], the squares of the norms of the
    basis functions.

    Local, double ALPHA[NP], BETA[NP], the basis function 
    recurrence coefficients.

    Local, double F[NP+1].
    F contains the basis function coefficients that form the
    representation of the solution U.  That is,
      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
    where "BASIS(I)(X)" means the I-th basis function
    evaluated at the point X.

    Local, int NP.
    The highest degree polynomial to use.

    Local, int NPRINT.
    The number of points at which the computed solution
    should be printed out at the end of the computation.

    Local, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Local, int QUAD_NUM, the order of the quadrature rule.

    Local, double QUAD_W[QUAD_NUM], the quadrature weights.

    Local, double QUAD_X[QUAD_NUM], the quadrature abscissas.
*/
{
# define NP 2
# define QUAD_NUM 10

  double a[NP+1];
  double alpha[NP];
  double beta[NP];
  double f[NP+1];
  int nprint = 10;
  int problem = 2;
  double quad_w[QUAD_NUM];
  double quad_x[QUAD_NUM];

  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_PMETHOD\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Solve the two-point boundary value problem\n" );
  printf ( "\n" );
  printf ( "  - d/dX (P dU/dX) + Q U  =  F\n" );
  printf ( "\n" );
  printf ( "  on the interval [-1,1], with\n" );
  printf ( "  U(-1) = U(1) = 0.\n" );
  printf ( "\n" );
  printf ( "  The P method is used, which represents U as\n" );
  printf ( "  a weighted sum of orthogonal polynomials.\n" );
  printf ( "\n" );
  printf ( "  Highest degree polynomial to use is %d\n", NP );
  printf ( "  Number of points to be used for output = %d\n", nprint );

  if ( problem == 1 )
  {
    printf ( "\n" );
    printf ( "  Problem #1:\n" );
    printf ( "  U=1-x^4,\n" );
    printf ( "  P=1,\n" );
    printf ( "  Q=1,\n" );
    printf ( "  F=1 + 12 * x^2 - x^4\n" );
  }
  else if ( problem == 2 )
  {
    printf ( "\n" );
    printf ( "  Problem #2:\n" );
    printf ( "  U=cos(0.5*pi*x),\n" );
    printf ( "  P=1,\n" );
    printf ( "  Q=0,\n" );
    printf ( "  F=0.25*pi*pi*cos(0.5*pi*x)\n" );
  }
/*
  Get quadrature abscissas and weights for interval [-1,1].
*/
  quad ( QUAD_NUM, quad_w, quad_x );
/*
  Compute the constants for the recurrence relationship
  that defines the basis functions.
*/
  alpbet ( a, alpha, beta, NP, problem, QUAD_NUM, quad_w, quad_x );
/*
  Test the orthogonality of the basis functions.
*/
  ortho ( a, alpha, beta, NP, problem, QUAD_NUM, quad_w, quad_x );
/*
  Solve for the solution of the problem, in terms of coefficients
  of the basis functions.
*/
  sol ( a, alpha, beta, f, NP, problem, QUAD_NUM, quad_w, quad_x );
/*
  Print out the solution, evaluated at each of the NPRINT points.
*/
  out ( alpha, beta, f, NP, nprint );
/*
  Compare the computed and exact solutions.
*/
  exact ( alpha, beta, f, NP, nprint, problem, QUAD_NUM, quad_w, quad_x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_PMETHOD\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
# undef NP
# undef QUAD_NUM
}
/******************************************************************************/

void alpbet ( double a[], double alpha[], double beta[], int np, int problem,
  int quad_num, double quad_w[], double quad_x[] )

/******************************************************************************/
/*
  Purpose:

    ALPBET calculates the coefficients in the recurrence relationship.

  Discussion:

    ALPHA and BETA are the coefficients in the three
    term recurrence relation for the orthogonal basis functions
    on [-1,1].

    The routine also calculates A, the square of the norm of each basis
    function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Output, double A[NP+1], the squares of the norms of the
    basis functions.

    Output, double ALPHA[NP], BETA[NP], the basis function 
    recurrence coefficients.

    Input, int NP.
    The highest degree polynomial to use.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Input, int QUAD_NUM, the order of the quadrature rule.

    Input, double QUAD_W(QUAD_NUM), the quadrature weights.

    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
*/
{
  int i;
  int iq;
  int k;
  double q;
  double qm1;
  double qm1x;
  double qm2;
  double qm2x;
  double qx;
  double s;
  double ss;
  double su;
  double sv;
  double t;
  double u;
  double v;
  double x;

  ss = 0.0;
  su = 0.0;

  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];

    s = 4.0 * pp ( x, problem ) * x * x
      + qq ( x, problem ) * ( 1.0 - x * x ) * ( 1.0 - x * x );

    u = 2.0 * pp ( x, problem ) * x * ( 3.0 * x * x - 1.0 )
      + x * qq ( x, problem ) * ( 1.0 - x * x ) * ( 1.0 - x * x );

    ss = ss + s * quad_w[iq];
    su = su + u * quad_w[iq];

  }

  a[0] = ss;
  alpha[0] = su / ss;
  beta[0] = 0.0;

  for ( i = 1; i <= np; i++ )
  {
    ss = 0.0;
    su = 0.0;
    sv = 0.0;

    for ( iq = 0; iq < quad_num; iq++ )
    {
      x = quad_x[iq];
/*
  Three term recurrence for Q.
*/
      qm1 = 0.0;
      q = 1.0;
      for ( k = 0; k <= i-1; k++ )
      {
        qm2 = qm1;
        qm1 = q;
        q = ( x - alpha[k] ) * qm1 - beta[k] * qm2;
      }
/*
  Three term recurrence for Q'.
*/
      qm1x = 0.0;
      qx = 0.0;
      for ( k = 0; k <= i-1; k++ )
      {
        qm2x = qm1x;
        qm1x = qx;
        qx = qm1 + ( x - alpha[k] ) * qm1x - beta[k] * qm2x;
      }

      t = 1.0 - x * x;
/*
  The basis function PHI = ( 1 - x^2 ) * q.

     s = pp * ( phi(i) )' * ( phi(i) )' + qq * phi(i) * phi(i)
*/
      s = pp ( x, problem ) * pow ( t * qx - 2.0 * x * q, 2 )
        + qq ( x, problem ) * t * t * q * q;
/*
     u = pp * ( x * phi(i) )' * phi(i)' + qq * x * phi(i) * phi(i)
*/
      u = pp ( x, problem )
        * ( x * t * qx + ( 1.0 - 3.0 * x * x ) * q )
        * ( t * qx - 2.0 * x * q ) + x * qq ( x, problem )
        * t * t * q * q;
/*
     v = pp * ( x * phi(i) )' * phi(i-1) + qq * x * phi(i) * phi(i-1)
*/
      v = pp ( x, problem )
        * ( x * t * qx + ( 1.0 - 3.0 * x * x ) * q )
        * ( t * qm1x - 2.0 * x * qm1 )
        + x * qq ( x, problem ) * t * t * q * qm1;
/*
  SS(i) = <   phi(i), phi(i)   > = integral ( S )
  SU(i) = < x phi(i), phi(i)   > = integral ( U )
  SV(i) = < x phi(i), phi(i-1) > = integral ( V )
*/
      ss = ss + s * quad_w[iq];
      su = su + u * quad_w[iq];
      sv = sv + v * quad_w[iq];
    }

    a[i] = ss;
/*
  ALPHA(i) = SU(i) / SS(i)
  BETA(i)  = SV(i) / SS(i-1)
*/
    if ( i < np )
    {
      alpha[i] = su / ss;
      beta[i] = sv / a[i-1];
    }
  }

  return;
}
/******************************************************************************/

void exact ( double alpha[], double beta[], double f[], int np, int nprint,
  int problem, int quad_num, double quad_w[], double quad_x[] )

/******************************************************************************/
/*
  Purpose:

    EXACT compares the computed and exact solutions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double ALPHA(NP), BETA(NP), the basis function 
    recurrence coefficients.

    Input, double F(0:NP).
    F contains the basis function coefficients that form the
    representation of the solution U.  That is,
      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
    where "BASIS(I)(X)" means the I-th basis function
    evaluated at the point X.

    Input, int NP.
    The highest degree polynomial to use.

    Input, int NPRINT.
    The number of points at which the computed solution
    should be printed out at the end of the computation.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Input, int QUAD_NUM, the order of the quadrature rule.

    Input, double QUAD_W(QUAD_NUM), the quadrature weights.

    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
*/
{
  double big_l2;
  double error;
  int i;
  int ip;
  int j;
  int k;
  int nsub = 10;
  double phii;
  double phiix;
  double ue;
  double up;
  double x;
  double xl;
  double xr;

  printf ( "\n" );
  printf ( "  Comparison of computed and exact solutions:\n" );
  printf ( "\n" );
  printf ( "    X        U computed    U exact     Difference\n" );
  printf ( "\n" );

  for ( i = 0; i <= nprint; i++ )
  {
    x = ( double ) ( 2 * i - nprint ) / ( double ) ( nprint );
    ue = uex ( x, problem );
    up = 0.0;
    for ( j = 0; j <= np; j++ )
    {
      phi ( alpha, beta, j, np, &phii, &phiix, x );
      up = up + phii * f[j];
    }
    printf ( "  %8g  %12g  %12g  %12g\n", x, up, ue, ue - up );
  }
/*
  Compute the big L2 error.
*/
  big_l2 = 0.0;

  for ( i = 1; i <= nsub; i++ )
  {
    xl = ( double ) ( 2 * i - nsub - 1 ) / ( double ) ( nsub );
    xr = ( double ) ( 2 * i - nsub     ) / ( double ) ( nsub );

    for ( j = 0; j < quad_num; j++ )
    {
      x = ( xl * ( 1.0 - quad_x[j] )
          + xr * ( 1.0 + quad_x[j] ) ) / 2.0;

      up = 0.0;
      for ( k = 0; k <= np; k++ )
      {
        phi ( alpha, beta, k, np, &phii, &phiix, x );
        up = up + phii * f[k];
      }

      big_l2 = big_l2 + pow ( up - uex ( x, problem ), 2 ) * quad_w[j]
        * ( xr - xl ) / 2.0;
    }
  }

  big_l2 = sqrt ( big_l2 );

  printf ( "\n" );
  printf ( "  Big L2 error = %g\n", big_l2 );

  return;
}
/******************************************************************************/

double ff ( double x, int problem )

/******************************************************************************/
/*
  Purpose:

    FF evaluates the right hand side function F(X) at any point X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double X, the evaluation point.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Output, double FF, the value of F(X).
*/
{
  double pi = 3.141592653589793;
  double value;
/*
  Test problem 1
*/
  if ( problem == 1 )
  {
    value = 1.0 + 12.0 * x * x - x * x * x * x;
  }
/*
  Test problem 2
*/
  else if ( problem == 2 )
  {
    value = 0.25 * pi * pi * cos ( 0.5 * pi * x );
  }

  return value;
}
/******************************************************************************/

void ortho ( double a[], double alpha[], double beta[], int np, int problem,
  int quad_num, double quad_w[], double quad_x[] )

/******************************************************************************/
/*
  Purpose:

    ORTHO tests the basis functions for orthogonality.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double A(0:NP), the squares of the norms of the
    basis functions.

    Input, double ALPHA(NP), BETA(NP), the basis function 
    recurrence coefficients.

    Input, int NP.
    The highest degree polynomial to use.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Input, int QUAD_NUM, the order of the quadrature rule.

    Input, double QUAD_W(QUAD_NUM), the quadrature weights.

    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
*/
{
  double *b;
  double bij;
  int i;
  int iq;
  int j;
  double phii;
  double phiix;
  double phij;
  double phijx;
  double x;
/*
  Zero out the B array, so we can start summing up the dot products.
*/
  b = ( double * ) malloc ( ( np + 1 ) * ( np + 1 ) * sizeof ( double ) );
 
  for ( j = 0; j <= np; j++ )
  {
    for ( i = 0; i <= np; i++ )
    {
       b[i+j*(np+1)] = 0.0;
    }
  }
/*
  Approximate the integral of the product of basis function
  I and basis function J over the interval [-1,1].

  We expect to get zero, except when I and J are equal,
  when we should get A(I).
*/
  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      for ( j = 0; j <= np; j++ )
      {
        phi ( alpha, beta, j, np, &phij, &phijx, x );

        bij = pp ( x, problem ) * phiix * phijx
            + qq ( x, problem ) * phii * phij;

        b[i+j*(np+1)] = b[i+j*(np+1)] + bij * quad_w[iq];
      }
    }
  }
/*
  Print out the results of the test.
*/
  printf ( "\n" );
  printf ( "  Basis function orthogonality test:\n" );
  printf ( "\n" );
  printf ( "   i   j     b(i,j)/a(i)\n" );
  printf ( "\n" );
  for ( i = 0; i <= np; i++ )
  {
    printf ( "\n" );
    for ( j = 0; j <= np; j++ )
    {
      printf ( "  %6d  %6d  %12g\n", i, j, b[i+j*(np+1)] / a[i] );
    }
  }

  free ( b );

  return;
}
/******************************************************************************/

void out ( double alpha[], double beta[], double f[], int np, int nprint )

/******************************************************************************/
/*
  Purpose:

    OUT prints the computed solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double ALPHA(NP), BETA(NP), the basis function 
    recurrence coefficients.

    Input, double F(0:NP).
    F contains the basis function coefficients that form the
    representation of the solution U.  That is,
      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
    where "BASIS(I)(X)" means the I-th basis function
    evaluated at the point X.

    Input, int NP.
    The highest degree polynomial to use.

    Input, int NPRINT.
    The number of points at which the computed solution
    should be printed out at the end of the computation.
*/
{
  int i;
  int ip;
  double phii;
  double phiix;
  double up;
  double x;

  printf ( "\n" );
  printf ( "  Representation of solution:\n" );
  printf ( "\n" );
  printf ( "  Basis function coefficients:\n" );
  printf ( "\n" );
  for ( i = 0; i <= np; i++ )
  {
    printf ( "  %8d  %12g\n", i, f[i] );
  }

  printf ( "\n" );
  printf ( "\n" );
  printf ( "       X     Approximate Solution\n" );
  printf ( "\n" );
  for ( ip = 0; ip <= nprint; ip++ )
  {
    x = ( double ) ( 2 * ip - nprint ) / ( double ) ( nprint );
    up = 0.0;
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      up = up + phii * f[i];
    }
    printf ( "  %12g  %12g\n", x, up );
  }

  printf ( "\n" );

  return;
}
/******************************************************************************/

void phi ( double alpha[], double beta[], int i, int np, double *phii,
  double *phiix, double x )

/******************************************************************************/
/*
  Purpose:

    PHI evaluates the I-th basis function at the point X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double ALPHA(NP), BETA(NP), the basis function 
    recurrence coefficients.

    Input, int I, the index of the basis function.

    Input, int NP.
    The highest degree polynomial to use.

    Output, double PHII, PHIIX, the value of the basis
    function and its derivative.

    Input, double X, the evaluation point.
*/
{
  int j;
  double q;
  double qm1;
  double qm1x;
  double qm2;
  double qm2x;
  double qx;
  double t;

  qm1 = 0.0;
  q = 1.0;
  qm1x = 0.0;
  qx = 0.0;

  for ( j = 1; j <= i; j++ )
  {
    qm2 = qm1;
    qm1 = q;
    qm2x = qm1x;
    qm1x = qx;
    t = x - alpha[j-1];
    q = t * qm1 - beta[j-1] * qm2;
    qx = qm1 + t * qm1x - beta[j-1] * qm2x;
  }

  t = 1.0 - x * x;
  *phii = t * q;
  *phiix = t * qx - 2.0 * x * q;

  return;
}
/******************************************************************************/

double pp ( double x, int problem )

/******************************************************************************/
/*
  Purpose:

    PP returns the value of the coefficient function P(X).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double X, the evaluation point.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Output, double PP, the value of P(X).
*/
{
  double value;
/*
  Test problem 1
*/
  if ( problem == 1 )
  {
    value = 1.0;
  }
/*
  Test problem 2
*/
  else if ( problem == 2 )
  {
    value = 1.0;
  }

  return value;
}
/******************************************************************************/

void quad ( int quad_num, double quad_w[], double quad_x[] )

/******************************************************************************/
/*
  Purpose:

    QUAD returns the abscissas and weights for gaussian quadrature on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, int QUAD_NUM, the order of the quadrature rule.

    Output, double QUAD_W(QUAD_NUM), the quadrature weights.

    Output, double QUAD_X(QUAD_NUM), the quadrature abscissas.
*/
{
/*
  Quadrature points on [-1,1]
*/
  quad_x[0] = -0.973906528517172;
  quad_x[1] = -0.865063366688985;
  quad_x[2] = -0.679409568299024;
  quad_x[3] = -0.433395394129247;
  quad_x[4] = -0.148874338981631;
  quad_x[5] =  0.148874338981631;
  quad_x[6] =  0.433395394129247;
  quad_x[7] =  0.679409568299024;
  quad_x[8] =  0.865063366688985;
  quad_x[9] = 0.973906528517172;
/*
  Weight factors
*/
  quad_w[0] =  0.066671344308688;
  quad_w[1] =  0.149451349150581;
  quad_w[2] =  0.219086362515982;
  quad_w[3] =  0.269266719309996;
  quad_w[4] =  0.295524224714753;
  quad_w[5] =  0.295524224714753;
  quad_w[6] =  0.269266719309996;
  quad_w[7] =  0.219086362515982;
  quad_w[8] =  0.149451349150581;
  quad_w[9] = 0.066671344308688;

  return;
}
/******************************************************************************/

double qq ( double x, int problem )

/******************************************************************************/
/*
  Purpose:

    QQ returns the value of the coefficient function Q(X).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double X, the evaluation point.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Output, double QQ, the value of Q(X).
*/
{
  double value;
/*
  Test problem 1
*/
  if ( problem == 1 )
  {
    value = 1.0;
  }
/*
  Test problem 2
*/
  else if ( problem == 2 )
  {
    value = 0.0;
  }

  return value;
}
/******************************************************************************/

void sol ( double a[], double alpha[], double beta[], double f[], int np,
  int problem, int quad_num, double quad_w[], double quad_x[] )

/******************************************************************************/
/*
  Purpose:

    SOL solves a linear system for the finite element coefficients.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double A(0:NP), the squares of the norms of the
    basis functions.

    Input, double ALPHA(NP), BETA(NP), the basis function 
    recurrence coefficients.

    Output, double F(0:NP).
    F contains the basis function coefficients that form the
    representation of the solution U.  That is,
      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
    where "BASIS(I)(X)" means the I-th basis function
    evaluated at the point X.

    Input, int NP.
    The highest degree polynomial to use.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Input, int QUAD_NUM, the order of the quadrature rule.

    Input, double QUAD_W(QUAD_NUM), the quadrature weights.

    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
*/
{
  int i;
  int iq;
  double phii;
  double phiix;
  double t;
  double x;

  for ( i = 0; i <= np; i++ )
  {
    f[i] = 0.0;
  }

  for ( iq = 0; iq < quad_num; iq++ )
  {
    x = quad_x[iq];
    t = ff ( x, problem ) * quad_w[iq];
    for ( i = 0; i <= np; i++ )
    {
      phi ( alpha, beta, i, np, &phii, &phiix, x );
      f[i] = f[i] + phii * t;
    }
  }

  for ( i = 0; i <= np; i++ )
  {
    f[i] = f[i] / a[i];
  }

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double uex ( double x, int problem )

/******************************************************************************/
/*
  Purpose:

    UEX returns the value of the exact solution at a point X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
    C version by John Burkardt.

  Parameters:

    Input, double X, the evaluation point.

    Input, int PROBLEM, indicates the problem being solved.
    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).

    Output, double UEX, the exact value of U(X).
*/
{
  double pi = 3.141592653589793;
  double value;
/*
  Test problem 1
*/
  if ( problem == 1 )
  {
    value = 1.0 - pow ( x, 4 );
  }
/*
  Test problem 2
*/
  else if ( problem == 2 )
  {
    value = cos ( 0.5 * pi * x );
  }

  return value;
}
