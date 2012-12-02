# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms655.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( int nt, int kind, double alpha, double beta );
double f ( double x, int i );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS655_PRB.

  Discussion:

    TOMS655_PRB calls a set of problems for DIVDIF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  int kind;
  int nt;

  timestamp ( );

  printf ( "\n" );
  printf ( "TOMS655_PRB\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Tests for the routines in the TOMS655 library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
/*
  Compute 15 points of an example of each rule.
*/
  for ( kind = 1; kind <= 8; kind++ )
  {
    nt = 15;
    if ( kind == 8 )
    {
      alpha = 1.0;
      beta = - alpha - 2 * nt - 2;
    }
    else
    {
      alpha = 0.0;
      beta = 0.0;
    }
    test10 ( nt, kind, alpha, beta );
  }

  printf ( "\n" );
  printf ( "TOMS655_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests CIQFS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double alpha;
  double beta;
  int i;
  int key;
  int kind;
  int lu;
  int *mlt;
  int *ndx;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double *t;
  double *wts;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test CIQFS.\n" );
/*
  Number of knots.
*/
  nt = 5;
/*
  Set the knots in the default interval [-1,+1].
*/
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
/*
  Set the knot multiplicities.
*/
  mlt = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
/*
  Set the size of the weights array.
*/
  nwts = 0;
  for ( i = 0; i < nt; i++ )
  {
    nwts = nwts + mlt[i];
  }
/*
  Because KEY = 1, NDX will be set up for us.
*/
  ndx = ( int * ) malloc ( nt * sizeof ( int ) );
/*
  KEY = 1 indicates that the WTS array should hold the weights
  in the usual order.
*/
  key = 1;
/*
  Request Legendre weight function.
*/
  kind = 1;
/*
  ALPHA, BETA not used in Legendre weight function but set anyway.
*/
  alpha = 0.0;
  beta  = 0.0;
/*
  LU controls printing.
  A positive value requests that we compute and print weights, and
  conduct a moments check.
*/
  lu = 6;
/*
  This call returns the WTS array.
*/
  wts = ciqfs ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, lu );

  free ( mlt );
  free ( ndx );
  free ( t );
  free ( wts );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CIQFS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int key;
  int kind;
  int lu;
  int *mlt;
  int *ndx;
  int nt;
  int nwts;
  double *t;
  double *wts;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test CIQF, CIQFS, CGQF and CGQFS\n" );
  printf ( "  with all classical weight functions.\n" );
/*
  Try all 8 weight functions.
*/
  for ( kind = 1; kind <= 8; kind++ )
  {
/*
  Number of knots.
*/
    nt = 5;
/*
  Set parameters ALPHA and BETA.
*/
    alpha = 0.5;
    if ( kind != 8 )
    {
      beta  = 2.0;
    }
    else
    {
      beta = - 16.0;
    }
/*
  Set A and B.
*/
    a = - 0.5;
    b = 2.0;
/*
  LU controls printing.
  A positive value requests that we compute and print weights, and
  conduct a moments check.
*/
    lu = 6;
/*
  Have CGQF compute the knots and weights.
*/
    t = ( double * ) malloc ( nt * sizeof ( double ) );
    wts = ( double * ) malloc ( nt * sizeof ( double ) );

    printf ( "\n" );
    printf ( "  Knots and weights of Gauss quadrature formula\n" );
    printf ( "  computed by CGQF.\n" );
    cgqf ( nt, kind, alpha, beta, a, b, lu, t, wts );
/*
  Now compute the weights for the same knots by CIQF.

  Set the knot multiplicities.
*/
    mlt = ( int * ) malloc ( nt * sizeof ( int ) );
    for ( i = 0; i < nt; i++ )
    {
      mlt[i] = 2;
    }
/*
  Set the size of the weights array.
*/
    nwts = 0;
    for ( i = 0; i < nt; i++ )
    {
      nwts = nwts + mlt[i];
    }
/*
  We need to deallocate and reallocate WTS because it is now of
  dimension NWTS rather than NT.
*/
    free ( wts );
    wts = ( double * ) malloc ( nwts * sizeof ( double ) );
/*
  Because KEY = 1, NDX will be set up for us.
*/
    ndx = ( int * ) malloc ( nt * sizeof ( int ) );
/*
  KEY = 1 indicates that the WTS array should hold the weights
  in the usual order.
*/
    key = 1;

    printf ( "\n" );
    printf ( "  Weights of Gauss quadrature formula computed from the\n" );
    printf ( "  knots by CIQF.\n" );

    wts = ciqf ( nt, t, mlt, nwts, ndx, key, kind, alpha, beta, a, b, lu );

    free ( mlt );
    free ( ndx );
    free ( t );
    free ( wts );
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CEIQFS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double alpha;
  double beta;
  int i;
  int kind;
  int lu;
  int *mlt;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Test CEIQFS.\n" );
/*
  Number of knots.
*/
  nt = 5;
/*
  Set the knots in the default interval [-1,+1].
*/
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
/*
  Set the knot multiplicities.
*/
  mlt = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
/*
  Set KIND to the Legendre weight function.
*/
  kind = 1;
/*
  ALPHA, BETA not used in Legendre weight function but set anyway.
*/
  alpha = 0.0;
  beta  = 0.0;
/*
  Call CEIQFS to set up the quadrature formula and evaluate it on F.
*/
  qfsum = ceiqfs ( nt, t, mlt, kind, alpha, beta, f );

  printf ( "\n" );
  printf ( "  Integral of sin(x) on -1, 1 by Fejer type rule\n" );
  printf ( "  with %d points of multiplicity 2.\n", nt );
  printf ( "  Quadrature formula:%24.16f\n", qfsum );

  qfsx = cos ( - 1.0 ) - cos ( 1.0 );
  printf ( "  Exact value       :%24.16f\n", qfsx );
  printf ( "  Error             :%13.3e\n", r8_abs ( qfsum - qfsx ) );

  free ( mlt );
  free ( t );

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CEIQF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int kind;
  int lu;
  int *mlt;
  int nt;
  int nwts;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test CEIQF.\n" );
/*
  Number of knots.
*/
  nt = 5;
/*
  Set the knots in the default interval [-1,+1].
*/
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
/*
  Set the knot multiplicities.
*/
  mlt = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 2;
  }
/*
  Set KIND to the Legendre weight function.
*/
  kind = 1;
/*
  ALPHA, BETA not used in Legendre weight function but set anyway.
*/
  alpha = 0.0;
  beta  = 0.0;
/*
  Set nonstandard interval A, B.
*/
  a = -0.5;
  b = 2.0;
/*
  Shift knots from [-1,1] to [A,B].
*/
  for ( i = 0; i < nt; i++ )
  {
    t[i] = ( ( b - a ) * t[i] + ( a + b ) ) / 2.0;
  }
/*
  Call CEIQF to set up the quadrature formula and evaluate it on F.
*/
  qfsum = ceiqf ( nt, t, mlt, kind, alpha, beta, a, b, f );

  printf ( "\n" );
  printf ( "  Integral of sin(x) from %f to %f\n", a, b );
  printf ( "  by Fejer type rule with %d points\n", nt );
  printf ( "  of multiplicity 2.\n" );
  printf ( "  Quadrature formula:%24.16f\n", qfsum );

  qfsx = cos ( a ) - cos ( b );
  printf ( "  Exact value       :%24.16f\n", qfsx );
  printf ( "  Error             :%13.3e\n", r8_abs ( qfsum - qfsx ) );

  free ( mlt );
  free ( t );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CLIQFS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double alpha;
  double beta;
  int i;
  int kind;
  int lu;
  int nt;
  double pi = 3.14159265358979323846264338327950;
  double *t;
  double *wts;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Test CLIQFS.\n" );
/*
  Number of knots.
*/
  nt = 5;
/*
  Set the knots in the default interval [-1,+1].
*/
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
/*
  Request Legendre weight function.
*/
  kind = 1;
/*
  ALPHA, BETA not used in Legendre weight function but set anyway.
*/
  alpha = 0.0;
  beta  = 0.0;
/*
  LU controls printing.
  A positive value requests that we compute and print weights, and
  conduct a moments check.
*/
  lu = 6;
/*
  This call returns the WTS array.
*/
  wts = cliqfs ( nt, t, kind, alpha, beta, lu );

  free ( t );
  free ( wts );
 
  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests CEIQF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int kind;
  int lu;
  int nt;
  double pi = 3.14159265358979323846264338327950;
  double qfsum;
  double qfsx;
  double *t;
  double *wts;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Test CLIQF and EIQFS.\n" );
/*
  Number of knots.
*/
  nt = 5;
/*
  Set the knots in the default interval [-1,+1].
*/
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = cos ( ( double ) ( 2 * i - 1 ) * pi / ( double ) ( 2 * nt ) );
  }
/*
  Set KIND to the Legendre weight function.
*/
  kind = 1;
/*
  ALPHA, BETA not used in Legendre weight function but set anyway.
*/
  alpha = 0.0;
  beta  = 0.0;
/*
  Set nonstandard interval A, B.
*/
  a = -0.5;
  b = 2.0;
/*
  Shift knots from [-1,1] to [A,B].
*/
  for ( i = 0; i < nt; i++ )
  {
    t[i] = ( ( b - a ) * t[i] + ( a + b ) ) / 2.0;
  }
/*
  LU controls printout.
*/
  lu = 6;
/*
  Call CLIQF to set up the quadrature formula.
*/
  wts = cliqf ( nt, t, kind, alpha, beta, a, b, lu );
/*
  Call EIQFS to evaluate the quadrature formula.
*/
  qfsum = eiqfs ( nt, t, wts, f );

  printf ( "\n" );
  printf ( "  Integral of sin(x) from %f to %f\n", a, b );
  printf ( "  by Fejer type rule with %d points\n", nt );
  printf ( "  of multiplicity 1.\n" );
  printf ( "  Quadrature formula:%24.16f\n", qfsum );

  qfsx = cos ( a ) - cos ( b );
  printf ( "  Exact value       :%24.16f\n", qfsx );
  printf ( "  Error             :%13.3e\n", r8_abs ( qfsum - qfsx ) );

  free ( t );
  free ( wts );

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests CEGQF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int kind;
  int nt;
  double qfsum;
  double qfsx;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Test CEGQF.\n" );
/*
  Number of knots.
*/
  nt = 12;
/*
  Request exponential weight function.
*/
  kind = 7;
/*
  Set ALPHA and BETA.
*/
  alpha = 1.0;
  beta  = 0.0;
/*
  Set interval [A,B].
*/
  a = -0.5;
  b = 2.0;
/*
  Call CEGQF to compute and evaluate the Gauss quadrature formula.
*/
  qfsum = cegqf ( nt, kind, alpha, beta, a, b, f );

  printf ( "\n" );
  printf ( "  Integral of x*sin(x) from %f to %f\n", a, b );
  printf ( "  by Gauss-exponential rule with %d points\n", nt );
  printf ( "  Quadrature formula:%24.16f\n", qfsum );

  qfsx = ( b - a ) * 0.5 * ( cos ( a ) - cos ( b ) ) 
    + sin ( b ) + sin ( a ) - 2.0 * sin ( ( a + b ) / 2.0 );

  printf ( "  Exact value       :%24.16f\n", qfsx );
  printf ( "  Error             :%13.3e\n", r8_abs ( qfsum - qfsx ) );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests CEGQFS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
*/
{
  double alpha;
  double beta;
  int kind;
  int nt;
  double qfsum;
  double qfsx;

  printf ( "  ----------------------------------------\n" );
  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Test CEGQFS.\n" );
/*
  Number of knots.
*/
  nt = 12;
/*
  Request exponential weight function.
*/
  kind = 7;
/*
  Set ALPHA and BETA.
*/
  alpha = 1.0;
  beta  = 0.0;
/*
  Call CEGQFS to compute and evaluate the Gauss quadrature formula.
*/
  qfsum = cegqfs ( nt, kind, alpha, beta, f );

  printf ( "\n" );
  printf ( "  Integral of x*sin(x) from -1 to +1\n" );
  printf ( "  by Gauss-exponential rule with %d points.\n", nt );
  printf ( "  Quadrature formula:%24.16f\n", qfsum );

  qfsx = cos ( -1.0 ) - cos ( +1.0 );

  printf ( "  Exact value       :%24.16f\n", qfsx );
  printf ( "  Error             :%13.3e\n", r8_abs ( qfsum - qfsx ) );

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 calls CGQFS to compute and print generalized Gauss-Hermite rules.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  int io;
  int kind;
  int nt;
  double *t;
  double *wts;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  Call CGQFS to compute generalized Hermite rules.\n" );

  nt = 15;
  kind = 6;
  alpha = 1.0;
  beta = 0.0;
  io = - 6;
  t = ( double * ) malloc ( nt * sizeof ( double ) );
  wts = ( double * ) malloc ( nt * sizeof ( double ) );

  printf ( "\n" );
  printf ( "  NT = %d\n", nt );
  printf ( "  ALPHA = %f\n", alpha );

  cgqfs ( nt, kind, alpha, beta, io, t, wts );

  free ( t );
  free ( wts );

  return;
}
/******************************************************************************/

void test10 ( int nt, int kind, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    TEST10 calls CDGQF to compute a quadrature formula.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    John Burkardt
*/
{
  int i;
  double *t;
  double *wts;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  Call CDGQF to compute a quadrature formula.\n" );
  printf ( "\n" );
  printf ( "  KIND = %d\n", kind );
  printf ( "  ALPHA = %f\n", alpha );
  printf ( "  BETA  = %f\n", beta );

  t = ( double * ) malloc ( nt * sizeof ( double ) );
  wts = ( double * ) malloc ( nt * sizeof ( double ) );

  cdgqf ( nt, kind, alpha, beta, t, wts );

  printf ( "\n" );
  printf ( " Index     Abscissas                 Weights\n" );
  printf ( "\n" );
  for ( i = 0; i < nt; i++ )
  {
    printf ( "  %4d  %24.16g  %24.16g\n", i, t[i], wts[i] );
  }

  free ( t );
  free ( wts );

  return;
}
/******************************************************************************/

double f ( double x, int i )

/******************************************************************************/
/*
  Purpose:

    F returns values of the integrand or its derivatives.

  Discussion:

    This function is an example of an integrand function.

    The package can generate quadrature formulas that use derivative 
    information as well as function values.  Therefore, this routine is
    set up to provide derivatives of any order as well as the function
    value.  In an actual application, the highest derivative needed
    is of order one less than the highest knot multiplicity.

    In other words, in the usual case where knots are not repeated,
    this routine only needs to return function values, not any derivatives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Parameters:

    Input, double X, the evaluation point.

    Input, int I, the order of the derivative of F to
    be evaluated.

    Output, double F, the value of the I-th derivative of F at X.
*/
{
  int l;
  double value;

  l = ( i % 4 );

  if ( l == 0 )
  {
    value = sin ( x );
  }
  else if ( l == 1 )
  {
    value = cos ( x );
  }
  else if ( l == 2 )
  {
    value = - sin ( x );
  }
  else if ( l == 3 )
  {
    value = - cos ( x );
  }

  return value;
}
