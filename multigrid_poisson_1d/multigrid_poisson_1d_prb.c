# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "multigrid_poisson_1d.h"

int main ( );

void test01_mono ( );
void test01_multi ( );
void test02_mono ( );
void test02_multi ( );
double exact1 ( double x );
double force1 ( double x );
double exact2 ( double x );
double force2 ( double x );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MULTIGRID_POISSON_1D_PRB.

  Discussion:

    MULTIGRID_POISSON_1D_PRB tests the MULTIGRID_POISSON_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "MULTIGRID_POISSON_1D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MULTIGRID_POISSON_1D library.\n" );

  test01_mono ( );
  test01_multi ( );
  test02_mono ( );
  test02_multi ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MULTIGRID_POISSON_1D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01_mono ( ) 

/******************************************************************************/
/*
  Purpose:

    TEST01_MONO tests MONOGRID_POISSON_1D on test case 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double ua;
  double ub;
  double *x;

  printf ( "\n" );
  printf ( "TEST01_MONO\n" );
  printf ( "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n" );
  printf ( "  using the Gauss-Seidel method.\n" );

  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;

  printf ( "\n" );
  printf ( "  -u''(x) = 1, for 0 < x < 1\n" );
  printf ( "  u(0) = u(1) = 0.\n" );
  printf ( "  Solution is u(x) = ( -x^2 + x ) / 2\n" );

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    x = r8vec_linspace_new ( n + 1, a, b );

    printf ( "\n" );
    printf ( "  Mesh index K = %d\n", k );
    printf ( "  Number of intervals N=2^K = %d\n", n );
    printf ( "  Number of nodes = 2^K+1 =   %d\n", n + 1 );

    monogrid_poisson_1d ( n, a, b, ua, ub, force1, exact1, &it_num, u );

    printf ( "\n" );
    printf ( "     I        X(I)      U(I)         U Exact(X(I))\n" );
    printf ( "\n" );
    for ( i = 0; i < n + 1; i++ )
    {
      printf ( "  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact1 ( x[i] ) );
    }

    printf ( "\n" );

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      difmax = r8_max ( difmax, fabs ( u[i] - exact1 ( x[i] ) ) );
    } 
    printf ( "  Maximum error = %g\n", difmax );
    printf ( "  Number of iterations = %d\n", it_num );

    free ( u );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test01_multi ( ) 

/******************************************************************************/
/*
  Purpose:

    TEST01_MULTI tests MULTIGRID_POISSON_1D on test case 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double ua;
  double ub;
  double *x;

  printf ( "\n" );
  printf ( "TEST01_MULTI\n" );
  printf ( "  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n" );
  printf ( "  using the multigrid method.\n" );

  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;

  printf ( "\n" );
  printf ( "  -u''(x) = 1, for 0 < x < 1\n" );
  printf ( "  u(0) = u(1) = 0.\n" );
  printf ( "  Solution is u(x) = ( -x^2 + x ) / 2\n" );

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    x = r8vec_linspace_new ( n + 1, a, b );

    printf ( "\n" );
    printf ( "  Mesh index K = %d\n", k );
    printf ( "  Number of intervals N=2^K = %d\n", n );
    printf ( "  Number of nodes = 2^K+1 =   %d\n", n + 1 );

    multigrid_poisson_1d ( n, a, b, ua, ub, force1, exact1, &it_num, u );

    printf ( "\n" );
    printf ( "     I        X(I)      U(I)         U Exact(X(I))\n" );
    printf ( "\n" );
    for ( i = 0; i < n + 1; i++ )
    {
      printf ( "  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact1 ( x[i] ) );
    }

    printf ( "\n" );

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      difmax = r8_max ( difmax, fabs ( u[i] - exact1 ( x[i] ) ) );
    } 
    printf ( "  Maximum error = %g\n", difmax );
    printf ( "  Number of iterations = %d\n", it_num );

    free ( u );
    free ( x );
  }
  return;
}
/******************************************************************************/

double exact1 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT1 evaluates the exact solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2011

  Author:

    John Burkardt

  Reference:

    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT1, the value of the exact solution at X.
*/
{
  double value;

  value = 0.5 * ( - x * x + x );

  return value;
}
/******************************************************************************/

double force1 ( double x )

/******************************************************************************/
/*
  Purpose:

    FORCE1 evaluates the forcing function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2011

  Author:

    John Burkardt

  Reference:

    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.

  Parameters:

    Input, double X, the evaluation point.

    Output, double FORCE1, the value of the forcing function at X.
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

void test02_mono ( ) 

/******************************************************************************/
/*
  Purpose:

    TEST02_MONO tests MONOGRID_POISSON_1D on test case 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double ua;
  double ub;
  double *x;

  printf ( "\n" );
  printf ( "TEST02_MONO\n" );
  printf ( "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n" );
  printf ( "  using the Gauss-Seidel method.\n" );

  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;

  printf ( "\n" );
  printf ( "  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n" );
  printf ( "  u(0) = u(1) = 0.\n" );
  printf ( "  Solution is u(x) = x * (x-1) * exp(x)\n" );

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    x = r8vec_linspace_new ( n + 1, a, b );

    printf ( "\n" );
    printf ( "  Mesh index K = %d\n", k );
    printf ( "  Number of intervals N=2^K = %d\n", n );
    printf ( "  Number of nodes = 2^K+1 =   %d\n", n + 1 );

    monogrid_poisson_1d ( n, a, b, ua, ub, force2, exact2, &it_num, u );

    printf ( "\n" );
    printf ( "     I        X(I)      U(I)         U Exact(X(I))\n" );
    printf ( "\n" );
    for ( i = 0; i < n + 1; i++ )
    {
      printf ( "  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact2 ( x[i] ) );
    }

    printf ( "\n" );

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      difmax = r8_max ( difmax, fabs ( u[i] - exact2 ( x[i] ) ) );
    } 
    printf ( "  Maximum error = %g\n", difmax );
    printf ( "  Number of iterations = %d\n", it_num );

    free ( u );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test02_multi ( ) 

/******************************************************************************/
/*
  Purpose:

    TEST02_MULTI tests MULTIGRID_POISSON_1D on test case 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double ua;
  double ub;
  double *x;

  printf ( "\n" );
  printf ( "TEST02_MULTI\n" );
  printf ( "  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n" );
  printf ( "  using the multigrid method.\n" );

  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;

  printf ( "\n" );
  printf ( "  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n" );
  printf ( "  u(0) = u(1) = 0.\n" );
  printf ( "  Solution is u(x) = x * (x-1) * exp(x)\n" );

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    x = r8vec_linspace_new ( n + 1, a, b );

    printf ( "\n" );
    printf ( "  Mesh index K = %d\n", k );
    printf ( "  Number of intervals N=2^K = %d\n", n );
    printf ( "  Number of nodes = 2^K+1 =   %d\n", n + 1 );

    multigrid_poisson_1d ( n, a, b, ua, ub, force2, exact2, &it_num, u );

    printf ( "\n" );
    printf ( "     I        X(I)      U(I)         U Exact(X(I))\n" );
    printf ( "\n" );
    for ( i = 0; i < n + 1; i++ )
    {
      printf ( "  %4d  %10f  %14g  %14g\n", i, x[i], u[i], exact2 ( x[i] ) );
    }

    printf ( "\n" );

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      difmax = r8_max ( difmax, fabs ( u[i] - exact2 ( x[i] ) ) );
    } 
    printf ( "  Maximum error = %g\n", difmax );
    printf ( "  Number of iterations = %d\n", it_num );

    free ( u );
    free ( x );
  }
  return;
}
/******************************************************************************/

double exact2 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT2 evaluates the exact solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2011

  Author:

    John Burkardt

  Reference:

    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT2, the value of the exact solution at X.
*/
{
  double value;

  value = x * ( x - 1.0 ) * exp ( x );

  return value;
}
/******************************************************************************/

double force2 ( double x )

/******************************************************************************/
/*
  Purpose:

    FORCE2 evaluates the forcing function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2011

  Author:

    John Burkardt

  Reference:

    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.

  Parameters:

    Input, double X, the evaluation point.

    Output, double FORCE2, the value of the forcing function at X.
*/
{
  double value;

  value = - x * ( x + 3.0 ) * exp ( x );

  return value;
}
