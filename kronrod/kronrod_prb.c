# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "kronrod.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
double f ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for KRONROD_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 April 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "KRONROD_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the KRONROD library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "KRONROD_PRB:\n" );
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

    TEST01 tests the code for the odd case N = 3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 August 2010

  Author:

    John Burkardt
*/
{
  double eps;
  int i;
  int i2;
  int n = 3;
  double s;
  double *w1;
  double *w2;
  double wg[3] = {
    0.555555555555555555556, 
    0.888888888888888888889, 
    0.555555555555555555556 };
  double wk[7] = {
    0.104656226026467265194, 
    0.268488089868333440729, 
    0.401397414775962222905, 
    0.450916538658474142345, 
    0.401397414775962222905, 
    0.268488089868333440729, 
    0.104656226026467265194 };
  double *x;
  double xg[3] = {
   -0.77459666924148337704, 
    0.0, 
    0.77459666924148337704 };
  double xk[7]= {
   -0.96049126870802028342, 
   -0.77459666924148337704, 
   -0.43424374934680255800, 
    0.0, 
    0.43424374934680255800, 
    0.77459666924148337704, 
    0.96049126870802028342 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Request KRONROD to compute the Gauss rule\n" );
  printf ( "  of order 3, and the Kronrod extension of\n" );
  printf ( "  order 3+4=7.\n" );
  printf ( "\n" );
  printf ( "  Compare to exact data.\n" );

  eps = 0.000001;
  w1 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  w2 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  x = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  kronrod ( n, eps, x, w1, w2 );

  printf ( "\n" );
  printf ( "  KRONROD returns 3 vectors of length %d\n", n + 1 );
  printf ( "\n" );
  printf ( "     I      X               WK              WG\n" );
  printf ( "\n" );
  for ( i = 1; i <= n + 1; i++ )
  {
    printf ( "  %4d  %14f  %14f  %14f\n", i, x[i-1], w1[i-1], w2[i-1] );
  }

  printf ( "\n" );
  printf ( "               Gauss Abscissas\n" );
  printf ( "            Exact           Computed\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
      s = +1.0;
    }
    printf ( "  %4d  %14f  %14f\n", i, xg[i-1], s * x[i2-1] );
  }
  printf ( "\n" );
  printf ( "               Gauss Weights\n" );
  printf ( "            Exact           Computed\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
    }
    printf ( "  %4d  %14f  %14f\n", i, wg[i-1], w2[i2-1] );
  }

  printf ( "\n" );
  printf ( "             Gauss Kronrod Abscissas\n" );
  printf ( "            Exact           Computed\n" );
  printf ( "\n" );
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
      s = +1.0;
    }
    printf ( "  %4d  %14f  %14f\n", i, xk[i-1], s * x[i2-1] );
  }
  printf ( "\n" );
  printf ( "             Gauss Kronrod Weights\n" );
  printf ( "            Exact           Computed\n" );
  printf ( "\n" );
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
    }
    printf ( "  %4d  %14f  %14f\n", i, wk[i-1], w1[i2-1] );
  }

  free ( w1 );
  free ( w2 );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests the code for the even case N = 4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 August 2010

  Author:

    John Burkardt
*/
{
  double eps;
  int i;
  int i2;
  int n = 4;
  double s;
  double *w1;
  double *w2;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Request KRONROD to compute the Gauss rule\n" );
  printf ( "  of order 4, and the Kronrod extension of\n" );
  printf ( "  order 4+5=9.\n" );

  eps = 0.000001;
  w1 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  w2 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  x = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  kronrod ( n, eps, x, w1, w2 );

  printf ( "\n" );
  printf ( "  KRONROD returns 3 vectors of length %d\n", n + 1 );
  printf ( "\n" );
  printf ( "     I      X               WK              WG\n" );
  printf ( "\n" );
  for ( i = 1; i <= n + 1; i++ )
  {
    printf ( "  %4d  %14f  %14f  %14f\n", i, x[i-1], w1[i-1], w2[i-1] );
  }

  free ( w1 );
  free ( w2 );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses the program to estimate an integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 April 2012

  Author:

    John Burkardt
*/
{
  double eps;
  double exact = 1.5643964440690497731;
  int i;
  double i1;
  double i2;
  int n;
  double *w1;
  double *w2;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Call Kronrod to estimate the integral of a function.\n" );
  printf ( "  Keep trying until the error is small.\n" );
/*
  EPS just tells KRONROD how carefully it must compute X, W1 and W2.
  It is NOT a statement about the accuracy of your integral estimate!
*/
  eps = 0.000001;
/*
  Start the process with a 1 point rule.
*/
  n = 1;

  for ( ; ; )
  {
/*
  Make space.
*/
    w1 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    w2 = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    x = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

    kronrod ( n, eps, x, w1, w2 );
/*
  Compute the estimates.
  There are two complications here:

  1) Both rules use all the points.  However, the lower order rule uses
     a zero weight for the points it doesn't need.

  2) The points X are all positive, and are listed in descending order.
     this means that 0 is always in the list, and always occurs as the
     last member.  Therefore, the integral estimates should use the
     function value at 0 once, and the function values at the other
     X values "twice", that is, once at X and once at -X.
*/
    i1 = w1[n] * f ( x[n] );
    i2 = w2[n] * f ( x[n] );

    for ( i = 0; i < n; i++ )
    {
      i1 = i1 + w1[i] * ( f ( - x[i] ) + f ( x[i] ) );
      i2 = i2 + w2[i] * ( f ( - x[i] ) + f ( x[i] ) );
    }

    free ( w1 );
    free ( w2 );
    free ( x );

    if ( fabs ( i1 - i2 ) < 0.0001 )
    {
      printf ( "\n" );
      printf ( "  Error tolerance satisfied with N = %d\n", n );
      printf ( "  Coarse integral estimate = %14.6g\n", i1 );
      printf ( "  Fine   integral estimate = %14.6g\n", i2 );
      printf ( "  Error estimate = %14.6g\n", fabs ( i2 - i1 ) );
      printf ( "  Actual error = %14.6g\n", fabs ( exact - i2 ) );
      break;
    }

    if ( 25 < n )
    {
      printf ( "\n" );
      printf ( "  Error tolerance failed even for n = %d\n", n );
      printf ( "  Canceling iteration, and accepting bad estimates!\n" );
      printf ( "  Coarse integral estimate = %14.6g\n", i1 );
      printf ( "  Fine   integral estimate = %14.6g\n", i2 );
      printf ( "  Error estimate = %14.6g\n", fabs ( i2 - i1 ) );
      printf ( "  Actual error = %14.6g\n", fabs ( exact - i2 ) );
      break;
    }
    n = 2 * n + 1;
  }

  return;
}
/******************************************************************************/

double f ( double x )

/******************************************************************************/
/*
  Purpose:

    F is a function whose integral from -1 to +1 is to be estimated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double F, the value.
*/
{
  double value;

  value = 1.0 / ( x * x + 1.005 );

  return value;
}
