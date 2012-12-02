# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "sine_transform.h"

int main ( );
void sine_transform_test01 ( );
void sine_transform_test02 ( );
void sine_transform_test03 ( );
void sine_transform_test04 ( );
void sine_transform_test05 ( );
double cosine_sum ( double x );
double poly5 ( double x );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SINE_TRANSFORM_TEST tests SINE_TRANSFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the SINE_TRANSFORM library.\n" );

  sine_transform_test01 ( );
  sine_transform_test02 ( );
  sine_transform_test03 ( );
  sine_transform_test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void sine_transform_test01 ( )

/******************************************************************************/
/*
  Purpose:

    SINE_TRANSFORM_TEST01 demonstrates that the transform is its own inverse.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  int n = 10;
  int seed;
  double *r;
  double *s;
  double *t;

  seed = 123456789;

  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST01:\n" );
  printf ( "  SINE_TRANSFORM_DATA does a sine transform of data\n" );
  printf ( "  defined by a vector.\n" );
  printf ( "\n" );
  printf ( "  Demonstrate that the transform is its own inverse.\n" );
  printf ( "  Let R be a random N vector.\n" );
  printf ( "  Let S be the transform of D.\n" );
  printf ( "  Let T be the transform of E.\n" );
  printf ( "  Then R and T will be equal.\n" );

  r = r8vec_uniform_01_new ( n, &seed );
  s = sine_transform_data ( n, r );
  t = sine_transform_data ( n, s );

  printf ( "\n" );
  printf ( "     I      R(I)        S(I)        T(I)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f  %10f  %10f\n", i, r[i], s[i], t[i] );
  }

  free ( r );
  free ( s );
  free ( t );

  return;
}
/******************************************************************************/

void sine_transform_test02 ( )

/******************************************************************************/
/*
  Purpose:

    SINE_TRANSFORM_TEST02 uses the functional form of the routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *f1;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 9;
  double *s;
  double *x;

  a = 1.0;
  b = 3.0;
/*
  Evenly spaced points between A and B, but omitting
  A and B themselves.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }

  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST02:\n" );
  printf ( "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n" );
  printf ( "  defined by a function F(X) evaluated at equally spaced\n" );
  printf ( "  points in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  Demonstrate that the transform is its own inverse.\n" );
  printf ( "  Let X(0:N+1) be N+2 equally spaced points in [A,B].\n" );
  printf ( "  Let S be the transform of F(X(1:N)).\n" );
  printf ( "  Let F1 be the linear interpolant of (A,F(A)), (B,F(B)).\n" );
  printf ( "  Let F2 be the transform of S.\n" );
  printf ( "  Then F(X(1:N)) = F1(X(1:N)) + F2(1:N).\n" );

  s = sine_transform_function ( n, a, b, poly5 );

  fa = poly5 ( a );
  fb = poly5 ( b );

  f1 = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    f1[i] = ( ( b - x[i]     ) * fa   
            + (     x[i] - a ) * fb ) 
            / ( b        - a );
  }

  f2 = sine_transform_data ( n, s );

  printf ( "\n" );
  printf ( "     I      X(I)      F(X(I))       S           F1          F2          F1+F2\n" );
  printf ( "\n" );
  printf ( "  %4d  %10f  %10f  %10f  %10f  %10f  %10f\n", 
    0, a, poly5 ( a ), 0.0, fa, 0.0, fa );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f  %10f  %10f  %10f  %10f  %10f\n", 
      i + 1, x[i], poly5 ( x[i] ), s[i], f1[i], f2[i], f1[i] + f2[i] );
  }

  printf ( "  %4d  %10f  %10f  %10f  %10f  %10f  %10f\n", 
    n + 1, b, poly5 ( b ), 0.0, fb, 0.0, fb );

  free ( f1 );
  free ( f2 );
  free ( s );
  free ( x );

  return;
}
/******************************************************************************/

void sine_transform_test03 ( )

/******************************************************************************/
/*
  Purpose:

    SINE_TRANSFORM_TEST03 evaluates the sine transform interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 9;
  int n2 = 1 + 2 * ( n + 1 );
  double *s;
  double *x;
  double *x2;

  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST03:\n" );
  printf ( "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n" );
  printf ( "  defined by a function F(X) evaluated at N equally spaced\n" );
  printf ( "  points in an interval [A,B].\n" );
  printf ( "  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.\n" );
  printf ( "\n" );
  printf ( "  The interpolant will be 0 at the 0th and (N+1)-th points.\n" );
  printf ( "  It equals the function at points 1 through N.\n" );
  printf ( "  In between, it can approximate smooth functions,\n" );
  printf ( "  and the approximation improves with N.\n" );
/*
  N determines the number of data points, indexed by 1 to N.  
  However, we essentially have N+2 data points, indexed 0 to N+1,
  with the data value being 0 at the first and last auxilliary points.
*/
  a = 1.0;
  b = 4.0;
/*
  Evenly spaced points between A and B, but omitting
  A and B themselves.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }
/*
  Determine the interpolant coefficients.
*/
  s = sine_transform_function ( n, a, b, poly5 );

  printf ( "\n" );
  printf ( "     I      X(I)      F(X(I))        S(I)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f  %10f  %10f\n", i, x[i], poly5 ( x[i] ), s[i] );
  }
/*
  Evaluate the interpolant.
*/
  fa = poly5 ( a );
  fb = poly5 ( b );
/*
  Evenly spaced points between A and B, including A and B,
  and twice the density of the previous set of points.
*/
  x2 = ( double * ) malloc ( n2 * sizeof ( double ) );

  for ( i = 0; i < n2; i++ )
  {
    x2[i] = ( ( double ) ( n2 - i - 1 ) * a   
            + ( double ) (      i     ) * b ) 
            / ( double ) ( n2     - 1 );
  }

  f2 = sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2 );

  printf ( "\n" );
  printf ( "     I      X            F(X)        FHAT(X)\n" );
  printf ( "\n" );

  for ( i = 0; i < n2; i++ )
  {
    printf ( "  %4d  %10f  %10f  %10f\n", i, x2[i], poly5 ( x2[i] ), f2[i] );
  }

  free ( f2 );
  free ( s );
  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

void sine_transform_test04 ( )

/******************************************************************************/
/*
  Purpose:

    SINE_TRANSFORM_TEST04 evaluates the sine transform interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 15;
  int n2 = 1 + 5 * ( n + 1 );
  double *s;
  double *x;
  double *x2;

  printf ( "\n" );
  printf ( "SINE_TRANSFORM_TEST04:\n" );
  printf ( "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n" );
  printf ( "  defined by a function F(X) evaluated at N equally spaced\n" );
  printf ( "  points in an interval [A,B].\n" );
  printf ( "  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.\n" );
  printf ( "\n" );
  printf ( "  The interpolant will be 0 at the 0th and (N+1)-th points.\n" );
  printf ( "  It equals the function at points 1 through N.\n" );
  printf ( "  In between, it can approximate smooth functions,\n" );
  printf ( "  and the approximation improves with N.\n" );
/*
  N determines the number of data points, indexed by 1 to N.  
  However, we essentially have N+2 data points, indexed 0 to N+1,
  with the data value being 0 at the first and last auxilliary points.
*/
  a = 0.0;
  b = 7.0;
/*
  Evenly spaced points between A and B, but omitting
  A and B themselves.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }
/*
  Determine the interpolant coefficients.
*/
  s = sine_transform_function ( n, a, b, cosine_sum );
/*
  Evaluate the interpolant.
*/
  fa = cosine_sum ( a );
  fb = cosine_sum ( b );
/*
  Evenly spaced points between A and B, including A and B,
  and twice the density of the previous set of points.
*/
  x2 = ( double * ) malloc ( n2 * sizeof ( double ) );

  for ( i = 0; i < n2; i++ )
  {
    x2[i] = ( ( double ) ( n2 - i - 1 ) * a   
            + ( double ) (      i     ) * b ) 
            / ( double ) ( n2     - 1 );
  }

  f2 = sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2 );

  printf ( "\n" );
  printf ( "  Expect exact agreement every 5th sample.\n" );
  printf ( "\n" );

  for ( i = 0; i < n2; i++ )
  {
    printf ( "  %4d  %10f  %10f  %10f\n", 
      i, x2[i], cosine_sum ( x2[i] ), f2[i] );
  }

  free ( f2 );
  free ( s );
  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

double cosine_sum ( double x )

/******************************************************************************/
/*
  Purpose:

    COSINE_SUM evaluates a function which is a cosine sum.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double COSINE_SUM, the value.
*/
{
  double value;

  value =   cos (       x ) 
    + 5.0 * cos ( 1.6 * x ) 
    - 2.0 * cos ( 2.0 * x ) 
    + 5.0 * cos ( 4.5 * x ) 
    + 7.0 * cos ( 9.0 * x );

  return value;
}
/******************************************************************************/

double poly5 ( double x )

/******************************************************************************/
/*
  Purpose:

    POLY5 evaluates a particular fifth-degree polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double POLY5, the value of the polynomial at X.
*/
{
  double value;

  value = ( x - 0.1 ) * 
          ( x - 0.2 ) * 
          ( x - 0.4 ) * 
          ( x - 2.1 ) * 
          ( x - 3.0 );

  return value;
}
