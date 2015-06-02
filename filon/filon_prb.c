# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "filon.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
double *zero_integrand ( int n, double x[] );
double *one_integrand ( int n, double x[] );
double *two_integrand ( int n, double x[] );
double *log_integrand ( int n, double x[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FILON_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FILON_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FILON library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FILON_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests FILON_TAB_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;

  a = 0.0;
  b = 2.0 * r8_pi;

  n = 11;
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  FILON_TAB_COS estimates the integral of.\n" );
  printf ( "  F(X) * COS ( T * X )\n" );
  printf ( "  Use integrands F(X)=1, X, X^2.\n" );
  printf ( "\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "       T                      Approximate             Exact\n" );

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      t = 1.0;
    }
    else if ( k == 2 )
    {
      t = 2.0;
    }
    else if ( k == 3 )
    {
      t = 10.0;
    }

    printf ( "\n" );

    for ( i = 1; i <= 3; i++ )
    {
      if ( i == 1 )
      {
        ftab = zero_integrand ( n, x );
      }
      else if ( i == 2 )
      {
        ftab = one_integrand ( n, x );
      }
      else if ( i == 3 )
      {
        ftab = two_integrand ( n, x );
      }

      result = filon_tab_cos ( n, ftab, a, b, t );

      if ( i == 1 )
      {
        exact = ( sin ( t * b ) - sin ( t * a ) ) / t;
      }
      else if ( i == 2 )
      {
        exact = ( ( cos ( t * b ) + t * b * sin ( t * b ) ) 
                - ( cos ( t * a ) + t * a * sin ( t * a ) ) ) / t / t;
      }
      else if ( i == 3 )
      {
        exact = ( ( 2.0 * t * b * cos ( t * b ) 
              + ( t * t * b * b - 2.0 ) * sin ( t * b ) ) 
                - ( 2.0 * t * a * cos ( t * a ) 
              + ( t * t * a * a - 2.0 ) * sin ( t * a ) ) ) / t / t / t;
      }
      printf ( "  %24.16g  %24.16g  %24.16g\n", t, result, exact );
    }
  }

  free ( ftab );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests FILON_TAB_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;
/*
  Example suggested by James Roedder.
*/
  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Integrate F(X) = log(1+X)*cos(T*X):\n" );
  printf ( "  Supply integrand as a table.\n" );
  printf ( "  T = 10, and N increases\n" );
  printf ( "\n" );
  printf ( "       N    Approximate             Exact                   Error\n" );
  printf ( "\n" );

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
/*
  Set the X values.
*/
    x = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * a   
             + ( double ) (     i     ) * b ) 
             / ( double ) ( n     - 1 );
    }

    ftab = log_integrand ( n, x );
 
    t = 10.0;

    result = filon_tab_cos ( n, ftab, a, b, t );

    exact = -0.008446594405;
    error = result - exact;

    printf ( "  %6d  %24.16g  %24.16g  %16.8g\n", n, result, exact, error );

    free ( ftab );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests FILON_FUN_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
/*
  Example suggested by James Roedder.
*/
  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Integrate F(X)=log(1+X)*cos(T*X):\n" );
  printf ( "  Supply integrand as a function.\n" );
  printf ( "  T = 10, and N increases\n" );
  printf ( "\n" );
  printf ( "       N    Approximate             Exact                   Error\n" );
  printf ( "\n" );

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
 
    t = 10.0;

    result = filon_fun_cos ( n, log_integrand, a, b, t );

    exact = -0.008446594405;
    error = result - exact;

    printf ( "  %6d  %24.16g  %24.16g  %16.8g\n", n, result, exact, error );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests FILON_TAB_SIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;

  a = 0.0;
  b = 2.0 * r8_pi;
  n = 11;
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  FILON_TAB_SIN estimates the integral of.\n" );
  printf ( "  F(X) * SIN ( T * X )\n" );
  printf ( "  Use integrands 1, X, X^2.\n" );
  printf ( "\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "       T                      Approximate             Exact\n" );
  printf ( "\n" );

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      t = 1.0;
    }
    else if ( k == 2 )
    {
      t = 2.0;
    }
    else if ( k == 3 )
    {
      t = 10.0;
    }

    printf ( "\n" );

    for ( i = 1; i <= 3; i++ )
    {
      if ( i == 1 )
      {
        ftab = zero_integrand ( n, x );
      }
      else if ( i == 2 )
      {
        ftab = one_integrand ( n, x );
      }
      else if ( i == 3 )
      {
        ftab = two_integrand ( n, x );
      }

      result = filon_tab_sin ( n, ftab, a, b, t );

      if ( i == 1 )
      {
        exact = ( - cos ( t * b ) + cos ( t * a ) ) / t;
      }
      else if ( i == 2 )
      {
        exact = ( ( sin ( t * b ) - t * b * cos ( t * b ) ) 
                - ( sin ( t * a ) - t * a * cos ( t * a ) ) ) / t / t;
      }
      else if ( i == 3 )
      {
        exact = ( ( 2.0 * t * b * sin ( t * b ) 
                + ( 2.0 - t * t * b * b ) * cos ( t * b ) ) 
                - ( 2.0 * t * a * sin ( t * a ) 
                + ( 2.0 - t * t * a * a ) * cos ( t * a ) ) ) / t / t / t;
      }
      printf ( "  %24.16g  %24.16g  %24.16g\n", t, result, exact );
    }
  }

  free ( ftab );
  free ( x );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests FILON_TAB_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  double *ftab;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
  double *x;
/*
  Example suggested by James Roedder.
*/
  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Integrate F(X)=log(1+X)*sin(T*X):\n" );
  printf ( "  Supply integrand as a table.\n" );
  printf ( "  T = 10, and N increases\n" );
  printf ( "\n" );
  printf ( "       N    Approximate             Exact                   Error\n" );
  printf ( "\n" );

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;
/*
  Set the X values.
*/
    x = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * a   
             + ( double ) (     i     ) * b ) 
             / ( double ) ( n     - 1 );
    }

    ftab = log_integrand ( n, x );

    t = 10.0;

    result = filon_tab_sin ( n, ftab, a, b, t );

    exact = -0.19762680771872;
    error = result - exact;

    printf ( "  %6d  %24.16g  %24.16g  %16.8g\n", n, result, exact, error );

    free ( ftab );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests FILON_FUN_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int k;
  int n;
  const double r8_pi = 3.141592653589793;
  double result;
  double t;
/*
  Example suggested by James Roedder.
*/
  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Integrate F(X)=log(1+X)*sin(T*X):\n" );
  printf ( "  Supply integrand as a function.\n" );
  printf ( "  T = 10, and N increases\n" );
  printf ( "\n" );
  printf ( "       N    Approximate             Exact                   Error\n" );
  printf ( "\n" );

  a = 0.0;
  b = 2.0 * r8_pi;

  for ( j = 1; j <= 6; j++ )
  {
    n = i4_power ( 2, j ) * 10 + 1;

    t = 10.0;

    result = filon_fun_sin ( n, log_integrand, a, b, t );

    exact = -0.19762680771872;
    error = result - exact;

    printf ( "  %6d  %24.16g  %24.16g  %16.8g\n", n, result, exact, error );
  }

  return;
}
/******************************************************************************/

double *zero_integrand ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    ZERO_INTEGRAND evaluates the integrand x^0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double X[N], the evaluation points.

    Output, double ZERO_INTEGRAND[N], the function values.
*/
{
  double *fx;
  int i;

  fx = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0;
  }
  return fx;
}
/******************************************************************************/

double *one_integrand ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    ONE_INTEGRAND evaluates the integrand X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double X[N], the evaluation points.

    Output, double ONE_INTEGRAND[N], the function values.
*/
{
  double *fx;
  int i;

  fx = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i];
  }
  return fx;
}
/******************************************************************************/

double *two_integrand ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    TWO_INTEGRAND evaluates the integrand X^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double X[N], the evaluation points.

    Output, double TWO_INTEGRAND[N], the function values.
*/
{
  double *fx;
  int i;

  fx = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] * x[i];
  }
  return fx;
}
/******************************************************************************/

double *log_integrand ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    LOG_INTEGRAND evaluates the logarithmic integrand.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double X[N], the evaluation points.

    Output, double LOG_INTEGRAND[N], the function values.
*/
{
  double *fx;
  int i;

  fx = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    fx[i] = log ( 1.0 + x[i] );
  }
  return fx;
}
