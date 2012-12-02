# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms446.h"

int main ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
double *functn ( double x );
double *functn_d ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS446_PRB.

  Discussion:

    TOMS446_PRB calls the TOMS446 tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TOMS446_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS446 library.\n" );

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
  printf ( "TOMS446_PRB\n" );
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

    TEST01 tests CHEBY, which computes Chebyshev series.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int nf = 5;
  int npl = 10;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test CHEBY, which computes the\n" );
  printf ( "  Chebyshev series for several functions.\n" );

  x = cheby ( nf, npl, functn );

  printf ( "\n" );
  printf ( "          Sin(x)          Cos(x)        Sin(2x)         Cos(2x)           X^5\n" );
  printf ( "\n" );

  for ( i = 0; i < npl; i++ )
  {
    for ( j = 0; j < nf; j++ )
    {
      printf ( "  %14f", x[i+j*npl] );
    }
    printf ( "\n" );
  }
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests MULTPLY, which multiplies two Chebyshev series.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int nf = 5;
  int npl = 10;
  double *x;
  double *x1;
  double *x2;
  double *x3;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test MLTPLY, which computes the\n" );
  printf ( "  product of two Chebyshev series.\n" );
  printf ( "\n" );
  printf ( "  Multiply series for SIN(X) and COS(X)\n" );
  printf ( "  and compare with series for 1/2*SIN(2X).\n" );

  x = cheby ( nf, npl, functn );

  x1 = ( double * ) malloc ( npl * sizeof ( double ) );
  x2 = ( double * ) malloc ( npl * sizeof ( double ) );

  for ( i = 0; i < npl; i++ )
  {
    x1[i] = x[i+0*npl];
    x2[i] = x[i+1*npl];
    x[i+2*npl] = 0.5 * x[i+2*npl];
  }
  x3 = mltply_new ( x1, x2, npl );

  printf ( "\n" );
  printf ( "          Sin(x)          Cos(x)       1/2*Sin(2x)         RESULT\n" );
  printf ( "\n" );

  for ( i = 0; i < npl; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i+0*npl], x[i+1*npl], x[i+2*npl], x3[i] );
  }
  free ( x );
  free ( x1 );
  free ( x2 );
  free ( x3 );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests ECHEB, which evaluates a Chebyshev series.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  double fval;
  double *fxj;
  int i;
  int j;
  int k;
  int nf = 5;
  int npl = 10;
  int nx;
  double *x;
  double *x2;
  double xval;

  nx = 6;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Test ECHEB, which evaluates a\n" );
  printf ( "  Chebyshev series.\n" );

  x = cheby ( nf, npl, functn );
  x2 = ( double * ) malloc ( npl * sizeof ( double ) );

  for ( j = 0; j < nf; j++ )
  {
    for ( i = 0; i < npl; i++ )
    {
      x2[i] = x[i+j*npl];
    }

    printf ( "\n" );
    if ( j == 0 )
    {
      printf ( "  Sin(x)\n" );
    }
    else if ( j == 1 )
    {
      printf ( "  Cos(x)\n" );
    }
    else if ( j == 2 )
    {
      printf ( "  Sin(2x)\n" );
    }
    else if ( j == 3 )
    {
      printf ( "  Cos(2x)\n" );
    }
    else if ( j == 4 )
    {
      printf ( "  x^5\n" );
    }

    printf ( "\n" );

    for ( k = 0; k < nx; k++ )
    {
      xval = 2.0 * ( double ) ( k ) / ( double ) ( nx - 1 ) - 1.0;

      fxj = functn ( xval );

      fval = echeb ( xval, x2, npl );

      printf ( "  %14g  %14g  %14g\n", xval, fxj[j], fval );
      free ( fxj );
    }
  }
  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests EDCHEB, which evaluates the derivative of a Chebyshev series.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  double fval;
  double *fxj;
  int i;
  int j;
  int k;
  int nf = 5;
  int npl = 10;
  int nx;
  double *x;
  double *x2;
  double xval;

  nx = 6;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test EDCHEB, which evaluates the\n" );
  printf ( "  derivative of a Chebyshev series.\n" );

  x = cheby ( nf, npl, functn );
  x2 = ( double * ) malloc ( npl * sizeof ( double ) );

  for ( j = 0; j < nf; j++ )
  {
    for ( i = 0; i < npl; i++ )
    {
      x2[i] = x[i+j*npl];
    }

    printf ( "\n" );
    if ( j == 0 )
    {
      printf ( "  Sin(x)\n" );
    }
    else if ( j == 1 )
    {
      printf ( "  Cos(x)\n" );
    }
    else if ( j == 2 )
    {
      printf ( "  Sin(2x)\n" );
    }
    else if ( j == 3 )
    {
      printf ( "  Cos(2x)\n" );
    }
    else if ( j == 4 )
    {
      printf ( "  x^5\n" );
    }

    printf ( "\n" );

    for ( k = 0; k < nx; k++ )
    {
      xval = 2.0 * ( double ) ( k ) / ( double ) ( nx - 1 ) - 1.0;

      fxj = functn_d ( xval );

      fval = edcheb ( xval, x2, npl );

      printf ( "  %14g  %14g  %14g\n", xval, fxj[j], fval );
      free ( fxj );
    }
  }
  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests DFRNT, which computes the Chebyshev series of a derivative.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int nf = 5;
  int npl = 10;
  double *x;
  double *x2;
  double *x3;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Test DFRNT, which computes the\n" );
  printf ( "  Chebyshev series for the derivative\n" );
  printf ( "  of several functions.\n" );

  x = cheby ( nf, npl, functn );
  x2 = ( double * ) malloc ( npl * sizeof ( double ) );

  for ( j = 0; j < nf; j++ )
  {
    for ( i = 0; i < npl; i++ )
    {
      x2[i] = x[i+j*npl];
    }
    x3 = dfrnt ( x2, npl );
    for ( i = 0; i < npl; i++ )
    {
      x[i+j*npl] = x3[i];
    }
    free ( x3 );
  }

  printf ( "\n" );
  printf ( "  Chebyshev series for d/dx of:\n" );
  printf ( "\n" );
  printf ( "        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5\n" );
  printf ( "\n" );

  for ( i = 0; i < npl; i++ )
  {
    for ( j = 0; j < nf; j++ )
    {
      printf ( "  %14g", x[i+j*npl] );
    }
    printf ( "\n" );
  }

  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests NTGRT, which computes the Chebyshev series of an indefinite integral.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int nf = 5;
  int npl = 10;
  double *x;
  double *x2;
  double *x3;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Test NTGRT, which computes the\n" );
  printf ( "  Chebyshev series for the indefinite\n" );
  printf ( "  integral of several functions.\n" );

  x = cheby ( nf, npl, functn );
  x2 = ( double * ) malloc ( npl * sizeof ( double ) );

  for ( j = 0; j < nf; j++ )
  {
    for ( i = 0; i < npl; i++ )
    {
      x2[i] = x[i+j*npl];
    }
    x3 = ntgrt ( x2, npl );
    for ( i = 0; i < npl; i++ )
    {
      x[i+j*npl] = x3[i];
    }
    free ( x3 );
  }

  printf ( "\n" );
  printf ( "  Chebyshev series for indefinite integral of:\n" );
  printf ( "\n" );
  printf ( "        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5\n" );
  printf ( "\n" );

  for ( i = 0; i < npl; i++ )
  {
    for ( j = 0; j < nf; j++ )
    {
      printf ( "  %14g", x[i+j*npl] );
    }
    printf ( "\n" );
  }
  free ( x );
  free ( x2 );

  return;
}
/******************************************************************************/

double *functn ( double x )

/******************************************************************************/
/*
  Purpose:

    FUNCTN evaluates several functions at X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double FXJ[5], the derivative values.
*/
{
  double *fxj;
  
  fxj = ( double * ) malloc ( 5 * sizeof ( double ) );

  fxj[0] = sin ( x );
  fxj[1] = cos ( x );
  fxj[2] = sin ( 2.0 * x );
  fxj[3] = cos ( 2.0 * x );
  fxj[4] = pow ( x, 5 );

  return fxj;
}
/******************************************************************************/

double *functn_d ( double x )

/******************************************************************************/
/*
  Purpose:

    FUNCTN_D evaluates the derivatives of several functions at X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double FXJ[5], the derivative values.
*/
{
  double *fxj;

  fxj = ( double * ) malloc ( 5 * sizeof ( double ) );

  fxj[0] =  cos ( x );
  fxj[1] = -sin ( x );
  fxj[2] =  2.0 * cos ( 2.0 * x );
  fxj[3] = -2.0 * sin ( 2.0 * x );
  fxj[4] =  5.0 * pow ( x, 4 );

  return fxj;
}
