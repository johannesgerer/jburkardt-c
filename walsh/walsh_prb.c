# include <stdlib.h>
# include <stdio.h>

# include "walsh.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    WALSH_PRB calls the WALSH test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 March 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "WALSH_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WALSH library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WALSH_PRB\n" );
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

    TEST01 tests FWT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2011

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  FWT computes a fast Walsh transform.\n" );

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = ( double * ) malloc ( n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );
    fwt ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    fwt ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    printf ( "\n" );
    printf ( "     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d  %10f  %10f  %10f\n", i, x[i], y[i], z[i] );
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  }
    
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests WALSH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2011

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  WALSH computes a fast Walsh transform.\n" );

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
       seed = 123456789;
       w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = ( double * ) malloc ( n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );
    walsh ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    walsh ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    printf ( "\n" );
    printf ( "     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d  %10f  %10f  %10f\n", i, x[i], y[i], z[i] );
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  }
    
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests HAAR, HAARIN and HNORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2011

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  HAAR computes a Haar transform.\n" );
  printf ( "  HNORM normalizes the transformed data.\n" );
  printf ( "  HAARIN computes an inverse Haar transform.\n" );

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = ( double * ) malloc ( n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }

    x = r8vec_copy_new ( n, w );

    haar ( n, w );

    y = r8vec_copy_new ( n, w );

    hnorm ( n, w );

    z = r8vec_copy_new ( n, w );

    haarin ( n, w );

    printf ( "\n" );
    printf ( "     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d  %10f  %10f  %10f  %10f\n", i, x[i], y[i], z[i], w[i] );
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests FFWT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2011

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n = 16;
  int seed;
  double *w;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  FFWT computes a fast Walsh transform.\n" );

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      seed = 123456789;
      w = r8vec_uniform_01_new ( n, &seed );
    }
    else
    {
      w = ( double * ) malloc ( n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        w[i] = ( double ) ( i + 1 );
      }
    }
    x = r8vec_copy_new ( n, w );
    ffwt ( n, w );
    y = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      y[i] = y[i] / ( double ) ( n );
    }
    ffwt ( n, w );
    z = r8vec_copy_new ( n, w );
    for ( i = 0; i < n; i++ )
    {
      z[i] = z[i] / ( double ) ( n );
    }

    printf ( "\n" );
    printf ( "     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d  %10f  %10f  %10f\n", i, x[i], y[i], z[i] );
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  } 
  return;
}
