# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "sftpack.h"

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
void test10 ( void );
void test11 ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SFTPACK_PRB.

  Discussion:

    SFTPACK_PRB calls the SFTPACK test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

   22 June 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SFTPACK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SFTPACK library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SFTPACK_PRB\n" );
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

    TEST01 tests R8VEC_SCT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int i;
  int n = 256;
  int seed;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For slow cosine transforms,\n" );
  printf ( "  R8VEC_SCT does a forward or backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
/*
  Compute the coefficients.
*/
  d = r8vec_sct ( n, c );

  r8vec_print_part ( n, d, 10, "  The cosine coefficients:" );
/*
  Now compute inverse transform of coefficients.  Should get back the
  original data.
*/
  e = r8vec_sct ( n, d );

  for ( i = 0; i < n; i++ )
  {
    e[i] = e[i] / ( double ) ( 2 * n );
  }

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  free ( c );
  free ( d );
  free ( e );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8VEC_SFTB and R8VEC_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double ahi = 5.0;
  double alo = 0.0;
  double azero;
  double *b;
  int i;
  int n = 36;
  int seed;
  double *x;
  double *z;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For real slow Fourier transforms,\n" );
  printf ( "  R8VEC_SFTF computes the forward transform.\n" );
  printf ( "  R8VEC_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data values, N = %d\n", n );

  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  a = ( double * ) malloc ( ( n / 2 ) * sizeof ( double ) );
  b = ( double * ) malloc ( ( n / 2 ) * sizeof ( double ) );

  r8vec_sftf ( n, x, &azero, a, b );

  printf ( "\n" );
  printf ( "  A (cosine) coefficients:\n" );
  printf ( "\n" );

  printf ( "  %4d  %14f\n", 0, azero );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    printf ( "  %4d  %14f\n", i, a[i] );
  }

  printf ( "\n" );
  printf ( "  B (sine) coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    printf ( "  %4d  %14f\n", i, b[i] );
  }
/*
  Now try to retrieve the data from the coefficients.
*/
  z = r8vec_sftb ( n, azero, a, b );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  free ( a );
  free ( b );
  free ( x );
  free ( z );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R8VEC_SHT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int n = 17;
  int seed;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For real slow Hartley transforms,\n" );
  printf ( "  R8VEC_SHT does a forward or backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
/*
  Compute the coefficients.
*/
  d = r8vec_sht ( n, c );

  r8vec_print_part ( n, d, 10, "  The Hartley coefficients:" );
/*
  Now compute inverse transform of coefficients.  Should get back the
  original data.
*/
  e = r8vec_sht ( n, d );

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  free ( c );
  free ( d );
  free ( e );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests R8VEC_SQCTB and R8VEC_SQCTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double ahi = 5.0;
  double alo = 0.0;
  int n = 256;
  int seed;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For real slow quarter wave cosine transforms,\n" );
  printf ( "  R8VEC_SQCTF does a forward transform;\n" );
  printf ( "  R8VEC_SQCTB does a backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the coefficients.
*/
  y = r8vec_sqctf ( n, x );

  r8vec_print_part ( n, y, 10, "  The cosine coefficients:" );
/*
  Now compute inverse transform of coefficients.  Should get back the
  original data.
*/
  z = r8vec_sqctb ( n, y );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  free ( x );
  free ( y );
  free ( z );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests R8VEC_SQSTB and R8VEC_SQSTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double ahi = 5.0;
  double alo = 0.0;
  int n = 256;
  int seed;
  double *x;
  double *y;
  double *z;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For real slow quarter wave sine transforms,\n" );
  printf ( "  R8VEC_SQSTF does a forward transform;\n" );
  printf ( "  R8VEC_SQSTB does a backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the coefficients.
*/
  y = r8vec_sqstf ( n, x );

  r8vec_print_part ( n, y, 10, "  The sine coefficients:" );
/*
  Now compute inverse transform of coefficients.  Should get back the
  original data.
*/
  z = r8vec_sqstb ( n, y );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  free ( x );
  free ( y );
  free ( z );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests R8VEC_SST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2010

  Author:

    John Burkardt
*/
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int i;
  int n = 256;
  int seed;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For slow sine transforms,\n" );
  printf ( "  R8VEC_SST does a forward or backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data items is N = %d\n", n );
/*
  Set the data values.
*/
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
/*
  Compute the coefficients;
*/
  d = r8vec_sst ( n, c );

  r8vec_print_part ( n, d, 10, "  The sine coefficients:" );
/*
  Now compute inverse transform of coefficients.  Should get back the
  original data.
*/
  e = r8vec_sst ( n, d );

  for ( i = 0; i < n; i++ )
  {
    e[i] = e[i] / ( double ) ( 2 * ( n + 1 ) );
  }

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  free ( c );
  free ( d );
  free ( e );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests C4VEC_SFTB and C4VEC_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 36;
  int seed;
  complex *x;
  complex *x2;
  complex *y;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  For complex slow Fourier transforms,\n" );
  printf ( "  C4VEC_SFTF computes the forward transform.\n" );
  printf ( "  C4VEC_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data values, N = %d\n", n );

  seed = 123456789;

  x = c4vec_uniform_01_new ( n, &seed );

  c4vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  y = c4vec_sftf ( n, x );

  c4vec_print_part ( n, y, 10, "  The Fourier coefficients:" );
/*
  Now try to retrieve the data from the coefficients.
*/
  x2 = c4vec_sftb ( n, y );

  c4vec_print_part ( n, x2, 10, "  The retrieved data:" );

  free ( x );
  free ( x2 );
  free ( y );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests C8VEC_SFTB and C8VEC_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 36;
  int seed;
  doublecomplex *x;
  doublecomplex *x2;
  doublecomplex *y;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  For complex slow Fourier transforms,\n" );
  printf ( "  C4VEC_SFTF computes the forward transform.\n" );
  printf ( "  C4VEC_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data values, N = %d\n", n );

  seed = 123456789;

  x = c8vec_uniform_01_new ( n, &seed );

  c8vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  y = c8vec_sftf ( n, x );

  c8vec_print_part ( n, y, 10, "  The Fourier coefficients:" );
/*
  Now try to retrieve the data from the coefficients.
*/
  x2 = c8vec_sftb ( n, y );

  c8vec_print_part ( n, x2, 10, "  The retrieved data:" );

  free ( x );
  free ( x2 );
  free ( y );

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests R4VEC_SFTB and R4VEC_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2010

  Author:

    John Burkardt
*/
{
  float *a;
  float ahi = 5.0;
  float alo = 0.0;
  float azero;
  float *b;
  int i;
  int n = 36;
  int seed;
  float *x;
  float *z;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  For real slow Fourier transforms,\n" );
  printf ( "  R4VEC_SFTF computes the forward transform.\n" );
  printf ( "  R4VEC_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The number of data values, N = %d\n", n );

  seed = 123456789;

  x = r4vec_uniform_new ( n, alo, ahi, &seed );

  r4vec_print_part ( n, x, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  a = ( float * ) malloc ( ( n / 2 ) * sizeof ( float ) );
  b = ( float * ) malloc ( ( n / 2 ) * sizeof ( float ) );

  r4vec_sftf ( n, x, &azero, a, b );

  printf ( "\n" );
  printf ( "  A (cosine) coefficients:\n" );
  printf ( "\n" );

  printf ( "  %4d  %14f\n", 0, azero );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    printf ( "  %4d  %14f\n", i, a[i] );
  }

  printf ( "\n" );
  printf ( "  B (sine) coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    printf ( "  %4d  %14f\n", i, b[i] );
  }
/*
  Now try to retrieve the data from the coefficients.
*/
  z = r4vec_sftb ( n, azero, a, b );

  r4vec_print_part ( n, z, 10, "  The retrieved data:" );

  free ( a );
  free ( b );
  free ( x );
  free ( z );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests C4MAT_SFTB and C4MAT_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 June 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n1 = 10;
  int n2 = 4;
  int seed;
  complex *x;
  complex *x2;
  complex *y;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  For complex slow Fourier transforms,\n" );
  printf ( "  C4MAT_SFTF computes the forward transform.\n" );
  printf ( "  C4MAT_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The data has dimensions N1 = %d by N2 = %d\n", n1, n2 );

  seed = 123456789;

  x = c4mat_uniform_01_new ( n1, n2, &seed );

  c4mat_print_some ( n1, n2, x, 1, 1, 10, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  y = c4mat_sftf ( n1, n2, x );

  c4mat_print_some ( n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:" );
/*
  Now try to retrieve the data from the coefficients.
*/
  x2 = c4mat_sftb ( n1, n2, y );

  c4mat_print_some ( n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:" );

  free ( x );
  free ( x2 );
  free ( y );

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests C8MAT_SFTB and C8MAT_SFTF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 June 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n1 = 10;
  int n2 = 4;
  int seed;
  doublecomplex *x;
  doublecomplex *x2;
  doublecomplex *y;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  For complex slow Fourier transforms,\n" );
  printf ( "  C8MAT_SFTF computes the forward transform.\n" );
  printf ( "  C8MAT_SFTB computes the backward transform.\n" );
  printf ( "\n" );
  printf ( "  The data has dimensions N1 = %d by N2 = %d\n", n1, n2 );

  seed = 123456789;

  x = c8mat_uniform_01_new ( n1, n2, &seed );

  c8mat_print_some ( n1, n2, x, 1, 1, 10, 10, "  The original data:" );
/*
  Compute the slow Fourier transform of the data.
*/
  y = c8mat_sftf ( n1, n2, x );

  c8mat_print_some ( n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:" );
/*
  Now try to retrieve the data from the coefficients.
*/
  x2 = c8mat_sftb ( n1, n2, y );

  c8mat_print_some ( n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:" );

  free ( x );
  free ( x2 );
  free ( y );

  return;
}
