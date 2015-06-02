# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "super_blas.h"

int main ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test08 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUPER_BLAS1_Z_PRB.

  Discussion:

    SUPER_BLAS1_Z_PRB tests the SUPER_BLAS library.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SUPER_BLAS1_Z_PRB:\n" );
  printf ( "  Z version\n" );
  printf ( "  Double precision complex arithmetic,\n" );
  printf ( "  Test the SUPER_BLAS library.\n" );
 
  test03 ( );
  test04 ( );
  test05 ( );
  test08 ( );

  test15 ( );
  test16 ( );
  test17 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SUPER_BLAS1_Z_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests ZAXPY.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  int inc1;
  int inc2;
  int ncopy;
  doublecomplex s;
  doublecomplex x[N] = {
    { 2.0, - 1.0 }, 
    {-4.0, - 2.0 }, 
    { 3.0, + 1.0 }, 
    { 2.0, + 2.0 }, 
    {-1.0, - 1.0 } };
  doublecomplex y[N] = {
    {-1.0, + 0.0 }, 
    { 0.0, - 3.0 },
    { 4.0, + 0.0 }, 
    {-3.0, + 4.0 }, 
    {-2.0, + 0.0 } };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  ZAXPY adds a multiple of one complex vector to another.\n" );
 
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }

  s.r = 0.50;
  s.i = - 1.00;

  printf ( "\n" );
  printf ( "  The scalar multiplier is: %f  %f\n", s.r, s.i );

  ncopy = N;
  inc1 = 1;
  inc2 = 1;

  zaxpy_ ( &ncopy, &s, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  A * X + Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }

  return;
# undef N
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests ZCOPY.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N1 5
# define N2 5
# define N 10

  doublecomplex a[N1*N2];
  int i;
  int inc1;
  int inc2;
  int j;
  int ncopy;
  doublecomplex x[N];
  doublecomplex y[N];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  ZCOPY copies one complex vector into another.\n" );
 
  for ( i = 0; i < N; i++ )
  {
    x[i].r = 10 * ( i + 1 );
    x[i].i =    + ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i].r = 20 * ( i + 1 );
    y[i].i =  2 * ( i + 1 );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1].r = 10 * ( i + 1 );
      a[i+j*N1].i =    + ( j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", a[i+j*N1].r, a[i+j*N1].i );
    }
    printf ( "\n" );
  }

  ncopy = 5;
  inc1 = 1;
  inc2 = 1;

  zcopy_ ( &ncopy, x, &inc1, y, &inc2 );
  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i].r = 20 * ( i + 1 );
    y[i].i =  2 * ( i + 1 );
  }

  ncopy = 3;
  inc1 = 2;
  inc2 = 3;

  zcopy_ ( &ncopy, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }

  ncopy = 5;
  inc1 = 1;
  inc2 = 1;

  zcopy_ ( &ncopy, x, &inc1, a, &inc2 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 1, A, 1 )\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", a[i+j*N1].r, a[i+j*N1].i );
    }
    printf ( "\n" );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1].r = 10 * ( i + 1 );
      a[i+j*N1].i =      ( j + 1 );
    }
  }

  ncopy = 5;
  inc1 = 2;
  inc2 = 5;

  zcopy_ ( &ncopy, x, &inc1, a, &inc2 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 2, A, 5 )\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", a[i+j*N1].r, a[i+j*N1].i );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests ZDOTC.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  int inc1;
  int inc2;
  int ncopy;
  doublecomplex x_norm;
  doublecomplex xy_dot;
  doublecomplex x[N] = {
    {  2.0, - 1.0 }, 
    { -4.0, - 2.0 }, 
    {  3.0, + 1.0 }, 
    {  2.0, + 2.0 }, 
    { -1.0, - 1.0 } };
  doublecomplex y[N] = {
    { -1.0, + 0.0 }, 
    {  0.0, - 3.0 }, 
    {  4.0, + 0.0 }, 
    { -3.0, + 4.0 }, 
    { -2.0, + 0.0 } };

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  ZDOTC computes the conjugated dot product of\n" );
  printf ( "  two complex vectors.\n" );
 
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }

  ncopy = N;
  inc1 = 1;
  inc2 = 1;

  zdotc_ ( &x_norm, &ncopy, x, &inc1, x, &inc2 );

  printf ( "\n" );
  printf ( "  The square of the norm of X, computed as\n" );
  printf ( "  ZDOTC(X,X) = (%f,  %f)\n", x_norm.r, x_norm.i );

  ncopy = N;
  inc1 = 1;
  inc2 = 1;

  zdotc_ ( &xy_dot, &ncopy, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }
  printf ( "\n" );
  printf ( "  The dot product X.Y* is (%f,  %f)\n", xy_dot.r, xy_dot.i );

  return;
# undef N
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests ZSCAL.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 6

  doublecomplex da;
  int i;
  int inc;
  int ncopy;
  doublecomplex x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i].r = 10 * ( i + 1 );
    x[i].i =      ( i + 1 );
  }

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  ZSCAL multiplies a complex scalar times a vector.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }

  da.r = 5.0;
  da.i = 0.0;

  ncopy = N;
  inc = 1;

  zscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  ZSCAL ( N, (%f, %f), X, 1 )\n", da.r, da.i );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i].r = 10 * ( i + 1 );
    x[i].i =      ( i + 1 );
  }

  ncopy = 3;
  da.r = -2.0;
  da.i = + 1.0;
  inc = 2;

  zscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  ZSCAL ( 3, (%f, %f), X, 2 )\n", da.r, da.i );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, x[i].r, x[i].i );
  }

  return;
# undef N
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests IZAMAX.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 5

  double cabs1;
  int i;
  int incx;
  int ncopy;
  doublecomplex x[N] = {
    {  2.0, - 1.0 }, 
    { -4.0, - 2.0 }, 
    {  3.0, + 1.0 }, 
    {  2.0, + 2.0 }, 
    { -1.0, - 1.0 } };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  IZAMAX returns the index of the entry of\n" );
  printf ( "  maximum magnitude in a complex vector.\n" );
 
  printf ( "\n" );
  printf ( "  The entries and CABS1 magnitudes:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    cabs1 = fabs ( x[i].r ) + fabs ( x[i].i );

    printf ( "  %6d  %6f  %6f    %6f\n", i, x[i].r, x[i].i, cabs1 );
  }

  ncopy = N;
  incx = 1;

  i = izamax_ ( &ncopy, x, &incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i );
  printf ( "\n" );
  printf ( "  Note that this is a 1-based index.\n" );
  printf ( "  Note that the L1 norm is used.\n" );

  return;
# undef N
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests DZASUM.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define MA 5
# define NA 4
# define NX 8

  doublecomplex a[MA*NA] = {
    { -3.0, + 4.0 }, 
    {  2.0, + 0.0 }, 
    {  3.0, - 4.0 }, 
    {  2.0, + 0.0 }, 
    {  2.0, - 1.0 }, 
    { -1.0, + 1.0 }, 
    {  0.0, + 5.0 }, 
    { -4.0, - 2.0 }, 
    { -4.0, + 1.0 }, 
    { -4.0, - 3.0 }, 
    {  0.0, - 2.0 }, 
    {  1.0, + 3.0 }, 
    { -3.0, + 3.0 }, 
    { -3.0, + 3.0 }, 
    { -1.0, - 2.0 }, 
    { -1.0, + 2.0 }, 
    {  2.0, - 4.0 }, 
    {  0.0, - 1.0 }, 
    {  0.0, - 1.0 }, 
    { -2.0, + 4.0 } };
  int i;
  int inc;
  int j;
  int ncopy;
  doublecomplex x[NX] = {
   {  2.0, - 1.0 }, 
   { -4.0, - 2.0 }, 
   {  3.0, + 1.0 }, 
   {  2.0, + 2.0 }, 
   { -1.0, - 1.0 }, 
   { -1.0, + 0.0 }, 
   {  0.0, - 3.0 }, 
   {  4.0, + 0.0 } };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  DZASUM adds the absolute values of elements\n" );
  printf ( "  of a complex vector.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < NX; i++ )
  {
    printf ( "  %6d  (%6.1f  %6.1f )\n", i, x[i].r, x[i].i );
  }

  ncopy = NX;
  inc = 1;
  printf ( "\n" );
  printf ( "  DZASUM ( NX,   X, 1    ) = %f\n",
    dzasum_ ( &ncopy,   x, &inc ) );

  ncopy = NX/2;
  inc = 2;
  printf ( "  DZASUM ( NX/2, X, 2    ) = %f\n",
    dzasum_ ( &ncopy, x, &inc ) );

  ncopy = 2;
  inc = NX/2;
  printf ( "  DZASUM ( 2,    X, NX/2 ) = %f\n",
    dzasum_ ( &ncopy,    x, &inc ) );

  printf ( "\n" );
  printf ( "  Demonstrate with a matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      printf ( "  (%6.1f  %6.1f )", 
        a[i+j*MA].r, a[i+j*MA].i );
    }
    printf ( "\n" );
  }

  ncopy = MA;
  inc = 1;
  printf ( "\n" );
  printf ( "  DZASUM ( MA, A[1,2], 1 )   = %f\n",
    dzasum_ ( &ncopy, a+0+1*MA, &inc ) );

  ncopy = NA;
  inc = MA;
  printf ( "  DZASUM ( NA, A[2,1], MA ) = %f\n",
    dzasum_ ( &ncopy, a+1+0*MA, &inc ) );

  return;
# undef MA
# undef NA
# undef NX
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests DZNRM2.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  int incx;
  int ncopy;
  double norm;
  doublecomplex x[N] = {
   { 2.0, - 1.0 }, 
   {-4.0, - 2.0 }, 
   { 3.0, + 1.0 }, 
   { 2.0, + 2.0 }, 
   {-1.0, - 1.0 } };

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  DZNRM2 returns the Euclidean norm of a complex vector.\n" );
 
  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", 
    i, x[i].r, x[i].i );
  }

  ncopy = N;
  incx = 1;

  norm = dznrm2_ ( &ncopy, x, &incx );

  printf ( "\n" ); 
  printf ( "  The L2 norm of X is %f\n", norm );

  return;
# undef N
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
