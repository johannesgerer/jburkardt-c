# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "super_blas.h"

int main ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test09 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUPER_BLAS1_C_PRB.

  Discussion:

    SUPER_BLAS1_C_PRB tests the SUPER_BLAS single precision complex routines.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SUPER_BLAS1_C_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Single precision complex arithmetic,\n" );
  printf ( "  Tests for the routines in the SUPER_BLAS library,\n" );
  printf ( "  the Level 1 Basic Linear Algebra Subprograms.\n" );
 
  test03 ( );
  test04 ( );
  test05 ( );
  test09 ( );

  test15 ( );
  test16 ( );
  test17 ( );

  printf ( "\n" );
  printf ( "SUPER_BLAS1_C_PRB:\n" );
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

    TEST03 tests CAXPY.

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
  complex s;
  complex x[N] = {
    { 2.0, - 1.0 }, 
    {-4.0, - 2.0 }, 
    { 3.0, + 1.0 }, 
    { 2.0, + 2.0 }, 
    {-1.0, - 1.0 } };
  complex y[N] = {
    {-1.0, + 0.0 }, 
    { 0.0, - 3.0 },
    { 4.0, + 0.0 }, 
    {-3.0, + 4.0 }, 
    {-2.0, + 0.0 } };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  CAXPY adds a multiple of one complex vector to another.\n" );
 
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

  caxpy_ ( &ncopy, &s, x, &inc1, y, &inc2 );

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

    TEST04 tests CCOPY.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N1 5
# define N2 5
# define N 10

  complex a[N1*N2];
  int i;
  int inc1;
  int inc2;
  int j;
  int ncopy;
  complex x[N];
  complex y[N];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CCOPY copies one complex vector into another.\n" );
 
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

  ccopy_ ( &ncopy, x, &inc1, y, &inc2 );
  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 1, Y, 1 )\n" );
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

  ccopy_ ( &ncopy, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  CCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, y[i].r, y[i].i );
  }

  ncopy = 5;
  inc1 = 1;
  inc2 = 1;

  ccopy_ ( &ncopy, x, &inc1, a, &inc2 );

  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 1, A, 1 )\n" );
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

  ccopy_ ( &ncopy, x, &inc1, a, &inc2 );

  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 2, A, 5 )\n" );
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

    TEST05 tests CDOTC.

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
  complex x_norm;
  complex xy_dot;
  complex x[N] = {
    {  2.0, - 1.0 }, 
    { -4.0, - 2.0 }, 
    {  3.0, + 1.0 }, 
    {  2.0, + 2.0 }, 
    { -1.0, - 1.0 } };
  complex y[N] = {
    { -1.0, + 0.0 }, 
    {  0.0, - 3.0 }, 
    {  4.0, + 0.0 }, 
    { -3.0, + 4.0 }, 
    { -2.0, + 0.0 } };

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CDOTC computes the conjugated dot product of\n" );
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

  cdotc_ ( &x_norm, &ncopy, x, &inc1, x, &inc2 );

  printf ( "\n" );
  printf ( "  The square of the norm of X, computed as\n" );
  printf ( "  CDOTC(X,X) = (%f,  %f)\n", x_norm.r, x_norm.i );

  ncopy = N;
  inc1 = 1;
  inc2 = 1;

  cdotc_ ( &xy_dot, &ncopy, x, &inc1, y, &inc2 );

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

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests CSCAL.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 6

  complex da;
  int i;
  int inc;
  int ncopy;
  complex x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i].r = 10 * ( i + 1 );
    x[i].i =      ( i + 1 );
  }

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  CSCAL multiplies a complex scalar times a vector.\n" );

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

  cscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  CSCAL ( N, (%f, %f), X, 1 )\n", da.r, da.i );
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

  cscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  CSCAL ( 3, (%f, %f), X, 2 )\n", da.r, da.i );
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

    TEST15 tests ICAMAX.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define N 5

  float cabs1;
  int i;
  int incx;
  int ncopy;
  complex x[N] = {
    {  2.0, - 1.0 }, 
    { -4.0, - 2.0 }, 
    {  3.0, + 1.0 }, 
    {  2.0, + 2.0 }, 
    { -1.0, - 1.0 } };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  ICAMAX returns the index of the entry of\n" );
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

  i = icamax_ ( &ncopy, x, &incx );

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

    TEST16 tests SCASUM.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
# define MA 5
# define NA 4
# define NX 8

  complex a[MA*NA] = {
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
  complex x[NX] = {
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
  printf ( "  SCASUM adds the absolute values of elements\n" );
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
  printf ( "  SCASUM ( NX,   X, 1    ) = %f\n",
    scasum_ ( &ncopy,   x, &inc ) );

  ncopy = NX/2;
  inc = 2;
  printf ( "  SCASUM ( NX/2, X, 2    ) = %f\n",
    scasum_ ( &ncopy, x, &inc ) );

  ncopy = 2;
  inc = NX/2;
  printf ( "  SCASUM ( 2,    X, NX/2 ) = %f\n",
    scasum_ ( &ncopy,    x, &inc ) );

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
  printf ( "  SCASUM ( MA, A[1,2], 1 )   = %f\n",
    scasum_ ( &ncopy, a+0+1*MA, &inc ) );

  ncopy = NA;
  inc = MA;
  printf ( "  SCASUM ( NA, A[2,1], MA ) = %f\n",
    scasum_ ( &ncopy, a+1+0*MA, &inc ) );

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

    TEST17 tests SCNRM2.

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
  float norm;
  complex x[N] = {
   { 2.0, - 1.0 }, 
   {-4.0, - 2.0 }, 
   { 3.0, + 1.0 }, 
   { 2.0, + 2.0 }, 
   {-1.0, - 1.0 } };

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  SCNRM2 returns the Euclidean norm of a complex vector.\n" );
 
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

  norm = scnrm2_ ( &ncopy, x, &incx );

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
