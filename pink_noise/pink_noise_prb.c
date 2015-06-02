# include <stdlib.h>
# include <stdio.h>

# include "pink_noise.h"

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

    MAIN is the main program for PINK_NOISE_PRB.

  Discussion:

    PINK_NOISE_PRB tests the PINK_NOISE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "PINK_NOISE_PRB:\n" );
  fprintf ( stdout, "  C version\n" );
  fprintf ( stdout, "  Test the PINK_NOISE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "PINK_NOISE_PRB:\n" );
  fprintf ( stdout, "  Normal end of execution.\n" );
  fprintf ( stdout, "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests WRAP2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2010

  Author:

    John Burkardt
*/
{
  int i;
  int m;
  int q;
  int q_in;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST01\n" );
  fprintf ( stdout, "  WRAP2 performs a circular wrap.\n" );
  fprintf ( stdout, "  Q is expected to range between 0 and M.\n" );
  fprintf ( stdout, "  WRAP2 takes an input value of Q, and either\n" );
  fprintf ( stdout, "  increments it by M+1 until in the range, or\n" );
  fprintf ( stdout, "  decrements it by M+1 until in the range,\n" );
  fprintf ( stdout, "  and returns the result as the function value.\n" );

  for ( m = 2; m <= 4; m++ )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "   M  Qin  Qout\n" );
    fprintf ( stdout, "\n" );
    for ( i = -5; i < 3 * m; i++ )
    {
      q = i;
      q_in = q;
      wrap2 ( m, &q );
      fprintf ( stdout, "  %2d  %2d  %2d\n", m, q_in, q );
    }
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CDELAY2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2010

  Author:

    John Burkardt
*/
{
  int i;
  int m;
  int q;
  int q_in;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST02\n" );
  fprintf ( stdout, "  CDELAY2 is a circular buffer implementation\n" );
  fprintf ( stdout, "  of an M-fold delay.  Q is a counter\n" );
  fprintf ( stdout, "  which is decremented by CDELAY2, but reset to M\n" );
  fprintf ( stdout, "  after it reaches 0.\n" );

  for ( m = 2; m <= 4; m++ )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "   I   M  Qin  Qout\n" );
    fprintf ( stdout, "\n" );
    q = m;
    for ( i = 1; i <= 3 * ( m + 1 ); i++ )
    {
      q_in = q;
      cdelay2 ( m, &q );
      fprintf ( stdout, "  %2d  %2d  %2d  %2d\n", i, m, q_in, q );
    }
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests RANH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2010

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int q;
  double u;
  double y;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST03\n" );
  fprintf ( stdout, "  RANH is a random hold function.\n" );
  fprintf ( stdout, "  Given a value U and a delay D, it returns the value\n" );
  fprintf ( stdout, "  U for D calls, then resets U.\n" );

  for ( d = 5; 1 <= d; d-- )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "   I   D   Q      U           Y\n" );
    fprintf ( stdout, "\n" );
    u = 0.5;
    q = 3;
    for ( i = 1; i <= 20; i++ )
    {
      y = ranh ( d, &u, &q );
      fprintf ( stdout, "  %2d  %2d  %2d  %10f  %10f\n", i, d, q, u, y );
    }
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests RAN1F.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2010

  Author:

    John Burkardt
*/
{
  int b;
  int i;
  int *q;
  int rep;
  double *u;
  double y;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST04\n" );
  fprintf ( stdout, "  RAN1F generates random values with an approximate\n" );
  fprintf ( stdout, "  1/F distribution.\n" );

  for ( b = 1; b < 32; b = b * 2 )
  {
    u = ( double * ) malloc ( b * sizeof ( double ) );
    q = ( int * ) malloc ( b * sizeof ( int ) );
    for ( rep = 1; rep <= 4; rep++ )
    {
      for ( i = 0; i < b; i++ )
      {
        u[i] = ( double ) rand ( ) / ( double ) ( RAND_MAX ) - 0.5;
      }
      for ( i = 0; i < b; i++ )
      {
        q[i] = 0;
      }
      fprintf ( stdout, "\n" );
      fprintf ( stdout, "   B   I      Y\n" );
      fprintf ( stdout, "\n" );

      for ( i = 1; i <= 20; i++ )
      {
        y = ran1f ( b, u, q );
        fprintf ( stdout, "  %2d  %2d  %10f\n", b, i, y );
      }
    }
    free ( q );
    free ( u );
  }
  return;
}
