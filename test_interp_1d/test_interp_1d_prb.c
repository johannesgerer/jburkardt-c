# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "test_interp_1d.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_INTERP_1D_TEST tests the TEST_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int nd;

  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_INTERP_1D_TEST\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );

  test01 ( );

  nd = 11;
  test02 ( nd );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "TEST_INTERP_1D_TEST\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 simply prints the title of each function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;
  char *title;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Print the title of each function.\n" );

  prob_num = p00_prob_num ( );
  
  printf ( "\n" );
  printf ( "  There are %d functions available:\n", prob_num );
  printf ( "  Index  Title\n" );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );
    printf ( "  %2d  %s\n", prob, title );
    free ( title );
  }
  return;
}
//****************************************************************************80

void test02 ( int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_INTERP_1D_TEST02 evaluates each test function at ND sample points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of sample points.
//
{
  double a;
  double b;
  double *f;
  int j;
  int prob;
  int prob_num;
  double *x;

  printf ( "\n" );
  printf ( "TEST_INTERP_1D_TEST02\n" );
  printf ( "  Use P00_F to sample each function.\n" );

  prob_num = p00_prob_num ( );

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( nd, a, b );

  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    f = p00_f ( prob, nd, x );
    printf ( "\n" );
    printf ( "  X, F(X) for problem %d\n", prob );
    printf ( "\n" );
    for ( j = 0; j < nd; j++ )
    {
      printf ( "  %2d  %10f  %10g\n", j, x[j], f[j] );
    }
    free ( f );
  }

  free ( x );

  return;
}
