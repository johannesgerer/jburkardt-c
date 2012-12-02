# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "test_interp_2d.h"
# include "r8lib.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN tests the TEST_INTERP_2D library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_INTERP_2D_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_INTERP_2D library.\n" );
  printf ( "  The R8LIB library is also required.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_INTERP_2D_PRB\n" );
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

    TEST01 simply prints the title of each grid and function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2012

  Author:

    John Burkardt
*/
{
  int f_num;
  int fi;
  char ft[100];
  int g_num;
  int gi;
  char gt[100];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For each grid and function, print the title.\n" );

  g_num = g00_num ( );

  printf ( "\n" );
  printf ( "  GRIDS:\n" );
  printf ( "  Index  Title\n" );
  printf ( "\n" );

  for ( gi = 1; gi <= g_num; gi++ )
  {
    g00_title ( gi, gt );
    printf ( "  %2d  %s\n", gi, gt );
  }

  f_num = f00_num ( );

  printf ( "\n" );
  printf ( "  FUNCTIONS:\n" );
  printf ( "  Index  Title\n" );
  printf ( "\n" );

  for ( fi = 1; fi <= f_num; fi++ )
  {
    f00_title ( fi, ft );
    printf ( "  %2d  %s\n", fi, ft );
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 samples each function using each grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2012

  Author:

    John Burkardt
*/
{
  double *f;
  double f_ave;
  double f_max;
  double f_min;
  int f_num;
  int fi;
  char ft[100];
  int g_num;
  int gi;
  int gn;
  char gt[100];
  double *gx;
  double *gy;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Sample each function over each grid.\n" );

  g_num = g00_num ( );
  f_num = f00_num ( );

  for ( fi = 1; fi <= f_num; fi++ )
  {
    f00_title ( fi, ft );
    printf ( "\n" );
    printf ( "  %2d  %s\n", fi, ft );
    printf ( "        Grid Title                     " );
    printf ( "Min(F)          Ave(F)           Max(F)\n" );
    printf ( "\n" );

    for ( gi = 1; gi <= g_num; gi++ )
    {
      g00_title ( gi, gt );
      gn = g00_size ( gi );

      gx = ( double * ) malloc ( gn * sizeof ( double ) );
      gy = ( double * ) malloc ( gn * sizeof ( double ) );

      g00_xy ( gi, gn, gx, gy );

      f = ( double * ) malloc ( gn * sizeof ( double ) );

      f00_f0 ( fi, gn, gx, gy, f );

      f_max = r8vec_max ( gn, f );
      f_min = r8vec_min ( gn, f );
      f_ave = r8vec_sum ( gn, f );
      f_ave = f_ave / ( double ) ( gn );

      printf ( "  %4d  %25s  %14g  %14g  %14g\n", gi, gt, f_min, f_ave, f_max );

      free ( f );
      free ( gx );
      free ( gy );
    }
  }
  return;
}
