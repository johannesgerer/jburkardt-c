# include <stdlib.h>
# include <stdio.h>

# include "test_interp.h"
# include "r8lib.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_INTERP_PRB.

  Discussion:

    TEST_INTERP_PRB calls the TEST_INTERP tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 February 2012

  Author:

    John Burkardt
*/
{
  timestamp (  );

  printf ( "\n" );
  printf ( "TEST_INTERP_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_INTERP library.\n" );
  printf ( "  This test also requires the R8LIB library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_INTERP_PRB\n" );
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

    TEST01 shows how P00_STORY can be called.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 February 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  P00_STORY prints the problem \"story\".\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem %d\n", prob );

    p00_story ( prob );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 prints the data for each problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 February 2012

  Author:

    John Burkardt
*/
{
  int data_num;
  int dim_num;
  double *p;
  int prob;
  int prob_num;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  P00_DATA_NUM returns N, the number of data points.\n" );
  printf ( "  P00_DIM_NUM returns M, the dimension of data.\n" );
  printf ( "  P00_DATA returns the actual (MxN) data.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem  %d\n", prob );

    data_num = p00_data_num ( prob );
    printf ( "  DATA_NUM %d\n", data_num );
    dim_num = p00_dim_num ( prob );
    printf ( "  DIM_NUM  %d\n", dim_num );

    p = p00_data ( prob, dim_num, data_num );

    r8mat_transpose_print ( dim_num, data_num, p, "  Data array:" );

    free ( p );
  }

  return;
}
