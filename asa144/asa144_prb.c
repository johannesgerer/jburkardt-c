# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa144.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA144_PRB.

  Discussion:

    ASA144_PRB tests the ASA144 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA144_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA144 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA144_PRB\n" );
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

    TEST01 tests RCONT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 November 2010

  Author:

    John Burkardt
*/
{
# define NROW 5
# define NCOL 5

  int i;
  int ifault;
  int key;
  int matrix[NROW*NCOL];
  int ncolt[NCOL] = { 2, 2, 2, 2, 1 };
  int nrowt[NROW] = { 3, 2, 2, 1, 1 };
  int nsubt[NCOL];
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  RCONT constructs a random matrix with\n" );
  printf ( "  given row and column sums.\n" );

  i4vec_print ( NROW, nrowt, "  The rowsum vector:" );
  i4vec_print ( NCOL, ncolt, "  The columnsum vector: " );

  key = 0;

  for ( test = 1; test <= test_num; test++ )
  {
    rcont ( NROW, NCOL, nrowt, ncolt, nsubt, matrix, &key, &ifault );

    if ( ifault != 0 )
    {
      printf ( "\n" );
      printf ( "  RCONT returned IFAULT = %d\n", ifault );
      return;
    }
    i4mat_print ( NROW, NCOL, matrix, "  The rowcolsum matrix:" );
  }

  return;
# undef NROW
# undef NCOL
}
