# include <stdlib.h>
# include <stdio.h>

# include "asa314.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA314_PRB.

  Discussion:

    ASA314_PRB tests the ASA314 library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2013

  Author:

    John Burkardt

  Reference:

    Roger Payne,
    Inversion of matrices with contents subject to modulo arithmetic,
    Applied Statistics,
    Volume 46, Number 2, 1997, pages 295-298.
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA314_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA314 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA314_PRB:\n" );
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

    TEST01 tests INVMOD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2013

  Author:

    John Burkardt

  Reference:

    Roger Payne,
    Inversion of matrices with contents subject to modulo arithmetic,
    Applied Statistics,
    Volume 46, Number 2, 1997, pages 295-298.
*/
{
  int cmod[3];
  int i;
  int ifault;
  int imat[3*3];
  int jmat[3*3] = { 1, 0, 0, 2, 1, 0, 1, 0, 1 };
  int mat[3*3] = { 1, 0, 0, 1, 1, 0, 2, 0, 1 };
  int nrow = 3;
  int rmod[3];

  for ( i = 0; i < nrow; i++ )
  {
    cmod[i] = 3;
  }
  for ( i = 0; i < nrow; i++ )
  {
    rmod[i] = 3;
  }

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  INVMOD computes the inverse of a matrix\n" );
  printf ( "  whose elements are subject to modulo arithmetic.\n" );

  i4mat_print ( nrow, nrow, mat, "  The matrix to be inverted:" );

  invmod ( mat, imat, rmod, cmod, nrow, &ifault );

  i4mat_print ( nrow, nrow, imat, "  The computed inverse:" );

  i4mat_print ( nrow, nrow, jmat, "  The correct inverse:" );

  return;
}

