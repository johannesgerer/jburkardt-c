# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "square_felippa_rule.h"

int main ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SQUARE_FELIPPA_RULE_PRB.

  Discussion:

    SQUARE_FELIPPA_RULE_PRB tests the SQUARE_FELIPPA_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2014

  Author:

    John Burkardt
*/
{
  int degree_max;

  timestamp ( );
  printf ( "\n" );
  printf ( "SQUARE_FELIPPA_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SQUARE_FELIPPA_RULE library.\n" );

  degree_max = 4;
  square_monomial_test ( degree_max );

  degree_max = 5;
  square_quad_test ( degree_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SQUARE_FELIPPA_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
