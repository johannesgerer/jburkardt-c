# include <stdlib.h>
# include <stdio.h>

# include "subset_sum_serial.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUBSET_SUM_SERIAL_PRB.

  Discussion:

    SUBSET_SUM_SERIAL_PRB tests the SUBSET_SUM_SERIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SUBSET_SUM_SERIAL_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the SUBSET_SUM_SERIAL library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SUBSET_SUM_SERIAL_PRB\n" );
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

    SUBSET_SUM_SERIAL_TEST01 tests the SUBSET_SUM_SERIAL program.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2013

  Author:

    John Burkardt
*/
{
  int *choice;
  int i;
  int n = 21;
  int target;
  int weight[21] = {
    518533, 1037066, 2074132, 1648264, 796528, 
   1593056,  686112, 1372224,  244448, 488896, 
    977792, 1955584, 1411168,  322336, 644672, 
   1289344,   78688,  157376,  314752, 629504, 
   1259008 };
  int w_sum;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the SUBSET_SUM_SERIAL function, which looks for a selection\n" );
  printf ( "  from a set of weights that adds up to a given target.\n" );
/*
  Define the problem data.
*/
  target = 2463098;
  printf ( "\n" );
  printf ( "  Target value:\n" );
  printf ( "  %d\n", target );
 
  printf ( "\n" );
  printf ( "   I      W(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %8d\n", i, weight[i] );
  }

  choice = subset_sum_serial ( n, weight, target );

  if ( choice[0] == -1 )
  {
    printf ( "\n" );
    printf ( "  No solution was found.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "   I*     W*\n" );
    printf ( "\n" );
    w_sum = 0;
    for ( i = 0; i < n; i++ )
    {
      if ( choice[i] == 1 )
      {
        w_sum = w_sum + weight[i];
        printf ( "  %2d  %8d\n", i, weight[i] );
      }
    }
    printf ( "\n" );
    printf ( "  Sum:    %d\n", w_sum );
    printf ( "  Target: %d\n", target );
  }

  free ( choice );

  return;
}

