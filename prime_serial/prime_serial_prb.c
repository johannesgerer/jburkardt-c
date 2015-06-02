# include <stdlib.h>
# include <stdio.h>

# include "prime_serial.h"

int main ( void );
void prime_number_sweep ( int n_lo, int n_hi, int n_factor );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PRIME_SERIAL_PRB.

  Discussion:

    PRIME_SERIAL_PRB tests the PRIME_SERIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 August 2009

  Author:

    John Burkardt
*/
{
  int n_factor;
  int n_hi;
  int n_lo;

  timestamp ( );
  printf ( "\n" );
  printf ( "PRIME_SERIAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PRIME_SERIAL library.\n" );

  n_lo = 1;
  n_hi = 131072;
  n_factor = 2;

  prime_number_sweep ( n_lo, n_hi, n_factor );

  n_lo = 5;
  n_hi = 500000;
  n_factor = 10;

  prime_number_sweep ( n_lo, n_hi, n_factor );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PRIME_SERIAL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void prime_number_sweep ( int n_lo, int n_hi, int n_factor )

/******************************************************************************/
/*
  Purpose:

   PRIME_NUMBER_SWEEP does repeated calls to PRIME_NUMBER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 August 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N_LO, the first value of N.

    Input, int N_HI, the last value of N.

    Input, int N_FACTOR, the factor by which to increase N after
    each iteration.

*/
{
  int i;
  int n;
  int primes;
  double ctime;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Call PRIME_NUMBER to count the primes from 1 to N.\n" );
  printf ( "\n" );
  printf ( "         N        Pi          Time\n" );
  printf ( "\n" );

  n = n_lo;

  while ( n <= n_hi )
  {
    ctime = cpu_time ( );

    primes = prime_number ( n );

    ctime = cpu_time ( ) - ctime;

    printf ( "  %8d  %8d  %14f\n", n, primes, ctime );
    n = n * n_factor;
  }
 
  return;
}


