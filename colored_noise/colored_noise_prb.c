# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "colored_noise.h"

int main ( void );
void test01 ( int n, double q_d, double alpha, int seed );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for COLORED_NOISE_PRB.

  Discussion:

    COLORED_NOISE_PRB tests the COLORED_NOISE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 June 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  int i;
  int n;
  double q_d;
  int seed_init;

  timestamp ( );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "COLORED_NOISE_PRB\n" );
  fprintf ( stdout, "  C version\n" );
  fprintf ( stdout, "  Test the COLORED_NOISE library.\n" );

  n = 128;
  q_d = 1.0;
  alpha = 0.00;
  seed_init = 123456789;

  for ( i = 0; i <= 8; i++ )
  {
    alpha = 0.25 * i;
    test01 ( n, q_d, alpha, seed_init );
  }
/*
  Terminate.
*/
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "COLORED_NOISE_PRB:\n" );
  fprintf ( stdout, "  Normal end of execution.\n" );
  fprintf ( stdout, "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int n, double q_d, double alpha, int seed_init )

/******************************************************************************/
/*
  Purpose:

    TEST01 calls F_ALPHA with particular parameters.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of the sequence to generate.

    Input, double Q_D, the variance of the sequence.

    Input, double ALPHA, the exponent of the power law.

    Input, int SEED_INIT, the initial seed for the random number generator.
*/
{
  int i;
  char output_filename[15];
  FILE *output_unit;
  int seed;
  double *x;

  sprintf ( output_filename, "alpha_%4.2f.txt", alpha );
/*
  Report parameters.
*/
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST01:\n" );
  fprintf ( stdout, "  Generating %d sample points.\n", n );
  fprintf ( stdout, "  1/F^ALPHA noise has ALPHA = %f\n", alpha );
  fprintf ( stdout, "  Variance is %f\n", q_d );
  fprintf ( stdout, "  Initial random number seed = %d\n", seed_init );

  seed = seed_init;

  x = f_alpha ( n, q_d, alpha, &seed );
/*
  Print no more than 10 entries of the data.
*/
  r8vec_print_part ( n, x, 10, "  Noise sample:" );
/*
  Write the data to a file.
*/
  output_unit = fopen ( output_filename, "wr" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( output_unit, "%f\n", x[i] );
  }
  fclose ( output_unit );

  fprintf ( stdout, "  Data written to file \"%s\"\n", output_filename );

  free ( x );

  return;
}
