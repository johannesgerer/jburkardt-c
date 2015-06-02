# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( );
double *correlation_damped_sine ( int n, double rho[], double rho0 );
void correlation_plot ( int n, double rho[], double c[], char *header, 
  char *title );
double *r8vec_linspace_new ( int n, double a, double b );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    DAMPED_SINE evaluates and plots the damped sine correlation function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 January 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int n = 101;
  double *rho;
  double rho0;

  printf ( "\n" );
  printf ( "DAMPED_SINE\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrating how a correlation function can be\n" );
  printf ( "  evaluated and plotted using GNUPLOT.\n" );

  rho0 = 1.0;
  rho = r8vec_linspace_new ( n, -12.0, 12.0 );
  c = correlation_damped_sine ( n, rho, rho0 );
  correlation_plot ( n, rho, c, "damped_sine", "Damped sine correlation" );

  free ( rho );
  free ( c );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DAMPED_SINE\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

double *correlation_damped_sine ( int n, double rho[], double rho0 )

/******************************************************************************/
/*
  Purpose:

    CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 January 2013

  Author:

    John Burkardt

  Reference:

    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.

  Parameters:

    Input, int N, the number of arguments.

    Input, double RHO[N], the arguments.

    Input, double RHO0, the correlation length.

    Output, double CORRELATION_DAMPED_SINE[N], the correlations.
*/
{
  double *c;
  int i;
  double rhohat;

  c = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    if ( rho[i] == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      rhohat = fabs ( rho[i] ) / rho0;
      c[i] = sin ( rhohat ) / rhohat;
    }
  }

  return c;
}
/******************************************************************************/

void correlation_plot ( int n, double rho[], double c[], char *header, 
  char *title )

/******************************************************************************/
/*
  Purpose:

    CORRELATION_PLOT makes a plot of a correlation function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of arguments.

    Input, double RHO[N], the arguments.

    Input, double C[N], the correlations.

    Input, char *HEADER, an identifier for the files.

    Input, char *TITLE, a title for the plot.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  int i;
  double rho0;

  strcpy ( data_filename, header );
  strcat ( data_filename, "_data.txt" );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( data_unit, "%14.6g  %14.6g\n", rho[i], c[i] );
  }
  fclose ( data_unit );
  printf ( "  Created data file \"%s\".\n", data_filename );

  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output \"%s.png\"\n", header );
  fprintf ( command_unit, "set xlabel 'Distance Rho'\n" );
  fprintf ( command_unit, "set ylabel 'Correlation C(Rho)'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot \"%s\" using 1:2 lw 3 linecolor rgb \"blue\"\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file \"%s\"\n", command_filename );

  return;
}
/******************************************************************************/

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return x;
}
