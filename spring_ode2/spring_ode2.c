# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPRING_ODE2.

  Discussion:

    This is a revision of the SPRING_ODE code.

    In this revision of the program, we want to use vectors (C arrays) to 
    store the data, and we want to write the data out to a file in a form 
    that Gnuplot (or other plotting programs) can use.

    Hooke's law for a spring observes that the restoring force is
    proportional to the displacement: F = - k x

    Newton's law relates the force to acceleration: F = m a

    Putting these together, we have

      m * d^2 x/dt^2 = - k * x

    We can add a damping force with coefficient c:

      m * d^2 x/dt^2 = - k * x - c * dx/dt

    If we write this as a pair of first order equations for (x,v), we have

          dx/dt = v
      m * dv/dt = - k * x - c * v

    and now we can approximate these values for small time steps.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 October 2013

  Author:

    John Burkardt

  Parameters:

    None
*/
{
  double c;
  char command_filename[] = "spring_ode2_commands.txt";
  FILE *command_unit;
  char data_filename[] = "spring_ode2_data.txt";
  FILE *data_unit;
  double dt;
  int i;
  double k;
  double m;
  int n = 101;
  double t[101];
  double t_final;
  double v[101];
  double x[101];

  timestamp ( );
  printf ( "\n" );
  printf ( "SPRING_ODE2\n" );
  printf ( "  C version\n" );
  printf ( "  Approximate the solution of a spring equation.\n" );
  printf ( "  Write data to a file for use by gnuplot.\n" );
/*
  Data
*/
  m = 1.0;
  k = 1.0;
  c = 0.3;
  t_final = 20.0;
  dt = t_final / ( double ) ( n - 1 );
/*
  Store the initial conditions in entry 0.
*/
  t[0] = 0.0;
  x[0] = 1.0;
  v[0] = 0.0;
/*
  Compute the approximate solution at equally spaced times 
  in entries 1 through N-1.
*/
  for ( i = 1; i < n; i++ )
  {
    t[i] = ( double ) ( i ) * t_final / ( double ) ( n - 1 );
    x[i] = x[i-1] + dt * v[i-1];
    v[i] = v[i-1] + ( dt / m ) * ( - k * x[i-1] - c * v[i-1] );
  }
/*
  Create the plot data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( data_unit, "  %14.6g  %14.6g  %14.6g\n", t[i], x[i], v[i] );
  }
  fclose ( data_unit );
  printf ( "  Created data file \"%s\".\n", data_filename );
/*
  Create the plot command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'xv_time.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- T --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- X(T), V(T) --->'\n" );
  fprintf ( command_unit, "set title 'Position and Velocity versus Time'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',", 
    data_filename );
  fprintf ( command_unit, "'' using 1:3 lw 3 linecolor rgb 'red'\n" );
  fprintf ( command_unit, "set output 'xv_phase.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- X(T) --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- V(T) --->'\n" );
  fprintf ( command_unit, "set title 'Position versus Velocity'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, 
    "plot '%s' using 2:3 lw 3 linecolor rgb 'green'\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPRING_ODE2:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
