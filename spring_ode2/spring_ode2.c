# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPRING_ODE2.
//
//  Discussion:
//
//    This is a revision of the SPRING_ODE code.
//
//    In this revision of the program, we want to use vectors (C arrays) to 
//    store the data, and we want to write the data out to a file in a form 
//    that Gnuplot (or other plotting programs) can use.
//
//    Hooke's law for a spring observes that the restoring force is
//    proportional to the displacement: F = - k x
//
//    Newton's law relates the force to acceleration: F = m a
//
//    Putting these together, we have
//
//      m * d^2 x/dt^2 = - k * x
//
//    We can add a damping force with coefficient c:
//
//      m * d^2 x/dt^2 = - k * x - c * dx/dt
//
//    If we write this as a pair of first order equations for (x,v), we have
//
//          dx/dt = v
//      m * dv/dt = - k * x - c * v
//
//    and now we can approximate these values for small time steps.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
  float c;
  float dt;
  int i;
  float k;
  float m;
  int n = 101;
  float t[101];
  float t_final;
  float v[101];
  float x[101];

  timestamp ( );
  fprintf ( stderr, "\n" );
  fprintf ( stderr, "SPRING_ODE2\n" );
  fprintf ( stderr, "  C version\n" );
  fprintf ( stderr, "  Approximate the solution of a spring equation.\n" );
  fprintf ( stderr, "  Write data to a file for use by gnuplot.\n" );
  fprintf ( stderr, "\n" );
//
//  Data
//
  m = 1.0;
  k = 1.0;
  c = 0.3;
  t_final = 20.0;
  dt = t_final / ( float ) ( n - 1 );
//
//  Store the initial conditions in entry 0.
//
  t[0] = 0.0;
  x[0] = 1.0;
  v[0] = 0.0;
//
//  Compute the approximate solution at equally spaced times in entries 1 through N-1.
//
  for ( i = 1; i < n; i++ )
  {
    t[i] = ( float ) ( i ) * t_final / ( float ) ( n - 1 );
    x[i] = x[i-1] + dt * v[i-1];
    v[i] = v[i-1] + ( dt / m ) * ( - k * x[i-1] - c * v[i-1] );
  }
//
//  Write the data to a file for plotting, possibly by Gnuplot.
//  Gnuplot expects T, X and V to be columns of output.
//
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "%f  %f  %f\n", t[i], x[i], v[i] );
  }
//
//  Terminate.
//
  fprintf ( stderr, "\n" );
  fprintf ( stderr, "SPRING_ODE2:\n" );
  fprintf ( stderr, "  Normal end of execution.\n" );
  fprintf ( stderr, "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stderr, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
