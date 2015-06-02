# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] );
void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_write ( char *output_filename, int n, double x[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SHALLOW_WATER_1D.

  Discussion:

    SHALLOW_WATER_1D approximates the 1D shallow water equations.

    This code can be considered a 1D version of Cleve Moler's shallow
    water equation solver.

    The version of the shallow water equations being solved here is in
    conservative form, and omits the Coriolis force.  The state variables
    are H (the height) and UH (the mass velocity).

    The equations have the form

      dH/dt + d UH/dx = 0

      d UH/dt + d ( U^2 H + 1/2 g H^2 )/dx = 0

    Here U is the ordinary velocity, U = UH/H, and g is the gravitational
    acceleration.

    The initial conditions are used to specify ( H, UH ) at an equally
    spaced set of points, and then the Lax-Wendroff method is used to advance
    the solution through a number of equally spaced points in time, with 
    boundary conditions supplying the first and last spatial values.


    Some input values will result in an unstable calculation that
    quickly blows up.  This is related to the Courant-Friedrichs-Lewy
    condition, which requires that DT be small enough, relative to DX and
    the velocity, that information cannot cross an entire cell.

    A "reasonable" set of input quantities is

      shallow_water_1d 41 100 1.0 0.2 9.8

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt

  Reference:

    Cleve Moler,
    "The Shallow Water Equations",
    Experiments with MATLAB.

  Parameters:

    Input, integer NX, the number of spatial nodes.

    Input, integer NT, the number of times steps.

    Input, real X_LENGTH, the length of the region.

    Input, real T_LENGTH, the time extent.

    Input, real G, the gravity constant.  G = 9.8 meters per second**2.

    Output, real H_ARRAY[NX*(NT+1)], the height for all space and time points.

    Output, real UH_ARRAY[NX*(NT+1], the mass velocity for all space and time points.

    Output, real X[NX], the X coordinates.

    Output, real T[NT+1], the T coordinates.
*/
{
  double dx;
  double dt;
  char filename_h[255];
  char filename_t[255];
  char filename_uh[255];
  char filename_x[255];
  double g;
  double *h;
  double *h_array;
  double *hm;
  int i;
  int it;
  int nt;
  int nx;
  double *t;
  double t_length;
  double *uh;
  double *uh_array;
  double *uhm;
  double *x;
  double x_length;

  timestamp ( );
  printf ( "\n" );
  printf ( "SHALLOW_WATER_1D\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
/*
  Get the quadrature file root name:
*/
  if ( argc <= 1 )
  {
    nx = 41;
  }
  else
  {
    nx = atoi ( argv[1] );
  }
  printf ( "  NX = %d\n", nx );

  if ( argc <= 2 )
  {
    nt = 100;
  }
  else
  {
    nt = atoi ( argv[2] );
  }
  printf ( "  NT = %d\n", nt );

  if ( argc <= 3 )
  {
    x_length = 1.0;
  }
  else
  {
    x_length = atof ( argv[3] );
  }
  printf ( "  X_LENGTH = %g\n", x_length );

  if ( argc <= 4 )
  {
    t_length = 0.2;
  }
  else
  {
    t_length = atof ( argv[4] );
  }
  printf ( "  T_LENGTH = %g\n", t_length );

  if ( argc <= 5 )
  {
    g = 9.8;
  }
  else
  {
    g = atof ( argv[5] );
  }
  printf ( "  G = %g\n", g );
/*
  Allocate space.
*/
  h = ( double * ) malloc ( nx * sizeof ( double ) );
  h_array = ( double * ) malloc ( nx * ( nt + 1 ) * sizeof ( double ) );
  hm = ( double * ) malloc ( ( nx - 1 ) * sizeof ( double ) );
  t = ( double * ) malloc ( ( nt + 1 ) * sizeof ( double ) );
  uh = ( double * ) malloc ( nx * sizeof ( double ) );
  uh_array = ( double * ) malloc ( nx * ( nt + 1 ) * sizeof ( double ) );
  uhm = ( double * ) malloc ( ( nx - 1 ) * sizeof ( double ) );
  x = ( double * ) malloc ( nx * sizeof ( double ) );
/*
  Define the locations of the nodes and time steps and the spacing.
*/
  x = r8vec_linspace_new ( nx, 0.0, x_length );
  t = r8vec_linspace_new ( nt + 1, 0.0, t_length );

  dx = x_length / ( double ) ( nx - 1 );
  dt = t_length / ( double ) ( nt );
/*
  Apply the initial conditions.
*/
  initial_conditions ( nx, nt, x, t[0], h, uh );
/*
  Apply the boundary conditions.
*/
  boundary_conditions ( nx, nt, x, t[0], h, uh );
/*
  Store the first time step into H_ARRAY and UH_ARRAY.
*/
  for ( i = 0; i < nx; i++ )
  {
    h_array[i+0*nx] = h[i];
    uh_array[i+0*nx] = uh[i];
  }
/*
  Take NT more time steps.
*/
  for ( it = 1; it <= nt; it++ )
  {
/*
  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
*/
    for ( i = 0; i < nx - 1; i++ )
    {
      hm[i] = ( h[i] + h[i+1] ) / 2.0 
        - ( dt / 2.0 ) * ( uh[i+1] - uh[i] ) / dx;
    }
    for ( i = 0; i < nx - 1; i++ )
    {
      uhm[i] = ( uh[i] + uh[i+1] ) / 2.0 
        - ( dt / 2.0 ) * ( 
          uh[i+1] * uh[i+1] / h[i+1] + 0.5 * g * h[i+1] * h[i+1]
        - uh[i] * uh[i]  / h[i] - 0.5 * g * h[i] * h[i] ) / dx;
    }
/*
  Take a full time step, evaluating the derivative at the half time step,
  to estimate the solution at the NX-2 nodes.
*/
    for ( i = 1; i < nx - 1; i++ )
    {
      h[i] = h[i] 
        - dt * ( uhm[i] - uhm[i-1] ) / dx;
    }
    for ( i = 1; i < nx; i++ )
    {
      uh[i] = uh[i] 
        - dt * ( 
          uhm[i] * uhm[i]  / hm[i] + 0.5 * g * hm[i] * hm[i]
        - uhm[i-1] * uhm[i-1]  / hm[i-1] - 0.5 * g * hm[i-1] * hm[i-1] ) / dx;
    }
/*
  Update the boundary conditions.
*/
    boundary_conditions ( nx, nt, x, t[it], h, uh );
/*
  Copy data into the big arrays.
*/
    for ( i = 0; i < nx; i++ )
    {
      h_array[i+it*nx] = h[i];
      uh_array[i+it*nx] = uh[i];
    } 
  }
/*
  Write data to files.
*/
  strcpy ( filename_x, "sw1d_x.txt" );
  strcpy ( filename_t, "sw1d_t.txt" );
  strcpy ( filename_h, "sw1d_h.txt" );
  strcpy ( filename_uh, "sw1d_uh.txt" );

  r8vec_write ( filename_x, nx, x );
  r8vec_write ( filename_t, nt + 1, t );
  r8mat_write ( filename_h, nx, nt + 1, h_array );
  r8mat_write ( filename_uh, nx, nt + 1, uh_array );

  printf ( "\n" );
  printf ( "  X  values saved in file \"%s\".\n", filename_x );
  printf ( "  T  values saved in file \"%s\".\n", filename_t );
  printf ( "  H  values saved in file \"%s\".\n", filename_h );
  printf ( "  UH values saved in file \"%s\".\n", filename_uh );
/*
  Free memory.
*/
  free ( h );
  free ( h_array );
  free ( hm );
  free ( t );
  free ( uh );
  free ( uh_array );
  free ( uhm );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SHALLOW_WATER_1D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void boundary_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )

/******************************************************************************/
/*
  Purpose:

    INITIAL_CONDITIONS sets the initial conditions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NX, the number of spatial nodes.

    Input, int NT, the number of times steps.

    Input, double X[NX], the coordinates of the nodes.

    Input, double T, the current time.

    Input/output, double H[NX], the height, with H(1) and H(NX) 
    adjusted for boundary conditions.

    Input/output, double UH[NX], the mass velocity, with UH(1) 
    and UH(NX) adjusted for boundary conditions.
*/
{
  int bc;

  bc = 1;
/*
  Periodic boundary conditions on H and UH.
*/
  if ( bc == 1 )
  {
    h[0]     = h[nx-2];
    h[nx-1]  = h[1];
    uh[0]    = uh[nx-2];
    uh[nx-1] = uh[1];
  }
/*
  Free boundary conditions on H and UH.
*/
  else if ( bc == 2 )
  {
    h[0]     = h[1];
    h[nx-1]  = h[nx-2];
    uh[0]    = uh[1];
    uh[nx-1] = uh[nx-2];
  }
/*
  Reflective boundary conditions on UH, free boundary conditions on H.
*/
  else if ( bc == 3 )
  {
    h[0]     =   h[1];
    h[nx-1]  =   h[nx-2];
    uh[0]    = - uh[1];
    uh[nx-1] = - uh[nx-2];
  }
  return;
}
/******************************************************************************/

void initial_conditions ( int nx, int nt, double x[], double t, double h[], 
  double uh[] )

/******************************************************************************/
/*
  Purpose:

    INITIAL_CONDITIONS sets the initial conditions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NX, the number of spatial nodes.

    Input, int NT, the number of times steps.

    Input, double X[NX], the coordinates of the nodes.

    Input, double T, the current time.

    Output, double H[NX], the initial height for all space.

    Output, double UH[NX], the initial mass velocity for all space.
*/
{
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < nx; i++ )
  {
    h[i] = 2.0 + sin ( 2.0 * pi * x[i] );
  }
  for ( i = 0; i < nx; i++ )
  {
    uh[i] = 0.0;
  }
  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

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
/******************************************************************************/

void r8vec_write ( char *output_filename, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_WRITE writes an R8VEC file.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int N, the number of points.

    Input, double X[N], the data.
*/
{
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    fprintf ( output, "  %24.16g\n", x[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
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
