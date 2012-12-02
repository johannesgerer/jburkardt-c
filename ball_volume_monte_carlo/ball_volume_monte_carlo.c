# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
double r8_abs ( double x );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double r8vec_sum ( int n, double a[] );
double *sphere_indicator ( int dim_num, int point_num, double x[] );
double sphere_unit_volume_nd ( int dim_num );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BALL_VOLUME_MONTE_CARLO.

  Discussion:

    DIM_NUM = 6 is a reasonable test.

    N_LOG2_MAX = 25 puts a strain on the system, since we generate that
    many temporary points at once.  To solve bigger problems, it would
    be better to compute the new points in batches whose maximum size
    is limited.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2012

  Author:

    John Burkardt
*/
{
  int dim_num;
  double estimate;
  double error;
  double exact;
  double *fx;
  int i;
  int j;
  int n;
  int n_more;
  int n_log2;
  int n_log2_max = 25;
  double quad;
  double quad_more;
  int seed;
  double volume;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "BALL_VOLUME_MONTE_CARLO:\n" );
  printf ( "  C version\n" );
  printf ( "  Estimate the volume of the unit sphere using\n" );
  printf ( "  a Monte Carlo procedure.\n" );
/*
  Get the quadrature file root name:
*/
  if ( 1 < argc )
  {
    dim_num = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "BALL_VOLUME_MONTE_CARLO:\n" );
    printf ( "  Enter the spatial dimension of the sphere\n" );

    scanf ( "%d", &dim_num );
  }
/*
  Get the random number seed, if supplied.
*/
  if ( 2 < argc )
  {
    seed = atoi ( argv[2] );
  }
  else
  {
    seed = 123456789;
    printf ( "\n" );
    printf ( "BALL_VOLUME_MONTE_CARLO:\n" );
    printf ( "  Using default seed for random number generator.\n" );
  }
/*
  Report user input.
*/
  printf ( "\n" );
  printf ( "  The spatial dimension is  %d\n", dim_num );
  printf ( "  The random number seed is %d\n", seed );
/*
  Begin computation.
*/
  printf ( "\n" );
  printf ( "    Log(N)         N      Estimate         Error\n" );
  printf ( "\n" );

  quad = 0.0;
  volume = pow ( 2.0, dim_num );

  for ( n_log2 = 0; n_log2 <= n_log2_max; n_log2++ )
  {
    if ( n_log2 == 0 )
    {
      quad = 0.0;
      n_more = 1;
      n = 0;
    }
    else if ( n_log2 == 1 )
    {
      n_more = 1;
    }
    else
    {
      n_more = 2 * n_more;
    }

    x = r8mat_uniform_01_new ( dim_num, n_more, &seed );
/*
  Rescale X from [0,1] to [-1,1].
*/
    for ( j = 0; j < n_more; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        x[i+j*dim_num] = 2.0 * x[i+j*dim_num] - 1.0;
      }
    }
    fx = sphere_indicator ( dim_num, n_more, x );

    quad_more = r8vec_sum ( n_more, fx );

    free ( fx );
    free ( x );
/*
  Incorporate the new data into the totals.
*/
    n = n + n_more;
    quad = quad + quad_more;

    estimate = volume * quad / ( double ) ( n );
    exact = sphere_unit_volume_nd ( dim_num );
    error = r8_abs ( exact - estimate );
    printf ( "  %8d  %8d  %16.10g  %16.2g\n", n_log2, n, estimate, error );
  }

  printf ( "\n" );
  printf ( "        oo        oo  %16.10g  %10.2g\n", exact, 0.0 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BALL_VOLUME_MONTE_CARLO:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with unit pseudorandom values.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

double r8vec_sum ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SUM returns the sum of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_SUM, the sum of the vector.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
}
/******************************************************************************/

double *sphere_indicator ( int dim_num, int point_num, double x[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE_INDICATOR evaluates the unit sphere indicator function.

  Discussion:

    F(X) = 1 if X is on or inside the unit sphere, and 0 elsewhere.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int POINT_NUM, the number of points to evaluate.

    Input, double X[DIM_NUM*POINT_NUM], the points.

    Output, double SPHERE_INDICATOR[POINT_NUM], the unit sphere indicator 
    function value.
*/
{
  int i;
  int j;
  double t;
  double *value;

  value = ( double * ) malloc ( point_num * sizeof ( double ) );

  for ( j = 0; j < point_num; j++ )
  {
    t = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      t = t + x[i+j*dim_num] * x[i+j*dim_num];
    }

    if ( t <= 1.0 )
    {
      value[j] = 1.0;
    }
    else
    {
      value[j] = 0.0;
    }
  }
  return value;
}
/******************************************************************************/

double sphere_unit_volume_nd ( int dim_num )

/******************************************************************************/
/*
  Purpose:

    SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.

  Discussion:

    The unit sphere in ND satisfies the equation:

      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1

     DIM_NUM  Volume

     1    2
     2    1        * PI
     3  ( 4 /   3) * PI
     4  ( 1 /   2) * PI^2
     5  ( 8 /  15) * PI^2
     6  ( 1 /   6) * PI^3
     7  (16 / 105) * PI^3
     8  ( 1 /  24) * PI^4
     9  (32 / 945) * PI^4
    10  ( 1 / 120) * PI^5

    For the unit sphere, Volume(N) = 2 * PI * Volume(N-2)/ N

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the dimension of the space.

    Output, double SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
*/
{
  int i;
  int m;
  double pi = 3.141592653589793;
  double volume;

  if ( dim_num % 2== 0 )
  {
    m = dim_num / 2;
    volume = 1.0;
    for ( i = 1; i <= m; i++ )
    {
      volume = volume * pi / ( ( double ) i );
    }
  }
  else
  {
    m = ( dim_num - 1 ) / 2;
    volume = pow ( pi, m ) * pow ( 2.0, dim_num );
    for ( i = m + 1; i <= 2 * m + 1; i++ )
    {
      volume = volume / ( ( double ) i );
    }
  }

  return volume;
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
