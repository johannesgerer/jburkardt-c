# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "simplex_monte_carlo.h"

/******************************************************************************/

double *monomial_value ( int m, int n, int e[], double x[] )

/******************************************************************************/
/*
  Purpose:

    MONOMIAL_VALUE evaluates a monomial.

  Discussion:

    This routine evaluates a monomial of the form

      product ( 1 <= i <= m ) x(i)^e(i)

    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points at which the
    monomial is to be evaluated.

    Input, int E[M], the exponents.

    Input, double X[M*N], the point coordinates.

    Output, double MONOMIAL_VALUE[N], the value of the monomial.
*/
{
  int i;
  int j;
  double *v;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
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

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

double simplex01_monomial_integral ( int m, int e[] )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX01_MONOMIAL_INTEGRAL: integrals in the unit simplex in M dimensions.

  Discussion:

    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int E[M], the exponents.  
    Each exponent must be nonnegative.

    Output, double SIMPLEX01_MONOMIAL_INTEGRAL, the integral.
*/
{
  int i;
  double integral;
  int j;
  int k;

  for ( i = 0; i < m; i++ )
  {
    if ( e[i] < 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "SIMPLEX01_MONOMIAL_INTEGRAL - Fatal error!\n" );
      fprintf ( stderr, "  All exponents must be nonnegative.\n" );
      exit ( 1 );
    }
  }

  k = 0;
  integral = 1.0;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 1; j <= e[i]; j++ )
    {
      k = k + 1;
      integral = integral * ( double ) ( j ) / ( double ) ( k );
    }
  }

  for ( i = 0; i < m; i++ )
  {
    k = k + 1;
    integral = integral / ( double ) ( k );
  }

  return integral;
}
/******************************************************************************/

double *simplex01_sample ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX01_SAMPLE samples the unit simplex in M dimensions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2014

  Author:

    John Burkardt

  Reference:

    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double SIMPLEX01_SAMPLE[M*N], the points.
*/
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;

  x = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    e = r8vec_uniform_01_new ( m + 1, seed );

    for ( i = 0; i < m + 1; i++ )
    {
      e[i] = - log ( e[i] );
    }
    e_sum = r8vec_sum ( m + 1, e );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = e[i] / e_sum;
    }
    free ( e );
  }

  return x;
}
/******************************************************************************/

double simplex01_volume ( int m )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX01_VOLUME returns the volume of the unit simplex in M dimensions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Output, double SIMPLEX01_VOLUME, the volume.
*/
{
  int i;
  double volume;

  volume = 1.0;
  for ( i = 1; i <= m; i++ )
  {
    volume = volume / ( double ) ( i );
  }

  return volume;
}
/******************************************************************************/

void timestamp ( )

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

