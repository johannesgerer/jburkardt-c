# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "line_monte_carlo.h"

/******************************************************************************/

double line01_length ( )

/******************************************************************************/
/*
  Purpose:

    LINE01_LENGTH: length of the unit line in 1D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2014

  Author:

    John Burkardt

  Parameters:

    Output, double LINE01_LENGTH, the length.
*/
{
  double length;

  length = 1.0;

  return length;
}
/******************************************************************************/

double line01_monomial_integral ( int e )

/******************************************************************************/
/*
  Purpose:

    LINE01_MONOMIAL_INTEGRAL: integrals on the unit line in 1D.

  Discussion:

    The integration region is 

      0 <= X <= 1.

    The monomial is F(X) = X^E.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2014

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.

  Parameters:

    Input, int E, the exponent.  E must be nonnegative.

    Output, double LINE01_MONOMIAL_INTEGRAL, the integral.
*/
{
  int i;
  double integral;

  if ( e == -1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CIRCLE01_MONOMIAL_INTEGRAL - Fatal error!\n" );
    fprintf ( stderr, "  E must not equal -1.\n" );
    exit ( 1 );
  }

  integral = 1.0 / ( double ) ( e + 1 );

  return integral;
}
/******************************************************************************/

double *line01_sample ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    LINE01_SAMPLE samples points on the unit line in 1D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2014

  Author:

    John Burkardt

  Reference:

    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.

    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity 
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.

  Parameters:

    Input, int N, the number of points.

    Input/output, int *SEED, a seed for the random 
    number generator.

    Output, double X[N], the points.
*/
{
  double *x;

  x = r8vec_uniform_01_new ( n, seed );

  return x;
}
/******************************************************************************/

double *monomial_value_1d ( int n, int e, double x[] )

/******************************************************************************/
/*
  Purpose:

    MONOMIAL_VALUE_1D evaluates a monomial in 1D.

  Discussion:

    This routine evaluates a monomial of the form

      x^e

    where the exponent is a nonnegative integer.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points at which the
    monomial is to be evaluated.

    Input, int E, the exponent.

    Input, double X[M*N], the point coordinates.

    Output, double MONOMIAL_VALUE_1D[N], the value of the monomial.
*/
{
  int j;
  double *v;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    v[j] = pow ( x[j], e );
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

