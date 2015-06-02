# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "poisson_simulation.h"

/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4vec_max ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MAX returns the value of the maximum element in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int IVEC_MAX, the value of the maximum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
/******************************************************************************/

double i4vec_mean ( int n, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MEAN returns the mean of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 May 1999

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int X[N], the vector whose mean is desired.

    Output, double I4VEC_MEAN, the mean, or average, of the vector entries.
*/
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + ( double ) x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
/******************************************************************************/

int i4vec_min ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MIN returns the minimum element in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int I4VEC_MIN, the value of the minimum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

double i4vec_variance ( int n, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_VARIANCE returns the variance of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 May 1999

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int X[N], the vector whose variance is desired.

    Output, double I4VEC_VARIANCE, the variance of the vector entries.
*/
{
  int i;
  double mean;
  double variance;

  if ( n < 2 )
  {
    variance = 0.0;
  }
  else
  {
    mean = i4vec_mean ( n, x );

    variance = 0.0;
    for ( i = 0; i < n; i++ )
    {
      variance = variance + pow ( ( double ) x[i] - mean, 2 );
    }

    variance = variance / ( double ) ( n - 1 );
  }

  return variance;
}
/******************************************************************************/

void poisson_fixed_events ( double lambda, int event_num, int *seed, 
  double t[], double w[] )

/******************************************************************************/
/*
  Purpose:

    POISSON_FIXED_EVENTS waits for a given number of Poisson events.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, double LAMBDA, the average number of events per 
    unit time.

    Input, int EVENT_NUM, the number of events to wait for.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double T[EVENT_NUM+1], the time at which a total 
    of 0, 1, 2, ... and EVENT_NUM events were observed.

    Output, double W[EVENT_NUM+1], the waiting time until the
    I-th event occurred.
*/
{
  int i;
  double *u;
/*
  Poisson waiting times follow an exponential distribution.
*/
  w[0] = 0.0;
  u = r8vec_uniform_01_new ( event_num, seed );
  for ( i = 1; i <= event_num; i++ )
  {
    w[i] = - log ( u[i-1] ) / lambda;
  }
/*
  The time til event I is the sum of the waiting times 0 through I.
*/
  r8vec_cum ( event_num + 1, w, t );

  free ( u );

  return;
}
/******************************************************************************/

int poisson_fixed_time ( double lambda, double time, int *seed )

/******************************************************************************/
/*
  Purpose:

    POISSON_FIXED_TIME counts the Poisson events in a fied time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, double LAMBDA, the average number of events 
    per unit time.

    Input, double TIME, the amount of time to observe.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, int POISSON_FIXED_TIME, the number of Poisson events observed.
*/
{
  double dt;
  int n;
  double t;
  double u;

  n = 0;
  t = 0.0;

  while ( t < time )
  {
    u = r8_uniform_01 ( seed );
    dt = - log ( u ) / lambda;
    n = n + 1;
    t = t + dt;
  }

  return n;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void r8vec_cum ( int n, double a[], double a_cum[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CUM computes the cumulutive sums of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Input:

      A = { 1.0, 2.0, 3.0, 4.0 }

    Output:

      A_CUM = { 1.0, 3.0, 6.0, 10.0 }

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector to be summed.

    Output, double A_CUM[N], the cumulative sums.
*/
{
  int i;

  a_cum[0] = a[0];

  for ( i = 1; i < n; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i];
  }

  return;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  double *r8vec_pointer;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_mean ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEAN returns the mean of a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector whose mean is desired.

    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
*/
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
/******************************************************************************/

double *r8vec_midspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIDSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    This function divides the interval [a,b] into n subintervals, and then
    returns the midpoints of those subintervals.

  Example:

    N = 5, A = 10, B = 20
    X = [ 11, 13, 15, 17, 19 ]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the endpoints of the interval.

    Output, double R8VEC_MIDSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  double *x;
  int i;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * a 
           + ( double ) (         2 * i + 1 ) * b ) 
           / ( double ) ( 2 * n );
  }

  return x;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
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

double r8vec_variance ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_VARIANCE returns the variance of an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector whose variance is desired.

    Output, double R8VEC_VARIANCE, the variance of the vector entries.
*/
{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
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
