# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
int i4_power ( int i, int j );
void i4mat_memory_test ( int n_log );
void i4vec_memory_test ( int n_log );
float r4_cpu_time ( );
float r4_real_time ( );
void r4mat_memory_test ( int n_log );
void r4vec_memory_test ( int n_log );
double r8_cpu_time ( );
double r8_real_time ( );
void r8mat_memory_test ( int n_log );
void r8vec_memory_test ( int n_log );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MEMORY_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt
*/
{
  int n;
  int n_log;
  int n_log_max;
  int n_log_min;

  timestamp ( );

  printf ( "\n" );
  printf ( "MEMORY_TEST\n" );
  printf ( "  C version\n" );
  printf ( "  Try to see how big vectors and matrices can be.\n" );

  if ( argc <= 1 )
  {
    n_log_min = 0;
    printf ( "\n" );
    printf ( "  Using default value of N_LOG_MIN = %d\n", n_log_min );
  }
  else
  {
    n_log_min = atoi ( argv[1] );
    printf ( "\n" );
    printf ( "  User value of N_LOG_MIN = %d\n", n_log_min );
  }
  if ( argc <= 2 )
  {
    n_log_max = 27;
    printf ( "\n" );
    printf ( "  Using default value of N_LOG_MAX = %d\n", n_log_max );
  }
  else
  {
    n_log_max = atoi ( argv[2] );
    printf ( "\n" );
    printf ( "  User value of N_LOG_MAX = %d\n", n_log_max );
  }
/*
  I4VEC test.
*/
  printf ( "\n" );
  printf ( "I4VEC Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    i4vec_memory_test ( n_log );
  }
/*
  R4VEC test.
*/
  printf ( "\n" );
  printf ( "R4VEC Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r4vec_memory_test ( n_log );
  }
/*
  R8VEC test.
*/
  printf ( "\n" );
  printf ( "R8VEC Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r8vec_memory_test ( n_log );
  }
/*
  I4MAT test.
*/
  printf ( "\n" );
  printf ( "I4MAT Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N            N1            N2     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    i4mat_memory_test ( n_log );
  }
/*
  R4MAT test.
*/
  printf ( "\n" );
  printf ( "R4MAT Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N            N1            N2     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r4mat_memory_test ( n_log );
  }
/*
  R8MAT test.
*/
  printf ( "\n" );
  printf ( "R8MAT Memory Test\n" );
  printf ( "\n" );
  printf ( "Log2(N)            N            N1            N2     Ave       CPU        Real\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r8mat_memory_test ( n_log );
  }

  printf ( "\n" );
  printf ( "MEMORY_TEST\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2004

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void i4mat_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    I4MAT_MEMORY_TEST declares and uses an I4MAT of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int *i4mat;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  i4mat = ( int * ) malloc ( n1 * n2 * sizeof ( int ) );

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( float ) random ( ) / ( float ) RAND_MAX;
      i4mat[i+j*n1] = ( int ) ( 3.0 * x );
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + ( float ) i4mat[i+j*n1];
    }
  }
  average = average / ( float ) n1 / ( float ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %12d  %12d  %4f  %10e  %10e\n",
    n1, n2, average, cpu_diff, real_diff );

  free ( i4mat );

  return;
}
/******************************************************************************/

void i4vec_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MEMORY_TEST declares and uses an I4VEC of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int *i4vec;
  int n;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  i4vec = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    x = ( float ) random ( ) / ( float ) RAND_MAX;
    i4vec[i] = ( int ) ( 3.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + ( float ) i4vec[i];
  }
  average = average / ( float ) n;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %4.2f  %10.2e  %10.2e\n", average, cpu_diff, real_diff );

  free ( i4vec );

  return;
}
/******************************************************************************/

float r4_cpu_time ( )

/******************************************************************************/
/*
  Purpose:

    R4_CPU_TIME reports the elapsed CPU time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2006

  Author:

    John Burkardt

  Parameters:

    Output, float R4_CPU_TIME, the current total elapsed CPU time in second.
*/
{
  float value;

  value = ( float ) clock ( ) / ( float ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

float r4_real_time ( )

/******************************************************************************/
/*
  Purpose:

    R4_REAL_TIME returns the current real time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2006

  Author:

    John Burkardt

  Parameters:

    Output, float R4_REAL_TIME, the real time in seconds.
*/
{
  time_t now;
  static time_t then = 0;
  float value;

  if ( then == 0 )
  {
    then = time ( NULL );
  }

  now = time ( NULL );

  value = ( float ) difftime ( now, then );

  return value;
}
/******************************************************************************/

void r4mat_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    R4MAT_MEMORY_TEST declares and uses an R4MAT of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  float *r4mat;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r4mat = ( float * ) malloc ( n1 * n2 * sizeof ( float ) );

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( float ) random ( ) / ( float ) RAND_MAX;
      r4mat[i+j*n1] = 2.0 * x;
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + r4mat[i+j*n1];
    }
  }
  average = average / ( float ) n1 / ( float ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %12d  %12d  %4.2f  %10.2e  %10.2e\n",
    n1, n2, average, cpu_diff, real_diff );

  free ( r4mat );

  return;
}
/******************************************************************************/

void r4vec_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    R4VEC_MEMORY_TEST declares and uses an R4VEC of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int n;
  float *r4vec;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r4vec = ( float * ) malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    x = ( float ) random ( ) / ( float ) RAND_MAX;
    r4vec[i] = ( 2.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + r4vec[i];
  }
  average = average / ( float ) n;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %4.2f  %10.2e  %10.2e\n", average, cpu_diff, real_diff );

  free ( r4vec );

  return;
}
/******************************************************************************/

double r8_cpu_time ( )

/******************************************************************************/
/*
  Purpose:

    R8_CPU_TIME reports the elapsed CPU time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2006

  Author:

    John Burkardt

  Parameters:

    Output, double R8_CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

double r8_real_time ( )

/******************************************************************************/
/*
  Purpose:

    R8_REAL_TIME returns the current real time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2006

  Author:

    John Burkardt

  Parameters:

    Output, double R8_REAL_TIME, the real time in seconds.
*/
{
  time_t now;
  static time_t then = 0;
  double value;

  if ( then == 0 )
  {
    then = time ( NULL );
  }

  now = time ( NULL );

  value = difftime ( now, then );

  return value;
}
/******************************************************************************/

void r8mat_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MEMORY_TEST declares and uses an R8MAT of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  double average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  double *r8mat;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  double x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r8mat = ( double * ) malloc ( n1 * n2 * sizeof ( double ) );

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( double ) random ( ) / ( double ) RAND_MAX;
      r8mat[i+j*n1] = 2.0 * x;
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + r8mat[i+j*n1];
    }
  }
  average = average / ( double ) n1 / ( double ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %12d  %12d  %4.2f  %10.2e  %10.2e\n",
    n1, n2, average, cpu_diff, real_diff );

  free ( r8mat );

  return;
}
/******************************************************************************/

void r8vec_memory_test ( int n_log )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEMORY_TEST declares and uses an R8VEC of size N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N_LOG, the logarithm base 2 of N.
*/
{
  double average;
  double cpu_diff;
  double cpu_time1;
  double cpu_time2;
  int i;
  int n;
  double *r8vec;
  double real_diff;
  double real_time1;
  double real_time2;
  unsigned int seed = 123456789;
  double x;

  n = i4_power ( 2, n_log );

  printf ( "  %4d  %12d", n_log, n );

  srandom ( seed );

  real_time1 = r8_real_time ( );
  cpu_time1 = r8_cpu_time ( );

  r8vec = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x = ( double ) random ( ) / ( double ) RAND_MAX;
    r8vec[i] = ( 2.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + r8vec[i];
  }
  average = average / ( double ) n;

  cpu_time2 = r8_cpu_time ( );
  real_time2 = r8_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  printf ( "  %4.2f  %10.2e  %10.2e\n", average, cpu_diff, real_diff );

  free ( r8vec );

  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    May 31 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2003

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
