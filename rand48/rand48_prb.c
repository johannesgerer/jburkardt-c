# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RAND48_PRB.

  Discussion:

    RAND48_PRB calls sample problems for the RAND48 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "RAND48_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the RAND48 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RAND48_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests MRAND48 + SRAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  long int seed;
  long int value;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  MRAND48 returns signed long integers in [-2^31,+2^31].\n" );
  printf ( "  SRAND48 is used to initialize the seed (but only 32 bits).\n" );

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests MRAND48 + SEED48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  MRAND48 returns signed long integers in [-2^31,+2^31].\n" );
  printf ( "  SEED48 is used to initialize the seed (all 48 bits).\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests JRAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  JRAND48 returns signed long integers in [-2^31,+2^31].\n" );
  printf ( "  The 48 bit seed is an explicit argument, 3 16 bit values.\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests LRAND48 + SRAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  long int seed;
  long int value;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  LRAND48 returns unsigned long integers in [0,+2^31].\n" );
  printf ( "  SRAND48 is used to initialize the seed (32 bits only).\n" );

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests LRAND48 + SEED48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  LRAND48 returns unsigned long integers in [0,+2^31].\n" );
  printf ( "  SEED48 is used to initialize the seed (all 48 bits).\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );


  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests NRAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  NRAND48 returns nonnegative long integers in [0,+2^31].\n" );
  printf ( "  The 48 bit seed is an explicit argument of 3 16 bit values.\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );


  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );


  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    printf ( "  %6d  %12d\n", i, value );
  }
  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests DRAND48 + SRAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  long int seed;
  double value;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  DRAND48 returns double precision floating point values in [0.0,1.0].\n" );
  printf ( "  SRAND48 is used to initialize the seed (32 bits only).\n" );

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 987654321L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 123456789L;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }
  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests DRAND48 + SEED48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  double value;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  DRAND48 returns double precision real values in [0.0,1.0].\n" );
  printf ( "  The 48 bit seed is an explicit argument of 3 16 bit values.\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    printf ( "  %6d  %12f\n", i, value );
  }
  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests ERAND48.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2008

  Author:

    John Burkardt
*/
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  double value;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  ERAND48 returns double precision real values in [0.0,1.0].\n" );
  printf ( "  The 48 bit seed is an explicit argument of 3 16 bit values.\n" );

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    printf ( "  %6d  %12f\n", i, value );
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "  The seed vector is %d, %d, %d\n", seedvec[0], seedvec[1], seedvec[2] );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    printf ( "  %6d  %12f\n", i, value );
  }
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
