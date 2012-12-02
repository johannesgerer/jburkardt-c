# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "ziggurat.h"

int main ( void );
double cpu_time ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( int sample_num );
void test06 ( int sample_num );
void test07 ( int sample_num );
void test08 ( int sample_num );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ZIGGURAT_PRB.

  Discussion:

    ZIGGURAT_PRB calls sample problems for the ZIGGURAT library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  int sample_num = 1000000;

  timestamp ( );

  printf ( "\n" );
  printf ( "ZIGGURAT_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the ZIGGURAT library.\n" );
/*
  Make sure that SEED controls the sequence, and can restart it.
*/
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Measure the time it takes to generate 10,000 variables.
*/
  test05 ( sample_num );
  test06 ( sample_num );
  test07 ( sample_num );
  test08 ( sample_num );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ZIGGURAT_PRB\n" );
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

    TEST01 tests SHR3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  unsigned long int seed;
  unsigned long int value;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SHR3 returns pseudorandom uniformly distributed\n" );
  printf ( "  unsigned long integers.\n" );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    printf ( "\n" );
    printf ( "  %6d  %12d  %12lu\n", 0, ( long int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = shr3 ( &seed );
      printf ( "  %6d  %12d  %12lu  %12d\n", i, ( long int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R4_UNI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  unsigned long int seed;
  float value;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R4_UNI returns pseudorandom uniformly distributed\n" );
  printf ( "  floats between 0 and 1.\n" );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    printf ( "\n" );
    printf ( "  %6d  %12d  %12lu\n", 0, ( long int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_uni ( &seed );

      printf ( "  %6d  %12d  %12lu  %14f\n", i, ( long int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R4_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  float fn[128];
  int i;
  int j;
  int kn[128];
  unsigned long int seed;
  float value;
  float wn[128];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  R4_NOR returns pseudorandom normally distributed\n" );
  printf ( "  real numbers between 0 and 1.\n" );

  r4_nor_setup ( kn, fn, wn );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    printf ( "\n" );
    printf ( "  %6d  %12d  %12lu\n", 0, ( long int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_nor ( &seed, kn, fn, wn );

      printf ( "  %6d  %12d  %12lu  %14f\n", i, ( long int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests R4_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  float fe[256];
  int i;
  int j;
  int ke[256];
  unsigned long seed;
  float value;
  float we[256];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  R4_EXP returns pseudorandom exponentially distributed\n" );
  printf ( "  real numbers between 0 and 1.\n" );

  r4_exp_setup ( ke, fe, we );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    printf ( "\n" );
    printf ( "  %6d  %12d  %12lu\n", 0, ( long int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_exp ( &seed, ke, fe, we );

      printf ( "  %6d  %12d  %12lu  %14f\n", i, ( long int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test05 ( int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST05 times SHR3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  double ctime;
  int sample;
  unsigned long int seed;
  unsigned long int value;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Measure the time it takes SHR3 to generate\n" );
  printf ( "  %d unsigned long int values.\n", sample_num );

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = shr3 ( &seed );
  }
  ctime = cpu_time ( ) - ctime;

  printf ( "\n" );
  printf ( "  %f seconds\n", ctime );

  return;
}
/******************************************************************************/

void test06 ( int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST06 times R4_UNI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  double ctime;
  int sample;
  unsigned long int seed;
  float value;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Measure the time it takes R4_UNI to generate\n" );
  printf ( "  %d uniformly random float values.\n", sample_num );

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_uni ( &seed );
  }
  ctime = cpu_time ( ) - ctime;

  printf ( "\n" );
  printf ( "  %f seconds\n", ctime );

  return;
}
/******************************************************************************/

void test07 ( int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST07 times R8_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  double ctime;
  float fn[128];
  int kn[128];
  int sample;
  unsigned long seed;
  float value;
  float wn[129];

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Measure the time it takes R8_NOR to generate\n" );
  printf ( "  %d normal random float values.\n", sample_num );

  r4_nor_setup ( kn, fn, wn );

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_nor ( &seed, kn, fn, wn );
  }
  ctime = cpu_time ( ) - ctime;

  printf ( "\n" );
  printf ( "  %f seconds\n", ctime );

  return;
}
/******************************************************************************/

void test08 ( int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST08 times R4_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt
*/
{
  double ctime;
  float fe[256];
  int ke[256];
  int sample;
  unsigned long int seed;
  float value;
  float we[256];

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Measure the time it takes R4_EXP to generate\n" );
  printf ( "  %d exponential random float values.\n", sample_num );

  r4_exp_setup ( ke, fe, we );

  seed = 123456789;


  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_exp ( &seed, ke, fe, we );
  }
  ctime = cpu_time ( ) - ctime;

  printf ( "\n" );
  printf ( "  %f seconds\n", ctime );

  return;
}
/******************************************************************************/

double cpu_time ( void )

/******************************************************************************/
/*
  Purpose:

    CPU_TIME returns the current reading on the CPU clock.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 December 2008

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
