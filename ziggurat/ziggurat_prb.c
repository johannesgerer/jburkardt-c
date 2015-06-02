# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <stdint.h>

# include "ziggurat.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( int sample_num );
void test06 ( int sample_num );
void test07 ( int sample_num );
void test08 ( int sample_num );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ZIGGURAT_PRB.

  Discussion:

    ZIGGURAT_PRB tests the ZIGGURAT library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 October 2013

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
  Sample 10 values of the unsigned integer 32 bit generators.
*/
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
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

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SHR3_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  uint32_t seed;
  uint32_t value;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SHR3_SEEDED returns pseudorandom uniformly distributed\n" );
  printf ( "  unsigned 32 bit integers.\n" );

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
    printf ( "  %6d  %12d  %12u\n", 0, ( int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = shr3_seeded ( &seed );
      printf ( "  %6d  %12d  %12u  %12u\n", i, ( int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R4_UNI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  uint32_t seed;
  float value;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R4_UNI returns pseudorandom uniformly distributed\n" );
  printf ( "  floats (single precision real numbers) between 0 and 1.\n" );

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
    printf ( "  %6d  %12d  %12u\n", 0, ( int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_uni ( &seed );

      printf ( "  %6d  %12d  %12u  %14f\n", i, ( int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R4_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2013

  Author:

    John Burkardt
*/
{
  float fn[128];
  int i;
  int j;
  uint32_t kn[128];
  uint32_t seed;
  float value;
  float wn[128];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  R4_NOR returns pseudorandom normally distributed\n" );
  printf ( "  floats (single precision real numbers) between 0 and 1.\n" );

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
    printf ( "  %6d  %12d  %12u\n", 0, ( int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_nor ( &seed, kn, fn, wn );

      printf ( "  %6d  %12d  %12u  %14f\n", i, ( int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests R4_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2013

  Author:

    John Burkardt
*/
{
  float fe[256];
  int i;
  int j;
  uint32_t ke[256];
  uint32_t seed;
  float value;
  float we[256];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  R4_EXP returns pseudorandom exponentially distributed\n" );
  printf ( "  floats (single precision real numbers) between 0 and 1.\n" );

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
    printf ( "  %6d  %12d  %12u\n", 0, ( int ) seed, seed );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_exp ( &seed, ke, fe, we );

      printf ( "  %6d  %12d  %12u  %14f\n", i, ( int ) seed, seed, value );
    }
  }

  return;
}
/******************************************************************************/

void test05 ( int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST05 times SHR3_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 2013

  Author:

    John Burkardt
*/
{
  double ctime;
  int sample;
  uint32_t seed;
  uint32_t value;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Measure the time it takes SHR3_SEEDED to generate\n" );
  printf ( "  %d unsigned 32 bit integers.\n", sample_num );

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = shr3_seeded ( &seed );
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

    04 October 2013

  Author:

    John Burkardt
*/
{
  double ctime;
  int sample;
  uint32_t seed;
  float value;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Measure the time it takes R4_UNI to generate\n" );
  printf ( "  %d uniformly random floats.\n", sample_num );

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

    TEST07 times R4_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2013

  Author:

    John Burkardt
*/
{
  double ctime;
  float fn[128];
  uint32_t kn[128];
  int sample;
  uint32_t seed;
  float value;
  float wn[129];

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Measure the time it takes R4_NOR to generate\n" );
  printf ( "  %d normal random floats.\n", sample_num );

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

    14 October 2013

  Author:

    John Burkardt
*/
{
  double ctime;
  float fe[256];
  uint32_t ke[256];
  int sample;
  uint32_t seed;
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

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests CONG_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2013

  Author:

    John Burkardt
*/
{
  int j;
  uint32_t jcong_new;
  uint32_t jcong_in;
  uint32_t jcong_old;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  CONG_SEEDED is a generator of pseudorandom uniformly\n" );
  printf ( "  distributed unsigned 32 bit integers.\n" );
  printf ( "\n" );
  printf ( "    Input Seed   Output Seed  Output Value\n" );
  printf ( "\n" );

  jcong_new = 234567891;

  for ( j = 1; j <= 10; j++ )
  {
    jcong_old = jcong_new;
    jcong_in = jcong_new;
    jcong_new = cong_seeded ( &jcong_in );
    printf ( "  %12u  %12u  %12u\n", jcong_old, jcong_in, jcong_new );
  }

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests KISS_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2013

  Author:

    John Burkardt
*/
{
  int j;
  uint32_t jcong_in;
  uint32_t jcong_old;
  uint32_t jsr_in;
  uint32_t jsr_old;
  uint32_t w_in;
  uint32_t w_old;
  uint32_t value;
  uint32_t z_in;
  uint32_t z_old;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  KISS_SEEDED is a generator of pseudorandom uniformly\n" );
  printf ( "  distributed unsigned 32 bit integers.\n" );
  printf ( "\n" );
  printf ( "              JCONG           JSR             W             Z         Value\n" );
  printf ( "\n" );

  jcong_in = 234567891;
  jsr_in = 123456789;
  w_in = 345678912;
  z_in = 456789123;

  for ( j = 1; j <= 10; j++ )
  {
    jcong_old = jcong_in;
    jsr_old = jsr_in;
    w_old = w_in;
    z_old = z_in;
    value = kiss_seeded ( &jcong_in, &jsr_in, &w_in, &z_in );
    printf ( "  In   %12u  %12u  %12u  %12u\n", jcong_old, jsr_old, w_old, z_old );
    printf ( "  Out  %12u  %12u  %12u  %12u  %12u\n", jcong_in, jsr_in, w_in, z_in, value );
  }

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests MWC_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2013

  Author:

    John Burkardt
*/
{
  int j;
  uint32_t w_in;
  uint32_t w_old;
  uint32_t value;
  uint32_t z_in;
  uint32_t z_old;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  MWC_SEEDED is a generator of pseudorandom uniformly\n" );
  printf ( "  distributed unsigned 32 bit integers.\n" );
  printf ( "\n" );
  printf ( "       Input W       Input Z      Output W      Output Z  Output Value\n" );
  printf ( "\n" );

  w_in = 345678912;
  z_in = 456789123;

  for ( j = 1; j <= 10; j++ )
  {
    w_old = w_in;
    z_old = z_in;
    value = mwc_seeded ( &w_in, &z_in );
    printf ( "  %12u  %12u  %12u  %12u  %12u\n", w_old, z_old, w_in, z_in, value );
  }

  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests SHR3_SEEDED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2013

  Author:

    John Burkardt
*/
{
  int j;
  uint32_t jsr_new;
  uint32_t jsr_in;
  uint32_t jsr_old;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  SHR3_SEEDED is a generator of pseudorandom uniformly\n" );
  printf ( "  distributed unsigned 32 bit integers.\n" );
  printf ( "\n" );
  printf ( "    Input Seed   Output Seed  Output Value\n" );
  printf ( "\n" );

  jsr_new = 123456789;

  for ( j = 1; j <= 10; j++ )
  {
    jsr_old = jsr_new;
    jsr_in = jsr_new;
    jsr_new = shr3_seeded ( &jsr_in );
    printf ( "  %12u  %12u  %12u\n", jsr_old, jsr_in, jsr_new );
  }

  return;
}
