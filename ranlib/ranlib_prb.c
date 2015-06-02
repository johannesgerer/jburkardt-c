# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "ranlib.h"
# include "rnglib.h"

int main ( );
void test_phrtsd ( char *phrase );
void test_bot ( );
void test_genbet ( char *phrase );
void test_ignbin ( char *phrase );
void test_genchi ( char *phrase );
void test_genexp ( char *phrase );
void test_genf ( char *phrase );
void test_gengam ( char *phrase );
void test_ignnbn ( char *phrase );
void test_gennch ( char *phrase );
void test_gennf ( char *phrase );
void test_gennor ( char *phrase );
void test_ignpoi ( char *phrase );
void test_genunf ( char *phrase );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RANLIB_PRB.

  Discussion:

    RANLIB_PRB tests the RANLIB library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  char phrase[] = "randomizer";

  timestamp ( );
  printf ( "\n" );
  printf ( "RANLIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RANLIB library.\n" );

  test_phrtsd ( phrase );

  test_bot ( );

  test_genbet ( phrase );
  test_ignbin ( phrase );
  test_genchi ( phrase );
  test_genexp ( phrase );
  test_genf ( phrase );
  test_gengam ( phrase );
  test_ignnbn ( phrase );
  test_gennch ( phrase );
  test_gennf ( phrase );
  test_gennor ( phrase );
  test_ignpoi ( phrase );
  test_genunf ( phrase );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RANLIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_phrtsd ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_PHRTSD tests PHRTSD, which generates two seeds from a phrase.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  int seed1;
  int seed2;

  printf ( "\n" );
  printf ( "TEST_PHRTSD\n" );
  printf ( "  Test PHRTST,\n" );
  printf ( "  which generates two seeds from a phrase.\n" );

  printf ( "\n" );
  printf ( "  Randomizing phrase is \"%s\"\n", phrase );

  phrtsd ( phrase, &seed1, &seed2 );

  printf ( "\n" );
  printf ( "  Seed1 = %d\n", seed1 );
  printf ( "  Seed2 = %d\n", seed2 );

  return;
}
/******************************************************************************/

void test_bot ( )

/******************************************************************************/
/*
  Purpose:

    TEST_BOT is a test program for the bottom level routines

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 April 2013

  Author:

    John Burkardt
*/
{
  int answer[10000];
  int genlst[5] = { 0, 4, 9, 19, 31 };
  int ians;
  int iblock;
  int igen;
  int itmp;
  int ix;
  int ixgen;
  int nbad;
  int seed1;
  int seed2;

  printf ( "\n" );
  printf ( "TEST_BOT\n" );
  printf ( "  Test the lower level random number generators.\n" );
  printf ( "\n" );
  printf ( "  Five of the 32 generators will be tested.\n" );
  printf ( "  We generate 100000 numbers, reset the block\n" );
  printf ( "  and do it again.  No disagreements should occur.\n" );
  printf ( "\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set up all generators.
*/
  seed1 = 12345;
  seed2 = 54321;
  set_initial_seed ( seed1, seed2 );
/*
  For a selected set of generators
*/
  nbad = 0;

  for ( ixgen = 0; ixgen < 5; ixgen++ )
  {
    igen = genlst[ixgen];
    cgn_set ( igen );
    printf ( "  Testing generator %d\n", igen );
/*
  Use 10 blocks, and generate 1000 numbers per block
*/
    init_generator ( 0 );

    for ( iblock = 0; iblock < 10; iblock++ )
    {
      for ( ians = 0; ians < 1000; ians++ )
      {
        ix = ians + iblock * 1000;
        answer[ix] = i4_uni ( );
      }
      init_generator ( 2 );
    }
/*
  Do it again and compare answers
  Use 10 blocks, and generate 1000 numbers.
*/
    init_generator ( 0 );

    for ( iblock = 0; iblock < 10; iblock++ )
    {
      for ( ians = 0; ians < 1000; ians++ )
      {
        ix = ians + iblock * 1000;
        itmp = i4_uni ( );

        if ( itmp != answer[ix] )
        {
          printf ( "\n" );
          printf ( "TEST_BOT - Warning!\n" );
          printf ( "  Data disagreement:\n" );
          printf ( "  Block = %d\n", iblock );
          printf ( "  N within block = %d\n", ians );
          printf ( "  Index in ANSWER = %d\n", ix );
          printf ( "  First value =  %d\n", answer[ix] );
          printf ( "  Second value = %d\n", itmp );

          nbad = nbad + 1;

          if ( 10 < nbad )
          {
            printf ( "\n" );
            printf ( "TEST_BOT - Warning!\n" );
            printf ( "  More than 10 mismatches!\n" );
            printf ( "  Tests terminated early.\n" );
            return;
          }
        }
      }
      init_generator ( 2 );
    }
  }
  return;
}
/******************************************************************************/

void test_genbet ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENBET tests GENBET, which generates Beta deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float a;
  float *array;
  float av;
  float avtr;
  float b;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "bet";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENBET\n" );
  printf ( "  Test GENBET,\n" );
  printf ( "  which generates Beta deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = 1.0;
  high = 10.0;
  b = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = genbet ( a, b );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = a;
  param[1] = b;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_ignbin ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_IGNBIN tests IGNBIN, which generates Binomial deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 10000;
  int nn;
  float param[2];
  char pdf[] = "bin";
  float pp;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_IGNBIN\n" );
  printf ( "  Test IGNBIN,\n" );
  printf ( "  which generates binomial deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 0.5;
  high = 20.0;
  nn = ( int ) genunf ( low, high );

  low = 0.0;
  high = 1.0;
  pp = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  NN = %d\n", nn );
  printf ( "  PP = %g\n", pp );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignbin ( nn, pp );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = ( float ) ( nn );
  param[1] = pp;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_genchi ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENCHI tests GENCHI, which generates Chi-Square deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float df;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[1];
  char pdf[] = "chi";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENCHI\n" );
  printf ( "  Test GENCHI,\n" );
  printf ( "  which generates Chi-square deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 1.0;
  high = 10.0;
  df = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  DF = %g\n", df );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = genchi ( df );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = df;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_genexp ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENEXP tests GENEXP, which generates exponential deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[2];
  char pdf[] = "exp";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENEXP\n" );
  printf ( "  Test GENEXP,\n" );
  printf ( "  which generates exponential deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low =  0.5;
  high = 10.0;
  mu = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = genexp ( mu );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = mu;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_genf ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENF tests GENF, which generates F deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float dfd;
  float dfn;
  float high;
  int i;
  float low;
  int n = 10000;
  float param[2];
  char pdf[] = "f";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENF\n" );
  printf ( "  Test GENF,\n" );
  printf ( "  which generates F deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 3.0;
  high = 10.0;
  dfn = genunf ( low, high );

  low = 5.0;
  high = 10.0;
  dfd = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  DFN =   %g\n", dfn );
  printf ( "  DFD =   %g\n", dfd );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = genf ( dfn, dfd );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = dfn;
  param[1] = dfd;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_gengam ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENGAM tests GENGAM, which generates Gamma deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float a;
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "gam";
  float r;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENGAM\n" );
  printf ( "  Test GENGAM,\n" );
  printf ( "  which generates Gamma deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = 1.0;
  high = 10.0;
  r = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  A = %g\n", a );
  printf ( "  R = %g\n", r );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = gengam ( a, r );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = a;
  param[1] = r;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_ignnbn ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_IGNNBN tests IGNNBN, which generates Negative Binomial deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 10000;
  int nn;
  float param[2];
  char pdf[] = "nbn";
  float pp;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_IGNNBN\n" );
  printf ( "  Test IGNNBN,\n" );
  printf ( "  which generates negative binomial deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 3.0;
  high = 20.0;
  nn = ( int ) genunf ( low, high );

  low = 0.0;
  high = 1.0;
  pp = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  NN = %d\n", nn );
  printf ( "  PP = %g\n", pp );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignnbn ( nn, pp );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = ( float ) ( nn );
  param[1] = pp;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_gennch ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENNCH tests GENNCH, which generates noncentral Chi-Square deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float df;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "nch";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;
  float xnonc;

  printf ( "\n" );
  printf ( "TEST_GENNCH\n" );
  printf ( "  Test GENNCH,\n" );
  printf ( "  which generates noncentral Chi-square deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 2.0;
  high = 10.0;
  df = genunf ( low, high );

  low = 0.0;
  high = 2.0;
  xnonc = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  DF =    %g\n", df );
  printf ( "  XNONC = %g\n", xnonc );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennch ( df, xnonc );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = df;
  param[1] = xnonc;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_gennf ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENNF tests GENNF, which generates noncentral F deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float dfd;
  float dfn;
  float high;
  int i;
  float low;
  int n = 10000;
  float param[3];
  char pdf[] = "nf";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;
  float xnonc;

  printf ( "\n" );
  printf ( "TEST_GENNF\n" );
  printf ( "  Test GENNF,\n" );
  printf ( "  which generates noncentral F deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 3.0;
  high = 10.0;
  dfn = genunf ( low, high );

  low = 5.0;
  high = 10.0;
  dfd = genunf ( low, high );

  low = 0.0;
  high = 2.0;
  xnonc = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  DFN =   %g\n", dfn );
  printf ( "  DFD =   %g\n", dfd );
  printf ( "  XNONC = %g\n", xnonc );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennf ( dfn, dfd, xnonc );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = dfn;
  param[1] = dfd;
  param[2] = xnonc;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_gennor ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENNOR tests GENNOR, which generates normal deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[2];
  char pdf[] = "nor";
  float sd;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENNOR\n" );
  printf ( "  Test GENNOR,\n" );
  printf ( "  which generates normal deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = -10.0;
  high = 10.0;
  mu = genunf ( low, high );

  low = 0.25;
  high = 4.0;
  sd = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( " '  MU =   %g\n", mu );
  printf ( " '  SD =   %g\n", sd );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennor ( mu, sd );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = mu;
  param[1] = sd;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_ignpoi ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_IGNPOI tests IGNPOI, which generates Poisson deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[1];
  char pdf[] = "poi";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_IGNPOI\n" );
  printf ( "  Test IGNPOI,\n" );
  printf ( "  which generates Poisson deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 0.5;
  high = 20.0;
  mu = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignpoi ( mu );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = mu;

  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
/******************************************************************************/

void test_genunf ( char *phrase )

/******************************************************************************/
/*
  Purpose:

    TEST_GENUNF tests GENUNF, which generates uniform deviates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
  float a;
  float *array;
  float av;
  float avtr;
  float b;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "unf";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST_GENUNF\n" );
  printf ( "  Test GENUNF,\n" );
  printf ( "  which generates uniform deviates.\n" );
/*
  Initialize the generators.
*/
  initialize ( );
/*
  Set the seeds based on the phrase.
*/
  phrtsd ( phrase, &seed1, &seed2 );
/*
  Initialize all generators.
*/
  set_initial_seed ( seed1, seed2 );
/*
  Select the parameters at random within a given range.
*/
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = a + 1.0;
  high = a + 10.0;
  b = genunf ( low, high );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
/*
  Generate N samples.
*/
  array = ( float * ) malloc ( n * sizeof ( float ) );
  for ( i = 0; i < n; i++ )
  {
    array[i] = genunf ( a, b );
  }
/*
  Compute statistics on the samples.
*/
  stats ( array, n, &av, &var, &xmin, &xmax );
/*
  Request expected value of statistics for this distribution.
*/
  param[0] = a;
  param[1] = b;
  trstat ( pdf, param, &avtr, &vartr );

  printf ( "\n" );
  printf ( "  Sample data range:          %14g  %14g\n", xmin, xmax );
  printf ( "  Sample mean, variance:      %14g  %14g\n", av,   var );
  printf ( "  Distribution mean, variance %14g  %14g\n", avtr, vartr );

  free ( array );

  return;
}
