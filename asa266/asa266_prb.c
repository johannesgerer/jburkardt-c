# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa266.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test085 ( );
void test09 ( );
void test10 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA266_PRB.

  Discussion:

    ASA266_PRB tests the ASA266 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA266_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA266 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test085 ( );
  test09 ( );
  test10 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA266_PRB:\n" );
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

    TEST01 tests ALNORM, NORMP, NPROB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  double ccdf1;
  double ccdf2;
  double ccdf3;
  double cdf1;
  double cdf2;
  double cdf3;
  int i;
  int ntest = 16;
  double pdf2;
  double pdf3;
  int upper;
  double x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ALNORM,\n" );
  printf ( "  NORMP, and\n" );
  printf ( "  NPROB are routines that compute the cumulative\n" );
  printf ( "  density function for the normal distribution.\n" );
  printf ( "\n" );
  printf ( "  X  CDF1  1-CDF1\n" );
  printf ( "     CDF2  1-CDF2  PDF2\n" );
  printf ( "     CDF3  1-CDF3  PDF3\n" );

  for ( i = 0; i < ntest; i++ )
  {
    x = 3.0 * ( double ) ( i ) / ( double ) ( ntest - 1 );

    upper = 0;
    cdf1 = alnorm ( x, upper );

    upper = 1;
    ccdf1 = alnorm ( x, upper );

    normp ( x, &cdf2, &ccdf2, &pdf2 );

    nprob ( x, &cdf3, &ccdf3, &pdf3 );

    printf ( "\n" );
    printf ( "%14.6g%14.6g%14.6g\n", x, cdf1, ccdf1 );
    printf ( "              %14.6g%14.6g%14.6g\n", cdf2, ccdf2, pdf2 );
    printf ( "              %14.6g%14.6g%14.6g\n", cdf3, ccdf3, pdf3 );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PPND, PPND16.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  int ifault;
  int ntest = 9;
  double x1;
  double x2;

  ifault = 0;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PPND,\n" );
  printf ( "  PPND16 compute the percentage \n" );
  printf ( "  points of the normal distribution.\n" );
  printf ( "\n" );
  printf ( "           CDF     PPND(CDF)   PPND16(CDF)\n" );
  printf ( "\n" );

  for ( i = 1; i <= ntest; i++ )
  {
    cdf = ( double ) ( i ) / ( double ) ( ntest + 1 );
    x1 = ppnd ( cdf, &ifault );
    x2 = ppnd16 ( cdf, &ifault );

    printf ( "%14.6g%14.6g%14.6g\n", cdf, x1, x2 );
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests DIGAMMA, R8_PSI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  int i;
  int ntest = 10;
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  digamma(X) = d ( Log ( Gamma ( X ) ) ) / dX.\n" );
  printf ( "\n" );
  printf ( "  DIGAMMA and\n" );
  printf ( "  R8_PSI compute the digamma function:\n" );
  printf ( "\n" );
  printf ( "             X       DIGAMMA        R8_PSI\n" );
  printf ( "\n" );

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    printf ( "%14.6g%14.6g%14.6g\n", x, digamma ( x ), r8_psi ( x ) );
  }
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests TRIGAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  int i;
  int ifault;
  int ntest = 10;
  double t;
  double x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  TRIGAMMA computes the trigamma function:\n" );
  printf ( "    trigamma(X) = d^2 ( Log ( Gamma ( X ) ) ) / dX^2.\n" );
  printf ( "\n" );
  printf ( "             X       TRIGAMMA\n" );
  printf ( "\n" );

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    t = trigamma ( x, &ifault );
    printf ( "%14.6g%14.6g\n", x, t );
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests ALNGAM, ALOGAM, R8_GAMMA_LOG, LNGAMMA;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  int i;
  int ifault;
  double log1;
  double log2;
  double log3;
  double log4;
  int ntest = 10;
  double x;

  ifault = 0;
  
  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  ALNGAM\n" );
  printf ( "  ALOGAM,\n" );
  printf ( "  R8_GAMMA_LOG, and\n" );
  printf ( "  LNGAMMA compute the logarithm of the gamma function.\n" );
  printf ( "\n" );
  printf ( "             X        ALNGAM        ALOGAM    R8_GAMMA_LOG     LNGAMMA\n" );
  printf ( "\n" );

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    log1 = alngam ( x, &ifault );
    log2 = alogam ( x, &ifault );
    log3 = r8_gamma_log ( x );
    log4 = lngamma ( x, &ifault );

    printf ( "%14.6g%14.6g%14.6g%14.6g%14.6g\n", x, log1, log2, log3, log4 );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests GAMAIN, GAMMDS, GAMMAD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  double g1;
  double g2;
  double g3;
  int i;
  int ifault;
  int j;
  int ntest = 10;
  double p;
  double x;

  ifault = 0;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  GAMAIN, \n" );
  printf ( "  GAMMDS and \n" );
  printf ( "  GAMMAD compute the incomplete Gamma integral.\n" );
  printf ( "\n" );
  printf ( "             X             P        GAMMDS        GAMMAD        GAMAIN\n" );
  printf ( "\n" );

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    printf ( "\n" );
    for ( j = 1; j <= ntest; j++ )
    {
      p = ( double ) ( j ) / ( double ) ( ntest );
      g1 = gammds ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g1 = -99.0;
      }

      g2 = gammad ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g2 = -99.0;
      }

      g3 = gamain ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g3 = - 99.0;
      }
      printf ( "%14.6g%14.6g%14.6g%14.6g%14.6g\n", x, p, g1, g2, g3 );
    }
  }
  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests PPCHI2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  double gg;
  int i;
  int ifault;
  int j;
  int nitest = 9;
  int njtest = 9;
  double v;
  double x1;

  ifault = 0;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  PPCHI2 computes the percentage points\n" );
  printf ( "  of the chi squared distribution.\n" );
  printf ( "\n" );
  printf ( "      CDF      PPCHI2(CDF)\n" );
  printf ( "\n" );

  for ( j = 1; j <= njtest; j++ )
  {
    v = ( double ) ( j );

    printf ( "\n" );
    printf ( "  For Chi^2 parameter value = %g\n", v );
    printf ( "\n" );

    for ( i = 1; i <= nitest; i++ )
    {
      cdf = ( double ) ( i ) / ( double ) ( nitest + 1 );
      gg = alngam ( v / 2.0, &ifault );
      x1 = ppchi2 ( cdf, v, gg, &ifault );
      printf ( "%14.6g%14.6g\n", cdf, x1 );
    }
  }
  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE.

  Discussion:

    Canned data is used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2013

  Author:

    John Burkardt
*/
{
# define ELEM_NUM 3
# define SAMPLE_NUM 23

  double *alpha;
  double alpha_sum;
  double aminus;
  double aplus;
  int elem_i;
  int elem_num = ELEM_NUM;
  double eps;
  double *g;
  int ifault;
  int init;
  double *mean;
  int niter;
  double rlogl;
  double s;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  double *v;
  double vari;
  double *variance;
  double x[SAMPLE_NUM*ELEM_NUM] = {
    0.178, 0.162, 0.083, 0.087, 0.078, 0.040, 0.049, 0.100, 0.075, 0.084,
    0.060, 0.089, 0.050, 0.073, 0.064, 0.085, 0.094, 0.014, 0.060, 0.031,
    0.025, 0.045, 0.0195,
    0.346, 0.307, 0.448, 0.474, 0.503, 0.456, 0.363, 0.317, 0.394, 0.445,
    0.435, 0.418, 0.485, 0.378, 0.562, 0.465, 0.388, 0.449, 0.544, 0.569,
    0.491, 0.613, 0.526,
    0.476, 0.531, 0.469, 0.439, 0.419, 0.504, 0.588, 0.583, 0.531, 0.471,
    0.505, 0.493, 0.465, 0.549, 0.374, 0.450, 0.518, 0.537, 0.396, 0.400,
    0.484, 0.342, 0.4545 };

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  For samples of a Dirichlet PDF,\n" );
  printf ( "  DIRICHLET_ESTIMATE estimates the parameters.\n" );
  printf ( "  DIRICHLET_MEAN finds the means;\n" );
  printf ( "  DIRICHLET_VARIANCE finds the variances;\n" );

  r8mat_print ( sample_num, elem_num, x, "  Sampled data:" );
/*
  Compute the observed averages.
*/
  mean = r8col_mean ( sample_num, elem_num, x );

  variance = r8col_variance ( sample_num, elem_num, x );

  printf ( "\n" );
  printf ( "  Observed means, variances are:\n" );
  printf ( "\n" );
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    printf ( "%6d%14.6g%14.6g\n", elem_i, mean[elem_i], variance[elem_i] );
  }

  init = 1;
  alpha = ( double * ) malloc ( elem_num * sizeof ( double ) );
  g = ( double * ) malloc ( elem_num * sizeof ( double ) );
  v = ( double *) malloc ( elem_num * elem_num * sizeof ( double ) );

  dirichlet_estimate ( elem_num, sample_num, x, sample_num, 
    init, alpha, &rlogl, v, g, &niter, &s, &eps, &ifault );

  if ( ifault != 0 )
  {
    printf ( "\n" );
    printf ( "WARNING!\n" );
    printf ( "  DIRICHLET_ESTIMATE error code:\n" );
    printf ( "  IFAULT = %d\n", ifault );
  }

  printf ( "\n" );
  printf ( "  Index, Estimate, Lower Limit, Upper Limit:\n" );
  printf ( "\n" );

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = alpha[elem_i] - 1.96 * sqrt ( vari );
    aplus = alpha[elem_i] + 1.96 * sqrt ( vari );
    printf ( "%6d%14.6g%14.6g%14.6g\n", elem_i, alpha[elem_i], aminus, aplus );
  }

  free ( mean );
  free ( variance );

  mean = dirichlet_mean ( elem_num, alpha );

  variance = dirichlet_variance ( elem_num, alpha );

  printf ( "\n" );
  printf ( "  Expected means, variances are:\n" );
  printf ( "\n" );
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    printf ( "%6d%14.6g%14.6g\n", elem_i, mean[elem_i], variance[elem_i] );
  }

  alpha_sum = r8vec_sum ( elem_num, alpha );

  printf ( "\n" );
  printf ( "  Alpha sum is %g\n", alpha_sum );
  printf ( "\n" );
  printf ( "  NORMALIZED VALUES:\n" );
  printf ( "  Index, Estimate, Lower Limit, Upper Limit:\n" );
  printf ( "\n" );

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = ( alpha[elem_i] - 1.96 * sqrt ( vari ) ) / alpha_sum;
    aplus = ( alpha[elem_i] + 1.96 * sqrt ( vari ) ) / alpha_sum;
    printf ( "%6d%14.6g%14.6g%14.6g\n", 
      elem_i, alpha[elem_i] / alpha_sum, aminus, aplus );
  }

  printf ( "\n" );
  printf ( "  Log likelikhood function = %g\n", rlogl );

  free ( alpha );
  free ( g );
  free ( mean );
  free ( v );
  free ( variance );

  return;
# undef ELEM_NUM
# undef SAMPLE_NUM
}
/******************************************************************************/

void test085 ( )

/******************************************************************************/
/*
  Purpose:

    TEST085 tests GAMMA_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int rep;
  int rep_num = 5;
  int seed;
  int test;
  int test_num = 5;
  double x;

  printf ( "\n" );
  printf ( "TEST085\n" );
  printf ( "  GAMMA_SAMPLE samples a Gamma distribution.\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8_uniform_ab ( 0.1, 2.0, &seed );
    b = r8_uniform_ab ( 0.1, 2.0, &seed );
    printf ( "\n" );
    printf ( "  A = %g, B = %g\n", a, b );
    for ( rep = 1; rep <= rep_num; rep++ )
    {
      x = gamma_sample ( a, b, &seed );
      printf ( "  %2d  %14.6g\n", rep, x );
    }
  }
  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE, DIRICHLET_SAMPLE.

  Discussion:

    Data is generated by sampling a distribution with known parameters.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2008

  Author:

    John Burkardt
*/
{
# define ELEM_NUM 3
# define SAMPLE_NUM 1000

  double alpha[ELEM_NUM] = { 3.22, 20.38, 21.68 };
  double alpha_sum;
  double aminus;
  double aplus;
  int elem_i;
  int elem_num = ELEM_NUM;
  double eps;
  double *g;
  int ifault;
  int init;
  double *mean;
  int niter;
  double rlogl;
  double s;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  int seed;
  double *v;
  double vari;
  double *variance;
  double *x_sample;
  double *x;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  For a Dirichlet distribution,\n" );
  printf ( "  DIRICHLET_SAMPLE samples;\n" );
  printf ( "  DIRICHLET_MEAN finds the means;\n" );
  printf ( "  DIRICHLET_VARIANCE finds the variances;\n" );
  printf ( "  DIRICHLET_ESTIMATE estimates the parameters.\n" );
/*
  Report.
*/
  r8vec_print ( elem_num, alpha, "  Distribution parameters:" );

  mean = dirichlet_mean ( elem_num, alpha );

  variance = dirichlet_variance ( elem_num, alpha );

  printf ( "\n" );
  printf ( "  Distribution means, variances are:\n" );
  printf ( "\n" );
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    printf ( "%6d%14.6g%14.6g\n", elem_i, mean[elem_i], variance[elem_i] );
  }
/*
  Sample the distribution.
*/
  x_sample = ( double * ) malloc ( sample_num * elem_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "  Number of samples is %d\n", sample_num );

  for ( sample_i = 0; sample_i < sample_num; sample_i++ )
  {
    x = dirichlet_sample ( elem_num, alpha, &seed );

    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      x_sample[sample_i+elem_i*sample_num] = x[elem_i];
    }
    free ( x );
  }
/*
  Print some results.
*/
  printf ( "\n" );
  printf ( "  First few samples:\n" );
  printf ( "\n" );

  for ( sample_i = 0; sample_i < i4_min ( sample_num, 10 ); sample_i++ )
  {
    printf ( "%6d", sample_i );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      printf ( "%14.6g", x_sample[sample_i+elem_i*sample_num] );
    }
    printf ( "\n" );
  }
/*
  Compute means, variances.
*/
  free ( mean );
  free ( variance );

  mean = r8col_mean ( sample_num, elem_num, x_sample );

  variance = r8col_variance ( sample_num, elem_num, x_sample );

  printf ( "\n" );
  printf ( "  Observed means, variances are:\n" );
  printf ( "\n" );
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    printf ( "%6d%14.6g%14.6g\n", elem_i, mean[elem_i], variance[elem_i] );
  }
/*
  Destroy the values of ALPHA.
*/
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    alpha[elem_i] = 0.0;
  }
/*
  Try to recover the values of ALPHA.
*/
  init = 1;
  v = ( double * ) malloc ( elem_num * elem_num * sizeof ( double ) );
  g = ( double * ) malloc ( elem_num * sizeof ( double ) );

  dirichlet_estimate ( elem_num, sample_num, x_sample, sample_num, 
    init, alpha, &rlogl, v, g, &niter, &s, &eps, &ifault );

  if ( ifault != 0 )
  {
    printf ( "\n" );
    printf ( "Warning!\n" );
    printf ( "  DIRICHLET_ESTIMATE error code:\n" );
    printf ( "  IFAULT = %d\n", ifault );
  }

  printf ( "\n" );
  printf ( "  Index, Estimate, Lower Limit, Upper Limit:\n" );
  printf ( "\n" );

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = alpha[elem_i] - 1.96 * sqrt ( vari );
    aplus = alpha[elem_i] + 1.96 * sqrt ( vari );
    printf ( "%6d%14.6g%14.6g%14.6g\n",
      elem_i, alpha[elem_i], aminus, aplus );
  }

  alpha_sum = r8vec_sum ( elem_num, alpha );

  printf ( "\n" );
  printf ( "  Alpha sum is %g\n", alpha_sum );
  printf ( "\n" );
  printf ( "  NORMALIZED VALUES:\n" );
  printf ( "  Index, Estimate, Lower Limit, Upper Limit:\n" );
  printf ( "\n" );

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = ( alpha[elem_i] - 1.96 * sqrt ( vari ) ) / alpha_sum;
    aplus = ( alpha[elem_i] + 1.96 * sqrt ( vari ) ) / alpha_sum;
    printf ( "%6d%14.6g%14.6g%14.6g\n",
      elem_i, alpha[elem_i] / alpha_sum, aminus, aplus );
  }

  printf ( "\n" );
  printf ( "  Log likelikhood function = %g\n", rlogl );

  free ( mean );
  free ( v );
  free ( variance );
  free ( x_sample );

  return;
# undef ELEM_NUM
# undef SAMPLE_NUM
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests DIRICHLET_MIX_SAMPLE, DIRICHLET_MIX_MEAN, DIRICHLET_MIX_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2013

  Author:

    John Burkardt
*/
{
# define COMP_NUM 3
# define COMP_MAX 3
# define ELEM_NUM 3
# define SAMPLE_NUM 200

  double a[ELEM_NUM];
  double alpha[COMP_MAX*ELEM_NUM] = {
    0.05, 0.85, 0.00,
    0.20, 0.10, 0.50,
    0.75, 0.05, 0.50 };
  int comp_max = COMP_MAX;
  int comp_num = COMP_NUM;
  int *comp_sample;
  int comp_i;
  double comp_weight[COMP_NUM] = { 3.0, 2.0, 1.0 };
  int elem_i;
  int elem_num = ELEM_NUM;
  double *mean;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  int seed;
  double *variance;
  double *x;
  double *x_sample;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  For a Dirichlet mixture distribution,\n" );
  printf ( "  DIRICHLET_MIX_SAMPLE samples;\n" );
  printf ( "  DIRICHLET_MIX_MEAN computes means;\n" );
  printf ( "  DIRICHLET_MIX_VARIANCE computes variances.\n" );
/*
  Report.
*/
  r8vec_print ( comp_num, comp_weight, "  Component weight:" );

  printf ( "\n" );
  printf ( "  Component  Parameters Means Variances\n" );
  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    printf ( "\n" );
    printf ( "%6d\n", comp_i );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      a[elem_i] = alpha[comp_i+elem_i*comp_max];
    }
    mean = dirichlet_mean ( elem_num, a );
    variance = dirichlet_variance ( elem_num, a );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      printf ( "%6d  %10.6f%10.6f%10.6f\n", elem_i, 
        alpha[comp_i+elem_i*comp_max], mean[elem_i], variance[elem_i] );
    }
    free ( mean );
    free ( variance );
  }

  mean = dirichlet_mix_mean ( comp_max, comp_num, elem_num, alpha, 
    comp_weight );

  r8vec_print ( elem_num, mean, "  Element means:" );
  free ( mean );
/*
  Sample the distribution.
*/
  comp_sample = ( int * ) malloc ( sample_num * sizeof ( int ) );
  x_sample = ( double * ) malloc ( elem_num * sample_num * sizeof ( double ) );
  printf ( "\n" );
  printf ( "  Number of samples is %d\n", sample_num );

  for ( sample_i = 0; sample_i < sample_num; sample_i++ )
  {
    x = dirichlet_mix_sample ( comp_max, comp_num, elem_num, alpha, 
      comp_weight, &seed, &comp_i );

    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      x_sample[elem_i+sample_i*elem_num] = x[elem_i];
    }

    comp_sample[sample_i] = comp_i;

    free ( x );
  }
/*
  Print some results.
*/
  printf ( "\n" );
  printf ( "  First few samples:\n" );
  printf ( "\n" );
  printf ( "  Sample  Component  X\n" );
  printf ( "\n" );

  for ( sample_i = 0; sample_i < i4_min ( sample_num, 10 ); sample_i++ )
  {
    printf ( "  %2d  %2d", sample_i, comp_sample[sample_i] );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      printf ( "  %10.6f", x_sample[elem_i+sample_i*elem_num] );
    }
    printf ( "\n" );
  }
  free ( comp_sample );
  free ( x_sample );
/*
  Compute the observed averages.
*/
  mean = r8col_mean ( sample_num, elem_num, x_sample );

  variance = r8col_variance ( sample_num, elem_num, x_sample );

  printf ( "\n" );
  printf ( "  Element  Observed mean, variance\n" );
  printf ( "\n" );
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    printf ( "%6d %10.6f%10.6f\n", elem_i, mean[elem_i], variance[elem_i] );
  }

  free ( mean );
  free ( variance );

  return;
# undef COMP_MAX
# undef COMP_NUM
# undef ELEM_NUM
# undef SAMPLE_NUM
}

