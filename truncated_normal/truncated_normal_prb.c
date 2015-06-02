# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "truncated_normal.h"

int main ( );

void i4_uniform_ab_test ( );

void normal_01_cdf_test ( );
void normal_01_cdf_inv_test ( );
void normal_01_mean_test ( );
void normal_01_moment_test ( );
void normal_01_pdf_test ( );
void normal_01_sample_test ( );
void normal_01_variance_test ( );

void normal_ms_cdf_test ( );
void normal_ms_cdf_inv_test ( );
void normal_ms_mean_test ( );
void normal_ms_moment_test ( );
void normal_ms_moment_central_test ( );
void normal_ms_pdf_test ( );
void normal_ms_sample_test ( );
void normal_ms_variance_test ( );

void r8_choose_test ( );
void r8_factorial2_test ( );
void r8_mop_test ( );
void r8_uniform_01_test ( );

void r8poly_print_test ( );
void r8poly_value_horner_test ( );

void r8vec_linspace_new_test ( );
void r8vec_print_test ( );

void truncated_normal_a_cdf_test ( );
void truncated_normal_a_cdf_inv_test ( );
void truncated_normal_a_mean_test ( );
void truncated_normal_a_moment_test ( );
void truncated_normal_a_pdf_test ( );
void truncated_normal_a_sample_test ( );
void truncated_normal_a_variance_test ( );

void truncated_normal_ab_cdf_test ( );
void truncated_normal_ab_cdf_inv_test ( );
void truncated_normal_ab_mean_test ( );
void truncated_normal_ab_moment_test ( );
void truncated_normal_ab_pdf_test ( );
void truncated_normal_ab_sample_test ( );
void truncated_normal_ab_variance_test ( );

void truncated_normal_b_cdf_test ( );
void truncated_normal_b_cdf_inv_test ( );
void truncated_normal_b_mean_test ( );
void truncated_normal_b_moment_test ( );
void truncated_normal_b_pdf_test ( );
void truncated_normal_b_sample_test ( );
void truncated_normal_b_variance_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRUNCATED_NORMAL_PRB.

  Discussion:

    TRUNCATED_NORMAL_PRB tests the TRUNCATED_NORMAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRUNCATED_NORMAL library.\n" );
/*
  Utilities.
*/
  i4_uniform_ab_test ( );

  r8_choose_test ( );
  r8_factorial2_test ( );
  r8_mop_test ( );
  r8_uniform_01_test ( );

  r8poly_print_test ( );
  r8poly_value_horner_test ( ); 

  r8vec_linspace_new_test ( );
  r8vec_print_test ( );
/*
  Library.
*/
  normal_01_cdf_test ( );
  normal_01_cdf_inv_test ( );
  normal_01_mean_test ( );
  normal_01_moment_test ( );
  normal_01_pdf_test ( );
  normal_01_sample_test ( );
  normal_01_variance_test ( );

  normal_ms_cdf_test ( );
  normal_ms_cdf_inv_test ( );
  normal_ms_mean_test ( );
  normal_ms_moment_test ( );
  normal_ms_moment_central_test ( );
  normal_ms_pdf_test ( );
  normal_ms_sample_test ( );
  normal_ms_variance_test ( );

  truncated_normal_a_cdf_test ( );
  truncated_normal_a_cdf_inv_test ( );
  truncated_normal_a_mean_test ( );
  truncated_normal_a_moment_test ( );
  truncated_normal_a_pdf_test ( );
  truncated_normal_a_sample_test ( );
  truncated_normal_a_variance_test ( );

  truncated_normal_ab_cdf_test ( );
  truncated_normal_ab_cdf_inv_test ( );
  truncated_normal_ab_mean_test ( );
  truncated_normal_ab_moment_test ( );
  truncated_normal_ab_pdf_test ( );
  truncated_normal_ab_sample_test ( );
  truncated_normal_ab_variance_test ( );

  truncated_normal_b_cdf_test ( );
  truncated_normal_b_cdf_inv_test ( );
  truncated_normal_b_mean_test ( );
  truncated_normal_b_moment_test ( );
  truncated_normal_b_pdf_test ( );
  truncated_normal_b_sample_test ( );
  truncated_normal_b_variance_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void i4_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_TEST tests I4_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "I4_UNIFORM_TEST\n" );
  printf ( "  I4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, &seed );
    printf ( "  %8d  %d\n", i, j );
  }

  return;
}
/******************************************************************************/

void normal_01_cdf_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_CDF_TEST tests NORMAL_01_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 February 2015

  Author:

    John Burkardt
*/
{
  double cdf1;
  double cdf2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_01_CDF_TEST\n" );
  printf ( "  NORMAL_01_CDF evaluates the Normal 01 CDF;\n" );
  printf ( "\n" );
  printf ( "       X              CDF                       CDF\n" );
  printf ( "                     (exact)                   (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &cdf1 );

    if ( n_data == 0 )
    {
      break;
    }

    cdf2 = normal_01_cdf ( x );

    printf ( "  %14.6g  %24.16g  %24.16g\n", x, cdf1, cdf2 );
  }

  return;
}
/******************************************************************************/

void normal_01_cdf_inv_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_CDF_INV_TEST tests NORMAL_01_CDF_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2015

  Author:

    John Burkardt
*/
{
  double cdf;
  int n_data;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "NORMAL_01_CDF_INV_TEST\n" );
  printf ( "  NORMAL_01_CDF_INV inverts the Normal 01 CDF;\n" );
  printf ( "\n" );
  printf ( "      CDF             X                         X\n" );
  printf ( "                     (exact)                   (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x1, &cdf );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = normal_01_cdf_inv ( cdf );

    printf ( "  %14.6g  %24.16g  %24.16g\n", cdf, x1, x2 );
  }

  return;
}
/******************************************************************************/

void normal_01_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_MEAN_TEST tests NORMAL_01_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double mean;
  int sample_num;
  int seed = 123456789;
  double *x;
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "NORMAL_01_MEAN_TEST\n" );
  printf ( "  NORMAL_01_MEAN computes the Normal 01 mean.\n" );

  mean = normal_01_mean ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );

  sample_num = 1000;
  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_01_sample ( &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void normal_01_moment_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_MOMENT_TEST tests NORMAL_01_MOMENT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2015

  Author:

    John Burkardt
*/
{
  double moment;
  int order;

  printf ( "\n" );
  printf ( "NORMAL_01_MOMENT_TEST\n" );
  printf ( "  NORMAL_01_MOMENT evaluates Normal 01 moments;\n" );
  printf ( "\n" );
  printf ( "   Order    Moment\n" );
  printf ( "\n" );

  for ( order = 0; order <= +10; order++ )
  {
    moment = normal_01_moment ( order );
    printf ( "  %6d  %24.16g\n", order, moment );
  }

  return;
}
/******************************************************************************/

void normal_01_pdf_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_PDF_TEST tests NORMAL_01_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double pdf;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_01_PDF_TEST\n" );
  printf ( "  NORMAL_01_PDF evaluates the Normal 01 PDF;\n" );
  printf ( "\n" );
  printf ( "       X              PDF\n" );
  printf ( "\n" );

  for ( i = - 20; i <= +20; i++ )
  {
    x = ( double ) ( i ) / 10.0;
    pdf = normal_01_pdf ( x );
    printf ( "  %14.6g  %24.16g\n", x, pdf );
  }

  return;
}
/******************************************************************************/

void normal_01_sample_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_SAMPLE_TEST tests NORMAL_01_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 February 2015

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_01_SAMPLE_TEST\n" );
  printf ( "  NORMAL_01_SAMPLE returns samples from the normal\n" );
  printf ( "  distribution with mean 0 and standard deviation 1.\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_01_sample ( &seed );
    printf ( "  %4d  %14.6g\n", i, x );
  }

  return;
}
/******************************************************************************/

void normal_01_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_VARIANCE_TEST tests NORMAL_01_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  int sample_num;
  int seed = 123456789;
  double variance;
  double *x;

  printf ( "\n" );
  printf ( "NORMAL_01_VARIANCE_TEST\n" );
  printf ( "  NORMAL_01_VARIANCE computes the Normal 01 variance;\n" );

  variance = normal_01_variance ( );

  printf ( "\n" );
  printf ( "  PDF variance = %g\n", variance );

  sample_num = 1000;
  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_01_sample ( &seed );
  }

  variance = r8vec_variance ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample variance = %g\n", variance );

  free ( x );

  return;
}
/******************************************************************************/

void normal_ms_cdf_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_CDF_TEST tests NORMAL_MS_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double mu;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_MS_CDF_TEST\n" );
  printf ( "  NORMAL_MS_CDF evaluates the Normal MS CDF;\n" );

  mu = 100.0;
  sigma = 15.0;

  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );
  printf ( "\n" );
  printf ( "       X              CDF\n" );
  printf ( "\n" );

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    cdf = normal_ms_cdf ( x, mu, sigma );
    printf ( "  %14.6g  %24.16g\n", x, cdf );
  }

  return;
}
/******************************************************************************/

void normal_ms_cdf_inv_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_CDF_INV_TEST tests NORMAL_MS_CDF_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double mu;
  double sigma;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "NORMAL_MS_CDF_INV_TEST\n" );
  printf ( "  NORMAL_MS_CDF_INV inverts the Normal MS CDF;\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  printf ( "\n" );
  printf ( "       X            CDF           CDF_INV\n" );
  printf ( "\n" );

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    cdf = normal_ms_cdf ( x, mu, sigma );
    x2 = normal_ms_cdf_inv ( cdf, mu, sigma );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void normal_ms_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_MEAN_TEST tests NORMAL_MS_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double mean;
  double mu;
  int sample_num;
  int seed = 123456789;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "NORMAL_MS_MEAN_TEST\n" );
  printf ( "  NORMAL_MS_MEAN computes the Normal MS mean.\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  mean = normal_ms_mean ( mu, sigma );

  printf ( "\n" );
  printf ( "  PDF mean = %g\n", mean );

  sample_num = 1000;
  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_ms_sample ( mu, sigma, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void normal_ms_moment_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_MOMENT_TEST tests NORMAL_MS_MOMENT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  double moment;
  double mu;
  int order;
  double sigma;

  printf ( "\n" );
  printf ( "NORMAL_MS_MOMENT_TEST\n" );
  printf ( "  NORMAL_MS_MOMENT evaluates Normal MS moments;\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  printf ( "\n" );
  printf ( "   Order    Moment\n" );
  printf ( "\n" );

  for ( order = 0; order <= +10; order++ )
  {
    moment = normal_ms_moment ( order, mu, sigma );
    printf ( "  %6d  %24.16g\n", order, moment );
  }

  return;
}
/******************************************************************************/

void normal_ms_moment_central_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_MOMENT_CENTRAL_TEST tests NORMAL_MS_MOMENT_CENTRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  double moment;
  double mu;
  int order;
  double sigma;

  printf ( "\n" );
  printf ( "NORMAL_MS_MOMENT_CENTRAL_TEST\n" );
  printf ( "  NORMAL_MS_MOMENT_CENTRAL evaluates Normal MS central moments;\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  printf ( "\n" );
  printf ( "   Order    Moment\n" );
  printf ( "\n" );

  for ( order = 0; order <= +10; order++ )
  {
    moment = normal_ms_moment_central ( order, mu, sigma );
    printf ( "  %6d  %24.16g\n", order, moment );
  }

  return;
}
/******************************************************************************/

void normal_ms_pdf_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_PDF_TEST tests NORMAL_MS_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double mu;
  double pdf;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_MS_PDF_TEST\n" );
  printf ( "  NORMAL_MS_PDF evaluates the Normal MS PDF;\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  printf ( "\n" );
  printf ( "       X              PDF\n" );
  printf ( "\n" );

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    pdf = normal_ms_pdf ( mu, sigma, x );
    printf ( "  %14.6g  %24.16g\n", x, pdf );
  }

  return;
}
/******************************************************************************/

void normal_ms_sample_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_SAMPLE_TEST tests NORMAL_MS_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double mu;
  int seed;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_MS_SAMPLE_TEST\n" );
  printf ( "  NORMAL_MS_SAMPLE returns samples from the Normal MS PDF.\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_ms_sample ( mu, sigma, &seed );
    printf ( "  %4d  %14.6g\n", i, x );
  }

  return;
}
/******************************************************************************/

void normal_ms_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_MS_VARIANCE_TEST tests NORMAL_MS_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double mu;
  int sample_num;
  int seed = 123456789;
  double sigma;
  double variance;
  double *x;

  printf ( "\n" );
  printf ( "NORMAL_MS_VARIANCE_TEST\n" );
  printf ( "  NORMAL_MS_VARIANCE computes the Normal MS variance;\n" );

  mu = 100.0;
  sigma = 15.0;
  printf ( "\n" );
  printf ( "  Parameter MU = %g\n", mu );
  printf ( "  Parameteter SIGMA = %g\n", sigma );

  variance = normal_ms_variance ( mu, sigma );

  printf ( "\n" );
  printf ( "  PDF variance = %g\n", variance );

  sample_num = 1000;
  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_ms_sample ( mu, sigma, &seed );
  }

  variance = r8vec_variance ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample variance = %g\n", variance );

  free ( x );

  return;
}
/******************************************************************************/

void r8_choose_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE_TEST tests R8_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "R8_CHOOSE_TEST\n" );
  printf ( "  R8_CHOOSE evaluates C(N,K).\n" );
  printf ( "\n" );
  printf ( "         N         K       CNK\n" );
 
  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );
      printf ( "  %8d  %8d  %14.6g\n", n, k, cnk );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_factorial2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL2_TEST tests R8_FACTORIAL2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2015

  Author:

    John Burkardt
*/
{
  double f1;
  double f2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL2_TEST\n" );
  printf ( "  R8_FACTORIAL2 evaluates the double factorial.\n" );
  printf ( "\n" );
  printf ( "    N                Exact                  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( &n_data, &n, &f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial2 ( n );

    printf ( "  %4d  %24.16g  %24.16g\n", n, f1, f2 );

  }
 
  return;
}
/******************************************************************************/

void r8_mop_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MOP_TEST tests R8_MOP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2014

  Author:

    John Burkardt
*/
{
  int i4;
  int i4_max;
  int i4_min;
  double r8;
  int seed = 123456789;
  int test;

  printf ( "\n" );
  printf ( "R8_MOP_TEST\n" );
  printf ( "  R8_MOP evaluates (-1.0)^I4 as an R8.\n" );
  printf ( "\n" );
  printf ( "    I4  R8_MOP(I4)\n" );
  printf ( "\n" );

  i4_min = -100;
  i4_max = +100;

  for ( test = 1; test <= 10; test++ )
  {
    i4 = i4_uniform_ab ( i4_min, i4_max, &seed );
    r8 = r8_mop ( i4 );
    printf ( "  %4d  %4.1f\n", i4, r8 );
  }

  return;
}
/******************************************************************************/

void r8_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
# define N 1000

  int i;
  double max;
  double mean;
  double min;
  int n;
  int seed = 123456789;
  double x[N];
  double variance;

  printf ( "\n" );
  printf ( "R8_UNIFORM_01_TEST\n" );
  printf ( "  R8_UNIFORM_01 samples a uniform random distribution in [0,1].\n" );
  printf ( "  distributed random numbers.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  for ( i = 0; i < N; i++ )
  {
    x[i] = r8_uniform_01 ( &seed );
  }

  printf ( "\n" );
  printf ( "  First few values:\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, x[i] );
  }
  min = r8vec_min ( N, x );
  max = r8vec_max ( N, x );
  mean = r8vec_mean ( N, x );
  variance = r8vec_variance ( N, x );

  printf ( "\n" );
  printf ( "  Number of samples was %d\n", N );
  printf ( "  Minimum value was %f\n", min );
  printf ( "  Maximum value was %f\n", max );
  printf ( "  Average value was %f\n", mean );
  printf ( "  Variance was      %f\n", variance );

  return;
# undef N
}
/******************************************************************************/

void r8poly_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_PRINT_TEST tests R8POLY_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 January 2015

  Author:

    John Burkardt
*/
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  printf ( "\n" );
  printf ( "R8POLY_PRINT_TEST\n" );
  printf ( "  R8POLY_PRINT prints an R8POLY.\n" );

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
/******************************************************************************/

void r8poly_value_horner_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2015

  Author:

    John Burkardt
*/
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  printf ( "\n" );
  printf ( "R8POLY_VALUE_HORNER_TEST\n" );
  printf ( "  R8POLY_VALUE_HORNER evaluates a polynomial at\n" );
  printf ( "  one point, using Horner's method.\n" );

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  printf ( "\n" );
  printf ( "   I    X    P(X)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    printf ( "  %2d  %8.4f  %14.6g\n", i, x[i], p );
  }

  free ( x );

  return;
}
/******************************************************************************/

void r8vec_linspace_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW_TEST tests R8VEC_LINSPACE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int n = 5;
  double *x;

  printf ( "\n" );
  printf ( "R8VEC_LINSPACE_NEW_TEST\n" );
  printf ( "  For a R8VEC:\n" );
  printf ( "  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;\n" );

  a = 10.0;
  b = 20.0;

  x = r8vec_linspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_linspace ( 5, 10, 20 )" );
  free ( x );

  return;
}
/******************************************************************************/

void r8vec_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_TEST tests R8VEC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  printf ( "\n" );
  printf ( "R8VEC_PRINT_TEST\n" );
  printf ( "  R8VEC_PRINT prints an R8VEC.\n" );

  r8vec_print ( n, a, "  The R8VEC:" );

  return;
}
/******************************************************************************/

void truncated_normal_a_cdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_A_CDF_TEST tests TRUNCATED_NORMAL_A_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_CDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_A_CDF evaluates\n" );
  printf ( "  the lower Truncated Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU       S         A         X        CDF1           CDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_cdf_values ( &n_data, &mu, &sigma, &a, &x, &cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_a_cdf ( x, mu, sigma, a );

    printf ( "  %8.1f  %8.1f  %8.1f %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, a, x, cdf1, cdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_a_cdf_inv_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_CDF_INV_TEST tests TRUNCATED_NORMAL_A_CDF_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;
  double x2;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_CDF_INV_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_A_CDF_INV inverts the CDF of\n" );
  printf ( "  the lower Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo)\n", a );
  printf ( "\n" );
  printf ( "       X            CDF           CDF_INV\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_a_sample ( mu, sigma, a, &seed );
    cdf = truncated_normal_a_cdf ( x, mu, sigma, a );
    x2 = truncated_normal_a_cdf_inv ( cdf, mu, sigma, a );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void truncated_normal_a_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_MEAN_TEST tests TRUNCATED_NORMAL_A_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_MEAN_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_A_MEAN computes the mean\n" );
  printf ( "  of the lower Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo)\n", a );

  mean = truncated_normal_a_mean ( mu, sigma, a );

  printf ( "\n" );
  printf ( "  PDF mean = %g\n", mean );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_a_sample ( mu, sigma, a, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void truncated_normal_a_moment_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_MOMENT_TEST tests TRUNCATED_NORMAL_A_MOMENT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double a_test[6] = {
    0.0, -10.0, 10.0, -10.0, +10.0, -10.0 };
  double moment;
  double mu;
  double mu_test[6] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  -5.0 };
  int order;
  double sigma;
  double sigma_test[6] = {
    1.0,  1.0,  1.0,  2.0,  2.0,  1.0 };
  int test;
  int test_num;

  test_num = 6;
 
  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_MOMENT_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_A_MOMENT evaluates the moments\n" );
  printf ( "  of the Lower Truncated Normal PDF:\n" );

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    a = a_test[test];
    printf ( "\n" );
    printf ( "  Test = %d, Mu = %g, Sigma = %g, A = %g\n", test, mu, sigma, a );
    printf ( " Order  Moment\n" );
    printf ( "\n" );
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_a_moment ( order, mu, sigma, a );
      printf ( "  %2d  %14.6g\n", order, moment );
    }
  }
  return;
}
/******************************************************************************/

void truncated_normal_a_pdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_A_PDF_TEST tests TRUNCATED_NORMAL_A_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_PDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_A_PDF evaluates\n" );
  printf ( "  the lower Truncated Normal Probability Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU       S         A         X        PDF1        PDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_pdf_values ( &n_data, &mu, &sigma, &a, &x, &pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_a_pdf ( x, mu, sigma, a );

    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, a, x, pdf1, pdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_a_sample_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_SAMPLE_TEST tests TRUNCATED_NORMAL_A_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_SAMPLE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_A_SAMPLE samples;\n" );
  printf ( "  the lower Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo)\n", a );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_a_sample ( mu, sigma, a, &seed );
    printf ( "  %4d  %14.6g\n", i, x );
  }

  return;
}
/******************************************************************************/

void truncated_normal_a_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_VARIANCE_TEST tests TRUNCATED_NORMAL_A_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;
 
  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_VARIANCE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_A_VARIANCE computes the variance\n" );
  printf ( "  of the lower Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo)\n", a );

  variance = truncated_normal_a_variance ( mu, sigma, a );

  printf ( "\n" );
  printf ( "  PDF variance = %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_a_sample ( mu, sigma, a, &seed );
  }

  variance = r8vec_variance ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample variance = %g\n", variance );

  free ( x );

  return;
}
/******************************************************************************/

void truncated_normal_ab_cdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_AB_CDF_TEST tests TRUNCATED_NORMAL_AB_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_CDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_AB_CDF evaluates\n" );
  printf ( "  the Truncated Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU       S         A         B         X        CDF1           CDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_cdf_values ( &n_data, &mu, &sigma, &a, &b, &x, &cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_ab_cdf ( x, mu, sigma, a, b );

    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, a, b, x, cdf1, cdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_ab_cdf_inv_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_CDF_INV_TEST tests TRUNCATED_NORMAL_AB_CDF_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;
  double x2;

  a = 50.0;
  b = 150;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_CDF_INV_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_AB_CDF_INV inverts the CDF of\n" );
  printf ( "  the Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g)\n", a, b );
  printf ( "\n" );
  printf ( "       X            CDF           CDF_INV\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_ab_sample ( mu, sigma, a, b, &seed );
    cdf = truncated_normal_ab_cdf ( x, mu, sigma, a, b );
    x2 = truncated_normal_ab_cdf_inv ( cdf, mu, sigma, a, b );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void truncated_normal_ab_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_MEAN_TEST tests TRUNCATED_NORMAL_AB_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_MEAN_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_AB_MEAN computes the mean\n" );
  printf ( "  of the Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g]\n", a, b );

  mean = truncated_normal_ab_mean ( mu, sigma, a, b );

  printf ( "\n" );
  printf ( "  PDF mean = %g\n", mean );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_ab_sample ( mu, sigma, a, b, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void truncated_normal_ab_moment_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_MOMENT_TEST tests TRUNCATED_NORMAL_AB_MOMENT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double a_test[9] = {
    -1.0, 0.0, -1.0, -1.0,  0.0, 0.5, -2.0, -4.0, 4.0 };
  double b;
  double b_test[9] = {
    1.0, 1.0,  0.0,  1.0,  2.0, 2.0,  2.0,  4.0, 7.0 };
  double moment;
  double mu;
  double mu_test[9] = {
    0.0, 0.0,  0.0,  0.0,  1.0, 0.0,  0.0,  0.0, 5.0 };
  int order;
  double sigma;
  double sigma_test[9] = {
    1.0, 1.0,  1.0,  2.0,  1.0, 1.0,  1.0,  1.0, 0.5 };
  int test;
  int test_num;

  test_num = 9;
 
  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_MOMENT_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_AB_MOMENT evaluates the moments\n" );
  printf ( "  of the Truncated Normal PDF:\n" );

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    a = a_test[test];
    b = b_test[test];
    printf ( "\n" );
    printf ( "  Test = %d, Mu = %g, Sigma = %g, A = %g, B = %g\n", test, mu, sigma, a, b );
    printf ( " Order  Moment\n" );
    printf ( "\n" );
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_ab_moment ( order, mu, sigma, a, b );
      printf ( "  %2d  %14.6g\n", order, moment );
    }
  }
  return;
}
/******************************************************************************/

void truncated_normal_ab_pdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_AB_PDF_TEST tests TRUNCATED_NORMAL_AB_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_PDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_AB_PDF evaluates the PDF of\n" );
  printf ( "  the Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "        MU       S         A         B         X        PDF1        PDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_pdf_values ( &n_data, &mu, &sigma, &a, &b, &x, &pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_ab_pdf ( x, mu, sigma, a, b );

    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, a, b, x, pdf1, pdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_ab_sample_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_SAMPLE_TEST tests TRUNCATED_NORMAL_AB_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_SAMPLE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_AB_SAMPLE samples;\n" );
  printf ( "  the Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g]\n", a, b );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_ab_sample ( mu, sigma, a, b, &seed );
    printf ( "  %4d  %14.6g\n", i, x );
  }

  return;
}
/******************************************************************************/

void truncated_normal_ab_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_VARIANCE_TEST tests TRUNCATED_NORMAL_AB_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;
 
  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_VARIANCE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_AB_VARIANCE computes the variance\n" );
  printf ( "  of the Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g]\n", a, b );

  variance = truncated_normal_ab_variance ( mu, sigma, a, b );

  printf ( "\n" );
  printf ( "  PDF variance = %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_ab_sample ( mu, sigma, a, b, &seed );
  }

  variance = r8vec_variance ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample variance = %g\n", variance );

  free ( x );

  return;
}
/******************************************************************************/

void truncated_normal_b_cdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_B_CDF_TEST tests TRUNCATED_NORMAL_B_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_CDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_B_CDF_TEST evaluates the CDF of\n" );
  printf ( "  the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "        MU       S         B         X        CDF1           CDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_cdf_values ( &n_data, &mu, &sigma, &b, &x, &cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_b_cdf ( x, mu, sigma, b );

    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, b, x, cdf1, cdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_b_cdf_inv_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_CDF_INV_TEST tests TRUNCATED_NORMAL_B_CDF_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double b;
  double cdf;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;
  double x2;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_CDF_INV_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_B_CDF_INV inverts the CDF of\n" );
  printf ( "  the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval (-oo,%g]\n", b );
  printf ( "\n" );
  printf ( "       X            CDF           CDF_INV\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_b_sample ( mu, sigma, b, &seed );
    cdf = truncated_normal_b_cdf ( x, mu, sigma, b );
    x2 = truncated_normal_b_cdf_inv ( cdf, mu, sigma, b );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void truncated_normal_b_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_MEAN_TEST tests TRUNCATED_NORMAL_B_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double b;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_MEAN_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_B_MEAN computes the mean\n" );
  printf ( "  of the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval (-oo,%g]\n", b );

  mean = truncated_normal_b_mean ( mu, sigma, b );

  printf ( "\n" );
  printf ( "  PDF mean = %g\n", mean );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_b_sample ( mu, sigma, b, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void truncated_normal_b_moment_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_MOMENT_TEST tests TRUNCATED_NORMAL_B_MOMENT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double b_test[6] = {
    0.0, 10.0, -10.0, 10.0, -10.0, 10.0 };
  double moment;
  double mu;
  double mu_test[6] = {
    0.0,  0.0,  0.0,  0.0,  0.0,  5.0 };
  int order;
  double sigma;
  double sigma_test[6] = {
    1.0,  1.0,  1.0,  2.0,  2.0,  1.0 };
  int test;
  int test_num;

  test_num = 6;
 
  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_MOMENT_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_B_MOMENT evaluates the moments\n" );
  printf ( "  of the Upper Truncated Normal PDF:\n" );

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    b = b_test[test];
    printf ( "\n" );
    printf ( "  Test = %d, Mu = %g, Sigma = %g, B = %g\n", test, mu, sigma, b );
    printf ( " Order  Moment\n" );
    printf ( "\n" );
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_b_moment ( order, mu, sigma, b );
      printf ( "  %2d  %14.6g\n", order, moment );
    }
  }
  return;
}
/******************************************************************************/

void truncated_normal_b_pdf_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_B_PDF_TEST tests TRUNCATED_NORMAL_B_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_PDF_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_B_PDF evaluates the PDF of\n" );
  printf ( "  the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "        MU       S         B         X        PDF1        PDF2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_pdf_values ( &n_data, &mu, &sigma, &b, &x, &pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_b_pdf ( x, mu, sigma, b );

    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16g  %24.16g\n", 
      mu, sigma, b, x, pdf1, pdf2 );
  }
  return;
}
/******************************************************************************/

void truncated_normal_b_sample_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_SAMPLE_TEST tests TRUNCATED_NORMAL_B_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double b;
  int i;
  double mu;
  int seed;
  double sigma;
  double x;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_SAMPLE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_B_SAMPLE samples;\n" );
  printf ( "  the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval (-oo,%g]\n", b );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = truncated_normal_b_sample ( mu, sigma, b, &seed );
    printf ( "  %4d  %14.6g\n", i, x );
  }

  return;
}
/******************************************************************************/

void truncated_normal_b_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_VARIANCE_TEST tests TRUNCATED_NORMAL_B_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double b;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;
 
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_VARIANCE_TEST\n" );
  printf ( "  TRUNCATED_NORMAL_B_VARIANCE computes the variance\n" );
  printf ( "  of the upper Truncated Normal Distribution.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", sigma );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval (-oo,%g]\n", b );

  variance = truncated_normal_b_variance ( mu, sigma, b );

  printf ( "\n" );
  printf ( "  PDF variance = %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_b_sample ( mu, sigma, b, &seed );
  }

  variance = r8vec_variance ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample variance = %g\n", variance );

  free ( x );

  return;
}
