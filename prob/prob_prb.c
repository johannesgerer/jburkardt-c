# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "prob.h"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );

void test010 ( );
void test0105 ( );
void test0106 ( );
void test011 ( );
void test012 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test016 ( );

void test020 ( );
void test021 ( );
void test022 ( );
void test023 ( );
void test0235 ( );
void test024 ( );
void test025 ( );
void test0251 ( );
void test0252 ( );
void test0253 ( );
void test0254 ( );
void test026 ( );
void test027 ( );
void test0275 ( );
void test0276 ( );
void test028 ( );
void test029 ( );

void test030 ( );
void test031 ( );
void test032 ( );
void test033 ( );
void test034 ( );
void test035 ( );
void test036 ( );
void test037 ( );
void test0375 ( );
void test038 ( );
void test039 ( );
void test0395 ( );

void test040 ( );
void test041 ( );
void test042 ( );
void test043 ( );
void test044 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test048 ( );
void test049 ( );

void test050 ( );
void test051 ( );
void test052 ( );
void test053 ( );
void test054 ( );
void test055 ( );
void test056 ( );
void test0563 ( );
void test0564 ( );
void test0565 ( );
void test0566 ( );
void test057 ( );
void test058 ( );
void test059 ( );

void test060 ( );
void test061 ( );
void test062 ( );
void test063 ( );
void test064 ( );
void test065 ( );
void test066 ( );
void test067 ( );
void test068 ( );
void test069 ( );

void test070 ( );
void test0705 ( );
void test071 ( );
void test072 ( );
void test073 ( );
void test074 ( );
void test0744 ( );
void test0745 ( );
void test075 ( );
void test076 ( );
void test077 ( );
void test078 ( );
void test079 ( );

void test080 ( );
void test081 ( );
void test082 ( );
void test083 ( );
void test084 ( );
void test085 ( );
void test086 ( );
void test087 ( );
void test088 ( );
void test089 ( );

void test090 ( );
void test091 ( );
void test092 ( );
void test093 ( );
void test094 ( );
void test095 ( );
void test096 ( );
void test0965 ( );
void test097 ( );
void test098 ( );
void test099 ( );

void test100 ( );
void test101 ( );
void test102 ( );
void test103 ( );
void test104 ( );
void test105 ( );
void test106 ( );
void test107 ( );
void test108 ( );
void test109 ( );

void test110 ( );
void test111 ( );
void test112 ( );
void test113 ( );
void test114 ( );
void test1145 ( );
void test1146 ( );
void test115 ( );
void test116 ( );
void test117 ( );
void test118 ( );
void test1184 ( );
void test1185 ( );
void test1186 ( );
void test1187 ( );
void test1188 ( );
void test1189 ( );
void test119 ( );

void test120 ( );
void test123 ( );
void test124 ( );
void test125 ( );
void test126 ( );
void test127 ( );
void test128 ( );
void test129 ( );

void test130 ( );
void test1304 ( );
void test1306 ( );
void test131 ( );
void test132 ( );
void test133 ( );
void test134 ( );
void test1341 ( );
void test1342 ( );
void test1344 ( );
void test135 ( );
void test136 ( );
void test137 ( );
void test138 ( );
void test139 ( );

void test140 ( );
void test141 ( );
void test142 ( );
void test1425 ( );
void test143 ( );
void test144 ( );
void test145 ( );
void test146 ( );
void test147 ( );
void test148 ( );
void test1485 ( );
void test1486 ( );
void test149 ( );

void test150 ( );
void test151 ( );
void test152 ( );
void test153 ( );
void test154 ( );
void test155 ( );
void test1555 ( );
void test156 ( );
void test157 ( );
void test158 ( );
void test159 ( );

void test160 ( );
void test161 ( );
void test162 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PROB_PRB.

  Discussion:

    PROB_PRB tests the PROB library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PROB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PROB library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );

  test010 ( );
  test0105 ( );
  test0106 ( );
  test011 ( );
  test012 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test016 ( );

  test020 ( );
  test021 ( );
  test022 ( );
  test023 ( );
  test0235 ( );
  test024 ( );
  test025 ( );
  test0251 ( );
  test0252 ( );
  test0253 ( );
  test0254 ( );
  test026 ( );
  test027 ( );
  test0275 ( );
  test0276 ( );
  test028 ( );
  test029 ( );

  test030 ( );
  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test035 ( );
  test036 ( );
  test037 ( );
  test0375 ( );
  test038 ( );
  test039 ( );
  test0395 ( );

  test040 ( );
  test041 ( );
  test042 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test051 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test056 ( );
  test0563 ( );
  test0564 ( );
  test0565 ( );
  test0566 ( );
  test057 ( );
  test058 ( );
  test059 ( );

  test060 ( );
  test061 ( );
  test062 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test067 ( );
  test068 ( );
  test069 ( );

  test070 ( );
  test0705 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test0744 ( );
  test0745 ( );
  test075 ( );
  test076 ( );
  test077 ( );
  test078 ( );
  test079 ( );

  test080 ( );
  test081 ( );
  test082 ( );
  test083 ( );
  test084 ( );
  test085 ( );
  test086 ( );
  test087 ( );
  test088 ( );
  test089 ( );

  test090 ( );
  test091 ( );
  test092 ( );
  test093 ( );
  test094 ( );
  test095 ( );
  test096 ( );
  test0965 ( );
  test097 ( );
  test098 ( );
  test099 ( );

  test100 ( );
  test101 ( );
  test102 ( );
  test103 ( );
  test104 ( );
  test105 ( );
  test106 ( );
  test107 ( );
  test108 ( );
  test109 ( );

  test110 ( );
  test111 ( );
  test112 ( );
  test113 ( );
  test114 ( );
  test1145 ( );
  test1146 ( );
  test115 ( );
  test116 ( );
  test117 ( );
  test118 ( );
  test1184 ( );
  test1185 ( );
  test1186 ( );
  test1187 ( );
  test1188 ( );
  test1189 ( );
  test119 ( );

  test120 ( );
  test123 ( );
  test124 ( );
  test125 ( );
  test126 ( );
  test127 ( );
  test128 ( );
  test129 ( );

  test130 ( );
  test1304 ( );
  test1306 ( );
  test131 ( );
  test132 ( );
  test133 ( );
  test134 ( );
  test1341 ( );
  test1342 ( );
  test1344 ( );
  test135 ( );
  test136 ( );
  test137 ( );
  test138 ( );
  test139 ( );

  test140 ( );
  test141 ( );
  test142 ( );
  test1425 ( );
  test143 ( );
  test144 ( );
  test145 ( );
  test146 ( );
  test147 ( );
  test148 ( );
  test1485 ( );
  test1486 ( );
  test149 ( );

  test150 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test154 ( );
  test155 ( );
  test1555 ( );
  test156 ( );
  test157 ( );
  test158 ( );
  test159 ( );

  test160 ( );
  test161 ( );
  test162 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PROB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests ANGLE_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int n;
  double x;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  For the ANGLE PDF:\n" );
  printf ( "  ANGLE_CDF evaluates the CDF;\n" );

  n = 5;
  x = 0.50E+00;

  cdf = angle_cdf ( x, n );

  printf ( "\n" );
  printf ( "  Parameter N =     %d\n", n );
  printf ( "  PDF argument X =  %g\n", x );
  printf ( "  CDF value =       %g\n", cdf );

  return;
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests ANGLE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 June 2013

  Author:

    John Burkardt
*/
{
  int n;
  double pdf;
  double x;

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  For the ANGLE PDF:\n" );
  printf ( "  ANGLE_PDF evaluates the PDF;\n" );

  n = 5;
  x = 0.50E+00;

  pdf = angle_pdf ( x, n );

  printf ( "\n" );
  printf ( "  Parameter N =    %d\n", n );
  printf ( "  PDF argument X = %g\n", x );
  printf ( "  PDF value =      %g\n", pdf );

  return;
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests ANGLE_MEAN;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 June 2013

  Author:

    John Burkardt
*/
{
  double mean;
  int n;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  For the ANGLE PDF:\n" );
  printf ( "  ANGLIT_MEAN computes the mean;\n" );

  n = 5;
  mean = angle_mean ( n );

  printf ( "\n" );
  printf ( "  Parameter N = %d\n", n );
  printf ( "  PDF mean =    %g\n", mean );

  return;
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
/*
  Purpose:

    TEST004 tests ANGLIT_CDF, ANGLIT_CDF_INV, ANGLIT_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST004\n" );
  printf ( "  For the Anglit PDF:\n" );
  printf ( "  ANGLIT_CDF evaluates the CDF;\n" );
  printf ( "  ANGLIT_CDF_INV inverts the CDF.\n" );
  printf ( "  ANGLIT_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = anglit_sample ( &seed );
    pdf = anglit_pdf ( x );
    cdf = anglit_cdf ( x );
    x2 = anglit_cdf_inv ( cdf );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests ANGLIT_MEAN, ANGLIT_SAMPLE, ANGLIT_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  For the Anglit PDF:\n" );
  printf ( "  ANGLIT_MEAN computes the mean;\n" );
  printf ( "  ANGLIT_SAMPLE samples;\n" );
  printf ( "  ANGLIT_VARIANCE computes the variance.\n" );

  mean     = anglit_mean ( );
  variance = anglit_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = anglit_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests ARCSIN_CDF, ARCSIN_CDF_INV, ARCSIN_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  For the Arcsin PDF:\n" );
  printf ( "  ARCSIN_CDF evaluates the CDF;\n" );
  printf ( "  ARCSIN_CDF_INV inverts the CDF.\n" );
  printf ( "  ARCSIN_PDF evaluates the PDF;\n" );

  a = 1.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !arcsin_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST006 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = arcsin_sample ( a, &seed );
    pdf = arcsin_pdf ( x, a );
    cdf = arcsin_cdf ( x, a );
    x2 = arcsin_cdf_inv ( cdf, a );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test007 ( )

/******************************************************************************/
/*
  Purpose:

    TEST007 tests ARCSIN_MEAN, ARCSIN_SAMPLE, ARCSIN_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST007\n" );
  printf ( "  For the Arcsin PDF:\n" );
  printf ( "  ARCSIN_MEAN computes the mean;\n" );
  printf ( "  ARCSIN_SAMPLE samples;\n" );
  printf ( "  ARCSIN_VARIANCE computes the variance.\n" );

  for ( i = 1; i <= 2; i++ )
  {
    if ( i == 1 )
    {
      a = 1.0;
    }
    else if ( i == 2 )
    {
      a = 16.0;
    }

    printf ( "\n" );
    printf ( "  PDF parameter A =             %g\n", a );

    if ( !arcsin_check ( a ) )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "TEST007 - Fatal error!\n" );
      fprintf ( stderr, "  The parameters are not legal.\n" );
      return;
    }

    mean = arcsin_mean ( a );
    variance = arcsin_variance ( a );

    printf ( "  PDF mean =                    %g\n", mean );
    printf ( "  PDF variance =                %g\n", variance );

    for ( j = 0; j < SAMPLE_NUM; j++ )
    {
      x[j] = arcsin_sample ( a, &seed );
    }

    mean = r8vec_mean ( SAMPLE_NUM, x );
    variance = r8vec_variance ( SAMPLE_NUM, x );
    xmax = r8vec_max ( SAMPLE_NUM, x );
    xmin = r8vec_min ( SAMPLE_NUM, x );

    printf ( "\n" );
    printf ( "  Sample size =     %d\n", SAMPLE_NUM );
    printf ( "  Sample mean =     %g\n", mean );
    printf ( "  Sample variance = %g\n", variance );
    printf ( "  Sample maximum =  %g\n", xmax );
    printf ( "  Sample minimum =  %g\n", xmin );
  }

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test008 ( )

/******************************************************************************/
/*
  Purpose:

    TEST008 tests BENFORD_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  int n;
  double pdf;

  printf ( "\n" );
  printf ( "TEST008\n" );
  printf ( "  For the Benford PDF:\n" );
  printf ( "  BENFORD_PDF evaluates the PDF.\n" );

  printf ( "\n" );
  printf ( "  N    PDF(N)\n" );
  printf ( "\n" );

  for ( n = 1; n <= 19; n++ )
  {
    pdf = benford_pdf ( n );
    printf ( "  %6d  %14g\n", n, pdf );
  }

  return;
}
/******************************************************************************/

void test009 ( )

/******************************************************************************/
/*
  Purpose:

    TEST009 tests BERNOULLI_CDF, BERNOULLI_CDF_INV, BERNOULLI_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST009\n" );
  printf ( "  For the Bernoulli PDF,\n" );
  printf ( "  BERNOULLI_CDF evaluates the CDF;\n" );
  printf ( "  BERNOULLI_CDF_INV inverts the CDF;\n" );
  printf ( "  BERNOULLI_PDF evaluates the PDF.\n" );

  a = 0.75E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !bernoulli_check ( a ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TEST009 - Fatal error!\n" );
    fprintf ( stderr, "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = bernoulli_sample ( a, &seed );
    pdf = bernoulli_pdf ( x, a );
    cdf = bernoulli_cdf ( x, a );
    x2 = bernoulli_cdf_inv ( cdf, a );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test010 ( )

/******************************************************************************/
/*
  Purpose:

    TEST010 tests BERNOULLI_MEAN, BERNOULLI_SAMPLE, BERNOULLI_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST010\n" );
  printf ( "  For the Bernoulli PDF:\n" );
  printf ( "  BERNOULLI_MEAN computes the mean;\n" );
  printf ( "  BERNOULLI_SAMPLE samples;\n" );
  printf ( "  BERNOULLI_VARIANCE computes the variance.\n" );

  a = 0.75E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !bernoulli_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST010 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = bernoulli_mean ( a );
  variance = bernoulli_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = bernoulli_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0105 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0105 tests BESSEL_I0_INC, BESSEL_I0_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0105:\n" );
  printf ( "  BESSEL_I0 evaluates the Bessel function of the\n" );
  printf ( "  first kind and order 0;\n" );
  printf ( "  BESSEL_I0_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      X       Exact F       BESSEL_I0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = bessel_i0 ( x );

    printf ( "  %8g  %16g  %16g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test0106 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0106 tests BESSEL_I1_INC, BESSEL_I1_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0106:\n" );
  printf ( "  BESSEL_I1 evaluates the Bessel function of the\n" );
  printf ( "  first kind and order 1;\n" );
  printf ( "  BESSEL_I1_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      X       Exact F       BESSEL_I1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = bessel_i1 ( x );

    printf ( "  %8g  %16g  %16g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test011 ( )

/******************************************************************************/
/*
  Purpose:

    TEST011 tests BETA, GAMMA;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double beta1;
  double beta2;

  printf ( "\n" );
  printf ( "TEST011\n" );
  printf ( "  BETA evaluates the Beta function;\n" );
  printf ( "  TGAMMA evaluates the Gamma function.\n" );

  a = 2.2;
  b = 3.7;

  beta1 = beta ( a, b );
  beta2 = tgamma ( a ) * tgamma ( b ) / tgamma ( a + b );

  printf ( "\n" );
  printf ( "  Argument A =                   %g\n", a );
  printf ( "  Argument B =                   %g\n", b );
  printf ( "  Beta(A,B) =                    %g\n", beta1 );
  printf ( "  (Expected value = 0.0454 )\n" );
  printf ( "\n" );
  printf ( "  Gamma(A)*Gamma(B)/Gamma(A+B) = %g\n", beta2 );

  return;
}
/******************************************************************************/

void test012 ( )

/******************************************************************************/
/*
  Purpose:

    TEST012 tests BETA_CDF, BETA_CDF_INV, BETA_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST012\n" );
  printf ( "  For the Beta PDF:\n" );
  printf ( "  BETA_CDF evaluates the CDF;\n" );
  printf ( "  BETA_CDF_INV inverts the CDF.\n" );
  printf ( "  BETA_PDF evaluates the PDF;\n" );

  a = 12.0;
  b = 12.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !beta_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST012 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "             A             B        X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = beta_sample ( a, b, &seed );
    pdf = beta_pdf ( x, a, b );
    cdf = beta_cdf ( x, a, b );
    x2 = beta_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g  %12g  %12g\n", a, b, x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test013 ( )

/******************************************************************************/
/*
  Purpose:

    TEST013 tests BETA_INC, BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST013:\n" );
  printf ( "  BETA_INC evaluates the normalized incomplete Beta\n" );
  printf ( "  function BETA_INC(A,B,X).\n" );
  printf ( "  BETA_INC_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "         A         B         X       Exact F       BETA_INC(A,B,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = beta_inc ( a, b, x );

    printf ( "  %8g  %8g  %8g  %16g  %16g\n", a, b, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test014 ( )

/******************************************************************************/
/*
  Purpose:

    TEST014 tests BETA_MEAN, BETA_SAMPLE, BETA_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST014\n" );
  printf ( "  For the Beta PDF:\n" );
  printf ( "  BETA_MEAN computes the mean;\n" );
  printf ( "  BETA_SAMPLE samples;\n" );
  printf ( "  BETA_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !beta_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST014 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = beta_mean ( a, b );
  variance = beta_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = beta_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test015 ( )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests BETA_BINOMIAL_CDF, BETA_BINOMIAL_CDF_INV, BETA_BINOMIAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST015\n" );
  printf ( "  For the Beta Binomial PDF:\n" );
  printf ( "  BETA_BINOMIAL_CDF evaluates the CDF;\n" );
  printf ( "  BETA_BINOMIAL_CDF_INV inverts the CDF.\n" );
  printf ( "  BETA_BINOMIAL_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;
  c = 4;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %d\n", c );

  if ( !beta_binomial_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST015 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = beta_binomial_sample ( a, b, c, &seed );
    pdf = beta_binomial_pdf ( x, a, b, c );
    cdf = beta_binomial_cdf ( x, a, b, c );
    x2 = beta_binomial_cdf_inv ( cdf, a, b, c );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test016 ( )

/******************************************************************************/
/*
  Purpose:

    TEST016 tests BETA_BINOMIAL_MEAN, BETA_BINOMIAL_SAMPLE, BETA_BINOMIAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST016\n" );
  printf ( "  For the Beta Binomial PDF:\n" );
  printf ( "  BETA_BINOMIAL_MEAN computes the mean;\n" );
  printf ( "  BETA_BINOMIAL_SAMPLE samples;\n" );
  printf ( "  BETA_BINOMIAL_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;
  c = 4;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %d\n", c );

  if ( !beta_binomial_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST016 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = beta_binomial_mean ( a, b, c );
  variance = beta_binomial_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = beta_binomial_sample ( a, b, c, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test020 ( )

/******************************************************************************/
/*
  Purpose:

    TEST020 tests BINOMIAL_CDF, BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  double b;
  double fx;
  double fx2;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "TEST020:\n" );
  printf ( "  BINOMIAL_CDF evaluates the cumulative distribution\n" );
  printf ( "  function for the discrete binomial probability\n" );
  printf ( "  density function.\n" );
  printf ( "  BINOMIAL_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  A is the number of trials;\n" );
  printf ( "  B is the probability of success on one trial;\n" );
  printf ( "  X is the number of successes;\n" );
  printf ( "  BINOMIAL_CDF is the probability of having up to X\n" );
  printf ( "  successes.\n" );
  printf ( "\n" );
  printf ( "      A     B         X   Exact F     BINOMIAL_CDF(A,B,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = binomial_cdf ( x, a, b );

    printf ( "  %8d  %8g  %8d  %16g  %16g\n", a, b, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test021 ( )

/******************************************************************************/
/*
  Purpose:

    TEST021 tests BINOMIAL_CDF, BINOMIAL_CDF_INV, BINOMIAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST021\n" );
  printf ( "  For the Binomial PDF:\n" );
  printf ( "  BINOMIAL_CDF evaluates the CDF;\n" );
  printf ( "  BINOMIAL_CDF_INV inverts the CDF.\n" );
  printf ( "  BINOMIAL_PDF evaluates the PDF;\n" );

  a = 5;
  b = 0.65;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !binomial_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST021 - Warning!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = binomial_sample ( a, b, &seed );
    pdf = binomial_pdf ( x, a, b );
    cdf = binomial_cdf ( x, a, b );
    x2 = binomial_cdf_inv ( cdf, a, b );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test022 ( )

/******************************************************************************/
/*
  Purpose:

    TEST022 tests BINOMIAL_COEF, BINOMIAL_COEF_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 September 2013

  Author:

    John Burkardt
*/
{
  int cnk1;
  double cnk2_log;
  double cnk2;
  int k;
  int n;

  printf ( "\n" );
  printf ( "TEST022\n" );
  printf ( "  BINOMIAL_COEF evaluates binomial coefficients.\n" );
  printf ( "  BINOMIAL_COEF_LOG evaluates the logarithm.\n" );

  printf ( "\n" );
  printf ( "    N     K       C(N,K)\n" );
  printf ( "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk1 = binomial_coef ( n, k );
      cnk2_log = binomial_coef_log ( n, k );
      cnk2 = exp ( cnk2_log );
      printf ( "  %6d  %6d  %6d  %14g\n", n, k, cnk1, cnk2 );
    }
  }
  return;
}
/******************************************************************************/

void test023 ( )

/******************************************************************************/
/*
  Purpose:

    TEST023 tests BINOMIAL_MEAN, BINOMIAL_SAMPLE, BINOMIAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST023\n" );
  printf ( "  For the Binomial PDF:\n" );
  printf ( "  BINOMIAL_MEAN computes the mean;\n" );
  printf ( "  BINOMIAL_SAMPLE samples;\n" );
  printf ( "  BINOMIAL_VARIANCE computes the variance;\n" );

  a = 5;
  b = 0.30;

  printf ( "\n" );
  printf ( "  PDF parameter A = %d\n", a );
  printf ( "  PDF parameter B = %g\n", b );

  if ( !binomial_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST023 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = binomial_mean ( a, b );
  variance = binomial_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = binomial_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0235 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0235 tests BIRTHDAY_CDF, BIRTHDAY_CDF_INV, BIRTHDAY_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2006

  Author:

    John Burkardt
*/
{
  double cdf;
  int n;
  int n2;
  double pdf;

  printf ( "\n" );
  printf ( "TEST0235\n" );
  printf ( "  For the Birthday PDF,\n" );
  printf ( "  BIRTHDAY_CDF evaluates the CDF;\n" );
  printf ( "  BIRTHDAY_CDF_INV inverts the CDF.\n" );
  printf ( "  BIRTHDAY_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       N            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( n = 1; n <= 30; n++ )
  {
    pdf = birthday_pdf ( n );

    cdf = birthday_cdf ( n );

    n2 = birthday_cdf_inv ( cdf );

    printf ( "  %8d  %14g  %14g  %8d\n", n, pdf, cdf, n2 );
  }
  return;
}
/******************************************************************************/

void test024 ( )

/******************************************************************************/
/*
  Purpose:

    TEST024 tests BRADFORD_CDF, BRADFORD_CDF_INV, BRADFORD_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST024\n" );
  printf ( "  For the Bradford PDF:\n" );
  printf ( "  BRADFORD_CDF evaluates the CDF;\n" );
  printf ( "  BRADFORD_CDF_INV inverts the CDF.\n" );
  printf ( "  BRADFORD_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !bradford_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST024 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = bradford_sample ( a, b, c, &seed );
    pdf = bradford_pdf ( x, a, b, c );
    cdf = bradford_cdf ( x, a, b, c );
    x2 = bradford_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test025 ( )

/******************************************************************************/
/*
  Purpose:

    TEST025 tests BRADFORD_MEAN, BRADFORD_SAMPLE, BRADFORD_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST025\n" );
  printf ( "  For the Bradford PDF:\n" );
  printf ( "  BRADFORD_MEAN computes the mean;\n" );
  printf ( "  BRADFORD_SAMPLE samples;\n" );
  printf ( "  BRADFORD_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !bradford_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST025 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = bradford_mean ( a, b, c );
  variance = bradford_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = bradford_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0251 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0251 tests BUFFON_LAPLACE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  int k;
  double l;
  double pdf;

  printf ( "\n" );
  printf ( "TEST0251\n" );
  printf ( "  BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF, the probability\n" );
  printf ( "  that, on a grid of cells of width A and height B,\n" );
  printf ( "  a needle of length L, dropped at random, will cross\n" );
  printf ( "  at least one grid line.\n" );
  printf ( "\n" );
  printf ( "      A         B         L        PDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    a = ( double ) ( i );
    for ( j = 1; j <= 5; j++ )
    {
      b = ( double ) ( j );
      for ( k = 0; k <= 5; k++ )
      {
        l = ( double ) ( k ) * r8_min ( a, b ) / 5.0;
        pdf = buffon_laplace_pdf ( a, b, l );
        printf ( "  %8g  %8g  %8g  %14g\n", a, b, l, pdf  );
      }
      printf ( "\n" );
    }
  }
  return;
}
/******************************************************************************/

void test0252 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0252 tests BUFFON_LAPLACE_SIMULATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2013

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  double a;
  double b;
  double err;
  int hits;
  double l;
  double pi = 3.141592653589793238462643;
  double pi_est;
  int test;
  int trial_num;
  int trial_num_test[TEST_NUM] = { 10, 100, 10000, 1000000 };

  a = 1.0;
  b = 1.0;
  l = 1.0;

  printf ( "\n" );
  printf ( "TEST0252\n" );
  printf ( "  BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle dropping\n" );
  printf ( "  experiment.  On a grid of cells of width A and height B,\n" );
  printf ( "  a needle of length L is dropped at random.  We count\n" );
  printf ( "  the number of times it crosses at least one grid line,\n" );
  printf ( "  and use this to estimate the value of PI.\n" );

  printf ( "\n" );
  printf ( "  Cell width A =    %g\n", a );
  printf ( "  Cell height B =   %g\n", b );
  printf ( "  Needle length L = %g\n", l );
  printf ( "\n" );
  printf ( "    Trials      Hits          Est(Pi)     Err\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    trial_num = trial_num_test[test];

    hits = buffon_laplace_simulate ( a, b, l, trial_num );

    if ( 0 < hits )
    {
      pi_est = ( 2.0 * l * ( a + b ) - l * l ) * ( double ) trial_num
        / ( a * b * ( double ) hits );
    }
    else
    {
      pi_est = r8_huge ( );
    }

    err = r8_abs ( pi_est - pi );

    printf ( "  %8d  %8d  %14g  %14g\n", trial_num, hits, pi_est, err );
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test0253 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0253 tests BUFFON_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  int j;
  int k;
  double l;
  double pdf;

  printf ( "\n" );
  printf ( "TEST0253\n" );
  printf ( "  BUFFON_PDF evaluates the Buffon PDF, the probability\n" );
  printf ( "  that, on a grid of cells of width A,\n" );
  printf ( "  a needle of length L, dropped at random, will cross\n" );
  printf ( "  at least one grid line.\n" );
  printf ( "\n" );
  printf ( "      A         L        PDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    a = ( double ) ( i );
    for ( k = 0; k <= 5; k++ )
    {
      l = ( double ) ( k ) * a / 5.0;
      pdf = buffon_pdf ( a, l );
      printf ( "  %8g  %8g  %14g\n", a, l, pdf );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void test0254 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0254 tests BUFFON_SIMULATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2013

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  double a;
  double err;
  int hits;
  double l;
  double pi = 3.141592653589793238462643;
  double pi_est;
  int test;
  int trial_num;
  int trial_num_test[TEST_NUM] = { 10, 100, 10000, 1000000 };

  a = 1.0;
  l = 1.0;

  printf ( "\n" );
  printf ( "TEST0254\n" );
  printf ( "  BUFFON_SIMULATE simulates a Buffon needle dropping\n" );
  printf ( "  experiment.  On a grid of cells of width A,\n" );
  printf ( "  a needle of length L is dropped at random.  We count\n" );
  printf ( "  the number of times it crosses at least one grid line,\n" );
  printf ( "  and use this to estimate the value of PI.\n" );

  printf ( "\n" );
  printf ( "  Cell width A =    %g\n", a );
  printf ( "  Needle length L = %g\n", l );
  printf ( "\n" );
  printf ( "    Trials      Hits          Est(Pi)     Err\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    trial_num = trial_num_test[test];

    hits = buffon_simulate ( a, l, trial_num );

    if ( 0 < hits )
    {
      pi_est = ( 2.0 * l * ( double ) trial_num ) / ( a * ( double ) hits );
    }
    else
    {
      pi_est = r8_huge ( );
    }

    err = r8_abs ( pi_est - pi );

    printf ( "  %8d  %8d  %14g  %14g\n", trial_num, hits, pi_est, err );
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test026 ( )

/******************************************************************************/
/*
  Purpose:

    TEST026 tests BURR_CDF, BURR_CDF_INV, BURR_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  double d;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST026\n" );
  printf ( "  For the Burr PDF:\n" );
  printf ( "  BURR_CDF evaluates the CDF;\n" );
  printf ( "  BURR_CDF_INV inverts the CDF.\n" );
  printf ( "  BURR_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;
  d = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );
  printf ( "  PDF parameter D =      %g\n", d );

  if ( !burr_check ( a, b, c, d ) )
  {
    printf ( "\n" );
    printf ( "TEST026 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = burr_sample ( a, b, c, d, &seed );
    pdf = burr_pdf ( x, a, b, c, d );
    cdf = burr_cdf ( x, a, b, c, d );
    x2 = burr_cdf_inv ( cdf, a, b, c, d );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test027 ( )

/******************************************************************************/
/*
  Purpose:

    TEST027 tests BURR_MEAN, BURR_SAMPLE, BURR_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  double d;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST027\n" );
  printf ( "  For the Burr PDF:\n" );
  printf ( "  BURR_MEAN computes the mean;\n" );
  printf ( "  BURR_SAMPLE samples;\n" );
  printf ( "  BURR_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;
  d = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );
  printf ( "  PDF parameter D =      %g\n", d );

  if ( !burr_check ( a, b, c, d ) )
  {
    printf ( "\n" );
    printf ( "TEST027 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = burr_mean ( a, b, c, d );
  variance = burr_variance ( a, b, c, d );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = burr_sample ( a, b, c, d, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0275 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0275 tests CARDIOID_CDF, CARDIOID_CDF_INV, CARDIOID_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2013

  Author:

    John Burkardt
*/
{
  double a = 0.0;
  double b = 0.25;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST0275\n" );
  printf ( "  For the Cardioid PDF:\n" );
  printf ( "  CARDIOID_CDF evaluates the CDF;\n" );
  printf ( "  CARDIOID_CDF_INV inverts the CDF.\n" );
  printf ( "  CARDIOID_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a );
  printf ( "  PDF parameter B = %g\n", b );

  if ( !cardioid_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = cardioid_sample ( a, b, &seed );
    pdf = cardioid_pdf ( x, a, b );
    cdf = cardioid_cdf ( x, a, b );
    x2 = cardioid_cdf_inv ( cdf, a, b );
    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test0276 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0276 tests CARDIOID_MEAN, CARDIOID_SAMPLE, CARDIOID_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a = 0.0;
  double b = 0.25;
  int i;
  int imax;
  int imin;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST0276\n" );
  printf ( "  For the Cardioid PDF:\n" );
  printf ( "  CARDIOID_MEAN computes the mean;\n" );
  printf ( "  CARDIOID_SAMPLE samples;\n" );
  printf ( "  CARDIOID_VARIANCE computes the variance.\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a );
  printf ( "  PDF parameter B = %g\n", b );

  if ( !cardioid_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = cardioid_mean ( a, b );
  variance = cardioid_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cardioid_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test028 ( )

/******************************************************************************/
/*
  Purpose:

    TEST028 tests CAUCHY_CDF, CAUCHY_CDF_INV, CAUCHY_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST028\n" );
  printf ( "  For the Cauchy PDF:\n" );
  printf ( "  CAUCHY_CDF evaluates the CDF;\n" );
  printf ( "  CAUCHY_CDF_INV inverts the CDF.\n" );
  printf ( "  CAUCHY_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !cauchy_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST028 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = cauchy_sample ( a, b, &seed );
    pdf = cauchy_pdf ( x, a, b );
    cdf = cauchy_cdf ( x, a, b );
    x2 = cauchy_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test029 ( )

/******************************************************************************/
/*
  Purpose:

    TEST029 tests CAUCHY_MEAN, CAUCHY_SAMPLE, CAUCHY_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST029\n" );
  printf ( "  For the Cauchy PDF:\n" );
  printf ( "  CAUCHY_MEAN computes the mean;\n" );
  printf ( "  CAUCHY_SAMPLE samples;\n" );
  printf ( "  CAUCHY_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !cauchy_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST029 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = cauchy_mean ( a, b );
  variance = cauchy_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cauchy_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean    );
  printf ( "  Sample variance = %g\n", variance);
  printf ( "  Sample maximum =  %g\n", xmax    );
  printf ( "  Sample minimum =  %g\n", xmin    );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test030 ( )

/******************************************************************************/
/*
  Purpose:

    TEST030 tests CHI_CDF, CHI_CDF_INV, CHI_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST030\n" );
  printf ( "  For the Chi PDF:\n" );
  printf ( "  CHI_CDF evaluates the CDF;\n" );
  printf ( "  CHI_CDF_INV inverts the CDF.\n" );
  printf ( "  CHI_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !chi_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST030 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = chi_sample ( a, b, c, &seed );
    pdf = chi_pdf ( x, a, b, c );
    cdf = chi_cdf ( x, a, b, c );
    x2 = chi_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests CHI_MEAN, CHI_SAMPLE, CHI_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST031\n" );
  printf ( "  For the Chi PDF:\n" );
  printf ( "  CHI_MEAN computes the mean;\n" );
  printf ( "  CHI_SAMPLE samples;\n" );
  printf ( "  CHI_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !chi_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST031 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = chi_mean ( a, b, c );
  variance = chi_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = chi_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double a2;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST032:\n" );
  printf ( "  CHI_SQUARE_CDF evaluates the cumulative\n" );
  printf ( "  distribution function for the chi-square central\n" );
  printf ( "  probability density function.\n" );
  printf ( "  CHI_SQUARE_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      A     X   Exact F     CHI_SQUARE_CDF(A,X)\n" );

  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    a2 = ( double ) a;

    fx2 = chi_square_cdf ( x, a2 );

    printf ( "  %8d  %8g  %16g  %16g\n", a, x, fx, fx2 );
  }

  return;
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_INV, CHI_SQUARE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST033\n" );
  printf ( "  For the Chi Square PDF:\n" );
  printf ( "  CHI_SQUARE_CDF evaluates the CDF;\n" );
  printf ( "  CHI_SQUARE_CDF_INV inverts the CDF.\n" );
  printf ( "  CHI_SQUARE_PDF evaluates the PDF;\n" );

  a = 4.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !chi_square_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST033 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = chi_square_sample ( a, &seed );
    pdf = chi_square_pdf ( x, a );
    cdf = chi_square_cdf ( x, a );
    x2 = chi_square_cdf_inv ( cdf, a );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test034 ( )

/******************************************************************************/
/*
  Purpose:

    TEST034 tests CHI_SQUARE_MEAN, CHI_SQUARE_SAMPLE, CHI_SQUARE_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST034\n" );
  printf ( "  For the Chi Square PDF:\n" );
  printf ( "  CHI_SQUARE_MEAN computes the mean;\n" );
  printf ( "  CHI_SQUARE_SAMPLE samples;\n" );
  printf ( "  CHI_SQUARE_VARIANCE computes the variance.\n" );

  a = 10.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !chi_square_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST034 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = chi_square_mean ( a );
  variance = chi_square_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = chi_square_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test035 ( )

/******************************************************************************/
/*
  Purpose:

    TEST035 tests CHI_SQUARE_NONCENTRAL_MEAN, CHI_SQUARE_NONCENTRAL_SAMPLE, CHI_SQUARE_NONCENTRAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  int j;
  double mean;
  int seed;
  int seed_init = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST035\n" );
  printf ( "  For the Chi Square Noncentral PDF:\n" );
  printf ( "  CHI_SQUARE_NONCENTRAL_MEAN computes the mean;\n" );
  printf ( "  CHI_SQUARE_NONCENTRAL_SAMPLE samples;\n" );
  printf ( "  CHI_SQUARE_NONCENTRAL_VARIANCE computes the variance;\n" );

  a = 3.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !chi_square_noncentral_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST035 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = chi_square_noncentral_mean ( a, b );
  variance = chi_square_noncentral_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  seed = seed_init;

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = chi_square_noncentral_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Initial seed =     %d\n", seed_init );
  printf ( "  Final seed =       %d\n", seed );

  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test036 ( )

/******************************************************************************/
/*
  Purpose:

    TEST036 tests CIRCLE_SAMPLE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  printf ( "\n" );
  printf ( "TEST036\n" );
  printf ( "  For the Circle PDF:\n" );
  printf ( "  CIRCLE_SAMPLE samples;\n" );

  a = 10.0;
  b = 4.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circle_sample ( a, b, c, &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    free ( y );
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %12g  %12g\n", mean[0], mean[1] );
  printf ( "  Sample variance = %12g  %12g\n", variance[0], variance[1] );
  printf ( "  Sample maximum =  %12g  %12g\n", xmax[0], xmax[1] );
  printf ( "  Sample minimum =  %12g  %12g\n", xmin[0], xmin[1] );

  free ( mean );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test037 ( )

/******************************************************************************/
/*
  Purpose:

    TEST037 tests CIRCULAR_NORMAL_01_*.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  printf ( "\n" );
  printf ( "TEST037\n" );
  printf ( "  For the Circular Normal 01 PDF:\n" );
  printf ( "  CIRCULAR_NORMAL_01_MEAN computes the mean;\n" );
  printf ( "  CIRCULAR_NORMAL_01_SAMPLE samples;\n" );
  printf ( "  CIRCULAR_NORMAL_01_VARIANCE computes the variance.\n" );

  mean = circular_normal_01_mean ( );
  variance = circular_normal_01_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %12g  %12g\n", mean[0], mean[1] );
  printf ( "  PDF variance = %12g  %12g\n", variance[0], variance[1] );

  free ( mean );
  free ( variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circular_normal_01_sample ( &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    free ( y );
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %12g  %12g\n", mean[0], mean[1] );
  printf ( "  Sample variance = %12g  %12g\n", variance[0], variance[1] );
  printf ( "  Sample maximum =  %12g  %12g\n", xmax[0], xmax[1] );
  printf ( "  Sample minimum =  %12g  %12g\n", xmin[0], xmin[1] );

  free ( mean );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0375 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0375 tests CIRCULAR_NORMAL_*.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a[2];
  double b;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  a[0] = 1.0;
  a[1] = 5.0;
  b = 0.75;

  printf ( "\n" );
  printf ( "TEST0375\n" );
  printf ( "  For the Circular Normal PDF:\n" );
  printf ( "  CIRCULAR_NORMAL_MEAN computes the mean;\n" );
  printf ( "  CIRCULAR_NORMAL_SAMPLE samples;\n" );
  printf ( "  CIRCULAR_NORMAL_VARIANCE computes the variance.\n" );

  mean = circular_normal_mean ( a, b );
  variance = circular_normal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %12g  %12g\n", mean[0], mean[1] );
  printf ( "  PDF variance = %12g  %12g\n", variance[0], variance[1] );

  free ( mean );
  free ( variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circular_normal_sample ( a, b, &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    free ( y );
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %12g  %12g\n", mean[0], mean[1] );
  printf ( "  Sample variance = %12g  %12g\n",variance[0], variance[1] );
  printf ( "  Sample maximum =  %12g  %12g\n",xmax[0], xmax[1] );
  printf ( "  Sample minimum =  %12g  %12g\n",xmin[0], xmin[1] );

  free ( mean );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test038 ( )

/******************************************************************************/
/*
  Purpose:

    TEST038 tests COSINE_CDF, COSINE_CDF_INV, COSINE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST038\n" );
  printf ( "  For the Cosine PDF:\n" );
  printf ( "  COSINE_CDF evaluates the CDF;\n" );
  printf ( "  COSINE_CDF_INV inverts the CDF.\n" );
  printf ( "  COSINE_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 1.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !cosine_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST038 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = cosine_sample ( a, b, &seed );
    pdf = cosine_pdf ( x, a, b );
    cdf = cosine_cdf ( x, a, b );
    x2 = cosine_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test039 ( )

/******************************************************************************/
/*
  Purpose:

    TEST039 tests COSINE_MEAN, COSINE_SAMPLE, COSINE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST039\n" );
  printf ( "  For the Cosine PDF:\n" );
  printf ( "  COSINE_MEAN computes the mean;\n" );
  printf ( "  COSINE_SAMPLE samples;\n" );
  printf ( "  COSINE_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 1.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !cosine_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST039 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = cosine_mean ( a, b );
  variance = cosine_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cosine_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n",SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0395 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0395 tests COUPON_COMPLETE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  int box_num;
  double cdf;
  double pdf;
  int type_num;

  printf ( "\n" );
  printf ( "TEST0395\n" );
  printf ( "  COUPON_COMPLETE_PDF evaluates the coupon collector's\n" );
  printf ( "  complete collection pdf.\n" );
  printf ( "\n" );

  for ( type_num = 2; type_num <= 4; type_num++ )
  {
    printf ( "\n" );
    printf ( "  Number of coupon types is %d\n", type_num );
    printf ( "\n" );
    printf ( "   BOX_NUM      PDF             CDF\n" );
    printf ( "\n" );
    cdf = 0.0;
    for ( box_num = 1; box_num <= 20; box_num++ )
    {
      pdf = coupon_complete_pdf ( type_num, box_num );
      cdf = cdf + pdf;
      printf ( "  %8d  %14g  %14g\n", box_num, pdf, cdf );
    }
  }

  return;
}
/******************************************************************************/

void test040 ( )

/******************************************************************************/
/*
  Purpose:

    TEST040 tests COUPON_SIMULATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
# define N_TRIAL 10
# define MAX_TYPE 25

  double average;
  int coupon[MAX_TYPE];
  double expect;
  int i;
  int n_coupon;
  int n_type;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST040:\n" );
  printf ( "  COUPON_SIMULATE simulates the coupon\n" );
  printf ( "  collector's problem.\n" );
  printf ( "\n" );

  for ( n_type = 5; n_type <= MAX_TYPE; n_type = n_type + 5 )
  {
    printf ( "\n" );
    printf ( "  Number of coupon types is %d\n", n_type );
    expect = ( double ) ( n_type ) * log ( ( double ) ( n_type ) );
    printf ( "  Expected wait is about %g\n", expect );
    printf ( "\n" );

    average = 0.0;
    for ( i = 1; i <= N_TRIAL; i++ )
    {
      coupon_simulate ( n_type, &seed, coupon, &n_coupon );

      printf ( "  %6d  %6d\n", i, n_coupon );

      average = average + ( double ) ( n_coupon );
    }

    average = average / ( double ) ( N_TRIAL );
    printf ( "\n" );
    printf ( "  Average wait was %g\n", average );
  }

  return;
# undef N_TRIAL
# undef MAX_TRIAL
}
/******************************************************************************/

void test041 ( )

/******************************************************************************/
/*
  Purpose:

    TEST041 tests DERANGED_CDF, DERANGED_CDF_INV, DERANGED_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST041\n" );
  printf ( "  For the Deranged PDF:\n" );
  printf ( "  DERANGED_CDF evaluates the CDF;\n" );
  printf ( "  DERANGED_CDF_INV inverts the CDF.\n" );
  printf ( "  DERANGED_PDF evaluates the PDF;\n" );

  a = 7;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %d\n", a   );

  if ( !deranged_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST041 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = deranged_sample ( a, &seed );
    pdf = deranged_pdf ( x, a );
    cdf = deranged_cdf ( x, a );
    x2 = deranged_cdf_inv ( cdf, a );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2  );
  }

  return;
}
/******************************************************************************/

void test042 ( )

/******************************************************************************/
/*
  Purpose:

    TEST042 tests DERANGED_CDF, DERANGED_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST042\n" );
  printf ( "  For the Deranged PDF:\n" );
  printf ( "  DERANGED_CDF evaluates the CDF;\n" );
  printf ( "  DERANGED_PDF evaluates the PDF;\n" );

  a = 7;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %d\n", a   );

  if ( !deranged_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST042 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 0; x <= a; x++ )
  {
    pdf = deranged_pdf ( x, a );
    cdf = deranged_cdf ( x, a );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test043 ( )

/******************************************************************************/
/*
  Purpose:

    TEST043 tests DERANGED_MEAN, DERANGED_SAMPLE, DERANGED_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST043\n" );
  printf ( "  For the Deranged PDF:\n" );
  printf ( "  DERANGED_MEAN computes the mean;\n" );
  printf ( "  DERANGED_SAMPLE samples;\n" );
  printf ( "  DERANGED_VARIANCE computes the variance.\n" );

  a = 7;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %d\n", a );

  if ( !deranged_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST042 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = deranged_mean ( a );
  variance = deranged_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = deranged_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test044 ( )

/******************************************************************************/
/*
  Purpose:

    TEST044 tests DIGAMMA, PSI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 September 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST044:\n" );
  printf ( "  DIGAMMA evaluates the DIGAMMA or PSI function.\n" );
  printf ( "  PSI_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      X   Exact F     DIGAMMA(X)\n" );

  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    if ( x <= 0.0 )
    {
      continue;
    }
    fx2 = digamma ( x );

    printf ( "  %8g  %16g  %16g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test045 ( )

/******************************************************************************/
/*
  Purpose:

    TEST045 tests DIPOLE_CDF, DIPOLE_CDF_INV, DIPOLE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 September 2013

  Author:

    John Burkardt
*/
{
# define PI 3.14159265358979323
# define TEST_NUM 3

  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int test_i;
  double test_a[TEST_NUM] = { 0.0, PI/4.0, PI/2.0 };
  double test_b[TEST_NUM] = { 1.0, 0.5,    0.0 };
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST045\n" );
  printf ( "  For the Cosine PDF:\n" );
  printf ( "  DIPOLE_CDF evaluates the CDF;\n" );
  printf ( "  DIPOLE_CDF_INV inverts the CDF.\n" );
  printf ( "  DIPOLE_PDF evaluates the PDF;\n" );

  for ( test_i = 0; test_i < TEST_NUM; test_i++ )
  {
    a = test_a[test_i];
    b = test_b[test_i];

    printf ( "\n" );
    printf ( "  PDF parameter A =      %g\n", a );
    printf ( "  PDF parameter B =      %g\n", b );

    if ( !dipole_check ( a, b ) )
    {
      printf ( "\n" );
      printf ( "TEST045 - Fatal error!\n" );
      printf ( "  The parameters are not legal.\n" );
      return;
    }

    printf ( "\n" );
    printf ( "       X            PDF           CDF            CDF_INV\n" );
    printf ( "\n" );

    for ( i = 1; i <= 10; i++ )
    {
      x = dipole_sample ( a, b, &seed );
      pdf = dipole_pdf ( x, a, b );
      cdf = dipole_cdf ( x, a, b );
      x2 = dipole_cdf_inv ( cdf, a, b );

      printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
    }

  }

  return;
# undef PI
# undef TEST_NUM
}
/******************************************************************************/

void test046 ( )

/******************************************************************************/
/*
  Purpose:

    TEST046 tests DIPOLE_SAMPLE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000
# define PI 3.14159265358979323
# define TEST_NUM 3

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double test_a[TEST_NUM] = { 0.0, PI/4.0, PI/2.0 };
  double test_b[TEST_NUM] = { 1.0, 0.5,    0.0 };
  int test_i;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST046\n" );
  printf ( "  For the Cosine PDF:\n" );
  printf ( "  DIPOLE_SAMPLE samples;\n" );

  for ( test_i = 0; test_i < TEST_NUM; test_i++ )
  {
    a = test_a[test_i];
    b = test_b[test_i];

    printf ( "\n" );
    printf ( "  PDF parameter A =      %g\n", a );
    printf ( "  PDF parameter B =      %g\n", b );

    if ( !dipole_check ( a, b ) )
    {
      printf ( "\n" );
      printf ( "TEST046 - Fatal error!\n" );
      printf ( "  The parameters are not legal.\n" );
      return;
    }

    for ( i = 0; i < SAMPLE_NUM; i++ )
    {
      x[i] = dipole_sample ( a, b, &seed );
    }

    mean = r8vec_mean ( SAMPLE_NUM, x );
    variance = r8vec_variance ( SAMPLE_NUM, x );
    xmax = r8vec_max ( SAMPLE_NUM, x );
    xmin = r8vec_min ( SAMPLE_NUM, x );

    printf ( "\n" );
    printf ( "  Sample size =     %d\n", SAMPLE_NUM );
    printf ( "  Sample mean =     %g\n", mean );
    printf ( "  Sample variance = %g\n", variance );
    printf ( "  Sample maximum =  %g\n", xmax );
    printf ( "  Sample minimum =  %g\n", xmin );
  }

  return;
# undef SAMPLE_NUM
# undef PI
# undef TEST_NUM
}
/******************************************************************************/

void test047 ( )

/******************************************************************************/
/*
  Purpose:

    TEST047 tests DIRICHLET_MEAN, DIRICHLET_SAMPLE, DIRICHLET_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define N 3
# define SAMPLE_NUM 1000

  double a[N] = { 0.250, 0.500, 1.250 };
  int i;
  int j;
  double *mean;
  double *m2;
  int seed = 123456789;
  double *variance;
  double x[N*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  printf ( "\n" );
  printf ( "TEST047\n" );
  printf ( "  For the Dirichlet PDF:\n" );
  printf ( "  DIRICHLET_MEAN computes the mean;\n" );
  printf ( "  DIRICHLET_SAMPLE samples;\n" );
  printf ( "  DIRICHLET_VARIANCE computes the variance;\n" );

  printf ( "\n" );
  printf ( "  Number of components N = %d\n", N );

  r8vec_print ( N, a, "  PDF parameter A:" );

  if ( !dirichlet_check ( N, a ) )
  {
    printf ( "\n" );
    printf ( "TEST047 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = dirichlet_mean ( N, a );
  variance = dirichlet_variance ( N, a );

  r8vec_print ( N, mean, "  PDF mean:" );
  r8vec_print ( N, variance, "  PDF variance:" );

  free ( mean );
  free ( variance );

  m2 = dirichlet_moment2 ( N, a );

  r8mat_print ( N, N, m2, "  Second moment matrix:" );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = dirichlet_sample ( N, a, &seed );
    for ( i = 0; i < N; i++ )
    {
      x[i+j*N] = y[i];
    }
    free ( y );
  }

  mean = r8row_mean ( N, SAMPLE_NUM, x );
  variance = r8row_variance ( N, SAMPLE_NUM, x );
  xmax = r8row_max ( N, SAMPLE_NUM, x );
  xmin = r8row_min ( N, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "\n" );
  printf ( "  Component Mean, Variance, Min, Max:\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %12g  %12g  %12g  %12g\n", 
      i, mean[i], variance[i], xmax[i], xmin[i] );
  }

  free ( mean );
  free ( m2 );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
# undef N
# undef SAMPLE_NUM
}
/******************************************************************************/

void test048 ( )

/******************************************************************************/
/*
  Purpose:

    TEST048 tests DIRICHLET_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N] = { 0.250, 0.500, 1.250 };
  int i;
  double pdf;
  double x[N] = { 0.500, 0.125, 0.375 };

  printf ( "\n" );
  printf ( "TEST048\n" );
  printf ( "  For the Dirichlet PDF:\n" );
  printf ( "  DIRICHLET_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "  Number of components N = %d\n", N );

  r8vec_print ( N, a, "  PDF parameter A:" );

  if ( !dirichlet_check ( N, a ) )
  {
    printf ( "\n" );
    printf ( "TEST048 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  r8vec_print ( N, x, "  PDF argument X:" );

  pdf = dirichlet_pdf ( x, N, a );

  printf ( "\n" );
  printf ( "  PDF value = %g\n", pdf );

  return;
# undef N
}
/******************************************************************************/

void test049 ( )

/******************************************************************************/
/*
  Purpose:

    TEST049 tests DIRICHLET_MIX_MEAN, DIRICHLET_MIX_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define COMP_NUM 2
# define ELEM_NUM 3
# define SAMPLE_NUM 1000

  double a[ELEM_NUM*COMP_NUM] = {
    0.250, 0.500, 1.250,
    1.500, 0.500, 2.000 };
  int comp;
  int comp_i;
  double comp_weight[COMP_NUM] = { 1.0, 2.0 };
  int elem_i;
  int i;
  int j;
  double *mean;
  double pdf;
  int seed = 123456789;
  double *variance;
  double x[ELEM_NUM*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  printf ( "\n" );
  printf ( "TEST049\n" );
  printf ( "  For the Dirichlet Mixture PDF:\n" );
  printf ( "  DIRICHLET_MIX_SAMPLE samples;\n" );
  printf ( "  DIRICHLET_MIX_MEAN computes the mean;\n" );

  printf ( "\n" );
  printf ( "  Number of elements ELEM_NUM =   %d\n", ELEM_NUM );
  printf ( "  Number of components COMP_NUM = %d\n", COMP_NUM );
  r8mat_print ( ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):" );
  r8vec_print ( COMP_NUM, comp_weight, "  Component weights" );

  if ( !dirichlet_mix_check ( COMP_NUM, ELEM_NUM, a, comp_weight ) )
  {
    printf ( "\n" );
    printf ( "TEST049 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = dirichlet_mix_mean ( COMP_NUM, ELEM_NUM, a, comp_weight );

  r8vec_print ( ELEM_NUM, mean, "  PDF mean" );

  free ( mean );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = dirichlet_mix_sample ( COMP_NUM, ELEM_NUM, a, comp_weight, &seed,
      &comp );

    for ( i = 0; i < ELEM_NUM; i++ )
    {
      x[i+j*ELEM_NUM] = y[i];
    }
    free ( y );
  }

  mean = r8row_mean ( ELEM_NUM, SAMPLE_NUM, x );
  variance = r8row_variance ( ELEM_NUM, SAMPLE_NUM, x );
  xmax = r8row_max ( ELEM_NUM, SAMPLE_NUM, x );
  xmin = r8row_min ( ELEM_NUM, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "\n" );
  printf ( "  Component Mean, Variance, Max, Min:\n" );
  printf ( "\n" );

  for ( i = 0; i < ELEM_NUM; i++ )
  {
    printf ( "  %6d  %12g  %12g  %12g  %12g\n", 
      i, mean[i], variance[i], xmax[i], xmin[i] );
  }

  free ( mean );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
}
/******************************************************************************/

void test050 ( )

/******************************************************************************/
/*
  Purpose:

    TEST050 tests DIRICHLET_MIX_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define COMP_NUM 2
# define ELEM_NUM 3

  double a[ELEM_NUM*COMP_NUM] = {
    0.250, 0.500, 1.250,
    1.500, 0.500, 2.000 };
  int comp_i;
  double comp_weight[COMP_NUM] = { 1.0, 2.0 };
  int elem_i;
  double pdf;
  double x[ELEM_NUM] = { 0.500, 0.125, 0.375 };

  printf ( "\n" );
  printf ( "TEST050\n" );
  printf ( "  For the Dirichlet mixture PDF:\n" );
  printf ( "  DIRICHLET_MIX_PDF evaluates the PDF.\n" );

  printf ( "\n" );
  printf ( "  Number of elements ELEM_NUM =   %d\n", ELEM_NUM );
  printf ( "  Number of components COMP_NUM = %d\n", COMP_NUM );
  r8mat_print ( ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):" );
  r8vec_print ( COMP_NUM, comp_weight, "  Component weights" );

  if ( !dirichlet_mix_check ( COMP_NUM, ELEM_NUM, a, comp_weight ) )
  {
    printf ( "\n" );
    printf ( "TEST050 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  r8vec_print ( ELEM_NUM, x, "  PDF argument X:" );

  pdf = dirichlet_mix_pdf ( x, COMP_NUM, ELEM_NUM, a, comp_weight );

  printf ( "\n" );
  printf ( "  PDF value =           %g\n", pdf );

  return;
# undef COMP_NUM
# undef ELEM_NUM
}
/******************************************************************************/

void test051 ( )

/******************************************************************************/
/*
  Purpose:

    TEST051 tests BETA_PDF, DIRICHLET_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define N 2

  double a;
  double aval;
  double avec[N];
  double b;
  double bval;
  int i;
  double pdf;
  double x;
  double xval;
  double xvec[N];

  xval = 0.25;
  aval = 2.50;
  bval = 3.50;

  printf ( "\n" );
  printf ( "TEST051\n" );
  printf ( "  BETA_PDF evaluates the Beta PDF.\n" );
  printf ( "  DIRICHLET_PDF evaluates the Dirichlet PDF.\n" );
  printf ( "\n" );
  printf ( "  For N = 2, Dirichlet = Beta.\n" );

  avec[0] = aval;
  avec[1] = bval;

  printf ( "\n" );
  printf ( "  Number of components N = %d", N );
  r8vec_print ( N, avec, "  PDF parameters A(1:N):" );

  if ( !dirichlet_check ( N, avec ) )
  {
    printf ( "\n" );
    printf ( "TEST051 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  xvec[0] = xval;
  xvec[1] = 1.0 - xval;

  r8vec_print ( N, xvec, "  PDF arguments X(1:N):" );

  pdf = dirichlet_pdf ( xvec, N, avec );

  printf ( "\n" );
  printf ( "  Dirichlet PDF value =  %g\n", pdf );

  x = xval;
  a = aval;
  b = bval;

  pdf = beta_pdf ( x, a, b );

  printf ( "  Beta PDF value =       %g\n", pdf );

  return;
# undef N
}
/******************************************************************************/

void test052 ( )

/******************************************************************************/
/*
  Purpose:

    TEST052 tests DISCRETE_CDF, DISCRETE_CDF_INV, DISCRETE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define A 6

  double b[A] = { 1.0, 2.0, 6.0, 2.0, 4.0, 1.0 };
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST052\n" );
  printf ( "  For the Discrete PDF:\n" );
  printf ( "  DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  DISCRETE_CDF_INV inverts the CDF.\n" );
  printf ( "  DISCRETE_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", A );
  r8vec_print ( A, b, "  PDF parameter B:" );

  if ( !discrete_check ( A, b ) )
  {
    printf ( "\n" );
    printf ( "TEST052 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = discrete_sample ( A, b, &seed );
    pdf = discrete_pdf ( x, A, b );
    cdf = discrete_cdf ( x, A, b );
    x2 = discrete_cdf_inv ( cdf, A, b );

    printf ( "  %8d  %12g  %12g  %8d\n", x, pdf, cdf, x2 );
  }

  return;
# undef A
}
/******************************************************************************/

void test053 ( )

/******************************************************************************/
/*
  Purpose:

    TEST053 tests DISCRETE_MEAN, DISCRETE_SAMPLE, DISCRETE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define A 6
# define SAMPLE_NUM 1000

  double b[A] = { 1.0, 2.0, 6.0, 2.0, 4.0, 1.0 };
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST053\n" );
  printf ( "  For the Discrete PDF:\n" );
  printf ( "  DISCRETE_MEAN computes the mean;\n" );
  printf ( "  DISCRETE_SAMPLE samples;\n" );
  printf ( "  DISCRETE_VARIANCE computes the variance;\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", A );
  r8vec_print ( A, b, "  PDF parameter B:" );

  if ( !discrete_check ( A, b ) )
  {
    printf ( "\n" );
    printf ( "TEST053 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = discrete_mean ( A, b );
  variance = discrete_variance ( A, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = discrete_sample ( A, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef A
# undef SAMPLE_NUM
}
/******************************************************************************/

void test054 ( )

/******************************************************************************/
/*
  Purpose:

    TEST054 tests EMPIRICAL_DISCRETE_CDF, EMPIRICAL_DISCRETE_CDF_INV, EMPIRICAL_DISCRETE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define A 6

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST054\n" );
  printf ( "  For the Empirical Discrete PDF:\n" );
  printf ( "  EMPIRICAL_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  EMPIRICAL_DISCRETE_CDF_INV inverts the CDF.\n" );
  printf ( "  EMPIRICAL_DISCRETE_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A = %d\n", A );
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST054 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = empirical_discrete_sample ( A, b, c, &seed );
    pdf = empirical_discrete_pdf ( x, A, b, c );
    cdf = empirical_discrete_cdf ( x, A, b, c );
    x2 = empirical_discrete_cdf_inv ( cdf, A, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
# undef A
}
/******************************************************************************/

void test055 ( )

/******************************************************************************/
/*
  Purpose:

    TEST055 tests EMPIRICAL_DISCRETE_MEAN, EMPIRICAL_DISCRETE_SAMPLE, EMPIRICAL_DISCRETE_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define A 6
# define SAMPLE_NUM 1000

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST055\n" );
  printf ( "  For the Empirical Discrete PDF:\n" );
  printf ( "  EMPIRICAL_DISCRETE_MEAN computes the mean;\n" );
  printf ( "  EMPIRICAL_DISCRETE_SAMPLE samples;\n" );
  printf ( "  EMPIRICAL_DISCRETE_VARIANCE computes the variance.\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A = %d\n", A );
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST055 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = empirical_discrete_mean ( A, b, c );
  variance = empirical_discrete_variance ( A, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = empirical_discrete_sample ( A, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef A
# undef SAMPLE_NUM
}
/******************************************************************************/

void test056 ( )

/******************************************************************************/
/*
  Purpose:

    TEST056 tests EMPIRICAL_DISCRETE_CDF, EMPIRICAL_DISCRETE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define A 6

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  double cdf;
  int i;
  double pdf;
  double x;

  printf ( "\n" );
  printf ( "TEST056\n" );
  printf ( "  For the Empirical Discrete PDF:\n" );
  printf ( "  EMPIRICAL_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  EMPIRICAL_DISCRETE_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "  PDF parameter A = %d\n", A );
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST056 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( i = -2; i <= 12; i++ )
  {
    x = ( double ) i;
    pdf = empirical_discrete_pdf ( x, A, b, c );
    cdf = empirical_discrete_cdf ( x, A, b, c );

    printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );
  }

  return;
# undef A
}
/******************************************************************************/

void test0563 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0563 tests ENGLISH_SENTENCE_LENGTH_CDF, ENGLISH_SENTENCE_LENGTH_CDF_INV and ENGLISH_SENTENCE_LENGTH_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST0563\n" );
  printf ( "  For the English Sentence Length PDF:\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_CDF evaluates the CDF;\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the CDF.\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = english_sentence_length_sample ( &seed );

    pdf = english_sentence_length_pdf ( x );

    cdf = english_sentence_length_cdf ( x );

    x2 = english_sentence_length_cdf_inv ( cdf );

    printf ( "  %8d  %12g  %12g  %8d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test0564 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0564 tests ENGLISH_SENTENCE_LENGTH_MEAN, ENGLISH_SENTENCE_LENGTH_SAMPLE and ENGLISH_SENTENCE_LENGTH_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST0564\n" );
  printf ( "  For the English Sentence Length PDF:\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_MEAN computes the mean;\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_SAMPLE samples;\n" );
  printf ( "  ENGLISH_SENTENCE_LENGTH_VARIANCE computes the variance.\n" );

  mean = english_sentence_length_mean ( );
  variance = english_sentence_length_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = english_sentence_length_sample ( &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0565 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0565 tests ENGLISH_WORD_LENGTH_CDF, ENGLISH_WORD_LENGTH_CDF_INV and ENGLISH_WORD_LENGTH_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST0565\n" );
  printf ( "  For the English Word Length PDF:\n" );
  printf ( "  ENGLISH_WORD_LENGTH_CDF evaluates the CDF;\n" );
  printf ( "  ENGLISH_WORD_LENGTH_CDF_INV inverts the CDF.\n" );
  printf ( "  ENGLISH_WORD_LENGTH_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = english_word_length_sample ( &seed );

    pdf = english_word_length_pdf ( x );

    cdf = english_word_length_cdf ( x );

    x2 = english_word_length_cdf_inv ( cdf );

    printf ( "  %8d  %12g  %12g  %8d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test0566 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0566 tests ENGLISH_WORD_LENGTH_MEAN, ENGLISH_WORD_LENGTH_SAMPLE and ENGLISH_WORD_LENGTH_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST0566\n" );
  printf ( "  For the English Word Length PDF:\n" );
  printf ( "  ENGLISH_WORD_LENGTH_MEAN computes the mean;\n" );
  printf ( "  ENGLISH_WORD_LENGTH_SAMPLE samples;\n" );
  printf ( "  ENGLISH_WORD_LENGTH_VARIANCE computes the variance.\n" );

  mean = english_word_length_mean ( );
  variance = english_word_length_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = english_word_length_sample ( &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax );
  printf ( "  Sample minimum =  %d\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test057 ( )

/******************************************************************************/
/*
  Purpose:

    TEST057 tests ERLANG_CDF, ERLANG_CDF_INV, ERLANG_PDF;

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
  int c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST057\n" );
  printf ( "  For the Erlang PDF:\n" );
  printf ( "  ERLANG_CDF evaluates the CDF;\n" );
  printf ( "  ERLANG_CDF_INV inverts the CDF.\n" );
  printf ( "  ERLANG_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %d\n", c );

  if ( !erlang_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST057 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = erlang_sample ( a, b, c, &seed );
    pdf = erlang_pdf ( x, a, b, c );
    cdf = erlang_cdf ( x, a, b, c );
    x2 = erlang_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test058 ( )

/******************************************************************************/
/*
  Purpose:

    TEST058 tests ERLANG_MEAN, ERLANG_SAMPLE, ERLANG_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST058\n" );
  printf ( "  For the Erlang PDF:\n" );
  printf ( "  ERLANG_MEAN computes the mean;\n" );
  printf ( "  ERLANG_SAMPLE samples;\n" );
  printf ( "  ERLANG_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %d\n", c );

  if ( !erlang_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST058 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = erlang_mean ( a, b, c );
  variance = erlang_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = erlang_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test059 ( )

/******************************************************************************/
/*
  Purpose:

    TEST059 tests ERROR_F and ERROR_F_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  double x;
  double y;
  double z;

  printf ( "\n" );
  printf ( "TEST059\n" );
  printf ( "  ERROR_F evaluates ERF(X).\n" );
  printf ( "  ERROR_F_INVERSE inverts ERF(X).\n" );
  printf ( "\n" );
  printf ( "X   -> Y = error_F(X) -> Z = error_f_inverse(Y)\n" );
  printf ( "\n" );

  seed = 123456789;

  x = 1.0;

  for ( i = 1; i <= 20; i++ )
  {
    x = normal_01_sample ( &seed );
    y = error_f ( x );
    z = error_f_inverse ( y );
    printf ( "  %12g  %12g  %12g\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test060 ( )

/******************************************************************************/
/*
  Purpose:

    TEST060 tests EXPONENTIAL_01_CDF, EXPONENTIAL_01_CDF_INV, EXPONENTIAL_01_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST060\n" );
  printf ( "  For the Exponential 01 PDF:\n" );
  printf ( "  EXPONENTIAL_01_CDF evaluates the CDF;\n" );
  printf ( "  EXPONENTIAL_01_CDF_INV inverts the CDF.\n" );
  printf ( "  EXPONENTIAL_01_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = exponential_01_sample ( &seed );
    pdf = exponential_01_pdf ( x );
    cdf = exponential_01_cdf ( x );
    x2 = exponential_01_cdf_inv ( cdf );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test061 ( )

/******************************************************************************/
/*
  Purpose:

    TEST061 tests EXPONENTIAL_01_MEAN, EXPONENTIAL_01_SAMPLE, EXPONENTIAL_01_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST061\n" );
  printf ( "  For the Exponential 01 PDF:\n" );
  printf ( "  EXPONENTIAL_01_MEAN computes the mean;\n" );
  printf ( "  EXPONENTIAL_01_SAMPLE samples;\n" );
  printf ( "  EXPONENTIAL_01_VARIANCE computes the variance.\n" );

  mean     = exponential_01_mean ( );
  variance = exponential_01_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = exponential_01_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test062 ( )

/******************************************************************************/
/*
  Purpose:

    TEST062 tests EXPONENTIAL_CDF, EXPONENTIAL_CDF_INV, EXPONENTIAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST062\n" );
  printf ( "  For the Exponential PDF:\n" );
  printf ( "  EXPONENTIAL_CDF evaluates the CDF;\n" );
  printf ( "  EXPONENTIAL_CDF_INV inverts the CDF.\n" );
  printf ( "  EXPONENTIAL_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !exponential_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST062 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = exponential_sample ( a, b, &seed );
    pdf = exponential_pdf ( x, a, b );
    cdf = exponential_cdf ( x, a, b );
    x2 = exponential_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2  );
  }

  return;
}
/******************************************************************************/

void test063 ( )

/******************************************************************************/
/*
  Purpose:

    TEST063 tests EXPONENTIAL_MEAN, EXPONENTIAL_SAMPLE, EXPONENTIAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST063\n" );
  printf ( "  For the Exponential PDF:\n" );
  printf ( "  EXPONENTIAL_MEAN computes the mean;\n" );
  printf ( "  EXPONENTIAL_SAMPLE samples;\n" );
  printf ( "  EXPONENTIAL_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !exponential_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST063 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = exponential_mean ( a, b );
  variance = exponential_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = exponential_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test064 ( )

/******************************************************************************/
/*
  Purpose:

    TEST064 tests EXTREME_VALUES_CDF, EXTREME_VALUES_CDF_INV, EXTREME_VALUES_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST064\n" );
  printf ( "  For the Extreme Values PDF:\n" );
  printf ( "  EXTREME_VALUES_CDF evaluates the CDF;\n" );
  printf ( "  EXTREME_VALUES_CDF_INV inverts the CDF.\n" );
  printf ( "  EXTREME_VALUES_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !extreme_values_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST064 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = extreme_values_sample ( a, b, &seed );
    pdf = extreme_values_pdf ( x, a, b );
    cdf = extreme_values_cdf ( x, a, b );
    x2 = extreme_values_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test065 ( )

/******************************************************************************/
/*
  Purpose:

    TEST065 tests EXTREME_VALUES_MEAN, EXTREME_VALUES_SAMPLE, EXTREME_VALUES_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST065\n" );
  printf ( "  For the Extreme Values PDF:\n" );
  printf ( "  EXTREME_VALUES_MEAN computes the mean;\n" );
  printf ( "  EXTREME_VALUES_SAMPLE samples;\n" );
  printf ( "  EXTREME_VALUES_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !extreme_values_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST065 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = extreme_values_mean ( a, b );
  variance = extreme_values_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = extreme_values_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test066 ( )

/******************************************************************************/
/*
  Purpose:

    TEST066 tests F_CDF, F_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST066:\n" );
  printf ( "  F_CDF evaluates the cumulative\n" );
  printf ( "  distribution function for the F\n" );
  printf ( "  probability density function.\n" );
  printf ( "  F_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      A     B     X   Exact F     F_CDF(A,B,X)\n" );

  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = f_cdf ( x, a, b );

    printf ( "  %8d  %8d  %8g  %8g  %8g\n", a, b, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test067 ( )

/******************************************************************************/
/*
  Purpose:

    TEST067 tests F_CDF, F_PDF, F_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  int m;
  int n;
  double pdf;
  int seed = 123456789;
  double x;

  printf ( "\n" );
  printf ( "TEST067\n" );
  printf ( "  For the F PDF:\n" );
  printf ( "  F_CDF evaluates the CDF;\n" );
  printf ( "  F_PDF evaluates the PDF;\n" );
  printf ( "  F_SAMPLE samples the PDF;\n" );

  m = 1;
  n = 1;

  printf ( "\n" );
  printf ( "  PDF parameter M = %d\n", m   );
  printf ( "  PDF parameter N = %d\n", n   );

  if ( !f_check ( m, n ) )
  {
    printf ( "\n" );
    printf ( "TEST067 - Fatal error!\n" );
    printf ( "  The parameter values are illegal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = f_sample ( m, n, &seed );
    pdf = f_pdf ( x, m, n );
    cdf = f_cdf ( x, m, n );

    printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test068 ( )

/******************************************************************************/
/*
  Purpose:

    TEST068 tests F_MEAN, F_SAMPLE, F_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST068\n" );
  printf ( "  For the F PDF:\n" );
  printf ( "  F_MEAN computes the mean;\n" );
  printf ( "  F_SAMPLE samples;\n" );
  printf ( "  F_VARIANCE computes the variance;\n" );

  m = 8;
  n = 6;

  printf ( "\n" );
  printf ( "  PDF parameter M = %d\n", m );
  printf ( "  PDF parameter N = %d\n", n );

  if ( !f_check ( m, n ) )
  {
    printf ( "\n" );
    printf ( "TEST068 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = f_mean ( m, n );
  variance = f_variance ( m, n );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = f_sample ( m, n, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test069 ( )

/******************************************************************************/
/*
  Purpose:

    TEST069 tests FACTORIAL_LOG, GAMMA_LOG_INT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double f;
  double g;
  int i;
  double x;

  printf ( "\n" );
  printf ( "TEST069\n" );
  printf ( "  FACTORIAL_LOG evaluates the log of the factorial function;\n" );
  printf ( "  GAMMA_LOG_INT evaluates the log for integer argument.\n" );

  printf ( "\n" );
  printf ( "  I, GAMMA_LOG_INT(I+1) FACTORIAL_LOG(I)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    g = gamma_log_int ( i+1 );

    f = factorial_log ( i );

    printf ( "  %6d  %12g  %12g\n", i, g, f );
  }

  return;
}
/******************************************************************************/

void test070 ( )

/******************************************************************************/
/*
  Purpose:

    TEST070 tests FACTORIAL_STIRLING, I4_FACTORIAL;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  int i;
  double value;

  printf ( "\n" );
  printf ( "TEST070\n" );
  printf ( "  FACTORIAL_STIRLING computes Stirling's\n" );
  printf ( "  approximate factorial function;\n" );
  printf ( "  I4_FACTORIAL evaluates the factorial function;\n" );
  printf ( "\n" );
  printf ( "  N      Stirling     N!\n" );
  printf ( "\n" );

  for ( i = 0; i <= 20; i++ )
  {
    value = factorial_stirling ( i );

    printf ( "  %6d  %12g  %20g\n", i, value, i4_factorial ( i ) );
  }

  return;
}
/******************************************************************************/

void test0705 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0705 tests FISHER_PDF and FISHER_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  int j;
  double kappa;
  double mu[3];
  int n = 10;
  double pdf;
  int seed;
  int test;
  int test_num = 3;
  double *x;

  printf ( "\n" );
  printf ( "TEST0705\n" );
  printf ( "  For the Fisher PDF:\n" );
  printf ( "  FISHER_PDF evaluates the PDF.\n" );
  printf ( "  FISHER_SAMPLE samples the PDF.\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      kappa = 0.0;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }
    else if ( test == 2 )
    {
      kappa = 0.5;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }
    else if ( test == 3 )
    {
      kappa = 10.0;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }

    printf ( "\n" );
    printf ( "  PDF parameters:\n" );
    printf ( "    Concentration parameter KAPPA = %g\n", kappa );
    printf ( "    Direction MU(1:3) = %g  %g  %g\n", mu[0], mu[1], mu[2] );

    printf ( "\n" );
    printf ( "      X                         PDF\n" );
    printf ( "\n" );

    seed = 123456789;
    x = fisher_sample ( kappa, mu, n, &seed );

    for ( j = 0; j < n; j++ )
    {
      pdf = fisher_pdf ( x+j*3, kappa, mu );

      printf ( "  %10g  %10g  %10g  %14g\n", x[0+j*3], x[1+j*3], x[2+j*3], pdf );
    }
    free ( x );
  }
  return;
}
/******************************************************************************/

void test071 ( )

/******************************************************************************/
/*
  Purpose:

    TEST071 tests FISK_CDF, FISK_CDF_INV, FISK_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST071\n" );
  printf ( "  For the Fisk PDF:\n" );
  printf ( "  FISK_CDF evaluates the CDF;\n" );
  printf ( "  FISK_CDF_INV inverts the CDF.\n" );
  printf ( "  FISK_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !fisk_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST071 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = fisk_sample ( a, b, c, &seed );
    pdf = fisk_pdf ( x, a, b, c );
    cdf = fisk_cdf ( x, a, b, c );
    x2 = fisk_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test072 ( )

/******************************************************************************/
/*
  Purpose:

    TEST072 tests FISK_MEAN, FISK_SAMPLE, FISK_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST072\n" );
  printf ( "  For the Fisk PDF:\n" );
  printf ( "  FISK_MEAN computes the mean;\n" );
  printf ( "  FISK_SAMPLE samples;\n" );
  printf ( "  FISK_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !fisk_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST072 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = fisk_mean ( a, b, c );
  variance = fisk_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = fisk_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test073 ( )

/******************************************************************************/
/*
  Purpose:

    TEST073 tests FOLDED_NORMAL_CDF, FOLDED_NORMAL_CDF_INV, FOLDED_NORMAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST073\n" );
  printf ( "  For the Folded Normal PDF:\n" );
  printf ( "  FOLDED_NORMAL_CDF evaluates the CDF;\n" );
  printf ( "  FOLDED_NORMAL_CDF_INV inverts the CDF.\n" );
  printf ( "  FOLDED_NORMAL_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !folded_normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST073 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = folded_normal_sample ( a, b, &seed );
    pdf = folded_normal_pdf ( x, a, b );
    cdf = folded_normal_cdf ( x, a, b );
    x2 = folded_normal_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test074 ( )

/******************************************************************************/
/*
  Purpose:

    TEST074 tests FOLDED_NORMAL_MEAN, FOLDED_NORMAL_SAMPLE, FOLDED_NORMAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST074\n" );
  printf ( "  For the Folded Normal PDF:\n" );
  printf ( "  FOLDED_NORMAL_MEAN computes the mean;\n" );
  printf ( "  FOLDED_NORMAL_SAMPLE samples;\n" );
  printf ( "  FOLDED_NORMAL_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !folded_normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST074 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = folded_normal_mean ( a, b );
  variance = folded_normal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = folded_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0744 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0744 tests FRECHET_CDF, FRECHET_CDF_INV, FRECHET_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double alpha;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST0744\n" );
  printf ( "  For the Frechet PDF:\n" );
  printf ( "  FRECHET_CDF evaluates the CDF;\n" );
  printf ( "  FRECHET_CDF_INV inverts the CDF.\n" );
  printf ( "  FRECHET_PDF evaluates the PDF;\n" );

  alpha = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter ALPHA =  %g\n", alpha );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = frechet_sample ( alpha, &seed );
    pdf = frechet_pdf ( x, alpha );
    cdf = frechet_cdf ( x, alpha );
    x2 = frechet_cdf_inv ( cdf, alpha );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test0745 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0745 tests FRECHET_MEAN, FRECHET_SAMPLE, FRECHET_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double alpha;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST0745\n" );
  printf ( "  For the Frechet PDF:\n" );
  printf ( "  FRECHET_MEAN computes the mean;\n" );
  printf ( "  FRECHET_SAMPLE samples;\n" );
  printf ( "  FRECHET_VARIANCE computes the variance;\n" );

  alpha = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter ALPHA =  %g\n", alpha );

  mean = frechet_mean ( alpha );
  variance = frechet_variance ( alpha );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = frechet_sample ( alpha, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test075 ( )

/******************************************************************************/
/*
  Purpose:

    TEST075 tests GAMMA, GAMMA_LOG, GAMMA_LOG_INT, I4_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2007

  Author:

    John Burkardt
*/
{
  double g1;
  double g2;
  double g3;
  double g4;
  int i;
  double x;

  printf ( "\n" );
  printf ( "TEST075\n" );
  printf ( "  GAMMA evaluates the Gamma function;\n" );
  printf ( "  GAMMA_LOG evaluates the log of the Gamma function;\n" );
  printf ( "  GAMMA_LOG_INT evaluates the log for integer argument;\n" );
  printf ( "  I4_FACTORIAL evaluates the factorial function.\n" );

  printf ( "\n" );
  printf ( "  X, GAMMA(X), Exp(GAMMA_LOG(X)), Exp(GAMMA_LOG_INT(X)) " );
  printf ( "I4_FACTORIAL(X+1)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = ( double ) i;
    g1 = tgamma ( x );
    g2 = exp ( gamma_log ( x ) );
    g3 = exp ( gamma_log_int ( i ) );
    g4 = i4_factorial ( i - 1 );
    printf ( "  %6g  %14g  %14g  %14g  %14g\n", x, g1, g2, g3, g4 );
  }

  return;
}
/******************************************************************************/

void test076 ( )

/******************************************************************************/
/*
  Purpose:

    TEST076 tests GAMMA_INC, GAMMA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST076:\n" );
  printf ( "  GAMMA_INC evaluates the normalized incomplete Gamma\n" );
  printf ( "  function GAMMA_INC(A,B,X).\n" );
  printf ( "  GAMMA_INC_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   A      X       Exact F       GAMMA_INC(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gamma_inc ( a, x );

    printf ( "  %8g  %8g  %16g  %16g\n", a, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test077 ( )

/******************************************************************************/
/*
  Purpose:

    TEST077 tests GAMMA_CDF, GAMMA_PDF, GAMMA_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  printf ( "\n" );
  printf ( "TEST077\n" );
  printf ( "  For the Gamma PDF:\n" );
  printf ( "  GAMMA_CDF evaluates the CDF;\n" );
  printf ( "  GAMMA_PDF evaluates the PDF;\n" );
  printf ( "  GAMMA_SAMPLE samples the PDF;\n" );

  a = 1.0;
  b = 1.5;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a );
  printf ( "  PDF parameter B = %g\n", b );
  printf ( "  PDF parameter B = %g\n", c );

  if ( !gamma_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST077 - Fatal error!\n" );
    printf ( "  The parameter values are illegal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = gamma_sample ( a, b, c, &seed );
    pdf = gamma_pdf ( x, a, b, c );
    cdf = gamma_cdf ( x, a, b, c );

    printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test078 ( )

/******************************************************************************/
/*
  Purpose:

    TEST078 tests GAMMA_MEAN, GAMMA_SAMPLE, GAMMA_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST078\n" );
  printf ( "  For the Gamma PDF:\n" );
  printf ( "  GAMMA_MEAN computes the mean;\n" );
  printf ( "  GAMMA_SAMPLE samples;\n" );
  printf ( "  GAMMA_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 3.0;
  c = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !gamma_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST078 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = gamma_mean ( a, b, c );
  variance = gamma_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gamma_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n",xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test079 ( )

/******************************************************************************/
/*
  Purpose:

    TEST079 tests GENLOGISTIC_CDF, GENLOGISTIC_CDF_INV, GENLOGISTIC_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST079\n" );
  printf ( "  For the Genlogistic PDF:\n" );
  printf ( "  GENLOGISTIC_CDF evaluates the CDF;\n" );
  printf ( "  GENLOGISTIC_CDF_INV inverts the CDF.\n" );
  printf ( "  GENLOGISTIC_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !genlogistic_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST079 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = genlogistic_sample ( a, b, c, &seed );
    pdf = genlogistic_pdf ( x, a, b, c );
    cdf = genlogistic_cdf ( x, a, b, c );
    x2 = genlogistic_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test080 ( )

/******************************************************************************/
/*
  Purpose:

    TEST080 tests GENLOGISTIC_MEAN, GENLOGISTIC_SAMPLE, GENLOGISTIC_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST080\n" );
  printf ( "  For the Genlogistic PDF:\n" );
  printf ( "  GENLOGISTIC_MEAN computes the mean;\n" );
  printf ( "  GENLOGISTIC_SAMPLE samples;\n" );
  printf ( "  GENLOGISTIC_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !genlogistic_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST080 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = genlogistic_mean ( a, b, c );
  variance = genlogistic_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = genlogistic_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test081 ( )

/******************************************************************************/
/*
  Purpose:

    TEST081 tests GEOMETRIC_CDF, GEOMETRIC_CDF_INV, GEOMETRIC_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST081\n" );
  printf ( "  For the Geometric PDF:\n" );
  printf ( "  GEOMETRIC_CDF evaluates the CDF;\n" );
  printf ( "  GEOMETRIC_CDF_INV inverts the CDF.\n" );
  printf ( "  GEOMETRIC_PDF evaluates the PDF;\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !geometric_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST081 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = geometric_sample ( a, &seed );
    pdf = geometric_pdf ( x, a );
    cdf = geometric_cdf ( x, a );
    x2 = geometric_cdf_inv ( cdf, a );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test082 ( )

/******************************************************************************/
/*
  Purpose:

    TEST082 tests GEOMETRIC_MEAN, GEOMETRIC_SAMPLE, GEOMETRIC_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST082\n" );
  printf ( "  For the Geometric PDF:\n" );
  printf ( "  GEOMETRIC_MEAN computes the mean;\n" );
  printf ( "  GEOMETRIC_SAMPLE samples;\n" );
  printf ( "  GEOMETRIC_VARIANCE computes the variance.\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a );

  if ( !geometric_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST082 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = geometric_mean ( a );
  variance = geometric_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = geometric_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test083 ( )

/******************************************************************************/
/*
  Purpose:

    TEST083 tests GEOMETRIC_CDF, GEOMETRIC_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST083\n" );
  printf ( "  For the Geometric PDF:\n" );
  printf ( "  GEOMETRIC_CDF evaluates the CDF;\n" );
  printf ( "  GEOMETRIC_PDF evaluates the PDF;\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !geometric_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST083 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 0; x <= 10; x++ )
  {
    pdf = geometric_pdf ( x, a );
    cdf = geometric_cdf ( x, a );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test084 ( )

/******************************************************************************/
/*
  Purpose:

    TEST084 tests GOMPERTZ_CDF, GOMPERTZ_CDF_INV, GOMPERTZ_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST084\n" );
  printf ( "  For the Gompertz PDF:\n" );
  printf ( "  GOMPERTZ_CDF evaluates the CDF;\n" );
  printf ( "  GOMPERTZ_CDF_INV inverts the CDF.\n" );
  printf ( "  GOMPERTZ_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !gompertz_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST084 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = gompertz_sample ( a, b, &seed );
    pdf = gompertz_pdf ( x, a, b );
    cdf = gompertz_cdf ( x, a, b );
    x2 = gompertz_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test085 ( )

/******************************************************************************/
/*
  Purpose:

    TEST085 tests GOMPERTZ_MEAN, GOMPERTZ_SAMPLE, GOMPERTZ_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST085\n" );
  printf ( "  For the Gompertz PDF:\n" );
  printf ( "  GOMPERTZ_MEAN computes the mean;\n" );
  printf ( "  GOMPERTZ_SAMPLE samples;\n" );
  printf ( "  GOMPERTZ_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !gompertz_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST085 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gompertz_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test086 ( )

/******************************************************************************/
/*
  Purpose:

    TEST086 tests GUMBEL_CDF, GUMBEL_CDF_INV, GUMBEL_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST086\n" );
  printf ( "  For the Gumbel PDF:\n" );
  printf ( "  GUMBEL_CDF evaluates the CDF;\n" );
  printf ( "  GUMBEL_CDF_INV inverts the CDF.\n" );
  printf ( "  GUMBEL_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = gumbel_sample ( &seed );
    pdf = gumbel_pdf ( x );
    cdf = gumbel_cdf ( x );
    x2 = gumbel_cdf_inv ( cdf );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test087 ( )

/******************************************************************************/
/*
  Purpose:

    TEST087 tests GUMBEL_MEAN, GUMBEL_SAMPLE, GUMBEL_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST087\n" );
  printf ( "  For the Gumbel PDF:\n" );
  printf ( "  GUMBEL_MEAN computes the mean;\n" );
  printf ( "  GUMBEL_SAMPLE samples;\n" );
  printf ( "  GUMBEL_VARIANCE computes the variance.\n" );

  mean     = gumbel_mean ( );
  variance = gumbel_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gumbel_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test088 ( )

/******************************************************************************/
/*
  Purpose:

    TEST088 tests HALF_NORMAL_CDF, HALF_NORMAL_CDF_INV, HALF_NORMAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST088\n" );
  printf ( "  For the Half Normal PDF:\n" );
  printf ( "  HALF_NORMAL_CDF evaluates the CDF;\n" );
  printf ( "  HALF_NORMAL_CDF_INV inverts the CDF.\n" );
  printf ( "  HALF_NORMAL_PDF evaluates the PDF;\n" );

  a = 0.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !half_normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST088 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = half_normal_sample ( a, b, &seed );
    pdf = half_normal_pdf ( x, a, b );
    cdf = half_normal_cdf ( x, a, b );
    x2 = half_normal_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test089 ( )

/******************************************************************************/
/*
  Purpose:

    TEST089 tests HALF_NORMAL_MEAN, HALF_NORMAL_SAMPLE, HALF_NORMAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST089\n" );
  printf ( "  For the Half Normal PDF:\n" );
  printf ( "  HALF_NORMAL_MEAN computes the mean;\n" );
  printf ( "  HALF_NORMAL_SAMPLE samples;\n" );
  printf ( "  HALF_NORMAL_VARIANCE computes the variance;\n" );

  a = 0.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !half_normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST089 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = half_normal_mean ( a, b );
  variance = half_normal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = half_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test090 ( )

/******************************************************************************/
/*
  Purpose:

  TEST090 tests HYPERGEOMETRIC_CDF, HYPERGEOMETRIC_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int l;
  int m;
  int n;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST090\n" );
  printf ( "  For the Hypergeometric PDF:\n" );
  printf ( "  HYPERGEOMETRIC_CDF evaluates the CDF.\n" );
  printf ( "  HYPERGEOMETRIC_PDF evaluates the PDF.\n" );

  n = 100;
  m = 70;
  l = 1000;

  printf ( "\n" );
  printf ( "  Total number of balls L =         %d\n", l );
  printf ( "  Number of white balls M =         %d\n", m );
  printf ( "  Number of balls taken N =         %d\n", n );

  if ( !hypergeometric_check ( n, m, l ) )
  {
    printf ( "\n" );
    printf ( "TEST090 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  x = 7;

  pdf = hypergeometric_pdf ( x, n, m, l );

  cdf = hypergeometric_cdf ( x, n, m, l );

  printf ( "  PDF argument X =                %d\n", x   );
  printf ( "  PDF value =                   = %g\n", pdf );
  printf ( "  CDF value =                   = %g\n", cdf );

  return;
}
/******************************************************************************/

void test091 ( )

/******************************************************************************/
/*
  Purpose:

    TEST091 tests HYPERGEOMETRIC_MEAN, HYPERGEOMETRIC_SAMPLE, HYPERGEOMETRIC_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  int j;
  int l;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST091\n" );
  printf ( "  For the Hypergeometric PDF:\n" );
  printf ( "  HYPERGEOMETRIC_MEAN computes the mean;\n" );
  printf ( "  HYPERGEOMETRIC_SAMPLE samples;\n" );
  printf ( "  HYPERGEOMETRIC_VARIANCE computes the variance.\n" );

  n = 100;
  m = 70;
  l = 1000;

  printf ( "\n" );
  printf ( "  Total number of balls L =         %d\n", l );
  printf ( "  Number of white balls M =         %d\n", m );
  printf ( "  Number of balls taken N =         %d\n", n );

  if ( !hypergeometric_check ( n, m, l ) )
  {
    printf ( "\n" );
    printf ( "TEST090 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = hypergeometric_mean ( n, m, l );
  variance = hypergeometric_variance ( n, m, l );

  printf ( "  PDF mean =                    %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  printf ( "\n" );
  printf ( "THIS CALL IS TAKING FOREVER!\n" );
  return;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = hypergeometric_sample ( n, m, l, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test092 ( )

/******************************************************************************/
/*
  Purpose:

    TEST092 tests R8_CEILING.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  int i;
  int ival;
  double rval;

  printf ( "\n" );
  printf ( "TEST092\n" );
  printf ( "  R8_CEILING rounds an R8 up.\n" );
  printf ( "\n" );
  printf ( "       X           R8_CEILING(X)\n" );
  printf ( "\n" );

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( double ) ( i ) / 5.0;
    ival = r8_ceiling ( rval );
    printf ( "  %14g  %6d\n", rval, ival );
  }

  return;
}
/******************************************************************************/

void test093 ( )

/******************************************************************************/
/*
  Purpose:

    TEST093 tests INVERSE_GAUSSIAN_CDF, INVERSE_GAUSSIAN_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  printf ( "\n" );
  printf ( "TEST093\n" );
  printf ( "  For the Inverse Gaussian PDF:\n" );
  printf ( "  INVERSE_GAUSSIAN_CDF evaluates the CDF;\n" );
  printf ( "  INVERSE_GAUSSIAN_PDF evaluates the PDF;\n" );

  a = 5.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !inverse_gaussian_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST093 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = inverse_gaussian_sample ( a, b, &seed );
    pdf = inverse_gaussian_pdf ( x, a, b );
    cdf = inverse_gaussian_cdf ( x, a, b );

    printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test094 ( )

/******************************************************************************/
/*
  Purpose:

    TEST094 tests INVERSE_GAUSSIAN_MEAN, INVERSE_GAUSSIAN_SAMPLE, INVERSE_GAUSSIAN_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST094\n" );
  printf ( "  For the Inverse Gaussian PDF:\n" );
  printf ( "  INVERSE_GAUSSIAN_MEAN computes the mean;\n" );
  printf ( "  INVERSE_GAUSSIAN_SAMPLE samples;\n" );
  printf ( "  INVERSE_GAUSSIAN_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !inverse_gaussian_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST094 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = inverse_gaussian_mean ( a, b );
  variance = inverse_gaussian_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = inverse_gaussian_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test095 ( )

/******************************************************************************/
/*
  Purpose:

    TEST095 tests LAPLACE_CDF, LAPLACE_CDF_INV, LAPLACE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST095\n" );
  printf ( "  For the Laplace PDF:\n" );
  printf ( "  LAPLACE_CDF evaluates the CDF;\n" );
  printf ( "  LAPLACE_CDF_INV inverts the CDF.\n" );
  printf ( "  LAPLACE_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !laplace_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST095 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = laplace_sample ( a, b, &seed );
    pdf = laplace_pdf ( x, a, b );
    cdf = laplace_cdf ( x, a, b );
    x2 = laplace_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test096 ( )

/******************************************************************************/
/*
  Purpose:

    TEST096 tests LAPLACE_MEAN, LAPLACE_SAMPLE, LAPLACE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST096\n" );
  printf ( "  For the Laplace PDF:\n" );
  printf ( "  LAPLACE_MEAN computes the mean;\n" );
  printf ( "  LAPLACE_SAMPLE samples;\n" );
  printf ( "  LAPLACE_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !laplace_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST096 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = laplace_mean ( a, b );
  variance = laplace_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = laplace_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test0965 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0965 tests LEVY_CDF, LEVY_CDF_INV, LEVY_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST0965\n" );
  printf ( "  For the Levy PDF:\n" );
  printf ( "  LEVY_CDF evaluates the CDF;\n" );
  printf ( "  LEVY_CDF_INV inverts the CDF.\n" );
  printf ( "  LEVY_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = levy_sample ( a, b, &seed );
    pdf = levy_pdf ( x, a, b );
    cdf = levy_cdf ( x, a, b );
    x2 = levy_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test097 ( )

/******************************************************************************/
/*
  Purpose:

    TEST097 tests LOGISTIC_CDF, LOGISTIC_CDF_INV, LOGISTIC_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST097\n" );
  printf ( "  For the Logistic PDF:\n" );
  printf ( "  LOGISTIC_CDF evaluates the CDF;\n" );
  printf ( "  LOGISTIC_CDF_INV inverts the CDF.\n" );
  printf ( "  LOGISTIC_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !logistic_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST097 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = logistic_sample ( a, b, &seed );
    pdf = logistic_pdf ( x, a, b );
    cdf = logistic_cdf ( x, a, b );
    x2 = logistic_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test098 ( )

/******************************************************************************/
/*
  Purpose:

    TEST098 tests LOGISTIC_MEAN, LOGISTIC_SAMPLE, LOGISTIC_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 February 2007

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST098\n" );
  printf ( "  For the Logistic PDF:\n" );
  printf ( "  LOGISTIC_MEAN computes the mean;\n" );
  printf ( "  LOGISTIC_SAMPLE samples;\n" );
  printf ( "  LOGISTIC_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !logistic_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST098 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = logistic_mean ( a, b );
  variance = logistic_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = logistic_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test099 ( )

/******************************************************************************/
/*
  Purpose:

    TEST099 tests LOG_NORMAL_CDF, LOG_NORMAL_CDF_INV, LOG_NORMAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST099\n" );
  printf ( "  For the Log Normal PDF:\n" );
  printf ( "  LOG_NORMAL_CDF evaluates the CDF;\n" );
  printf ( "  LOG_NORMAL_CDF_INV inverts the CDF.\n" );
  printf ( "  LOG_NORMAL_PDF evaluates the PDF;\n" );

  a = 10.0;
  b = 2.25;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !log_normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST099 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = log_normal_sample ( a, b, &seed );
    pdf = log_normal_pdf ( x, a, b );
    cdf = log_normal_cdf ( x, a, b );
    x2 = log_normal_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test100 ( )

/******************************************************************************/
/*
  Purpose:

    TEST100 tests LOG_NORMAL_MEAN, LOG_NORMAL_SAMPLE, LOG_NORMAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST100\n" );
  printf ( "  For the LogNormal PDF:\n" );
  printf ( "  LOG_NORMAL_MEAN computes the mean;\n" );
  printf ( "  LOG_NORMAL_SAMPLE samples;\n" );
  printf ( "  LOG_NORMAL_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !normal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST100 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = log_normal_mean ( a, b);
  variance = log_normal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = log_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test101 ( )

/******************************************************************************/
/*
  Purpose:

    TEST101 tests LOG_SERIES_CDF, LOG_SERIES_CDF_INV, LOG_SERIES_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST101\n" );
  printf ( "  For the Log Series PDF:\n" );
  printf ( "  LOG_SERIES_CDF evaluates the CDF;\n" );
  printf ( "  LOG_SERIES_CDF_INV inverts the CDF.\n" );
  printf ( "  LOG_SERIES_PDF evaluates the PDF;\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !log_series_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST101 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = log_series_sample ( a, &seed );
    pdf = log_series_pdf ( x, a );
    cdf = log_series_cdf ( x, a );
    x2 = log_series_cdf_inv ( cdf, a );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test102 ( )

/******************************************************************************/
/*
  Purpose:

    TEST102 tests LOG_SERIES_CDF, LOG_SERIES_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST102\n" );
  printf ( "  For the Log Series PDF:\n" );
  printf ( "  LOG_SERIES_CDF evaluates the CDF;\n" );
  printf ( "  LOG_SERIES_PDF evaluates the PDF;\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !log_series_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST101 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 1; x <= 10; x++ )
  {
    pdf = log_series_pdf ( x, a );
    cdf = log_series_cdf ( x, a );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test103 ( )

/******************************************************************************/
/*
  Purpose:

    TEST103 tests LOG_SERIES_MEAN, LOG_SERIES_SAMPLE, LOG_SERIES_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST103\n" );
  printf ( "  For the Log Series PDF:\n" );
  printf ( "  LOG_SERIES_MEAN computes the mean;\n" );
  printf ( "  LOG_SERIES_SAMPLE samples;\n" );
  printf ( "  LOG_SERIES_VARIANCE computes the variance.\n" );

  a = 0.25E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a        );

  if ( !log_series_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST103 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = log_series_mean ( a );
  variance = log_series_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = log_series_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test104 ( )

/******************************************************************************/
/*
  Purpose:

    TEST104 tests LOG_UNIFORM_CDF, LOG_UNIFORM_CDF_INV, LOG_UNIFORM_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST104\n" );
  printf ( "  For the Log Uniform PDF:\n" );
  printf ( "  LOG_UNIFORM_CDF evaluates the CDF;\n" );
  printf ( "  LOG_UNIFORM_CDF_INV inverts the CDF.\n" );
  printf ( "  LOG_UNIFORM_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 20.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !log_uniform_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST104 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = log_uniform_sample ( a, b, &seed );
    pdf = log_uniform_pdf ( x, a, b );
    cdf = log_uniform_cdf ( x, a, b );
    x2 = log_uniform_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test105 ( )

/******************************************************************************/
/*
  Purpose:

    TEST105 tests LOG_UNIFORM_MEAN, LOG_UNIFORM_SAMPLE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST105\n" );
  printf ( "  For the Log Uniform PDF:\n" );
  printf ( "  LOG_UNIFORM_MEAN computes the mean;\n" );
  printf ( "  LOG_UNIFORM_SAMPLE samples;\n" );

  a = 2.0;
  b = 20.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !log_uniform_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST105 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = log_uniform_mean ( a, b);

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = log_uniform_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test106 ( )

/******************************************************************************/
/*
  Purpose:

    TEST106 tests LORENTZ_CDF, LORENTZ_CDF_INV, LORENTZ_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST106\n" );
  printf ( "  For the Lorentz PDF:\n" );
  printf ( "  LORENTZ_CDF evaluates the CDF;\n" );
  printf ( "  LORENTZ_CDF_INV inverts the CDF.\n" );
  printf ( "  LORENTZ_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = lorentz_sample ( &seed );
    pdf = lorentz_pdf ( x );
    cdf = lorentz_cdf ( x );
    x2 = lorentz_cdf_inv ( cdf );

    printf ( "  %12g   %12g  %12g  %12g\n", x, pdf, cdf, x2  );
  }

  return;
}
/******************************************************************************/

void test107 ( )

/******************************************************************************/
/*
  Purpose:

    TEST107 tests LORENTZ_MEAN, LORENTZ_SAMPLE, LORENTZ_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST107\n" );
  printf ( "  For the Lorentz PDF:\n" );
  printf ( "  LORENTZ_MEAN computes the mean;\n" );
  printf ( "  LORENTZ_SAMPLE samples;\n" );
  printf ( "  LORENTZ_VARIANCE computes the variance.\n" );

  mean     = lorentz_mean ( );
  variance = lorentz_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = lorentz_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test108 ( )

/******************************************************************************/
/*
  Purpose:

    TEST108 tests MAXWELL_CDF, MAXWELL_CDF_INV, MAXWELL_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST108\n" );
  printf ( "  For the Maxwell PDF:\n" );
  printf ( "  MAXWELL_CDF evaluates the CDF;\n" );
  printf ( "  MAXWELL_CDF_INV inverts the CDF.\n" );
  printf ( "  MAXWELL_PDF evaluates the PDF;\n" );

  a = 2.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !maxwell_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST108 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }
  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = maxwell_sample ( a, &seed );
    pdf = maxwell_pdf ( x, a );
    cdf = maxwell_cdf ( x, a );
    x2 = maxwell_cdf_inv ( cdf, a );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test109 ( )

/******************************************************************************/
/*
  Purpose:

    TEST109 tests MAXWELL_MEAN, MAXWELL_SAMPLE, MAXWELL_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST109\n" );
  printf ( "  For the Maxwell PDF:\n" );
  printf ( "  MAXWELL_MEAN computes the mean;\n" );
  printf ( "  MAXWELL_SAMPLE samples;\n" );
  printf ( "  MAXWELL_VARIANCE computes the variance.\n" );

  a = 2.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a        );

  if ( !maxwell_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST109 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = maxwell_mean ( a );
  variance = maxwell_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = maxwell_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test110 ( )

/******************************************************************************/
/*
  Purpose:

    TEST110 tests MULTINOMIAL_COEF1, MULTINOMIAL_COEF2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define MAXFACTOR 5

  int factor[MAXFACTOR];
  int i;
  int j;
  int n;
  int ncomb1;
  int ncomb2;
  int nfactor;

  printf ( "\n" );
  printf ( "TEST110\n" );
  printf ( "  MULTINOMIAL_COEF1 computes multinomial\n" );
  printf ( "  coefficients using the Gamma function;\n" );
  printf ( "  MULTINOMIAL_COEF2 computes multinomial\n" );
  printf ( "  coefficients directly.\n" );

  printf ( "\n" );
  printf ( "  Line 10 of the BINOMIAL table:\n" );
  printf ( "\n" );

  n = 10;
  nfactor = 2;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;
    factor[1] = n - i;

    ncomb1 = multinomial_coef1 ( nfactor, factor );

    ncomb2 = multinomial_coef2 ( nfactor, factor );

    printf ( "  %2d  %2d  %5d  %5d\n",
      factor[0], factor[1], ncomb1, ncomb2 );
  }

  printf ( "\n" );
  printf ( "  Level 5 of the TRINOMIAL coefficients:\n" );

  n = 5;
  nfactor = 3;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;

    printf ( "\n" );

    for ( j = 0; j <= n - factor[0]; j++ )
    {
      factor[1] = j;
      factor[2] = n - factor[0] - factor[1];

      ncomb1 = multinomial_coef1 ( nfactor, factor );

      ncomb2 = multinomial_coef2 ( nfactor, factor );

      printf ( "  %2d  %2d  %2d  %5d  %5d\n",
        factor[0], factor[1], factor[2], ncomb1, ncomb2 );
    }
  }

  return;
# undef MAXFACTOR
}
/******************************************************************************/

void test111 ( )

/******************************************************************************/
/*
  Purpose:

    TEST111 tests MULTINOMIAL_MEAN, MULTINOMIAL_SAMPLE, MULTINOMIAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define B 3
# define SAMPLE_NUM 1000

  int a;
  double c[B] = { 0.125, 0.500, 0.375 };
  int i;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  int x[B*SAMPLE_NUM];
  int *xmax;
  int *xmin;
  int *y;

  printf ( "\n" );
  printf ( "TEST111\n" );
  printf ( "  For the Multinomial PDF:\n" );
  printf ( "  MULTINOMIAL_MEAN computes the mean;\n" );
  printf ( "  MULTINOMIAL_SAMPLE samples;\n" );
  printf ( "  MULTINOMIAL_VARIANCE computes the variance;\n" );

  a = 5;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %d\n", B );
  r8vec_print ( B, c, "  PDF parameter C:" );

  if ( !multinomial_check ( a, B, c ) )
  {
    printf ( "\n" );
    printf ( "TEST111 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = multinomial_mean ( a, B, c );
  variance = multinomial_variance ( a, B, c );
  r8vec_print ( B, mean, "  PDF mean:" );
  r8vec_print ( B, variance, "  PDF variance:" );

  free ( mean );
  free ( variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = multinomial_sample ( a, B, c, &seed );
    for ( i = 0; i < B; i++ )
    {
      x[i+j*B] = y[i];
    }
    free ( y );
  }

  mean = i4row_mean ( B, SAMPLE_NUM, x );
  variance = i4row_variance ( B, SAMPLE_NUM, x );
  xmax = i4row_max ( B, SAMPLE_NUM, x );
  xmin = i4row_min ( B, SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "\n" );
  printf ( "  Component Mean, Variance, Min, Max:\n" );
  printf ( "\n" );

  for ( i = 0; i < B; i++ )
  {
    printf ( "  %6d  %12g  %12g  %12d  %12d\n", 
      i+1, mean[i], variance[i], xmin[i], xmax[i] );
  }

  free ( mean );
  free ( variance );
  free ( xmax );
  free ( xmin );

  return;
# undef B
# undef SAMPLE_NUM
}
/******************************************************************************/

void test112 ( )

/******************************************************************************/
/*
  Purpose:

    TEST112 tests MULTINOMIAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define B 3

  int a;
  double c[B] = { 0.1, 0.5, 0.4 };
  int i;
  double pdf;
  int x[B] = { 0, 2, 3 };

  printf ( "\n" );
  printf ( "TEST112\n" );
  printf ( "  For the Multinomial PDF:\n" );
  printf ( "  MULTINOMIAL_PDF evaluates the PDF;\n" );

  a = 5;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %d\n", B );
  r8vec_print ( B, c, "  PDF parameter C:" );

  if ( !multinomial_check ( a, B, c ) )
  {
    printf ( "\n" );
    printf ( "TEST112 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  i4vec_print ( B, x, "  PDF argument X:" );

  pdf = multinomial_pdf ( x, a, B, c );

  printf ( "\n" );
  printf ( "  PDF value = %g\n", pdf);

  return;
# undef B
# undef SAMPLE_NUM
}
/******************************************************************************/

void test113 ( )

/******************************************************************************/
/*
  Purpose:

    TEST113 tests NAKAGAMI_CDF, NAKAGAMI_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  double pdf;
  double x;

  printf ( "\n" );
  printf ( "TEST113\n" );
  printf ( "  For the Nakagami PDF:\n" );
  printf ( "  NAKAGAMI_CDF evaluates the CDF;\n" );
  printf ( "  NAKAGAMI_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !nakagami_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST113 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  x = 1.25;
  pdf = nakagami_pdf ( x, a, b, c );
  cdf = nakagami_cdf ( x, a, b, c );

  printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );

  return;
}
/******************************************************************************/

void test114 ( )

/******************************************************************************/
/*
  Purpose:

    TEST114 tests NAKAGAMI_MEAN, NAKAGAMI_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double mean;
  double variance;
  double x;

  printf ( "\n" );
  printf ( "TEST114\n" );
  printf ( "  For the Nakagami PDF:\n" );
  printf ( "  NAKAGAMI_MEAN evaluates the mean;\n" );
  printf ( "  NAKAGAMI_VARIANCE evaluates the variance;\n" );

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !nakagami_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST114 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = nakagami_mean ( a, b, c );
  variance = nakagami_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =      %g\n", mean     );
  printf ( "  PDF variance =  %g\n", variance );

  return;
}
/******************************************************************************/

void test1145 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1145 tests NEGATIVE_BINOMIAL_CDF, NEGATIVE_BINOMIAL_CDF_INV, NEGATIVE_BINOMIAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST1145\n" );
  printf ( "  For the Negative Binomial PDF:\n" );
  printf ( "  NEGATIVE_BINOMIAL_CDF evaluates the CDF;\n" );
  printf ( "  NEGATIVE_BINOMIAL_CDF_INV inverts the CDF.\n" );
  printf ( "  NEGATIVE_BINOMIAL_PDF evaluates the PDF;\n" );

  a = 2;
  b = 0.25;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !negative_binomial_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST1145 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = negative_binomial_sample ( a, b, &seed );
    pdf = negative_binomial_pdf ( x, a, b );
    cdf = negative_binomial_cdf ( x, a, b );
    x2 = negative_binomial_cdf_inv ( cdf, a, b );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test1146 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1146 tests NEGATIVE_BINOMIAL_MEAN, NEGATIVE_BINOMIAL_SAMPLE, NEGATIVE_BINOMIAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST1146\n" );
  printf ( "  For the Negative Binomial PDF:\n" );
  printf ( "  NEGATIVE_BINOMIAL_MEAN computes the mean;\n" );
  printf ( "  NEGATIVE_BINOMIAL_SAMPLE samples;\n" );
  printf ( "  NEGATIVE_BINOMIAL_VARIANCE computes the variance;\n" );

  a = 2;
  b = 0.75;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !negative_binomial_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST1146 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = negative_binomial_mean ( a, b );
  variance = negative_binomial_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = negative_binomial_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test115 ( )

/******************************************************************************/
/*
  Purpose:

    TEST115 tests NORMAL_01_CDF, NORMAL_01_CDF_INV, NORMAL_01_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST115\n" );
  printf ( "  For the Normal 01 PDF:\n" );
  printf ( "  NORMAL_01_CDF evaluates the CDF;\n" );
  printf ( "  NORMAL_01_CDF_INV inverts the CDF.\n" );
  printf ( "  NORMAL_01_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_01_sample ( &seed );
    pdf = normal_01_pdf ( x );
    cdf = normal_01_cdf ( x );
    x2 = normal_01_cdf_inv ( cdf );

    printf ( "  %24.16g  %12.6g  %12.6g  %24.16g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test116 ( )

/******************************************************************************/
/*
  Purpose:

    TEST116 tests NORMAL_01_MEAN, NORMAL_01_SAMPLE, NORMAL_01_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST116\n" );
  printf ( "  For the Normal 01 PDF:\n" );
  printf ( "  NORMAL_01_MEAN computes the mean;\n" );
  printf ( "  NORMAL_01_SAMPLE samples;\n" );
  printf ( "  NORMAL_01_VARIANCE computes the variance;\n" );

  mean = normal_01_mean ( );
  variance = normal_01_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = normal_01_sample ( &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test117 ( )

/******************************************************************************/
/*
  Purpose:

    TEST117 tests NORMAL_CDF, NORMAL_CDF_INV, NORMAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST117\n" );
  printf ( "  For the Normal PDF:\n" );
  printf ( "  NORMAL_CDF evaluates the CDF;\n" );
  printf ( "  NORMAL_CDF_INV inverts the CDF.\n" );
  printf ( "  NORMAL_PDF evaluates the PDF;\n" );

  a = 100.0;
  b = 15.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !normal_check ( a, b ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TEST117 - Fatal error!\n" );
    fprintf ( stderr, "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_sample ( a, b, &seed );
    pdf = normal_pdf ( x, a, b );
    cdf = normal_cdf ( x, a, b );
    x2 = normal_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test118 ( )

/******************************************************************************/
/*
  Purpose:

    TEST118 tests NORMAL_MEAN, NORMAL_SAMPLE, NORMAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST118\n" );
  printf ( "  For the Normal PDF:\n" );
  printf ( "  NORMAL_MEAN computes the mean;\n" );
  printf ( "  NORMAL_SAMPLE samples;\n" );
  printf ( "  NORMAL_VARIANCE computes the variance;\n" );

  a = 100.0;
  b = 15.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !normal_check ( a, b ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TEST118 - Fatal error!\n" );
    fprintf ( stderr, "  The parameters are not legal.\n" );
    return;
  }

  mean = normal_mean ( a, b );
  variance = normal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test1184 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1184 tests NORMAL_TRUNCATED_AB_CDF, NORMAL_TRUNCATED_AB_CDF_INV, NORMAL_TRUNCATED_AB_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double mu;
  double pdf;
  double s;
  int seed;
  double x;
  double x2;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1184\n" );
  printf ( "  For the Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_AB_CDF evaluates the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_AB_CDF_INV inverts the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_AB_PDF evaluates the PDF.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g]\n", a, b );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_truncated_ab_sample ( mu, s, a, b, &seed );

    pdf = normal_truncated_ab_pdf ( x, mu, s, a, b );

    cdf = normal_truncated_ab_cdf ( x, mu, s, a, b );

    x2 = normal_truncated_ab_cdf_inv ( cdf, mu, s, a, b );

    printf ( "  %14.6g%14.6g%14.6g%14.6g\n", x, pdf, cdf, x2 );
  }
  return;
}
/******************************************************************************/

void test1185 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1185 tests NORMAL_TRUNCATED_AB_MEAN, NORMAL_TRUNCATED_AB_SAMPLE, NORMAL_TRUNCATED_AB_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  double mean;
  double mu;
  double s;
  int sample_num = 1000;
  int seed;
  double variance;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1185\n" );
  printf ( "  For the Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_AB_MEAN computes the mean;\n" );
  printf ( "  NORMAL_TRUNCATED_AB_SAMPLE samples;\n" );
  printf ( "  NORMAL_TRUNCATED_AB_VARIANCE computes the variance.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,%g]\n", a, b );

  mean = normal_truncated_ab_mean ( mu, s, a, b );

  variance = normal_truncated_ab_variance ( mu, s, a, b );

  printf ( "\n" );
  printf ( "  PDF mean      =               %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_truncated_ab_sample ( mu, s, a, b, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  variance = r8vec_variance ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void test1186 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1186 tests NORMAL_TRUNCATED_A_CDF, NORMAL_TRUNCATED_A_CDF_INV, NORMAL_TRUNCATED_A_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double mu;
  double pdf;
  double s;
  int seed;
  double x;
  double x2;

  a = 50.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1186\n" );
  printf ( "  For the Lower Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_A_CDF evaluates the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_A_CDF_INV inverts the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_A_PDF evaluates the PDF.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo]\n", a );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_truncated_a_sample ( mu, s, a, &seed );

    pdf = normal_truncated_a_pdf ( x, mu, s, a );

    cdf = normal_truncated_a_cdf ( x, mu, s, a );

    x2 = normal_truncated_a_cdf_inv ( cdf, mu, s, a );

    printf ( "  %14.6g%14.6g%14.6g%14.6g\n", x, pdf, cdf, x2 );
  }
  return;
}
/******************************************************************************/

void test1187 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1187 tests NORMAL_TRUNCATED_A_MEAN, NORMAL_TRUNCATED_A_SAMPLE, NORMAL_TRUNCATED_A_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double mean;
  double mu;
  double s;
  int sample_num = 1000;
  int seed;
  double variance;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1187\n" );
  printf ( "  For the Lower Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_A_MEAN computes the mean;\n" );
  printf ( "  NORMAL_TRUNCATED_A_SAMPLE samples;\n" );
  printf ( "  NORMAL_TRUNCATED_A_VARIANCE computes the variance.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [%g,+oo]\n", a );

  mean = normal_truncated_a_mean ( mu, s, a );

  variance = normal_truncated_a_variance ( mu, s, a );

  printf ( "\n" );
  printf ( "  PDF mean      =               %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_truncated_a_sample ( mu, s, a, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  variance = r8vec_variance ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void test1188 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1188 tests NORMAL_TRUNCATED_B_CDF, NORMAL_TRUNCATED_B_CDF_INV, NORMAL_TRUNCATED_B_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double b;
  double cdf;
  int i;
  double mu;
  double pdf;
  double s;
  int seed;
  double x;
  double x2;

  b = 150.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1188\n" );
  printf ( "  For the Upper Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_B_CDF evaluates the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_B_CDF_INV inverts the CDF.\n" );
  printf ( "  NORMAL_TRUNCATED_B_PDF evaluates the PDF.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [-oo,%g]\n", b );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_truncated_b_sample ( mu, s, b, &seed );

    pdf = normal_truncated_b_pdf ( x, mu, s, b );

    cdf = normal_truncated_b_cdf ( x, mu, s, b );

    x2 = normal_truncated_b_cdf_inv ( cdf, mu, s, b );

    printf ( "  %14.6g%14.6g%14.6g%14.6g\n", x, pdf, cdf, x2 );
  }
  return;
}
/******************************************************************************/

void test1189 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1189 tests NORMAL_TRUNCATED_B_MEAN, NORMAL_TRUNCATED_B_SAMPLE, NORMAL_TRUNCATED_B_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2013

  Author:

    John Burkardt
*/
{
  double b;
  int i;
  double mean;
  double mu;
  double s;
  int sample_num = 1000;
  int seed;
  double variance;
  double *x;
  double xmax;
  double xmin;

  b = 150.0;
  mu = 100.0;
  s = 25.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST1189\n" );
  printf ( "  For the Upper Truncated Normal PDF:\n" );
  printf ( "  NORMAL_TRUNCATED_B_MEAN computes the mean;\n" );
  printf ( "  NORMAL_TRUNCATED_B_SAMPLE samples;\n" );
  printf ( "  NORMAL_TRUNCATED_B_VARIANCE computes the variance.\n" );
  printf ( "\n" );
  printf ( "  The parent normal distribution has\n" );
  printf ( "    mean =               %g\n", mu );
  printf ( "    standard deviation = %g\n", s );
  printf ( "  The parent distribution is truncated to\n" );
  printf ( "  the interval [-oo,%g]\n", b );

  mean = normal_truncated_b_mean ( mu, s, b );

  variance = normal_truncated_b_variance ( mu, s, b );

  printf ( "\n" );
  printf ( "  PDF mean      =               %g\n", mean );
  printf ( "  PDF variance =                %g\n", variance );

  x = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_truncated_b_sample ( mu, s, b, &seed );
  }

  mean = r8vec_mean ( sample_num, x );
  variance = r8vec_variance ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", sample_num );
  printf ( "  Sample mean =     %g\n", mean );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax );
  printf ( "  Sample minimum =  %g\n", xmin );

  free ( x );

  return;
}
/******************************************************************************/

void test119 ( )

/******************************************************************************/
/*
  Purpose:

    TEST119 tests PARETO_CDF, PARETO_CDF_INV, PARETO_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST119\n" );
  printf ( "  For the Pareto PDF:\n" );
  printf ( "  PARETO_CDF evaluates the CDF;\n" );
  printf ( "  PARETO_CDF_INV inverts the CDF.\n" );
  printf ( "  PARETO_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !pareto_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST119 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = pareto_sample ( a, b, &seed );
    pdf = pareto_pdf ( x, a, b );
    cdf = pareto_cdf ( x, a, b );
    x2 = pareto_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test120 ( )

/******************************************************************************/
/*
  Purpose:

    TEST120 tests PARETO_MEAN, PARETO_SAMPLE, PARETO_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST120\n" );
  printf ( "  For the Pareto PDF:\n" );
  printf ( "  PARETO_MEAN computes the mean;\n" );
  printf ( "  PARETO_SAMPLE samples;\n" );
  printf ( "  PARETO_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !pareto_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST120 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = pareto_mean ( a, b );
  variance = pareto_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = pareto_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test123 ( )

/******************************************************************************/
/*
  Purpose:

    TEST123 tests PEARSON_05_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double pdf;
  double x;

  printf ( "\n" );
  printf ( "TEST123\n" );
  printf ( "  For the Pearson 05 PDF:\n" );
  printf ( "  PEARSON_05_PDF evaluates the PDF.\n" );

  x = 5.0;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a   );
  printf ( "  PDF parameter B = %g\n", b   );
  printf ( "  PDF parameter C = %g\n", c   );

  if ( !pearson_05_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST123 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  pdf = pearson_05_pdf ( x, a, b, c );

  printf ( "\n" );
  printf ( "  PDF argument X =  %g\n", x   );
  printf ( "  PDF value =       %g\n", pdf );

  return;
}
/******************************************************************************/

void test124 ( )

/******************************************************************************/
/*
  Purpose:

    TEST124 tests PLANCK_PDF, PLANCK_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  printf ( "\n" );
  printf ( "TEST124\n" );
  printf ( "  For the Planck PDF:\n" );
  printf ( "  PLANCK_PDF evaluates the PDF.\n" );
  printf ( "  PLANCK_SAMPLE samples the PDF.\n" );

  a = 2.0E+00;
  b = 3.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a );
  printf ( "  PDF parameter B = %g\n", b );

  if ( !planck_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST124 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = planck_sample ( a, b, &seed );

    pdf = planck_pdf ( x, a, b );

    printf ( "  %12g  %12g\n", x, pdf );
  }

  return;
}
/******************************************************************************/

void test125 ( )

/******************************************************************************/
/*
  Purpose:

    TEST125 tests PLANCK_MEAN, PLANCK_SAMPLE, PLANCK_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST125\n" );
  printf ( "  For the Planck PDF:\n" );
  printf ( "  PLANCK_MEAN computes the mean;\n" );
  printf ( "  PLANCK_SAMPLE samples;\n" );
  printf ( "  PLANCK_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !planck_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST125 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = planck_mean ( a, b );
  variance = planck_variance ( a, b );

  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = planck_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test126 ( )

/******************************************************************************/
/*
  Purpose:

    TEST126 tests POISSON_CDF, POISSON_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "TEST126:\n" );
  printf ( "  POISSON_CDF evaluates the cumulative distribution\n" );
  printf ( "  function for the discrete Poisson probability\n" );
  printf ( "  density function.\n" );
  printf ( "  POISSON_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  A is the expected mean number of successes per unit time;\n" );
  printf ( "  X is the number of successes;\n" );
  printf ( "  POISSON_CDF is the probability of having up to X\n" );
  printf ( "  successes in unit time.\n" );
  printf ( "\n" );
  printf ( "   A          X   Exact F     POISSON_CDF(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = poisson_cdf ( x, a );

    printf ( "  %8d  %8d  %16g  %12g\n", a, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test127 ( )

/******************************************************************************/
/*
  Purpose:

    TEST127 tests POISSON_CDF, POISSON_CDF_INV, POISSON_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST127\n" );
  printf ( "  For the Poisson PDF:\n" );
  printf ( "  POISSON_CDF evaluates the CDF;\n" );
  printf ( "  POISSON_CDF_INV inverts the CDF.\n" );
  printf ( "  POISSON_PDF evaluates the PDF;\n" );

  a = 10.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !poisson_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST127 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = poisson_sample ( a, &seed );
    pdf = poisson_pdf ( x, a );
    cdf = poisson_cdf ( x, a );
    x2 = poisson_cdf_inv ( cdf, a );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test128 ( )

/******************************************************************************/
/*
  Purpose:

    TEST128 tests POISSON_MEAN, POISSON_SAMPLE, POISSON_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST128\n" );
  printf ( "  For the Poisson PDF,\n" );
  printf ( "  POISSON_SAMPLE samples the Poisson PDF.\n" );
  printf ( "  POISSON_SAMPLE samples the Poisson PDF.\n" );
  printf ( "  POISSON_SAMPLE samples the Poisson PDF.\n" );

  a = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );

  if ( !poisson_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST128 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = poisson_mean ( a );
  variance = poisson_variance ( a );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = poisson_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test129 ( )

/******************************************************************************/
/*
  Purpose:

    TEST129 tests POWER_CDF, POWER_CDF_INV, POWER_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST129\n" );
  printf ( "  For the Power PDF:\n" );
  printf ( "  POWER_CDF evaluates the CDF;\n" );
  printf ( "  POWER_CDF_INV inverts the CDF.\n" );
  printf ( "  POWER_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !power_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST129 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = power_sample ( a, b, &seed );
    pdf = power_pdf ( x, a, b );
    cdf = power_cdf ( x, a, b );
    x2 = power_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test130 ( )

/******************************************************************************/
/*
  Purpose:

    TEST130 tests POWER_MEAN, POWER_SAMPLE, POWER_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST130\n" );
  printf ( "  For the Power PDF:\n" );
  printf ( "  POWER_MEAN computes the mean;\n" );
  printf ( "  POWER_SAMPLE samples;\n" );
  printf ( "  POWER_VARIANCE computes the variance;\n" );

  a = 2.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !power_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST130 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = power_mean ( a, b );
  variance = power_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = power_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test1304 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1304 tests QUASIGEOMETRIC_CDF, *_CDF_INV, *_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST1304\n" );
  printf ( "  For the Quasigeometric PDF:\n" );
  printf ( "  QUASIGEOMETRIC_CDF evaluates the CDF;\n" );
  printf ( "  QUASIGEOMETRIC_CDF_INV inverts the CDF.\n" );
  printf ( "  QUASIGEOMETRIC_PDF evaluates the PDF;\n" );

  a = 0.4825;
  b = 0.5893;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a   );
  printf ( "  PDF parameter B = %g\n", b   );

  if ( !quasigeometric_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST1304 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = quasigeometric_sample ( a, b, &seed );
    pdf = quasigeometric_pdf ( x, a, b );
    cdf = quasigeometric_cdf ( x, a, b );
    x2 = quasigeometric_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test1306 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1306 tests QUASIGEOMETRIC_MEAN, *_SAMPLE, *_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int *x;
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST1306\n" );
  printf ( "  For the Quasigeometric PDF:\n" );
  printf ( "  QUASIGEOMETRIC_MEAN computes the mean;\n" );
  printf ( "  QUASIGEOMETRIC_SAMPLE samples;\n" );
  printf ( "  QUASIGEOMETRIC_VARIANCE computes the variance.\n" );

  a = 0.4825;
  b = 0.5893;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a );
  printf ( "  PDF parameter B = %g\n", b );

  if ( !quasigeometric_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST1306 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = quasigeometric_mean ( a, b );
  variance = quasigeometric_variance ( a, b );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  x = ( int * ) malloc ( SAMPLE_NUM * sizeof ( int ) );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = quasigeometric_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  free ( x );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test131 ( )

/******************************************************************************/
/*
  Purpose:

    TEST131 tests RAYLEIGH_CDF, RAYLEIGH_CDF_INV, RAYLEIGH_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST131\n" );
  printf ( "  For the Rayleigh PDF:\n" );
  printf ( "  RAYLEIGH_CDF evaluates the CDF;\n" );
  printf ( "  RAYLEIGH_CDF_INV inverts the CDF.\n" );
  printf ( "  RAYLEIGH_PDF evaluates the PDF;\n" );

  a = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !rayleigh_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST131 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = rayleigh_sample ( a, &seed );
    pdf = rayleigh_pdf ( x, a );
    cdf = rayleigh_cdf ( x, a );
    x2 = rayleigh_cdf_inv ( cdf, a );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test132 ( )

/******************************************************************************/
/*
  Purpose:

    TEST132 tests RAYLEIGH_MEAN, RAYLEIGH_SAMPLE, RAYLEIGH_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST132\n" );
  printf ( "  For the Rayleigh PDF:\n" );
  printf ( "  RAYLEIGH_MEAN computes the mean;\n" );
  printf ( "  RAYLEIGH_SAMPLE samples;\n" );
  printf ( "  RAYLEIGH_VARIANCE computes the variance.\n" );

  a = 2.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a        );

  if ( !rayleigh_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST132 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = rayleigh_mean ( a );
  variance = rayleigh_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = rayleigh_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test133 ( )

/******************************************************************************/
/*
  Purpose:

    TEST133 tests RECIPROCAL_CDF, RECIPROCAL_CDF_INV, RECIPROCAL_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST133\n" );
  printf ( "  For the Reciprocal PDF:\n" );
  printf ( "  RECIPROCAL_CDF evaluates the CDF;\n" );
  printf ( "  RECIPROCAL_CDF_INV inverts the CDF.\n" );
  printf ( "  RECIPROCAL_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !reciprocal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST133 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = reciprocal_sample ( a, b, &seed );
    pdf = reciprocal_pdf ( x, a, b );
    cdf = reciprocal_cdf ( x, a, b );
    x2 = reciprocal_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test134 ( )

/******************************************************************************/
/*
  Purpose:

    TEST134 tests RECIPROCAL_MEAN, RECIPROCAL_SAMPLE, RECIPROCAL_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST134\n" );
  printf ( "  For the Reciprocal PDF:\n" );
  printf ( "  RECIPROCAL_MEAN computes the mean;\n" );
  printf ( "  RECIPROCAL_SAMPLE samples;\n" );
  printf ( "  RECIPROCAL_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !reciprocal_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST134 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = reciprocal_mean ( a, b );
  variance = reciprocal_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = reciprocal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test1341 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1341 checks RIBESL against BESSEL_IX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define NB_MAX 10

  double alpha;
  double alpha_frac;
  double b[NB_MAX];
  double fx;
  double fx2;
  int ize;
  int n_data;
  int nb;
  int ncalc;
  double x;

  printf ( "\n" );
  printf ( "TEST1341:\n" );
  printf ( "  RIBESL computes values of Bessel functions\n" );
  printf ( "  of NONINTEGER order.\n" );
  printf ( "  BESSEL_IX_VALUES returns selected values of the\n" );
  printf ( "  Bessel function In for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      ALPHA         X             FX                        FX2\n" );
  printf ( "                                  (table)                   (RIBESL)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_ix_values ( &n_data, &alpha, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    ize = 1;

    nb = ( int ) alpha + 1;

    if ( NB_MAX < nb )
    {
      printf ( "  [Skipping calculation, NB_MAX too small.]\n" );
      continue;
    }

    alpha_frac = alpha - ( double ) ( ( int ) alpha );

    ncalc = ribesl ( x, alpha_frac, nb, ize, b );

    fx2 = b[nb-1];

    printf ( "  %12g  %12g  %12g  %12g\n", alpha, x, fx, fx2 );
  }

  return;
# undef NB_MAX
}
/******************************************************************************/

void test1342 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1342 tests RUNS_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  double pdf;
  double pdf_total;
  int r;

  printf ( "\n" );
  printf ( "TEST1342\n" );
  printf ( "  For the RUNS PDF:\n" );
  printf ( "  RUNS_PDF evaluates the PDF;\n" );
  printf ( "\n" );
  printf ( "  M is the number of symbols of one kind,\n" );
  printf ( "  N is the number of symbols of the other kind,\n" );
  printf ( "  R is the number of runs (sequences of one symbol)\n" );
  printf ( "\n" );
  printf ( "         M         N         R      PDF\n" );
  printf ( "\n" );

  m = 6;

  for ( n = 0; n <= 9; n++ )
  {
    printf ( "\n" );
    pdf_total = 0.0;

    for ( r = 1; r <= 2 * i4_min ( m, n ) + 2; r++ )
    {
      pdf = runs_pdf ( m, n, r );

      printf ( "  %8d  %8d  %8d  %14g\n", m, n, r, pdf );

      pdf_total = pdf_total + pdf;
    }

    printf ( "  %8d                      %14g\n", m, pdf_total );

  }

  return;
}
/******************************************************************************/

void test1344 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1344 tests RUNS_MEAN, RUNS_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST1344\n" );
  printf ( "  For the RUNS PDF:\n" );
  printf ( "  RUNS_MEAN computes the mean;\n" );
  printf ( "  RUNS_VARIANCE computes the variance\n" );

  m = 10;
  n = 5;

  printf ( "\n" );
  printf ( "  PDF parameter M = %d\n", m );
  printf ( "  PDF parameter N = %d\n", n );

  mean = runs_mean ( m, n );
  variance = runs_variance ( m, n );

  printf ( "  PDF mean =        %g\n", mean );
  printf ( "  PDF variance =    %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = runs_sample ( m, n, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test135 ( )

/******************************************************************************/
/*
  Purpose:

    TEST135 tests SECH_CDF, SECH_CDF_INV, SECH_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST135\n" );
  printf ( "  For the Sech PDF:\n" );
  printf ( "  SECH_CDF evaluates the CDF;\n" );
  printf ( "  SECH_CDF_INV inverts the CDF.\n" );
  printf ( "  SECH_PDF evaluates the PDF;\n" );

  a = 3.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !sech_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST135 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = sech_sample ( a, b, &seed );
    pdf = sech_pdf ( x, a, b );
    cdf = sech_cdf ( x, a, b );
    x2 = sech_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test136 ( )

/******************************************************************************/
/*
  Purpose:

    TEST136 tests SECH_MEAN, SECH_SAMPLE, SECH_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST136\n" );
  printf ( "  For the Sech PDF:\n" );
  printf ( "  SECH_MEAN computes the mean;\n" );
  printf ( "  SECH_SAMPLE samples;\n" );
  printf ( "  SECH_VARIANCE computes the variance;\n" );

  a = 3.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !sech_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST136 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = sech_mean ( a, b );
  variance = sech_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = sech_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test137 ( )

/******************************************************************************/
/*
  Purpose:

    TEST137 tests SEMICIRCULAR_CDF, SEMICIRCULAR_CDF_INV, SEMICIRCULAR_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST137\n" );
  printf ( "  For the Semicircular PDF:\n" );
  printf ( "  SEMICIRCULAR_CDF evaluates the CDF;\n" );
  printf ( "  SEMICIRCULAR_CDF_INV inverts the CDF.\n" );
  printf ( "  SEMICIRCULAR_PDF evaluates the PDF;\n" );

  a = 3.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !semicircular_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST137 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = semicircular_sample ( a, b, &seed );
    pdf = semicircular_pdf ( x, a, b );
    cdf = semicircular_cdf ( x, a, b );
    x2 = semicircular_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test138 ( )

/******************************************************************************/
/*
  Purpose:

    TEST138 tests SEMICIRCULAR_MEAN, SEMICIRCULAR_SAMPLE and SEMICIRCULAR_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST138\n" );
  printf ( "  For the Semicircular PDF:\n" );
  printf ( "  SEMICIRCULAR_MEAN computes the mean;\n" );
  printf ( "  SEMICIRCULAR_SAMPLE samples;\n" );
  printf ( "  SEMICIRCULAR_VARIANCE computes the variance;\n" );

  a = 3.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !semicircular_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST138 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = semicircular_mean ( a, b );
  variance = semicircular_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = semicircular_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test139 ( )

/******************************************************************************/
/*
  Purpose:

    TEST139 tests STUDENT_CDF and STUDENT_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST139:\n" );
  printf ( "  STUDENT_CDF evaluates the cumulative density function\n" );
  printf ( "  for the Student's T PDF.\n" );
  printf ( "  STUDENT_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   A      B      C     X       Exact F       STUDENT_CDF(A,B,C,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    a = 0.0;
    b = 1.0;

    fx2 = student_cdf ( x, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g  %12g  %12g\n", a, b, c, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test140 ( )

/******************************************************************************/
/*
  Purpose:

    TEST140 tests STUDENT_CDF, STUDENT_PDF and STUDENT_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  printf ( "\n" );
  printf ( "TEST140\n" );
  printf ( "  For the Student PDF:\n" );
  printf ( "  STUDENT_CDF evaluates the CDF;\n" );
  printf ( "  STUDENT_PDF evaluates the PDF;\n" );
  printf ( "  STUDENT_SAMPLE samples the PDF;\n" );

  a = 0.5;
  b = 2.0;
  c = 6.0;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a   );
  printf ( "  PDF parameter B = %g\n", b   );
  printf ( "  PDF parameter C = %g\n", c   );

  if ( !student_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST140 - Fatal error!\n" );
    printf ( "  The parameter values are illegal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = student_sample ( a, b, c, &seed );
    pdf = student_pdf ( x, a, b, c );
    cdf = student_cdf ( x, a, b, c );

    printf ( "  %12g  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test141 ( )

/******************************************************************************/
/*
  Purpose:

    TEST141 tests STUDENT_MEAN, STUDENT_SAMPLE and STUDENT_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST141\n" );
  printf ( "  For the Student PDF:\n" );
  printf ( "  STUDENT_MEAN evaluates the mean;\n" );
  printf ( "  STUDENT_SAMPLE samples the PDF;\n" );
  printf ( "  STUDENT_VARIANCE computes the variance;\n" );

  a = 0.5;
  b = 2.0;
  c = 6.0;

  printf ( "\n" );
  printf ( "  PDF parameter A = %g\n", a   );
  printf ( "  PDF parameter B = %g\n", b   );
  printf ( "  PDF parameter C = %g\n", c   );

  if ( !student_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST141 - Fatal error!\n" );
    printf ( "  The parameter values are illegal.\n" );
    return;
  }

  mean = student_mean ( a, b, c );
  variance = student_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = student_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );


  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test142 ( )

/******************************************************************************/
/*
  Purpose:

    TEST142 tests STUDENT_NONCENTRAL_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double cdf;
  int idf;
  double x;

  printf ( "\n" );
  printf ( "TEST142\n" );
  printf ( "  For the Noncentral Student PDF:\n" );
  printf ( "  STUDENT_NONCENTRAL_CDF evaluates the CDF;\n" );

  x = 0.50;
  idf = 10;
  b = 1.0;

  cdf = student_noncentral_cdf ( x, idf, b );

  printf ( "\n" );
  printf ( "  PDF argument X =              %g\n", x   );
  printf ( "  PDF parameter IDF =           %d\n", idf );
  printf ( "  PDF parameter B =             %g\n", b   );
  printf ( "  CDF value =                   %g\n", cdf );

  return;
}
/******************************************************************************/

void test1425 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1425 tests TFN, OWEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double h;
  int n_data;
  double t;
  double t2;

  printf ( "\n" );
  printf ( "TEST1425\n" );
  printf ( "  TFN evaluates Owen's T function;\n" );
  printf ( "  OWEN_VALUES returns some exact values;\n" );

  printf ( "\n" );
  printf ( "      H             A           T(H,A)          Exact\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t );

    if ( n_data <= 0 )
    {
      break;
    }

    t2 = tfn ( h, a );

    printf ( "  %12g  %12g  %12g  %12g\n", h, a, t2, t );
  }

  return;
}
/******************************************************************************/

void test143 ( )

/******************************************************************************/
/*
  Purpose:

    TEST143 tests TRIANGLE_CDF, TRIANGLE_CDF_INV, TRIANGLE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST143\n" );
  printf ( "  For the Triangle PDF:\n" );
  printf ( "  TRIANGLE_CDF evaluates the CDF;\n" );
  printf ( "  TRIANGLE_CDF_INV inverts the CDF.\n" );
  printf ( "  TRIANGLE_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 3.0;
  c = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !triangle_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST143 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = triangle_sample ( a, b, c, &seed );
    pdf = triangle_pdf ( x, a, b, c );
    cdf = triangle_cdf ( x, a, b, c );
    x2 = triangle_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test144 ( )

/******************************************************************************/
/*
  Purpose:

    TEST144 tests TRIANGLE_MEAN, TRIANGLE_SAMPLE, TRIANGLE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST144\n" );
  printf ( "  For the Triangle PDF:\n" );
  printf ( "  TRIANGLE_MEAN computes the mean;\n" );
  printf ( "  TRIANGLE_SAMPLE samples;\n" );
  printf ( "  TRIANGLE_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 3.0;
  c = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =        %g\n", a );
  printf ( "  PDF parameter B =        %g\n", b );
  printf ( "  PDF parameter C =        %g\n", c );

  if ( !triangle_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST144 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = triangle_mean ( a, b, c );
  variance = triangle_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = triangle_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test145 ( )

/******************************************************************************/
/*
  Purpose:

    TEST145 tests TRIANGULAR_CDF, TRIANGULAR_CDF_INV, TRIANGULAR_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST145\n" );
  printf ( "  For the Triangular PDF:\n" );
  printf ( "  TRIANGULAR_CDF evaluates the CDF;\n" );
  printf ( "  TRIANGULAR_CDF_INV inverts the CDF.\n" );
  printf ( "  TRIANGULAR_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !triangular_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST145 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = triangular_sample ( a, b, &seed );
    pdf = triangular_pdf ( x, a, b );
    cdf = triangular_cdf ( x, a, b );
    x2 = triangular_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test146 ( )

/******************************************************************************/
/*
  Purpose:

    TEST146 tests TRIANGULAR_MEAN, TRIANGULAR_SAMPLE, TRIANGULAR_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST146\n" );
  printf ( "  For the Triangular PDF:\n" );
  printf ( "  TRIANGULAR_MEAN computes the mean;\n" );
  printf ( "  TRIANGULAR_SAMPLE samples;\n" );
  printf ( "  TRIANGULAR_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !triangular_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST146 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = triangular_mean ( a, b );
  variance = triangular_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = triangular_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test147 ( )

/******************************************************************************/
/*
  Purpose:

    TEST147 tests UNIFORM_01_ORDER_SAMPLE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  int n;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST147\n" );
  printf ( "  For the Uniform 01 Order PDF:\n" );
  printf ( "  UNIFORM_ORDER_SAMPLE samples.\n" );

  n = 10;
  x = uniform_01_order_sample ( n, &seed );

  r8vec_print ( n, x, "  Ordered sample:" );

  free ( x );

  return;
}
/******************************************************************************/

void test148 ( )

/******************************************************************************/
/*
  Purpose:

    TEST148 tests UNIFORM_NSPHERE_SAMPLE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n;
  int seed = 123456789;
  double *x;

  n = 3;

  printf ( "\n" );
  printf ( "TEST148\n" );
  printf ( "  For the Uniform PDF on the N-Sphere:\n" );
  printf ( "  UNIFORM_NSPHERE_SAMPLE samples.\n" );

  printf ( "\n" );
  printf ( "  Dimension N of sphere = %g\n", n );
  printf ( "\n" );
  printf ( "  Points on the sphere:\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_nsphere_sample ( n, &seed );
    printf ( "  %6d", i );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%12g", x[j] );
    }
    printf ( "\n" );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test1485 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1485 tests UNIFORM_01_CDF, UNIFORM_01_CDF_INV, UNIFORM_01_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST1485\n" );
  printf ( "  For the Uniform 01 PDF:\n" );
  printf ( "  UNIFORM_01_CDF evaluates the CDF;\n" );
  printf ( "  UNIFORM_01_CDF_INV inverts the CDF.\n" );
  printf ( "  UNIFORM_01_PDF evaluates the PDF;\n" );

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_01_sample ( &seed );
    pdf = uniform_01_pdf ( x );
    cdf = uniform_01_cdf ( x );
    x2 = uniform_01_cdf_inv ( cdf );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test1486 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1486 tests UNIFORM_01_MEAN, UNIFORM_01_SAMPLE, UNIFORM_01_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST1486\n" );
  printf ( "  For the Uniform 01 PDF:\n" );
  printf ( "  UNIFORM_01_MEAN computes the mean;\n" );
  printf ( "  UNIFORM_01_SAMPLE samples;\n" );
  printf ( "  UNIFORM_01_VARIANCE computes the variance.\n" );

  mean     = uniform_01_mean ( );
  variance = uniform_01_variance ( );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_01_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test149 ( )

/******************************************************************************/
/*
  Purpose:

    TEST149 tests UNIFORM_CDF, UNIFORM_CDF_INV, UNIFORM_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST149\n" );
  printf ( "  For the Uniform PDF:\n" );
  printf ( "  UNIFORM_CDF evaluates the CDF;\n" );
  printf ( "  UNIFORM_CDF_INV inverts the CDF.\n" );
  printf ( "  UNIFORM_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !uniform_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST149 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_sample ( a, b, &seed );
    pdf = uniform_pdf ( x, a, b );
    cdf = uniform_cdf ( x, a, b );
    x2 = uniform_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test150 ( )

/******************************************************************************/
/*
  Purpose:

    TEST150 tests UNIFORM_MEAN, UNIFORM_SAMPLE, UNIFORM_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST150\n" );
  printf ( "  For the Uniform PDF:\n" );
  printf ( "  UNIFORM_MEAN computes the mean;\n" );
  printf ( "  UNIFORM_SAMPLE samples;\n" );
  printf ( "  UNIFORM_VARIANCE computes the variance;\n" );

  a = 1.0;
  b = 10.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !uniform_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST150 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = uniform_mean ( a, b );
  variance = uniform_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test151 ( )

/******************************************************************************/
/*
  Purpose:

    TEST151 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_CDF_INV, UNIFORM_DISCRETE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST151\n" );
  printf ( "  For the Uniform Discrete PDF:\n" );
  printf ( "  UNIFORM_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  UNIFORM_DISCRETE_CDF_INV inverts the CDF.\n" );
  printf ( "  UNIFORM_DISCRETE_PDF evaluates the PDF;\n" );

  a = 1;
  b = 6;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %d\n", a );
  printf ( "  PDF parameter B =      %d\n", b );

  if ( !uniform_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST151 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_discrete_sample ( a, b, &seed );
    pdf = uniform_discrete_pdf ( x, a, b );
    cdf = uniform_discrete_cdf ( x, a, b );
    x2 = uniform_discrete_cdf_inv ( cdf, a, b );

    printf ( "  %12d  %12g  %12g  %12d\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test152 ( )

/******************************************************************************/
/*
  Purpose:

    TEST152 tests UNIFORM_DISCRETE_MEAN, UNIFORM_DISCRETE_SAMPLE, UNIFORM_DISCRETE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  int a;
  int b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST152\n" );
  printf ( "  For the Uniform Discrete PDF:\n" );
  printf ( "  UNIFORM_DISCRETE_MEAN computes the mean;\n" );
  printf ( "  UNIFORM_DISCRETE_SAMPLE samples;\n" );
  printf ( "  UNIFORM_DISCRETE_VARIANCE computes the variance;\n" );

  a = 1;
  b = 6;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !uniform_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST152 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = uniform_discrete_mean ( a, b );
  variance = uniform_discrete_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_discrete_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test153 ( )

/******************************************************************************/
/*
  Purpose:

    TEST153 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST153\n" );
  printf ( "  For the Uniform Discrete PDF:\n" );
  printf ( "  UNIFORM_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  UNIFORM_DISCRETE_PDF evaluates the PDF;\n" );

  a = 1;
  b = 6;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %d\n", a   );
  printf ( "  PDF parameter B =             %d\n", b   );

  if ( !uniform_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST153 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 0; x <= 10; x++ )
  {
    pdf = uniform_discrete_pdf ( x, a, b );
    cdf = uniform_discrete_cdf ( x, a, b );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test154 ( )

/******************************************************************************/
/*
  Purpose:

    TEST154 tests VON_MISES_CDF, VON_MISES_CDF_INV, VON_MISES_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST154\n" );
  printf ( "  For the Von Mises PDF:\n" );
  printf ( "  VON_MISES_CDF evaluates the CDF;\n" );
  printf ( "  VON_MISES_CDF_INV inverts the CDF.\n" );
  printf ( "  VON_MISES_PDF evaluates the PDF;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !von_mises_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST154 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = von_mises_sample ( a, b, &seed );
    pdf = von_mises_pdf ( x, a, b );
    cdf = von_mises_cdf ( x, a, b );
    x2 = von_mises_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test155 ( )

/******************************************************************************/
/*
  Purpose:

    TEST155 tests VON_MISES_MEAN, VON_MISES_SAMPLE, VON_MISES_CIRCULAR_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST155\n" );
  printf ( "  For the Von Mises PDF:\n" );
  printf ( "  VON_MISES_MEAN computes the mean;\n" );
  printf ( "  VON_MISES_SAMPLE samples;\n" );
  printf ( "  VON_MISES_CIRCULAR_VARIANCE computes the circular variance;\n" );

  a = 1.0;
  b = 2.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !von_mises_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST155 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = von_mises_mean ( a, b );
  variance = von_mises_circular_variance ( a, b );

  printf ( "\n" );
  printf ( "  PDF mean =              %g\n", mean     );
  printf ( "  PDF circular variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = von_mises_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_circular_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =              %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =              %g\n", mean     );
  printf ( "  Sample circular variance = %g\n", variance );
  printf ( "  Sample maximum =           %g\n", xmax     );
  printf ( "  Sample minimum =           %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test1555 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1555 tests VON_MISES_CDF, VON_MISES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1555:\n" );
  printf ( "  VON_MISES_CDF evaluates the von Mises CDF.\n" );
  printf ( "  VON_MISES_CDF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  A is the dominant angle;\n" );
  printf ( "  B is a measure of spread;\n" );
  printf ( "  X is the angle;\n" );
  printf ( "\n" );
  printf ( "      A     B         X   Exact F     Computed F\n" );

  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    von_mises_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = von_mises_cdf ( x, a, b );

    printf ( "  %12g  %12g  %12g  %12g  %12g\n", a, b, x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test156 ( )

/******************************************************************************/
/*
  Purpose:

    TEST156 tests WEIBULL_CDF, WEIBULL_CDF_INV, WEIBULL_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST156\n" );
  printf ( "  For the Weibull PDF:\n" );
  printf ( "  WEIBULL_CDF evaluates the CDF;\n" );
  printf ( "  WEIBULL_CDF_INV inverts the CDF.\n" );
  printf ( "  WEIBULL_PDF evaluates the PDF;\n" );

  a = 2.0;
  b = 3.0;
  c = 4.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !weibull_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST156 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = weibull_sample ( a, b, c, &seed );
    pdf = weibull_pdf ( x, a, b, c );
    cdf = weibull_cdf ( x, a, b, c );
    x2 = weibull_cdf_inv ( cdf, a, b, c );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test157 ( )

/******************************************************************************/
/*
  Purpose:

    TEST157 tests WEIBULL_MEAN, WEIBULL_SAMPLE, WEIBULL_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST157\n" );
  printf ( "  For the Weibull PDF:\n" );
  printf ( "  WEIBULL_MEAN computes the mean;\n" );
  printf ( "  WEIBULL_SAMPLE samples;\n" );
  printf ( "  WEIBULL_VARIANCE computes the variance.\n" );

  a = 2.0;
  b = 3.0;
  c = 4.0;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );
  printf ( "  PDF parameter C =      %g\n", c );

  if ( !weibull_check ( a, b, c ) )
  {
    printf ( "\n" );
    printf ( "TEST157 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = weibull_mean ( a, b, c );
  variance = weibull_variance ( a, b, c );

  printf ( "\n" );
  printf ( "  PDF mean =     %g\n", mean     );
  printf ( "  PDF variance = %g\n", variance );

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = weibull_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;

# undef SAMPLE_NUM
}
/******************************************************************************/

void test158 ( )

/******************************************************************************/
/*
  Purpose:

    TEST158 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_CDF_INV, WEIBULL_DISCRETE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  printf ( "\n" );
  printf ( "TEST158\n" );
  printf ( "  For the Weibull Discrete PDF:\n" );
  printf ( "  WEIBULL_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  WEIBULL_DISCRETE_CDF_INV inverts the CDF.\n" );
  printf ( "  WEIBULL_DISCRETE_PDF evaluates the PDF;\n" );

  a = 0.5;
  b = 1.5;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !weibull_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST158 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF            CDF_INV\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    x = weibull_discrete_sample ( a, b, &seed );
    pdf = weibull_discrete_pdf ( x, a, b );
    cdf = weibull_discrete_cdf ( x, a, b );
    x2 = weibull_discrete_cdf_inv ( cdf, a, b );

    printf ( "  %12g  %12g  %12g  %12g\n", x, pdf, cdf, x2 );
  }

  return;
}
/******************************************************************************/

void test159 ( )

/******************************************************************************/
/*
  Purpose:

    TEST159 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_PDF;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST159\n" );
  printf ( "  For the Weibull Discrete PDF:\n" );
  printf ( "  WEIBULL_DISCRETE_CDF evaluates the CDF;\n" );
  printf ( "  WEIBULL_DISCRETE_PDF evaluates the PDF;\n" );

  a = 0.5;
  b = 1.5;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !weibull_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST159 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 0; x <= 10; x++ )
  {
    pdf = weibull_discrete_pdf ( x, a, b );
    cdf = weibull_discrete_cdf ( x, a, b );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test160 ( )

/******************************************************************************/
/*
  Purpose:

    TEST160 tests WEIBULL_DISCRETE_MEAN, WEIBULL_DISCRETE_SAMPLE, WEIBULL_DISCRETE_VARIANCE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST160\n" );
  printf ( "  For the Weibull Discrete PDF:\n" );
  printf ( "  WEIBULL_DISCRETE_SAMPLE samples;\n" );

  a = 0.5;
  b = 1.5;

  printf ( "\n" );
  printf ( "  PDF parameter A =      %g\n", a );
  printf ( "  PDF parameter B =      %g\n", b );

  if ( !weibull_discrete_check ( a, b ) )
  {
    printf ( "\n" );
    printf ( "TEST160 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = weibull_discrete_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %g\n", xmax     );
  printf ( "  Sample minimum =  %g\n", xmin     );

  return;
# undef SAMPLE_NUM
}
/******************************************************************************/

void test161 ( )

/******************************************************************************/
/*
  Purpose:

    TEST161 tests ZIPF_CDF, ZIPF_PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double cdf;
  double pdf;
  int x;

  printf ( "\n" );
  printf ( "TEST161\n" );
  printf ( "  For the ZIPF PDF:\n" );
  printf ( "  ZIPF_CDF evaluates the CDF;\n" );
  printf ( "  ZIPF_PDF evaluates the PDF;\n" );

  a = 2.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a   );

  if ( !zipf_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST161 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "       X            PDF           CDF\n" );
  printf ( "\n" );

  for ( x = 1; x <= 20; x++ )
  {
    pdf = zipf_pdf ( x, a );
    cdf = zipf_cdf ( x, a );

    printf ( "  %12d  %12g  %12g\n", x, pdf, cdf );
  }

  return;
}
/******************************************************************************/

void test162 ( )

/******************************************************************************/
/*
  Purpose:

    TEST162 tests ZIPF_MEAN, ZIPF_SAMPLE, ZIPF_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2013

  Author:

    John Burkardt
*/
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  printf ( "\n" );
  printf ( "TEST162\n" );
  printf ( "  For the Zipf PDF:\n" );
  printf ( "  ZIPF_MEAN computes the mean;\n" );
  printf ( "  ZIPF_SAMPLE samples;\n" );
  printf ( "  ZIPF_VARIANCE computes the variance.\n" );

  a = 4.0E+00;

  printf ( "\n" );
  printf ( "  PDF parameter A =             %g\n", a        );

  if ( !zipf_check ( a ) )
  {
    printf ( "\n" );
    printf ( "TEST162 - Fatal error!\n" );
    printf ( "  The parameters are not legal.\n" );
    return;
  }

  mean = zipf_mean ( a );
  variance = zipf_variance ( a );

  printf ( "  PDF mean =                    %g\n", mean     );
  printf ( "  PDF variance =                %g\n", variance );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = zipf_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  printf ( "\n" );
  printf ( "  Sample size =     %d\n", SAMPLE_NUM  );
  printf ( "  Sample mean =     %g\n", mean     );
  printf ( "  Sample variance = %g\n", variance );
  printf ( "  Sample maximum =  %d\n", xmax     );
  printf ( "  Sample minimum =  %d\n", xmin     );

  return;
# undef SAMPLE_NUM
}
