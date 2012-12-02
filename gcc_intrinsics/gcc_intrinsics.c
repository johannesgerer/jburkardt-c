# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
void bessel_j0_values ( int *n_data, double *x, double *fx );
void erf_values ( int *n_data, double *x, double *fx );
void erfc_values ( int *n_data, double *x, double *fx );
void gamma_values ( int *n_data, double *x, double *fx );
double r8_abs ( double x );
void test_erf ( void );
void test_erfc ( void );
void test_gamma ( void );
void test_isnan ( void );
void test_j0 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for GCC_INTRINSICS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "GCC_INTRINSICS:\n" );
  printf ( "  Test the GCC intrinsic library.\n" );

  test_erf ( );
  test_erfc ( );
  test_gamma ( );
  test_isnan ( );
  test_j0 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "GCC_INTRINSICS:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_erf ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ERF tests ERF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 April 2008

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST_ERF:\n" );
  printf ( "  Test ERF, which evaluates the error function ERF(X).\n" );
  printf ( "\n" );
  printf ( "      X               ERF(X)          ERF(X)          DIFF\n" );
  printf ( "                     (tabulated)     (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = erf ( x );

    printf ( "  %14e  %14e  %14e  %14e\n", x, fx, fx2, r8_abs ( fx - fx2 ) );
  }
  return;
}
/******************************************************************************/

void test_erfc ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ERFC tests ERFC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 April 2008

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST_ERFC:\n" );
  printf ( "  Test ERFC, which evaluates the complementary error function.\n" );
  printf ( "\n" );
  printf ( "      X               ERFC(X)         ERFC(X)         DIFF\n" );
  printf ( "                     (tabulated)     (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = erfc ( x );

    printf ( "  %14e  %14e  %14e  %14e\n", x, fx, fx2, r8_abs ( fx - fx2 ) );
  }
  return;
}
/******************************************************************************/

void test_gamma ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_GAMMA tests GAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 April 2008

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST_GAMMA:\n" );
  printf ( "   Test GAMMA, which evaluates the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X              GAMMA(X)         GAMMA(X)        DIFF\n" );
  printf ( "                    (tabulated)      (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = gamma ( x );

    printf ( "  %14e  %14e  %14e  %14e\n", x, fx, fx2, r8_abs ( fx - fx2 ) );
  }
  return;
}
/******************************************************************************/

void test_isnan ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ISNAN tests ISNAN.

  Discussion:

    ISNAN ( X ) is a logical function which is TRUE if the number X
    is a "NaN", or "Not a Number".  This is a special value set to
    indicate that an illegal argument was input to a function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2010

  Author:

    John Burkardt
*/
{
  int i;
  int n = 4;
  float x[4] = { - 1.0, 0.0, + 1.0, + 2.0 };

  printf ( "\n" );
  printf ( "TEST_ISNAN\n" );
  printf ( "  ISNAN(X) is TRUE if X is \"Not a Number\".\n" );
  printf ( "\n" );
  printf ( "   Function  -1  0 +1 +2\n" );
  printf ( "\n" );
  
  printf ( "  ACOS(X)   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( acos ( x[i] ) ) );
  }
  printf ( "\n" );

  printf ( "  ASIN(X)   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( asin ( x[i] ) ) );
  }
  printf ( "\n" );

  printf ( "  ATAN(X)   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( atan ( x[i] ) ) );
  }
  printf ( "\n" );

  printf ( "   LOG(X)   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( log ( x[i] ) ) );
  }
  printf ( "\n" );

  printf ( "  SQRT(X)   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( sqrt ( x[i] ) ) );
  }
  printf ( "\n" );

  printf ( "    1 / X   " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %d", isnan ( 1.0 / x[i] ) );
  }
  printf ( "\n" );

  return;
}
/******************************************************************************/

void test_j0 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_J0 tests J0.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST_J0:\n" );
  printf ( "  Test J0, which evaluates the Bessel J0(X) function.\n" );
  printf ( "\n" );
  printf ( "      X                 J0(X)           J0(X)         DIFF\n" );
  printf ( "                     (tabulated)     (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = j0 ( x );

    printf ( "  %14e  %14e  %14e  %14e\n", x, fx, fx2, r8_abs ( fx - fx2 ) );
  }
  return;
}
/******************************************************************************/

void bessel_j0_values ( int *n_data, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J0_VALUES returns some values of the J0 Bessel function.

  Discussion:

    In Mathematica, the function can be evaluated by:

      BesselJ[0,x]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 August 2004

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  static double fx_vec[N_MAX] = { 
     -0.1775967713143383E+00,  
     -0.3971498098638474E+00,  
     -0.2600519549019334E+00,  
      0.2238907791412357E+00,  
      0.7651976865579666E+00,  
      0.1000000000000000E+01,  
      0.7651976865579666E+00,  
      0.2238907791412357E+00,  
     -0.2600519549019334E+00,  
     -0.3971498098638474E+00,  
     -0.1775967713143383E+00,  
      0.1506452572509969E+00,  
      0.3000792705195556E+00,  
      0.1716508071375539E+00,  
     -0.9033361118287613E-01,  
     -0.2459357644513483E+00,  
     -0.1711903004071961E+00,  
      0.4768931079683354E-01,  
      0.2069261023770678E+00,  
      0.1710734761104587E+00,  
     -0.1422447282678077E-01 };

  static double x_vec[N_MAX] = { 
     -5.0E+00,  
     -4.0E+00,  
     -3.0E+00,  
     -2.0E+00,  
     -1.0E+00,  
      0.0E+00,  
      1.0E+00,  
      2.0E+00,  
      3.0E+00,  
      4.0E+00,  
      5.0E+00,  
      6.0E+00,  
      7.0E+00,  
      8.0E+00,  
      9.0E+00,  
     10.0E+00,  
     11.0E+00,  
     12.0E+00,  
     13.0E+00,  
     14.0E+00,  
     15.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void erf_values ( int *n_data, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    ERF_VALUES returns some values of the ERF or "error" function.

  Discussion:

    The error function is defined by:

      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT

    In Mathematica, the function can be evaluated by:

      Erf[x]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 August 2004

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  double fx_vec[N_MAX] = { 
     0.0000000000000000E+00,  
     0.1124629160182849E+00,  
     0.2227025892104785E+00,  
     0.3286267594591274E+00,  
     0.4283923550466685E+00,  
     0.5204998778130465E+00,  
     0.6038560908479259E+00,  
     0.6778011938374185E+00,  
     0.7421009647076605E+00,  
     0.7969082124228321E+00,  
     0.8427007929497149E+00,  
     0.8802050695740817E+00,  
     0.9103139782296354E+00,  
     0.9340079449406524E+00,  
     0.9522851197626488E+00,  
     0.9661051464753107E+00,  
     0.9763483833446440E+00,  
     0.9837904585907746E+00,  
     0.9890905016357307E+00,  
     0.9927904292352575E+00,  
     0.9953222650189527E+00 }; 

  double x_vec[N_MAX] = { 
     0.0E+00,   
     0.1E+00,   
     0.2E+00,   
     0.3E+00,   
     0.4E+00,   
     0.5E+00,   
     0.6E+00,   
     0.7E+00,   
     0.8E+00,   
     0.9E+00,   
     1.0E+00,   
     1.1E+00,   
     1.2E+00,   
     1.3E+00,   
     1.4E+00,   
     1.5E+00,   
     1.6E+00,   
     1.7E+00,   
     1.8E+00,   
     1.9E+00,   
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void erfc_values ( int *n_data, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    ERFC_VALUES returns some values of the ERFC or "complementary error" function.

  Discussion:

    The complementary error function is defined by:

      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT

    In Mathematica, the function can be evaluated by:

      Erfc[x]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2007

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  double fx_vec[N_MAX] = { 
    1.000000000000000E+00, 
    0.7772974107895215E+00, 
    0.5716076449533315E+00, 
    0.3961439091520741E+00, 
    0.2578990352923395E+00, 
    0.1572992070502851E+00, 
    0.08968602177036462E+00, 
    0.04771488023735119E+00, 
    0.02365161665535599E+00, 
    0.01090949836426929E+00, 
    0.004677734981047266E+00,
    0.001862846297981891E+00, 
    0.0006885138966450786E+00, 
    0.0002360344165293492E+00, 
    0.00007501319466545902E+00, 
    0.00002209049699858544E+00, 
    6.025761151762095E-06, 
    1.521993362862285E-06, 
    3.558629930076853E-07, 
    7.700392745696413E-08, 
    1.541725790028002E-08 }; 

  double x_vec[N_MAX] = { 
    0.0E+00,  
    0.2E+00, 
    0.4E+00,  
    0.6E+00, 
    0.8E+00,  
    1.0E+00,  
    1.2E+00, 
    1.4E+00, 
    1.6E+00, 
    1.8E+00, 
    2.0E+00,
    2.2E+00, 
    2.4E+00, 
    2.6E+00, 
    2.8E+00, 
    3.0E+00, 
    3.2E+00, 
    3.4E+00, 
    3.6E+00, 
    3.8E+00, 
    4.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void gamma_values ( int *n_data, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    GAMMA_VALUES returns some values of the Gamma function.

  Discussion:

    The Gamma function is defined as:

      Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT

    It satisfies the recursion:

      Gamma(X+1) = X * Gamma(X)

    Gamma is undefined for nonpositive integral X.
    Gamma(0.5) = sqrt(PI)
    For N a positive integer, Gamma(N+1) = N!, the standard factorial.

    In Mathematica, the function can be evaluated by:

      Gamma[x]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2007

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 25

  double fx_vec[N_MAX] = { 
     -0.3544907701811032E+01,  
     -0.1005871979644108E+03,  
      0.9943258511915060E+02,  
      0.9513507698668732E+01,  
      0.4590843711998803E+01,  
      0.2218159543757688E+01,  
      0.1772453850905516E+01,  
      0.1489192248812817E+01,  
      0.1164229713725303E+01,  
      0.1000000000000000E+01,  
      0.9513507698668732E+00,  
      0.9181687423997606E+00,  
      0.8974706963062772E+00,  
      0.8872638175030753E+00,  
      0.8862269254527580E+00,  
      0.8935153492876903E+00,  
      0.9086387328532904E+00,  
      0.9313837709802427E+00,  
      0.9617658319073874E+00,  
      0.1000000000000000E+01,  
      0.2000000000000000E+01,  
      0.6000000000000000E+01,  
      0.3628800000000000E+06,  
      0.1216451004088320E+18,  
      0.8841761993739702E+31 };

  double x_vec[N_MAX] = { 
     -0.50E+00,  
     -0.01E+00,  
      0.01E+00,  
      0.10E+00,  
      0.20E+00,  
      0.40E+00,  
      0.50E+00,  
      0.60E+00,  
      0.80E+00,  
      1.00E+00,  
      1.10E+00,  
      1.20E+00,  
      1.30E+00,  
      1.40E+00,  
      1.50E+00,  
      1.60E+00,  
      1.70E+00,  
      1.80E+00,  
      1.90E+00,  
      2.00E+00,  
      3.00E+00,  
      4.00E+00,  
     10.00E+00,  
     20.00E+00,  
     30.00E+00 }; 

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
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
