# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "brent.h"

int main ( void );

void test_zero_all ( void );
void test_zero_rc_all ( void );
void test_local_min_all ( void );
void test_local_min_rc_all ( void );
void test_glomin_all ( void );

void test_zero_one ( double a, double b, double machep, double t, 
  double f ( double x ), char *title );
void test_zero_rc_one ( double a, double b, double machep, double t, 
  double f ( double x ), char *title );
void test_local_min_one ( double a, double b, double eps, double t, 
  double f ( double x ), char *title );
void test_local_min_rc_one ( double a, double b, double t, 
  double f ( double x ), char *title );
void test_glomin_one ( double a, double b, double c, double m, double machep, 
  double e, double t, double f ( double x ), char *title );

double f_01 ( double x );
double f_02 ( double x );
double f_03 ( double x );
double f_04 ( double x );
double f_05 ( double x );
double g_01 ( double x );
double g_02 ( double x );
double g_03 ( double x );
double g_04 ( double x );
double g_05 ( double x );
double h_01 ( double x );
double h_02 ( double x );
double h_03 ( double x );
double h_04 ( double x );
double h_05 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BRENT_PRB.

  Discussion:

    BRENT_PRB tests the BRENT library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BRENT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BRENT library.\n" );

  test_zero_all ( );
  test_zero_rc_all ( );
  test_local_min_all ( );
  test_local_min_rc_all ( );
  test_glomin_all ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BRENT_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_zero_all ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ZERO_ALL tests ZERO on all test functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double machep;
  double t;
  double z;

  printf ( "\n" );
  printf ( "TEST_ZERO_ALL\n" );
  printf ( "  Test the Brent ZERO routine, which seeks\n" );
  printf ( "  a root of a function F(X)\n" );
  printf ( "  in an interval [A,B].\n" );

  machep = r8_epsilon ( );
  t = machep;

  a = 1.0;
  b = 2.0;

  test_zero_one ( a, b, machep, t, f_01, 
    "f_01(x) = sin ( x ) - x / 2" );

  a = 0.0;
  b = 1.0;

  test_zero_one ( a, b, machep, t, f_02, 
    "f_02(x) = 2 * x - exp ( - x )" );

  a = -1.0;
  b =  0.5;

  test_zero_one ( a, b, machep, t, f_03, 
    "f_03(x) = x * exp ( - x )" );

  a =  0.0001;
  b =  20.0;

  test_zero_one ( a, b, machep, t, f_04, 
    "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )" );

  a = -5.0;
  b =  2.0;

  test_zero_one ( a, b, machep, t, f_05, 
    "f_05(x) = (x+3) * (x-1) * (x-1)" );

  return;
}
/******************************************************************************/

void test_zero_rc_all ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ZERO_RC_ALL tests ZERO_RC on all test functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double machep;
  double t;
  double z;

  printf ( "\n" );
  printf ( "TEST_ZERO_RC_ALL\n" );
  printf ( "  Test the ZERO_RC routine, which seeks\n" );
  printf ( "  a root of a function F(X)\n" );
  printf ( "  in an interval [A,B].\n" );

  machep = r8_epsilon ( );
  t = machep;

  a = 1.0;
  b = 2.0;

  test_zero_rc_one ( a, b, machep, t, f_01, 
    "f_01(x) = sin ( x ) - x / 2" );

  a = 0.0;
  b = 1.0;

  test_zero_rc_one ( a, b, machep, t, f_02, 
    "f_02(x) = 2 * x - exp ( - x )" );

  a = -1.0;
  b =  0.5;

  test_zero_rc_one ( a, b, machep, t, f_03, 
    "f_03(x) = x * exp ( - x )" );

  a =  0.0001;
  b =  20.0;

  test_zero_rc_one ( a, b, machep, t, f_04, 
    "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )" );

  a = -5.0;
  b =  2.0;

  test_zero_rc_one ( a, b, machep, t, f_05, 
    "f_05(x) = (x+3) * (x-1) * (x-1)" );

  return;
}
/******************************************************************************/

void test_local_min_all ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LOCAL_MIN_ALL tests LOCAL_MIN on all test functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double eps;
  double t;
  double z;

  printf ( "\n" );
  printf ( "TEST_LOCAL_MIN_ALL\n" );
  printf ( "  Test the Brent LOCAL_MIN routine, which seeks\n" );
  printf ( "  a local minimizer of a function F(X)\n" );
  printf ( "  in an interval [A,B].\n" );

  eps = 10.0 * sqrt ( r8_epsilon ( ) );
  t = eps;

  a = 0.0;
  b = 3.141592653589793;

  test_local_min_one ( a, b, eps, t, g_01, 
    "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1" );

  a = 0.0;
  b = 1.0;

  test_local_min_one ( a, b, eps, t, g_02, 
    "g_02(x) = x * x + exp ( - x )" );

  a = -2.0;
  b =  2.0;

  test_local_min_one ( a, b, eps, t, g_03, 
    "g_03(x) = x^4 + 2x^2 + x + 3" );

  a =  0.0001;
  b =  1.0;

  test_local_min_one ( a, b, eps, t, g_04, 
    "g_04(x) = exp ( x ) + 1 / ( 100 x )" );

  a =  0.0002;
  b = 2.0;

  test_local_min_one ( a, b, eps, t, g_05, 
    "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)" );

  return;
}
/******************************************************************************/

void test_local_min_rc_all ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LOCAL_MIN_RC_ALL tests LOCAL_MIN_RC on all test functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double t;
  double z;

  printf ( "\n" );
  printf ( "TEST_LOCAL_MIN_RC_ALL\n" );
  printf ( "  Test LOCAL_MIN_RC, a reverse communication version of\n" );
  printf ( "  the Brent LOCAL_MIN routine, which seeks\n" );
  printf ( "  a local minimizer of a function F(X)\n" );
  printf ( "  in an interval [A,B].\n" );

  t = 10.0 * sqrt ( r8_epsilon ( ) );

  a = 0.0;
  b = 3.141592653589793;

  test_local_min_rc_one ( a, b, t, g_01, 
    "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1" );

  a = 0.0;
  b = 1.0;

  test_local_min_rc_one ( a, b, t, g_02, 
    "g_02(x) = x * x + exp ( - x )" );

  a = -2.0;
  b =  2.0;

  test_local_min_rc_one ( a, b, t, g_03, 
    "g_03(x) = x^4 + 2x^2 + x + 3" );

  a =  0.0001;
  b =  1.0;

  test_local_min_rc_one ( a, b, t, g_04, 
    "g_04(x) = exp ( x ) + 1 / ( 100 x )" );

  a =  0.0002;
  b = 2.0;

  test_local_min_rc_one ( a, b, t, g_05, 
    "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)" );

  return;
}
/******************************************************************************/

void test_glomin_all ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_GLOMIN_ALL tests GLOMIN on all test functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double e;
  double m;
  double machep;
  double t;
  double z;

  printf ( "\n" );
  printf ( "TEST_GLOMIN_ALL\n" );
  printf ( "  Test the Brent GLOMIN routine, which seeks\n" );
  printf ( "  a global minimizer of a function F(X)\n" );
  printf ( "  in an interval [A,B],\n" );
  printf ( "  given some upper bound M \n" );
  printf ( "  for the second derivative of F.\n" );

  machep = r8_epsilon ( );
  e = sqrt ( machep );
  t = sqrt ( machep );

  a = 7.0;
  b = 9.0;
  c = ( a + b ) / 2.0;
  m = 0.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_01, 
    "h_01(x) = 2 - x" );

  a = 7.0;
  b = 9.0;
  c = ( a + b ) / 2.0;
  m = 100.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_01, 
    "h_01(x) = 2 - x" );

  a = -1.0;
  b = +2.0;
  c = ( a + b ) / 2.0;
  m = 2.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_02, 
    "h_02(x) = x * x" );

  a = -1.0;
  b = +2.0;
  c = ( a + b ) / 2.0;
  m = 2.1;

  test_glomin_one ( a, b, c, m, machep, e, t, h_02, 
    "h_02(x) = x * x" );

  a = -0.5;
  b =  +2.0;
  c = ( a + b ) / 2.0;
  m = 14.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_03, 
    "h_03(x) = x^3 + x^2" );

  a = -0.5;
  b =  +2.0;
  c = ( a + b ) / 2.0;
  m = 28.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_03, 
    "h_03(x) = x^3 + x^2" );

  a = -10.0;
  b = +10.0;
  c = ( a + b ) / 2.0;
  m = 72.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_04, 
    "h_04(x) = ( x + sin(x) ) * exp(-x*x)" );

  a = -10.0;
  b = +10.0;
  c = ( a + b ) / 2.0;
  m = 72.0;

  test_glomin_one ( a, b, c, m, machep, e, t, h_05, 
    "h_05(x) = ( x - sin(x) ) * exp(-x*x)" );

  return;
}
/******************************************************************************/

void test_zero_one ( double a, double b, double machep, double t, 
  double f ( double x ), char *title )

/******************************************************************************/
/*
  Purpose:

    TEST_ZERO_ONE tests ZERO on one test function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the two endpoints of the change of sign
    interval.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double T, a positive error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function which evaluates the function whose zero is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
  double fa;
  double fb;
  double fz;
  double z;
  
  z = zero ( a, b, machep, t, f );
  fz = f ( z );
  fa = f ( a );
  fb = f ( b );

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "      A                 Z             B\n" );
  printf ( "    F(A)              F(Z)          F(B)\n" );
  printf ( "\n" );
  printf ( "  %14f  %14f  %14f\n", a, z, b );
  printf ( "  %14e  %14e  %14e\n", fa, fz, fb );

  return;
}
/******************************************************************************/

void test_zero_rc_one ( double a, double b, double machep, double t, 
  double f ( double x ), char *title )

/******************************************************************************/
/*
  Purpose:

    TEST_ZERO_RC_ONE tests ZERO_RC on one test function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the two endpoints of the change of sign
    interval.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double T, a positive error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function which evaluates the function whose zero is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
  double arg;
  int status;
  double value;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "    STATUS      X               F(X)\n" );
  printf ( "\n" );

  status = 0;

  for ( ; ; ) 
  {
    zero_rc ( a, b, t, &arg, &status, value );

    if ( status < 0 )
    {
      printf ( "\n" );
      printf ( "  ZERO_RC returned an error flag!\n" );
      break;
    }

    value = f ( arg );

    printf ( "  %8d  %14e  %14e\n", status, arg, value );

    if ( status == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void test_local_min_one ( double a, double b, double eps, double t, 
  double f ( double x ), char *title )

/******************************************************************************/
/*
  Purpose:

    TEST_LOCAL_MIN_ONE tests LOCAL_MIN on one test function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the endpoints of the interval.

    Input, double EPS, a positive relative error tolerance.

    Input, double T, a positive absolute error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function, whose local minimum is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
  double fa;
  double fb;
  double fx;
  double x;

  fx = local_min ( a, b, eps, t, f, &x );
  fa = f ( a );
  fb = f ( b );

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "      A                 Z             B\n" );
  printf ( "    F(A)              F(Z)          F(B)\n" );
  printf ( "\n" );
  printf ( "  %14f  %14f  %14f\n",  a,  x,  b );
  printf ( "  %14e  %14e  %14e\n", fa, fx, fb );

  return;
}
/******************************************************************************/

void test_local_min_rc_one ( double a, double b, double t, 
  double f ( double x ), char *title )

/******************************************************************************/
/*
  Purpose:

    TEST_LOCAL_MIN_RC_ONE tests LOCAL_MIN_RC on one test function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the endpoints of the interval.

    Input, double T, a positive absolute error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function, whose local minimum is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
  double a2;
  double arg;
  double b2;
  int status;
  int step;
  double value;
  double x;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "  Step      X                          F(X)\n" );
  printf ( "\n" );
  step = 0;

  arg = a;
  value = f ( arg );
  printf ( "  %4d  %24.16e  %24.16e\n", step, arg, value );

  arg = b;
  value = f ( arg );
  printf ( "  %4d  %24.16e  %24.16e\n", step, arg, value );

  a2 = a;
  b2 = b;
  status = 0;

  for ( ; ; )
  {
    arg = local_min_rc ( &a2, &b2, &status, value );
 
    if ( status < 0 )
    {
      printf ( "\n" );
      printf ( "TEST_LOCAL_MIN_RC_ONE - Fatal error!\n" );
      printf ( "  LOCAL_MIN_RC returned negative status.\n" );
      break;
    }

    value = f ( arg );

    step = step + 1;
    printf ( "  %4d  %24.16e  %24.16e\n", step, arg, value );

    if ( status == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void test_glomin_one ( double a, double b, double c, double m, double machep, 
  double e, double t, double f ( double x ), char *title )

/******************************************************************************/
/*
  Purpose:

    TEST_GLOMIN_ONE tests GLOMIN on one test function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the endpoints of the interval.

    Input, double C, an initial guess for the global
    minimizer.  If no good guess is known, C = A or B is acceptable.

    Input, double M, the bound on the second derivative.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double E, a positive tolerance, a bound for the
    absolute error in the evaluation of F(X) for any X in [A,B].

    Input, double T, a positive absolute error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function whose global minimum is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
  double fa;
  double fb;
  double fx;
  double x;

  fx = glomin ( a, b, c, m, machep, e, t, f, &x );
  fa = f ( a );
  fb = f ( b );

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "      A                 X             B\n" );
  printf ( "    F(A)              F(X)          F(B)\n" );
  printf ( "\n" );
  printf ( "  %14f  %14f  %14f\n",  a,  x,  b );
  printf ( "  %14e  %14e  %14e\n", fa, fx, fb );

  return;
}
/******************************************************************************/

double f_01 ( double x )

/******************************************************************************/
/*
  Purpose:

    F_01 evaluates sin ( x ) - x / 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double F_01, the value of the function at X.
*/
{
  double value;

  value = sin ( x ) - 0.5 * x;

  return value;
}
/******************************************************************************/

double f_02 ( double x )

/******************************************************************************/
/*
  Purpose:

    F_02 evaluates 2*x-exp(-x).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double F_02, the value of the function at X.
*/
{
  double value;

  value = 2.0 * x - exp ( - x );

  return value;
}
/******************************************************************************/

double f_03 ( double x )

/******************************************************************************/
/*
  Purpose:

    F_03 evaluates x*exp(-x).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double F_03, the value of the function at X.
*/
{
  double value;

  value = x * exp ( - x );

  return value;
}
/******************************************************************************/

double f_04 ( double x )

/******************************************************************************/
/*
  Purpose:

    F_04 evaluates exp(x) - 1 / (100*x*x).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double F_04, the value of the function at X.
*/
{
  double value;

  value = exp ( x ) - 1.0 / 100.0 / x / x;

  return value;
}
/******************************************************************************/

double f_05 ( double x )

/******************************************************************************/
/*
  Purpose:

    F_05 evaluates (x+3)*(x-1)*(x-1).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double F_05, the value of the function at X.
*/
{
  double value;

  value = ( x + 3.0 ) * ( x - 1.0 ) * ( x - 1.0 );

  return value;
}
/******************************************************************************/

double g_01 ( double x )

/******************************************************************************/
/*
  Purpose:

    G_01 evaluates (x-2)^2 + 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double G_01, the value of the function at X.
*/
{
  double value;

  value = ( x - 2.0 ) * ( x - 2.0 ) + 1.0;

  return value;
}
/******************************************************************************/

double g_02 ( double x )

/******************************************************************************/
/*
  Purpose:

    G_02 evaluates x^2 + exp ( - x ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double G_02, the value of the function at X.
*/
{
  double value;

  value = x * x + exp ( - x );

  return value;
}
/******************************************************************************/

double g_03 ( double x )

/******************************************************************************/
/*
  Purpose:

    G_03 evaluates x^4+2x^2+x+3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double G_03, the value of the function at X.
*/
{
  double value;

  value = ( ( x * x + 2.0 ) * x + 1.0 ) * x + 3.0;

  return value;
}
/******************************************************************************/

double g_04 ( double x )

/******************************************************************************/
/*
  Purpose:

    G_04 evaluates exp(x)+1/(100X)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double G_04, the value of the function at X.
*/
{
  double value;

  value = exp ( x ) + 0.01 / x;

  return value;
}
/******************************************************************************/

double g_05 ( double x )

/******************************************************************************/
/*
  Purpose:

    G_05 evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double G_05, the value of the function at X.
*/
{
  double value;

  value = exp ( x ) - 2.0 * x + 0.01 / x - 0.000001 / x / x;

  return value;
}
/******************************************************************************/

double h_01 ( double x )

/******************************************************************************/
/*
  Purpose:

    H_01 evaluates 2 - x.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double H_01, the value of the function at X.
*/
{
  double value;

  value = 2.0 - x;

  return value;
}
/******************************************************************************/

double h_02 ( double x )

/******************************************************************************/
/*
  Purpose:

    H_02 evaluates x^2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double H_02, the value of the function at X.
*/
{
  double value;

  value = x * x;

  return value;
}
/******************************************************************************/

double h_03 ( double x )

/******************************************************************************/
/*
  Purpose:

    H_03 evaluates x^3+x^2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double H_03, the value of the function at X.
*/
{
  double value;

  value = x * x * ( x + 1.0 );

  return value;
}
/******************************************************************************/

double h_04 ( double x )

/******************************************************************************/
/*
  Purpose:

    H_04 evaluates ( x + sin ( x ) ) * exp ( - x * x ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double H_04, the value of the function at X.
*/
{
  double value;

  value = ( x + sin ( x ) ) * exp ( - x * x );

  return value;
}
/******************************************************************************/

double h_05 ( double x )

/******************************************************************************/
/*
  Purpose:

    H_05 evaluates ( x - sin ( x ) ) * exp ( - x * x ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which F is to be evaluated.

    Output, double H_05, the value of the function at X.
*/
{
  double value;

  value = ( x - sin ( x ) ) * exp ( - x * x );

  return value;
}
