# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <math.h>
# include <complex.h>

# include "polpak.h"

int main ( void );

void test001 ( void );
void test002 ( void );
void test003 ( void );
void test004 ( void );
void test005 ( void );
void test006 ( void );
void test007 ( void );
void test008 ( void );
void test009 ( void );

void test010 ( void );
void test011 ( void );
void test0115 ( void );
void test013 ( void );
void test012 ( void );
void test014 ( void );
void test015 ( void );
void test016 ( void );
void test017 ( void );
void test0175 ( void );
void test018 ( void );
void test0185 ( void );
void test019 ( void );

void test020 ( void );
void test021 ( void );
void test0215 ( void );
void test0216 ( void );
void test0217 ( void );
void test0218 ( void );
void test024 ( void );
void test0243 ( void );
void test025 ( void );
void test0255 ( void );
void test026 ( void );
void test0265 ( void );
void test028 ( void );
void test027 ( void );
void test029 ( void );

void test031 ( void );
void test032 ( void );
void test033 ( void );
void test034 ( void );
void test037 ( void );
void test052 ( void );
void test038 ( void );

void test040 ( void );
void test039 ( void );
void test041 ( void );
void test042 ( void );
void test0427 ( void );
void test043 ( void );
void test044 ( void );
void test045 ( void );
void test046 ( void );
void test047 ( void );
void test048 ( void );
void test049 ( void );

void test050 ( void );
void test0505 ( void );
void test051 ( void );
void test054 ( void );
void test055 ( void );
void test057 ( void );
void test058 ( void );
void test059 ( void );
void test0595 ( void );

void test060 ( void );
void test061 ( void );
void test0615 ( void );
void test0365 ( void );
void test062 ( void );
void test0623 ( void );
void test0625 ( void );
void test063 ( void );
void test0635 ( void );
void test064 ( void );
void test065 ( void );
void test066 ( void );
void test0665 ( void );
void test0667 ( void );
void test067 ( void );
void test0675 ( void );
void test0676 ( void );
void test06765 ( void );
void test022 ( void );
void test068 ( void );
void test0685 ( void );
void test06855 ( void );
void test036 ( void );
void test0425 ( void );
void test06856 ( void );
void test069 ( void );
void test0695 ( void );
void test0696 ( void );
void test0697 ( void );

void test070 ( void );
void test071 ( void );
void test072 ( void );
void test073 ( void );
void test074 ( void );
void test075 ( void );
void test076 ( void );
void test077 ( void );
void test0773 ( void );
void test0775 ( void );
void test078 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    POLPAK_PRB calls the POLPAK tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLPAK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the POLPAK library.\n" );

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
  test011 ( );
  test0115 ( );
  test013 ( );
  test012 ( );
  test014 ( );
  test015 ( );
  test016 ( );
  test017 ( );
  test0175 ( );
  test018 ( );
  test0185 ( );
  test019 ( );

  test020 ( );
  test021 ( );
  test0215 ( );
  test0216 ( );
  test0217 ( );
  test0218 ( );
  test024 ( );
  test0243 ( );
  test025 ( );
  test0255 ( );
  test026 ( );
  test0265 ( );
  test028 ( );
  test027 ( );
  test029 ( );

  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test037 ( );
  test052 ( );
  test038 ( );

  test040 ( );
  test039 ( );
  test041 ( );
  test042 ( );
  test0427 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test0505 ( );
  test051 ( );
  test054 ( );
  test055 ( );
  test057 ( );
  test058 ( );
  test059 ( );
  test0595 ( );

  test060 ( );
  test061 ( );
  test0615 ( );
  test0365 ( );
  test062 ( );
  test0623 ( );
  test0625 ( );
  test063 ( );
  test0635 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test0665 ( );
  test0667 ( );
  test067 ( );
  test0675 ( );
  test0676 ( );
  test06765 ( );
  test022 ( );
  test068 ( );
  test0685 ( );
  test06855 ( );
  test036 ( );
  test0425 ( );
  test06856 ( );
  test069 ( );
  test0695 ( );
  test0696 ( );
  test0697 ( );

  test070 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test075 ( );
  test076 ( );
  test077 ( );
  test0773 ( );
  test0775 ( );
  test078 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLPAK_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests AGM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 February 2010

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  AGM computes the arithmetic geometric mean.\n" );
  printf ( "\n" );
  printf ( "      A           B          " );
  printf ( "   AGM                       AGM                   Diff" );
  printf ( "                             " );
  printf ( "  (Tabulated)                AGM(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( &n_data, &a, &b, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = agm ( a, b );

    printf ( "  %10.6f  %10.6f  %24.16f  %24.16f  %10.6e\n",
      a, b, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test002 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests AGUD and GUD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2010

  Author:

    John Burkardt
*/
{
  double gamma;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  AGUD computes the inverse Gudermannian;\n" );
  printf ( "  GUD computes the Gudermannian.\n" );
  printf ( "\n" );
  printf ( "         X     GUD(X)     AGUD(GUD(X))\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    gamma = gud ( x );
    x2 = agud ( gamma );

    printf ( "  %10f  %10f  %10f\n", x, gamma, x2 );
  }

  return;
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests ALIGN_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
# define M_MAX 10
# define N_MAX 10

  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  ALIGN_ENUM counts the number of possible\n" );
  printf ( "  alignments of two biological sequences.\n" );

  printf ( "\n" );
  printf ( "  Alignment enumeration table:\n" );
  printf ( "\n" );

  printf ( "      " );
  for ( j = 0; j <= 5; j++ )
  {
    printf ( "%8d  ", j );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( i = 0; i <= M_MAX; i++ )
  {
    printf ( "   %2d  ", i );
    for ( j = 0; j <= 5; j++ )
    {
      printf ( "%8d  ", align_enum ( i, j ) );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "      " );
  for ( j = 6; j <= N_MAX; j++ )
  {
    printf ( "%8d  ", j );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( i = 0; i <= M_MAX; i++ )
  {
    printf ( "  %2d  ", i );
    for ( j = 6; j <= N_MAX; j++ )
    {
      printf ( "%8d  ", align_enum ( i, j ) );
    }
    printf ( "\n" );
  }
  return;
# undef M_MAX
# undef N_MAX
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
/*
  Purpose:

    TEST004 tests ARC_COSINE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST004\n" );
  printf ( "  ARC_COSINE computes the inverse cosine,\n" );
  printf ( "  and chops input arguments that are out of bounds.\n" );
  printf ( "\n" );
  printf ( "         X     ARC_COSINE(X)     COS(ARC_COSINE(X))\n" );
  printf ( "\n" );

  for ( i = -5; i <= 5; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = arc_cosine ( x );
    x2 = cos ( a );

    printf ( "  %10g  %10g  %10g\n", x, a, x2 );
  }

  return;
}
/******************************************************************************/

void test005 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests ATAN4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 August 2011

  Author:

    John Burkardt
*/
{
# define TEST_NUM 8

  int test;
  int test_num = TEST_NUM;
  double x;
  double x_test[TEST_NUM] = {
     1.0,  1.0, 0.0, -1.0, 
    -1.0, -1.0, 0.0,  1.0 };
  double y;
  double y_test[TEST_NUM] = {
    0.0,  1.0,  1.0,  1.0, 
    0.0, -1.0, -1.0, -1.0 };

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  ATAN4 computes the arc-tangent given Y and X;\n" );
  printf ( "  ATAN2 is the system version of this routine.\n" );
  printf ( "\n" );
  printf ( "       X             Y          ATAN2(Y,X)   ATAN4(Y,X)\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    x = x_test[test];
    y = y_test[test];
    printf ( "  %14g  %14g  %14g  %14g\n", 
      x, y, atan2 ( y, x ), atan4 ( y, x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests BELL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  BELL computes Bell numbers.\n" );
  printf ( "  BELL_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    bell ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void test007 ( )

/******************************************************************************/
/*
  Purpose:

    TEST007 tests BENFORD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "TEST007\n" );
  printf ( "  BENFORD(I) is the Benford probability of the\n" );
  printf ( "  initial digit sequence I.\n" );
  printf ( "\n" );
  printf ( "     I  BENFORD(I)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 9; i++ )
  {
    printf ( "  %4d  %10.4f\n", i, benford ( i ) );
  }

  return;
}
/******************************************************************************/

void test008 ( )

/******************************************************************************/
/*
  Purpose:

    TEST008 tests BERNOULLI_NUMBER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST008\n" );
  printf ( "  BERNOULLI_NUMBER computes Bernoulli numbers;\n" );
  printf ( "  BERNOULLI_NUMBER_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number ( n, c1 );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1[n] );
  }

  return;
}
/******************************************************************************/

void test009 ( )

/******************************************************************************/
/*
  Purpose:

    TEST009 tests BERNOULLI_NUMBER2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST009\n" );
  printf ( "  BERNOULLI_NUMBER2 computes Bernoulli numbers;\n" );
  printf ( "  BERNOULLI_NUMBER_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number2 ( n, c1 );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1[n] );
  }

  return;
}
/******************************************************************************/

void test010 ( )

/******************************************************************************/
/*
  Purpose:

    TEST010 tests BERNOULLI_NUMBER3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST010\n" );
  printf ( "  BERNOULLI_NUMBER3 computes Bernoulli numbers;\n" );
  printf ( "  BERNOULLI_NUMBER_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER3\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    c1 = bernoulli_number3 ( n );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1 );
  }

  return;
}
/******************************************************************************/

void test011 ( )

/******************************************************************************/
/*
  Purpose:

    TEST011 tests BERNOULLI_POLY;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  printf ( "\n" );
  printf ( "TEST011\n" );
  printf ( "  BERNOULLI_POLY evaluates Bernoulli polynomials;\n" );
  printf ( "\n" );
  printf ( "  X = %g\n", x );
  printf ( "\n" );
  printf ( "  I          BX\n" );
  printf ( "\n" );

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly ( i, x );

    printf ( "  %6d  %10g\n", i, bx );
  }

  return;
}
/******************************************************************************/

void test0115 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0115 tests BERNOULLI_POLY2;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  printf ( "\n" );
  printf ( "TEST0115\n" );
  printf ( "  BERNOULLI_POLY2 evaluates Bernoulli polynomials;\n" );
  printf ( "\n" );
  printf ( "  X = %g\n", x );
  printf ( "\n" );
  printf ( "  I          BX\n" );
  printf ( "\n" );

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly2 ( i, x );

    printf ( "  %6d  %10g\n", i, bx );
  }

  return;
}
/******************************************************************************/

void test013 ( )

/******************************************************************************/
/*
  Purpose:

    TEST013 tests BERNSTEIN_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double b;
  double bvec[11];
  int k;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST013:\n" );
  printf ( "  BERNSTEIN_POLY evaluates the Bernstein polynomials.\n" );
  printf ( "  BERNSTEIN_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "   N   K   X   Exact   B(N,K)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }

    bernstein_poly ( n, x, bvec );

    printf ( "  %4d  %4d  %7g  %14g  %14g\n", n, k, x, b, bvec[k] );
  }

  return;
}
/******************************************************************************/

void test012 ( )

/******************************************************************************/
/*
  Purpose:

    TEST012 tests BETA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double fxy;
  double fxy2;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST012:\n" );
  printf ( "  BETA evaluates the Beta function.\n" );
  printf ( "  BETA_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X      Y        Exact F       BETA(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = beta ( x, y );

    printf ( "  %10f  %10f  %10g  %10g\n", x, y, fxy, fxy2 );
  }

  return;
}
/******************************************************************************/

void test014 ( )

/******************************************************************************/
/*
  Purpose:

    TEST014 tests BPAB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double a;
  double b;
  double bern[N+1];
  int i;
  double x;

  printf ( "\n" );
  printf ( "TEST014\n" );
  printf ( "  BPAB evaluates Bernstein polynomials.\n" );
  printf ( "\n" );

  x = 0.3;
  a = 0.0;
  b = 1.0;

  bpab ( N, x, a, b, bern );

  printf ( "  The Bernstein polynomials of degree %d\n", N );
  printf ( "  based on the interval from %g\n", a );
  printf ( "  to %g\n", b );
  printf ( "  evaluated at X = %g\n", x );
  printf ( "\n" );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %14g\n", i, bern[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test015 ( )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests CARDAN and CARDAN_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  double *cx2;
  int i;
  int n;
  double s;
  double x;

  s = 1.0;

  printf ( "\n" );
  printf ( "TEST015\n" );
  printf ( "  CARDAN_POLY_COEF returns the coefficients of a\n" );
  printf ( "  Cardan polynomial.\n" );
  printf ( "  CARDAN evaluates a Cardan polynomial directly.\n" );
  printf ( "\n" );
  printf ( "  We use the parameter S = %g\n", s );
  printf ( "\n" );
  printf ( "  Table of polynomial coefficients:\n" );
  printf ( "\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );
    printf ( "  %2d  ", n );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%5d  ", c[i] );
    }
    printf ( "\n" );
  }

  s = 0.5;
  x = 0.25;

  printf ( "\n" );
  printf ( "  Compare CARDAN_POLY_COEF + R8POLY_VAL_HORNER\n" );
  printf ( "  versus CARDAN alone.\n" );
  printf ( "\n" );
  printf ( "  Evaluate polynomials at X = %g\n", x );
  printf ( "  We use the parameter S = %g\n", s );
  printf ( "\n" );
  printf ( "  Order       Horner          Direct\n" );
  printf ( "\n" );

  cx2 = cardan ( n, x, s );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );

    cx1 = r8poly_value ( n + 1, c, x );

    printf ( "  %2d  %14g  %14g\n", n, cx1, cx2[n] );
  }
  free ( cx2 );

  return;
# undef N_MAX
}
/******************************************************************************/

void test016 ( )

/******************************************************************************/
/*
  Purpose:

    TEST016 tests CATALAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST016\n" );
  printf ( "  CATALAN computes Catalan numbers.\n" );
  printf ( "  CATALAN_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    catalan ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void test017 ( )

/******************************************************************************/
/*
  Purpose:

    TEST017 tests CATALAN_ROW_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  bool next;

  printf ( "\n" );
  printf ( "TEST017\n" );
  printf ( "  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.\n" );
  printf ( "\n" );
  printf ( "  First, compute row 7:\n" );
  printf ( "\n" );

  next = false;
  n = 7;
  catalan_row_next ( next, n, c );

  printf ( "%4d  ", n );
  for ( i = 0; i <= n; i++ )
  {
    printf ( "%8d  ", c[i] );
  }
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  Now compute rows consecutively, one at a time:\n" );
  printf ( "\n" );

  next = false;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = true;

    printf ( "%4d  ", i );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%8d  ", c[i] );
    }
    printf ( "\n" );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test0175 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0175 tests CHARLIER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5
# define N 5

  double a;
  double a_test[TEST_NUM] = { 0.25, 0.5, 1.0, 2.0, 10.0 };
  int i;
  int j;
  int n;
  int test;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "TEST0175:\n" );
  printf ( "  CHARLIER evaluates Charlier polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      A         X        P(N,A,X)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    a = a_test[test];

    printf ( "\n" );

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      charlier ( n, a, x, value );

      printf ( "\n" );
      for ( i = 0; i <= 5; i++ )
      {

        printf ( "  %6d  %8g  %8g  %14g\n", i, a, x, value[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test018 ( )

/******************************************************************************/
/*
  Purpose:

    TEST018 tests CHEBY_T_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST018:\n" );
  printf ( "  CHEBY_T_POLY evaluates the Chebyshev T polynomial.\n" );
  printf ( "  CHEBY_T_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       T(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_t_poly ( 1, n, x_vec );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test0185 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0185 tests CHEBY_T_POLY_ZERO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  printf ( "\n" );
  printf ( "TEST0185:\n" );
  printf ( "  CHEBY_T_POLY_ZERO returns zeroes of T(N,X).\n" );
  printf ( "\n" );
  printf ( "       N      X        T(N,X)\n" );
  printf ( "\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_t_poly_zero ( n );
    fx = cheby_t_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[i+n*n] );
    }
    printf ( "\n" );
    free ( fx );
    free ( z );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test019 ( )

/******************************************************************************/
/*
  Purpose:

    TEST019 tests CHEBY_T_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 April 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 5;

  printf ( "\n" );
  printf ( "TEST019\n" );
  printf ( "  CHEBY_T_POLY_COEF determines the  polynomial coefficients\n" );
  printf ( "  of the Chebyshev polynomial T(n,x).\n" );

  c = cheby_t_poly_coef ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  T(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(n+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(n+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(n+1)], j );
        }
      }
    }
  }
 
  free ( c );

  return;
}
/******************************************************************************/

void test020 ( )

/******************************************************************************/
/*
  Purpose:

    TEST020 tests CHEBY_U_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST020:\n" );
  printf ( "  CHEBY_U_POLY evaluates the Chebyshev U polynomial.\n" );
  printf ( "  CHEBY_U_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       U(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    cheby_u_poly ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test021 ( )

/******************************************************************************/
/*
  Purpose:

    TEST021 tests CHEBY_U_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST021\n" );
  printf ( "  CHEBY_U_POLY_COEF determines the polynomial coefficients\n" );
  printf ( "  of the Chebyshev polynomial U(n,x).\n" );

  cheby_u_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  U(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(N+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(N+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
        }
      }
    }
  }
 
  return;
# undef N
}
/******************************************************************************/

void test0215 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0215 tests CHEBY_U_POLY_ZERO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 4

  double fx[N_MAX+1];
  int i;
  int n;
  double *z;

  printf ( "\n" );
  printf ( "TEST0215:\n" );
  printf ( "  CHEBY_U_POLY_ZERO returns zeroes of U(N,X).\n" );
  printf ( "\n" );
  printf ( "       N      X        U(N,X)\n" );
  printf ( "\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_u_poly_zero ( n );

    for ( i = 0; i < n; i++ )
    {
      cheby_u_poly ( n, z[i], fx );

      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[n] );
    }
    printf ( "\n" );
    free ( z );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test0216 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0216 tests CHEBYSHEV_DISCRETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5
# define N 5

  int i;
  int j;
  int m;
  int n;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "TEST0216:\n" );
  printf ( "  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      M         X        T(N,M,X)\n" );

  m = 5;
  n = N;

  for ( j = 0; j <= 5; j++ )
  {
    x = ( double ) ( j ) / 2.0;

    chebyshev_discrete ( n, m, x, value );

    printf ( "\n" );
    for ( i = 0; i <= 5; i++ )
    {
      printf ( "  %6d  %6d  %8g  %14g\n", i, m, x, value[i] );
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test0217 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0217 tests COLLATZ_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  int count;
  int count2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0217:\n" );
  printf ( "  COLLATZ_COUNT(N) counts the length of the\n" );
  printf ( "  Collatz sequence beginning with N.\n" );
  printf ( "\n" );
  printf ( "       N       COUNT(N)     COUNT(N)\n" );
  printf ( "              (computed)    (table)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }

    count2 = collatz_count ( n );

    printf ( "  %8d  %8d  %8d\n", n, count, count2 );
  }

  return;
}
/******************************************************************************/

void test0218 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0218 tests COLLATZ_COUNT_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
  int i_max;
  int j_max;
  int n;

  printf ( "\n" );
  printf ( "TEST0218:\n" );
  printf ( "  COLLATZ_COUNT_MAX(N) returns the length of the\n" );
  printf ( "  longest Collatz sequence from 1 to N.\n" );
  printf ( "\n" );
  printf ( "         N     I_MAX     J_MAX\n" );
  printf ( "\n" );

  n = 10;

  while ( n <= 100000 )
  {
    collatz_count_max ( n, &i_max, &j_max );

    printf ( "  %8d  %8d  %8d\n", n, i_max, j_max );

    n = n * 10;
  }

  return;
}
/******************************************************************************/

void test024 ( )

/******************************************************************************/
/*
  Purpose:

    TEST024 tests COMB_ROW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int c[N+1];
  int i;
  int j;
  bool next;

  printf ( "\n" );
  printf ( "TEST024\n" );
  printf ( "  COMB_ROW computes a row of Pascal's triangle.\n" );
  printf ( "\n" );

  next = false;

  for ( i = 0; i <= N; i++ )
  {
    comb_row ( next, i, c );
    next = true;
    printf ( "  %2d  ", i );
    for ( j = 0; j <= i; j++ )
    {
      printf ( "%5d", c[j] );
    }
    printf ( "\n" );
  }

  return;
# undef N
}
/******************************************************************************/

void test0243 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0243 tests COS_POWER_INT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0243:\n" );
  printf ( "  COS_POWER_INT computes the integral of the N-th power\n" );
  printf ( "  of the cosine function.\n" );
  printf ( "  COS_POWER_INT_VALUES returns selected values.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        Exact    Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = cos_power_int ( a, b, n );

    printf ( "  %8g  %8g  %6d  %12g  %12g\n", a, b, n, fx, fx2 );
  }
  return;
}
/******************************************************************************/

void test025 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST025 tests ERROR_F.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 August 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST025:\n" );
  printf ( "  ERROR_F evaluates the error function.\n" );
  printf ( "  ERF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X      Exact F       ERF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = error_f ( x );

    printf ( "  %8f  %14f  %14f\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test0255 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0255 tests ERROR_F_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2010

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST0255:\n" );
  printf ( "  ERROR_F_INVERSE inverts the error function.\n" );
  printf ( "  ERF_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "    FX           X1           X2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x1, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = error_f_inverse ( fx );

    printf ( "  %8f  %14f  %14f\n", fx, x1, x2 );
  }

  return;
}
/******************************************************************************/

void test026 ( )

/******************************************************************************/
/*
  Purpose:

    TEST026 tests EULER_NUMBER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int c1;
  int c2[13];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST026\n" );
  printf ( "  EULER_NUMBER computes Euler numbers.\n" );
  printf ( "  EULER_NUMBER_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact   EULER_NUMBER\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    euler_number ( n, c2 );

    printf ( "  %4d  %12d  %12d\n", n, c1, c2[n] );

  }
 
  return;
}
/******************************************************************************/

void test0265 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0265 tests EULER_NUMBER2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int c1;
  int c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0265\n" );
  printf ( "  EULER_NUMBER2 computes Euler numbers.\n" );
  printf ( "  EULER_NUMBER_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact   EULER_NUMBER2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = euler_number2 ( n );

    printf ( "  %4d  %12d  %12d\n", n, c1, c2 );

  }
 
  return;
}
/******************************************************************************/

void test028 ( )

/******************************************************************************/
/*
  Purpose:

    TEST028 tests EULER_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double f;
  int i;
  int n = 15;
  double x;

  x = 0.5;
 
  printf ( "\n" );
  printf ( "TEST028\n" );
  printf ( "  EULER_POLY evaluates Euler polynomials.\n" );
  printf ( "\n" );
  printf ( "  N         X              F(X)\n" );
  printf ( "\n" );
   
  for ( i = 0; i <= n; i++ )
  {
    f = euler_poly ( i, x );

    printf ( "  %2d  %14g  %14g\n", i, x, f );
  }
 
  return;
}
/******************************************************************************/

void test027 ( )

/******************************************************************************/
/*
  Purpose:

    TEST027 tests EULERIAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
# define N 7

  int e[N*N];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST027\n" );
  printf ( "  EULERIAN evaluates Eulerian numbers.\n" );
  printf ( "\n" );
 
  eulerian ( N, e );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      printf ( "%6d  ", e[i+j*N] );
    }
    printf ( "\n" );
  }
 
  return;
# undef N
}
/******************************************************************************/

void test029 ( )

/******************************************************************************/
/*
  Purpose:

    TEST029 tests F_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;

  printf ( "\n" );
  printf ( "TEST029\n" );
  printf ( "  F_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  F function.\n" );
  printf ( "\n" );
  printf ( "     N   F(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    f = f_hofstadter ( i );

    printf ( "  %6d  %6d\n", i, f );
  }

  return;
}
/******************************************************************************/

void test031 ( )

/******************************************************************************/
/*
  Purpose:

    TEST031 tests FIBONACCI_DIRECT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;
  int n = 20;

  printf ( "\n" );
  printf ( "TEST031\n" );
  printf ( "  FIBONACCI_DIRECT evalutes a Fibonacci number directly.\n" );
  printf ( "\n" );
  
  for ( i = 1; i <= n; i++ )
  {
    f = fibonacci_direct ( i );

    printf ( "  %6d  %10d\n", i, f );
  }
 
  return;
}
/******************************************************************************/

void test032 ( )

/******************************************************************************/
/*
  Purpose:

    TEST032 tests FIBONACCI_FLOOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;
  int n;

  printf ( "\n" );
  printf ( "TEST032\n" );
  printf ( "  FIBONACCI_FLOOR computes the largest Fibonacci number\n" );
  printf ( "  less than or equal to a given positive integer.\n" );
  printf ( "\n" );
  printf ( "     N  Fibonacci  Index\n" );
  printf ( "\n" );

  for ( n = 1; n <= 20; n++ )
  {
    fibonacci_floor ( n, &f, &i );

    printf ( "  %6d  %6d  %6d\n", n, f, i );
  }
 
  return;
}
/******************************************************************************/

void test033 ( )

/******************************************************************************/
/*
  Purpose:

    TEST033 tests FIBONACCI_RECURSIVE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
# define N 20

  int f[N];
  int i;

  printf ( "\n" );
  printf ( "TEST033\n" );
  printf ( "  FIBONACCI_RECURSIVE computes the Fibonacci sequence.\n" );
  printf ( "\n" );
 
  fibonacci_recursive ( N, f );
 
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %10d\n", i, f[i] );
  }
 
  return;
# undef N
}
/******************************************************************************/

void test034 ( )

/******************************************************************************/
/*
  Purpose:

    TEST034 tests G_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "TEST034\n" );
  printf ( "  G_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  G function.\n" );
  printf ( "\n" );
  printf ( "     N   G(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, g_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void test037 ( )

/******************************************************************************/
/*
  Purpose:

    TEST037 tests GEGENBAUER_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST037\n" );
  printf ( "  GEGENBAUER_POLY evaluates the Gegenbauer polynomials.\n" );
  printf ( "  GEGENBAUER_POLY_VALUES returns some exact values of\n" );
  printf ( "  the Gegenbauer polynomials.\n" );
  printf ( "\n" );
  printf ( "        N       A       X       GPV      GEGENBAUER\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {

    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

    gegenbauer_poly ( n, a, x, c );
    fx2 = c[n];

    printf ( "  %6d  %10g  %10g  %14g  %14g\n", n, a, x, fx, fx2 );

    free ( c );
  }
 
  return;
}
/******************************************************************************/

void test052 ( )

/******************************************************************************/
/*
  Purpose:

    TEST052 tests GEN_LAGUERRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
# define N 10
# define N_TEST 6

  double alpha;
  double alpha_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double c[N+1];
  int i;
  int j;
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  printf ( "\n" );
  printf ( "TEST052\n" );
  printf ( "  GEN_LAGUERRE_POLY evaluates the generalized Laguerre\n" );
  printf ( "  functions.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    alpha = alpha_test[i];

    printf ( "\n" );
    printf ( "  Table of L(N,ALPHA,X) for\n" );
    printf ( "\n" );
    printf ( "    N(max) = %d\n", N );
    printf ( "    ALPHA =  %g\n", alpha );
    printf ( "    X =      %g\n", x );
    printf ( "\n" );
  
    gen_laguerre_poly ( N, alpha, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      printf ( "  %6d  %14g\n", j, c[j] );
    }
  }
 
  return;
# undef N
# undef N_TEST
}
/******************************************************************************/

void test038 ( )

/******************************************************************************/
/*
  Purpose:

    TEST038 tests GUD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST038:\n" );
  printf ( "  GUD evaluates the Gudermannian function.\n" );
  printf ( "  GUD_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X      Exact F       GUD(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gud ( x );

    printf ( "  %10.6g  %10.6g  %10.6g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test040 ( )

/******************************************************************************/
/*
  Purpose:

    TEST040 tests H_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "TEST040\n" );
  printf ( "  H_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  H function.\n" );

  printf ( "\n" );
  printf ( "     N   H(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, h_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void test039 ( )

/******************************************************************************/
/*
  Purpose:

    TEST039 tests HAIL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "TEST039\n" );
  printf ( "  HAIL(I) computes the length of the hail sequence\n" );
  printf ( "  for I, also known as the 3*N+1 sequence.\n" );
  printf ( "\n" );
  printf ( "  I,  HAIL(I)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    printf ( "  %4d  %6d\n", i,  hail ( i ) );
  }
 
  return;
}
/******************************************************************************/

void test041 ( )

/******************************************************************************/
/*
  Purpose:

    TEST041 tests HERMITE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST041:\n" );
  printf ( "  HERMITE_POLY evaluates the Hermite polynomial.\n" );
  printf ( "  HERMITE_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       H(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    hermite_poly ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test042 ( )

/******************************************************************************/
/*
  Purpose:

    TEST042 tests HERMITE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST042\n" );
  printf ( "  HERMITE_POLY_COEF determines Hermite polynomial coefficients.\n" );

  hermite_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  H(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }

  return;
# undef N
}
/******************************************************************************/

void test0427 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0427 tests I4_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "TEST0427\n" );
  printf ( "  I4_CHOOSE evaluates C(N,K).\n" );
  printf ( "\n" );
  printf ( "   N     K    CNK\n" );
  printf ( "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      printf ( "  %6d  %6d  %6d\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void test043 ( )

/******************************************************************************/
/*
  Purpose:

    TEST043 tests I4_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  int fn;
  int fn2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST043:\n" );
  printf ( "  I4_FACTORIAL evaluates the factorial function.\n" );
  printf ( "  I4_FACTORIAL_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       I4_FACTORIAL(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial ( n );

    printf ( "  %4d  %12d  %12d\n", n, fn, fn2 );

  }

  return;
}
/******************************************************************************/

void test044 ( )

/******************************************************************************/
/*
  Purpose:

    TEST044 tests I4_FACTORIAL2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int fn;
  int fn2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST044:\n" );
  printf ( "  I4_FACTORIAL2 evaluates the double factorial function.\n" );
  printf ( "\n" );
  printf ( "   N   Exact  I4_FACTORIAL2(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial2 ( n );

    printf ( "  %4d  %8d  %8d\n", n, fn, fn2 );
  }

  return;
}
/******************************************************************************/

void test045 ( )

/******************************************************************************/
/*
  Purpose:

    TEST045 tests PARTITION_COUNT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST045:\n" );
  printf ( "  For the number of partitions of an integer,\n" );
  printf ( "  PARTITION_COUNT_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N       Exact F\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %10d  %10d\n", n, c );
  }

  return;
}
/******************************************************************************/

void test046 ( )

/******************************************************************************/
/*
  Purpose:

    TEST046 tests I4_PARTITION_DISTINCT_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int c2;
  int n;
  int n_data;
  int n_max = 20;

  printf ( "\n" );
  printf ( "TEST046:\n" );
  printf ( "  For the number of partitions of an integer\n" );
  printf ( "  into distinct parts,\n" );
  printf ( "  I4_PARTITION_DISTINCT_COUNT computes any value.\n" );
  printf ( "  PARTITION_DISTINCT_COUNT_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N       Exact F    Q(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n_max < n )
    {
      continue;
    }

    c2 = i4_partition_distinct_count ( n );

    printf ( "  %10d  %10d  %10d\n", n, c, c2 );
  }

  return;
}
/******************************************************************************/

void test047 ( )

/******************************************************************************/
/*
  Purpose:

    TEST047 tests I4_POCHHAMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "TEST047\n" );
  printf ( "  I4_POCHHAMMER evaluates the integer Pochhammer function.\n" );
  printf ( "\n" );
  printf ( "   I   J   I4_Pochhammer(I,J)\n" );
  printf ( "\n" );

  i = 3;
  j = 3;
  k = i4_pochhammer ( i, j );

  printf ( "  %4d  %4d  %4d\n", i, j, k );

  i = 3;
  j = 4;
  k = i4_pochhammer ( i, j );

  printf ( "  %4d  %4d  %4d\n", i, j, k );

  i = 3;
  j = 5;
  k = i4_pochhammer ( i, j );

  printf ( "  %4d  %4d  %4d\n", i, j, k );

  return;
}
/******************************************************************************/

void test048 ( )

/******************************************************************************/
/*
  Purpose:

    TEST048 tests I4_IS_TRIANGULAR, I4_TO_TRIANGLE and TRIANGLE_TO_I4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int i2;
  int j;
  int k;
  int k2;
  bool l;

  printf ( "\n" );
  printf ( "TEST048\n" );
  printf ( "  I4_TO_TRIANGLE converts a linear index to a\n" );
  printf ( "  triangular one.\n" );
  printf ( "  TRIANGLE_TO_I4 converts a triangular index to a\n" );
  printf ( "  linear one.\n" );
  printf ( "  I4_IS_TRIANGULAR returns 0 or 1 depending on\n" );
  printf ( "  whether I is triangular.\n" );
  printf ( "\n" );
  printf ( "   I  =>   J   K  =>   I   0/1\n" );
  printf ( "\n" );

  for ( i = 0; i <= 20; i++ )
  {
    i4_to_triangle ( i, &j, &k );

    i2 = triangle_to_i4 ( j, k );

    l = i4_is_triangular ( i );

    printf ( "  %4d  %4d  %4d  %4d  %1d\n", i, j, k, i2, l );
  }
 
  return;
}
/******************************************************************************/

void test049 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST049 tests JACOBI_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST049:\n" );
  printf ( "  JACOBI_POLY computes values of the Jacobi polynomial..\n" );
  printf ( "  JACOBI_POLY_VALUES returns values of the Jacobi polynomial\n" );
  printf ( "\n" );
  printf ( "       N       A       B      X       JPV      JACOBI\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_poly_values ( &n_data, &n, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = jacobi_poly ( n, a, b, x );
    fx2 = c[n];

    printf ( "  %8d  %8f  %8f  %10.4f  %14.6f  %14.6f\n", n, a, b, x, fx, fx2 );

    free ( c );
  }

  return;
}
/******************************************************************************/

void test050 ( )

/******************************************************************************/
/*
  Purpose:

    TEST050 tests JACOBI_SYMBOL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_TEST 4

  int i;
  int p;
  int ptest[N_TEST] = { 3, 9, 10, 12 };
  int q;

  printf ( "\n" );
  printf ( "TEST050\n" );
  printf ( "  JACOBI_SYMBOL computes the Jacobi symbol\n" );
  printf ( "  (Q/P), which records if Q is a quadratic\n" );
  printf ( "  residue modulo the number P.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    printf ( "\n" );
    printf ( "Jacobi Symbols for P = %d\n", p );
    printf ( "\n" );
    for ( q = 0; q <= p; q++ )
    {
      printf ( "  %8d  %8d  %8d\n", p, q, jacobi_symbol ( q, p ) );
    }
  }

  return;
# undef N_TEST
}
/******************************************************************************/

void test0505 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0505 tests KRAWTCHOUK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2
# define N 5

  int i;
  int j;
  int m;
  int n;
  double p;
  double p_test[TEST_NUM] = { 0.25, 0.5 };
  int test;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "TEST0505:\n" );
  printf ( "  KRAWTCHOUK evaluates Krawtchouk polynomials.\n" );
  printf ( "\n" );
  printf ( "        N         P         X          M      K(N,P,X,M)\n" );
  printf ( "\n" );

  m = 5;
  n = N;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];

    printf ( "\n" );

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      krawtchouk ( n, p, x, m, value );

      printf ( "\n" );
      for ( i = 0; i <= 5; i++ )
      {

        printf ( "  %8d  %8g  %8g  %8d  %14g\n", i, p, x, m, value[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test051 ( )

/******************************************************************************/
/*
  Purpose:

    TEST051 tests LAGUERRE_ASSOCIATED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 6
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  int m;
  int m_test[N_TEST] = { 0, 0, 1, 2, 3, 1 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  printf ( "\n" );
  printf ( "TEST051\n" );
  printf ( "  LAGUERRE_ASSOCIATED evaluates the associated Laguerre\n" );
  printf ( "  polynomials.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    m = m_test[i];
    x = x_test[i];

    printf ( "\n" );
    printf ( "  Table of L(N,M,X) for\n" );
    printf ( "\n" );
    printf ( "  N(max) = %d\n", N );
    printf ( "  M      = %d\n", m );
    printf ( "  X =      %g\n", x );
    printf ( "\n" );
 
    laguerre_associated ( N, m, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      printf ( "  %6d  %14g\n", j, c[j] );
    }
  }

  return;
# undef N
# undef N_TEST
}
/******************************************************************************/

void test054 ( )

/******************************************************************************/
/*
  Purpose:

    TEST054 tests LAGUERRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST054:\n" );
  printf ( "  LAGUERRE_POLY evaluates the Laguerre polynomial.\n" );
  printf ( "  LAGUERRE_POLYNOMIAL_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       L(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    laguerre_poly ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test055 ( )

/******************************************************************************/
/*
  Purpose:

    TEST055 tests LAGUERRE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  double fact;
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST055\n" );
  printf ( "  LAGUERRE_POLY_COEF determines Laguerre \n" );
  printf ( "  polynomial coefficients.\n" );

  laguerre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  L(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
 
  for ( i = 0; i <= N; i++ )
  {
    fact = r8_factorial ( i );
    printf ( "\n" );
    printf ( "  Factorially scaled L(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
  return;
# undef N
}
/******************************************************************************/

void test057 ( )

/******************************************************************************/
/*
  Purpose:

    TEST057 tests LEGENDRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fp2[N_MAX+1];
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST057:\n" );
  printf ( "  LEGENDRE_POLY evaluates the Legendre PN function.\n" );
  printf ( "  LEGENDRE_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       P(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_poly ( n, x, fx2, fp2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test058 ( )

/******************************************************************************/
/*
  Purpose:

    TEST058 tests LEGENDRE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST58\n" );
  printf ( "  LEGENDRE_POLY_COEF determines the Legendre P \n" );
  printf ( "  polynomial coefficients.\n" );

  legendre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  P(%8d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
 
  return;
# undef N
}
/******************************************************************************/

void test059 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST059 tests LEGENDRE_ASSOCIATED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST059:\n" );
  printf ( "  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.\n" );
  printf ( "  LEGENDRE_ASSOCIATED_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      N       M    X     Exact F     PNM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated ( n, m, x, fx2 );

    printf ( "  %8d  %8d  %8f  %14f  %14f\n", n, m, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test0595 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0595 tests LEGENDRE_ASSOCIATED_NORMALIZED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0595:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED evaluates \n" );
  printf ( "  normalized associated Legendre functions.\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      N       M    X     Exact F     PNM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated_normalized ( n, m, x, fx2 );

    printf ( "  %8d  %8d  %8f  %14f  %14f\n", n, m, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test060 ( )

/******************************************************************************/
/*
  Purpose:

    TEST060 tests LEGENDRE_FUNCTION_Q.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST060:\n" );
  printf ( "  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.\n" );
  printf ( "  LEGENDRE_FUNCTION_Q_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       Q(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_function_q ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test061 ( )

/******************************************************************************/
/*
  Purpose:

    TEST061 tests LEGENDRE_SYMBOL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N_TEST 4

  int i;
  int l;
  int p;
  int ptest[N_TEST] = { 7, 11, 13, 17 };
  int q;

  printf ( "\n" );
  printf ( "TEST061\n" );
  printf ( "  LEGENDRE_SYMBOL computes the Legendre\n" );
  printf ( "  symbol (Q/P) which records whether Q is \n" );
  printf ( "  a quadratic residue modulo the prime P.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    printf ( "\n" );
    printf ( "  Legendre Symbols for P = %d\n", p );
    printf ( "\n" );
    for ( q = 0; q <= p; q++ )
    {
      printf ( "  %8d  %8d  %8d\n", p, q, legendre_symbol ( q, p ) );
    }
  }

  return;
# undef N_TEST
}
/******************************************************************************/

void test0615 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0615 tests LERCH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int s;
  double z;

  printf ( "\n" );
  printf ( "TEST0615:\n" );
  printf ( "  LERCH evaluates the Lerch function.\n" );
  printf ( "  LERCH_VALUES returns some tabulated values.\n" );
  printf ( "\n" );
  printf ( "       Z       S       A         Lerch           Lerch\n" );
  printf ( "                             Tabulated        Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lerch ( z, s, a );

    printf ( "  %8g  %4d  %8g  %14g  %14g\n", z, s, a, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test0365 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0365 tests LGAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0365:\n" );
  printf ( "  LGAMMA is a C math library function which evaluates\n" );
  printf ( "  the logarithm of the Gamma function.\n" );
  printf ( "  GAMMA_LOG_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       LGAMMA(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lgamma ( x );

    printf ( "  %8g  %10g  %10g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void test062 ( )

/******************************************************************************/
/*
  Purpose:

    TEST062 tests LOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int a[N+1];
  int i;

  printf ( "\n" );
  printf ( "TEST062\n" );
  printf ( "  LOCK counts the combinations on a button lock.\n" );
  printf ( "\n" );
  printf ( "     I      LOCK(I)\n" );
  printf ( "\n" );

  lock ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %10d\n", i, a[i] );
  }
 
  return;
# undef N
}
/******************************************************************************/

void test0623 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0623 tests MEIXNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N 5
# define TEST_NUM 3

  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.0 };
  double c;
  double c_test[TEST_NUM] = { 0.125, 0.25, 0.5 };
  int i;
  int j;
  int n;
  int test;
  double v[N+1];
  double x;

  printf ( "\n" );
  printf ( "TEST0623:\n" );
  printf ( "  MEIXNER evaluates Meixner polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      BETA         C         X        M(N,BETA,C,X)\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    beta = beta_test[test];
    c = c_test[test];

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      meixner ( n, beta, c, x, v );

      printf ( "\n" );

      for ( i = 0; i <= n; i++ )
      {
        printf ( "  %8d  %8g  %8g  %8g  %14g\n", i, beta, c, x, v[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test0625 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0625 tests MERTENS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0625\n" );
  printf ( "  MERTENS computes the Mertens function.\n" );
  printf ( "  MERTENS_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      N   Exact   MERTENS(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
     mertens_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %10d  %10d\n", n, c, mertens ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test063 ( )

/******************************************************************************/
/*
  Purpose:

    TEST063 tests MOEBIUS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST063\n" );
  printf ( "  MOEBIUS computes the Moebius function.\n" );
  printf ( "  MOEBIUS_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "      N   Exact   MOEBIUS(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
     moebius_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %8d  %10d  %10d\n", n, c, moebius ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test0635 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0635 tests MOTZKIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int a[N+1];
  int i;

  printf ( "\n" );
  printf ( "TEST0635\n" );
  printf ( "  MOTZKIN computes the Motzkin numbers A(0:N).\n" );
  printf ( "  A(N) counts the paths from (0,0) to (N,0).\n" );
  printf ( "\n" );
  printf ( "  I,  A(I)\n" );
  printf ( "\n" );

  motzkin ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %10d\n", i, a[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test064 ( )

/******************************************************************************/
/*
  Purpose:

    TEST064 tests OMEGA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST064\n" );
  printf ( "  OMEGA computes the OMEGA function.\n" );
  printf ( "  OMEGA_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "          N   Exact   OMEGA(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %10d  %10d\n", n, c, omega ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test065 ( )

/******************************************************************************/
/*
  Purpose:

    TEST065 tests PENTAGON_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST065\n" );
  printf ( "  PENTAGON_NUM computes the pentagonal numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, pentagon_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test066 ( )

/******************************************************************************/
/*
  Purpose:

    TEST066 tests PHI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST066\n" );
  printf ( "  PHI computes the PHI function.\n" );
  printf ( "  PHI_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N   Exact   PHI(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, phi ( n ) );

  }
 
  return;
}
/******************************************************************************/

void test0665 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0665 tests POLY_BERNOULLI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int b;
  int k;
  int n;

  printf ( "\n" );
  printf ( "TEST0665\n" );
  printf ( "  POLY_BERNOULLI computes the poly-Bernoulli numbers\n" );
  printf ( "  of negative index, B_n^(-k)\n" );
  printf ( "\n" );
  printf ( "   N   K    B_N^(-K)\n" );
  printf ( "\n" );

  for ( k = 0; k <= 6; k++ )
  {
    printf ( "\n" );
    for ( n = 0; n <= 6; n++ )
    {
      b = poly_bernoulli ( n, k );

      printf ( "  %2d  %2d  %12d\n", n, k, b );
    }
  }

  return;
}
/******************************************************************************/

void test0667 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0667 tests POLY_COEF_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int degree;
  int dim;
  int n;

  printf ( "\n" );
  printf ( "TEST0667\n" );
  printf ( "  POLY_COEF_COUNT counts the number of coefficients\n" );
  printf ( "  in a polynomial of degree DEGREE and dimension DIM.\n" );
  printf ( "\n" );
  printf ( " Dimension    Degree     Count\n" );

  for ( dim = 1; dim <= 10; dim = dim + 3 )
  {
    printf ( "\n" );
    for ( degree = 0; degree <= 5; degree++ )
    {
      printf ( "  %8d  %8d  %8d\n", dim, degree, poly_coef_count ( dim, degree ) );
    }
  }

  return;
}
/******************************************************************************/

void test067 ( )

/******************************************************************************/
/*
  Purpose:

    TEST067 tests PYRAMID_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST067\n" );
  printf ( "  PYRAMID_NUM computes the pyramidal numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, pyramid_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test0675 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0675 tests R8_ACOSH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST0675\n" );
  printf ( "  R8_ACOSH computes the inverse hyperbolic cosine\n" );
  printf ( "  of a given value.\n" );
  printf ( "\n" );
  printf ( "         X  R8_ACOSH(X)      COSH(R8_ACOSH(X))\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = r8_acosh ( x );
    x2 = cosh ( a );

    printf ( "  %10g  %10g  %10g\n", x, a, x2 );
  }

  return;
}
/******************************************************************************/

void test0676 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0676 tests R8_ASINH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST0676\n" );
  printf ( "  R8_ASINH computes the inverse hyperbolic sine\n" );
  printf ( "  of a given value.\n" );
  printf ( "\n" );
  printf ( "         X     R8_ASINH(X)     SINH(R8_ASINH(X))\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = r8_asinh ( x );
    x2 = sinh ( a );

    printf ( "  %10f  %10f  %10f\n", x, a, x2 );
  }

  return;
}
/******************************************************************************/

void test06765 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06765 tests R8_ATANH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST06765\n" );
  printf ( "  R8_ATANH computes the inverse hyperbolic tangent\n" );
  printf ( "  of a given value.\n" );
  printf ( "\n" );
  printf ( "         X     R8_ATANH(X)     TANH(R8_ATANH(X))\n" );
  printf ( "\n" );

  for ( i = -2; i <= 9; i++ )
  {
    x = ( ( double ) i ) / 10.0;
    a = r8_atanh ( x );
    x2 = tanh ( a );

    printf ( "  %14g  %14g  %14g\n", x, a, x2 ); 
  }

  return;
}
/******************************************************************************/

void test022 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST022 tests R8_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
  double cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "TEST022\n" );
  printf ( "  R8_CHOOSE evaluates C(N,K) using real arithmetic.\n" );
  printf ( "\n" );
  printf ( "       N       K          CNK\n" );
  printf ( "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );

      printf ( "  %6d  %6d  %10g\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void test068 ( )

/******************************************************************************/
/*
  Purpose:

    TEST068 tests R8_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  double fn;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "TEST068:\n" );
  printf ( "  R8_FACTORIAL evaluates the factorial function.\n" );
  printf ( "  R8_FACTORIAL_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N       Exact F       R8_FACTORIAL(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %14g  %14g\n", n, fn, r8_factorial ( n ) );
  }

  return;
}
/******************************************************************************/

void test0685 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0685 tests R8_FACTORIAL_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double fn;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "TEST0685:\n" );
  printf ( "  R8_FACTORIAL_LOG evaluates the logarithm of the\n" );
  printf ( "  factorial function.\n" );
  printf ( "  R8_FACTORIAL_LOG_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N	   Exact F	 R8_FACTORIAL_LOG(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %5d  %14g  %14g\n", n, fn, r8_factorial_log ( n ) );
  }

  return;
}
/******************************************************************************/

void test06855 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06855 tests R8_GAMMA, GAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST06855:\n" );
  printf ( "  R8_GAMMA evaluates the Gamma function.\n" );
  printf ( "  GAMMA is a C MATH routine for the Gamma function.\n" );
  printf ( "  GAMMA_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "         X                  Gamma(X)         "  );
  printf ( "         Gamma(X)         "  );
  printf ( "         Gamma(X)\n" );
  printf ( "                         (Tabulated)         "  );
  printf ( "       (R8_GAMMA)         " );
  printf ( "       (GAMMA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_gamma ( x );

    fx3 = gamma ( x );

    printf ( "  %8g  %24.16g  %24.16g  %24.16g\n", fx, fx2, fx3 );
  }

  return;
}
/******************************************************************************/

void test036 ( )

/******************************************************************************/
/*
  Purpose:

    TEST036 tests R8_GAMMA_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST036:\n" );
  printf ( "  R8_GAMMA_LOG evaluates the logarithm of the Gamma function.\n" );
  printf ( "  GAMMA_LOG_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       R8_GAMMA_LOG(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_gamma_log ( x );

    printf ( "  %8g  %10g  %10g\n", x, fx, fx2 );
  }
  return;
}
/******************************************************************************/

void test0425 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0425 tests R8_HYPER_2F1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

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
  printf ( " TEST0425:\n" );
  printf ( "   R8_HYPER_2F1 evaluates the hypergeometric function 2F1.\n" );
  printf ( "\n" );
  printf ( "      A       B       C       X      " );
  printf ( " 2F1                       2F1                     DIFF\n" );
  printf ( "                                     " );
  printf ( "(tabulated)               (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_hyper_2f1 ( a, b, c, x );

    printf ( "  %6f  %6f  %6f  %6f  %24.16g  %24.16g  %10.4g\n",
      a, b, c, x, fx, fx2, r8_abs ( fx - fx2 )  );
  }
  return;
}
/******************************************************************************/

void test06856 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06856 tests R8_PSI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST06856:\n" );
  printf ( "  R8_PSI evaluates the Psi function.\n" );
  printf ( "  PSI_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "         X                  Psi(X)           " );
  printf ( "         Psi(X)          DIFF\n" );
  printf ( "                         (Tabulated)         " );
  printf ( "       (R8_PSI)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_psi ( x );

    printf ( "  %8.2g  %24.16g  %24.16g  %10.4g\n", x, fx, fx2, r8_abs ( fx - fx2 ) );

  }

  return;
}
/******************************************************************************/

void test069 ( )

/******************************************************************************/
/*
  Purpose:

    TEST069 tests SIGMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST069\n" );
  printf ( "  SIGMA computes the SIGMA function.\n" );
  printf ( "  SIGMA_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N   Exact   SIGMA(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    sigma_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, sigma ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test0695 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0695 tests SIN_POWER_INT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0695:\n" );
  printf ( "  SIN_POWER_INT computes the integral of the N-th power\n" );
  printf ( "  of the sine function.\n" );
  printf ( "  SIN_POWER_INT_VALUES returns selected values.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        Exact    Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = sin_power_int ( a, b, n );

    printf ( "  %8g  %8g  %6d  %12g  %12g\n", a, b, n, fx, fx2 );
  }
  return;
}
/******************************************************************************/

void test0696 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0696 tests SLICE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 August 2011

  Author:

    John Burkardt
*/
{
# define DIM_MAX 5
# define SLICE_MAX 8

  int dim_max = DIM_MAX;
  int dim_num;
  int p[DIM_MAX*SLICE_MAX];
  int piece_num;
  int slice_max = SLICE_MAX;
  int slice_num;

  printf ( "\n" );
  printf ( "TEST0696:\n" );
  printf ( "  SLICE determines the maximum number of pieces created\n" );
  printf ( "  by SLICE_NUM slices in a DIM_NUM space.\n" );

  for ( dim_num = 1; dim_num <= dim_max; dim_num++ )
  {
    for ( slice_num = 1; slice_num <= slice_max; slice_num++ )
    {
      piece_num = slice ( dim_num, slice_num );
      p[dim_num-1+(slice_num-1)*dim_max] = piece_num;
    }
  }

  i4mat_print ( dim_max, slice_max, p, "  Slice Array:" );

  return;
# undef DIM_MAX
# undef SLICE_MAX
}
/******************************************************************************/

void test0697 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0697 tests SPHERICAL_HARMONIC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double c[N_MAX+1];
  int l;
  int m;
  int n_data;
  double phi;
  double s[N_MAX+1];
  double theta;
  double yi;
  double yi2;
  double yr;
  double yr2;

  printf ( "\n" );
  printf ( "TEST0697:\n" );
  printf ( "  SPHERICAL_HARMONIC evaluates spherical harmonic functions.\n" );
  printf ( "  SPHERICAL_HARMONIC_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "         N         M    THETA      PHI            YR            YI\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }

    spherical_harmonic ( l, m, theta, phi, c, s );

    yr2 = c[l];
    yi2 = s[l];

    printf ( "  %8d  %8d  %8g  %8g  %14g  %14g\n", l, m, theta, phi, yr, yi );

    printf ( "                                          %14g  %14g\n", yr2, yi2 );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test070 ( )

/******************************************************************************/
/*
  Purpose:

    TEST070 tests STIRLING1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m = 8;
  int n = 8;
  int *s1;

  printf ( "\n" );
  printf ( "TEST070\n" );
  printf ( "  STIRLING1: Stirling numbers of first kind.\n" );
  printf ( "    Get rows 1 through %d\n", m );
  printf ( "\n" );
 
  s1 = stirling1 ( m, n );
 
  for ( i = 0; i < m; i++ )
  {
    printf ( "%6d  ", i + 1 );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%6d  ", s1[i+j*m] );
    }
    printf ( "\n" );
  }

  free ( s1 );
 
  return;
}
/******************************************************************************/

void test071 ( )

/******************************************************************************/
/*
  Purpose:

    TEST071 tests STIRLING2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m = 8;
  int n = 8;
  int *s2;

  printf ( "\n" );
  printf ( "TEST071\n" );
  printf ( "  STIRLING2: Stirling numbers of second kind.\n" );
  printf ( "  Get rows 1 through %d\n", m );
  printf ( "\n" );
 
  s2 = stirling2 ( m, n );
 
  for ( i = 0; i < m; i++ )
  {
    printf ( "%6d  ", i + 1 );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%6d  ", s2[i+j*m] );
    }
    printf ( "\n" );
  }
 
  free ( s2 );

  return;
}
/******************************************************************************/

void test072 ( )

/******************************************************************************/
/*
  Purpose:

    TEST072 tests TAU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST072\n" );
  printf ( "  TAU computes the Tau function.\n" );
  printf ( "  TAU_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    tau_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, tau ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test073 ( )

/******************************************************************************/
/*
  Purpose:

    TEST073 tests TETRAHEDRON_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST073\n" );
  printf ( "  TETRAHEDRON_NUM computes the tetrahedron numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, tetrahedron_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test074 ( )

/******************************************************************************/
/*
  Purpose:

    TEST074 tests TRIANGLE_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST074\n" );
  printf ( "  TRIANGLE_NUM computes the triangular numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, triangle_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void test075 ( )

/******************************************************************************/
/*
  Purpose:

    TEST075 tests V_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int v;

  printf ( "\n" );
  printf ( "TEST075\n" );
  printf ( "  V_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  V function.\n" );
  printf ( "\n" );
  printf ( "     N   V(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, v_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void test076 ( )

/******************************************************************************/
/*
  Purpose:

    TEST076 tests VIBONACCI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
# define N 20
  int i;
  int j;
  int seed;
  int v1[N];
  int v2[N];
  int v3[N];

  printf ( "\n" );
  printf ( "TEST076\n" );
  printf ( "  VIBONACCI computes a Vibonacci sequence.\n" );
  printf ( "\n" );
  printf ( "  We compute the series 3 times.\n" );
  printf ( "\n" );
  printf ( "     I      V1      V2      V3\n" );
  printf ( "\n" );

  seed = 123456789;

  vibonacci ( N, &seed, v1 );
  vibonacci ( N, &seed, v2 );
  vibonacci ( N, &seed, v3 );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6d  %6d  %6d\n", i, v1[i], v2[i], v3[i] );
  } 

  return;
# undef N
}
/******************************************************************************/

void test077 ( )

/******************************************************************************/
/*
  Purpose:

    TEST077 tests ZECKENDORF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
# define M_MAX 20

  int i;
  int i_list[M_MAX];
  int j;
  int f_list[M_MAX];
  int f_sum;
  int m;
  int n;

  printf ( "\n" );
  printf ( "TEST077\n" );
  printf ( "  ZECKENDORF computes the Zeckendorf decomposition of\n" );
  printf ( "  an integer N into nonconsecutive Fibonacci numbers.\n" );
  printf ( "\n" );
  printf ( "   N Sum M Parts\n" );
  printf ( "\n" );

  for ( n = 1; n <= 100; n++ )
  {
    zeckendorf ( n, M_MAX, &m, i_list, f_list );

    printf ( "%4d  ", n );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%4d  ", f_list[j] );
    }
    printf ( "\n" );

  }

  return;
# undef M_MAX
}
/******************************************************************************/

void test0773 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0773 tests ZERNIKE_POLY and ZERNIKE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int m;
  int n;
  double rho;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "TEST0773\n" );
  printf ( "  ZERNIKE_POLY_COEF returns the coefficients of a\n" );
  printf ( "  Zernike polynomial.\n" );
  printf ( "  ZERNIKE_POLY evaluates a Zernike polynomial directly.\n" );
  printf ( "\n" );
  printf ( "  Table of polynomial coefficients:\n" );
  printf ( "\n" );
  printf ( "   N   M\n" );
  printf ( "\n" );

  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      printf ( "  %2d  %2d", n, m );
      for ( i = 0; i <= n; i++ )
      {
        printf ( "  %7g", c[i] );
      }
      printf ( "\n" );
      free ( c );
    }
  }

  rho = 0.987654321;

  printf ( "\n" );
  printf ( "  Z1: Compute polynomial coefficients,\n" );
  printf ( "  then evaluate by Horner's method;\n" );
  printf ( "  Z2: Evaluate directly by recursion.\n" );
  printf ( "\n" );
  printf ( "   N   M       Z1              Z2\n" );
  printf ( "\n" );

  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      z1 = r8poly_value ( n + 1, c, rho );

      z2 = zernike_poly ( m, n, rho );
      printf ( "  %2d  %2d  %16g  %16g\n", n, m, z1, z2 );

      free ( c );
    }
  }

  return;
}
/******************************************************************************/

void test0775 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0775 tests ZERNIKE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int m;
  int n;

  printf ( "\n" );
  printf ( "TEST0775\n" );
  printf ( "  ZERNIKE_POLY_COEF determines the Zernike\n" );
  printf ( "  polynomial coefficients.\n" );

  n = 5;

  for ( m = 0; m <= n; m++ )
  {
    c = zernike_poly_coef ( m, n );
    r8poly_print ( n, c, "  Zernike polynomial" );
    free ( c );
  }

  return;
}
/******************************************************************************/

void test078 ( )

/******************************************************************************/
/*
  Purpose:

    TEST078 tests ZETA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;
  double n_real;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "TEST078\n" );
  printf ( "  ZETA computes the Zeta function.\n" );
  printf ( "  ZETA_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "       N            exact Zeta         computed Zeta\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &z1 );

    if ( n_data == 0 )
    {
      break;
    }

    n_real = ( double ) n;

    z2 = zeta ( n_real );

    printf ( "  %6d  %20.14g  %20.14g\n", n, z1, z2 );
  }

  return;
}
