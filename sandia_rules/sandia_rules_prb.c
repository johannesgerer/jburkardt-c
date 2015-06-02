# include "sandia_rules.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

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
void test10 ( void );
void test11 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void test18 ( void );
void test19 ( void );
void test20 ( void );
void test21 ( void );
void test22 ( void );
void test23 ( int r );
void test24 ( void );
void test25 ( void );
void test26 ( void );
void test27 ( void );
void test28 ( void );
void test29 ( void );
void test30 ( void );

void test01_np ( void );
void test02_np ( void );
void test03_np ( void );
void test04_np ( void );
void test05_np ( void );
void test06_np ( void );
void test07_np ( void );
void test08_np ( void );
void test09_np ( void );
void test10_np ( void );
void test11_np ( void );
void test12_np ( void );
void test13_np ( void );
void test14_np ( void );
void test15_np ( void );
void test16_np ( void );
void test17_np ( void );
void test18_np ( void );
void test19_np ( void );
void test20_np ( void );
void test21_np ( void );
void test22_np ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SANDIA_RULES_PRB.

  Discussion:

    SANDIA_RULES_PRB tests the SANDIA_RULES library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  int r;

  timestamp ( );
  printf ( "\n" );
  printf ( "SANDIA_RULES_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SANDIA_RULES library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
/*
  TEST23 takes an input argument of R, a rule index.
*/
  r = 1;
  test23 ( r );
  r = 3;
  test23 ( r );
  r = 4;
  test23 ( r );
  r = 11;
  test23 ( r );

  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
/*
  Repeat tests, but now call "NP" versions of routines.
*/
  test01_np ( );
  test02_np ( );
  test03_np ( );
  test04_np ( );
  test05_np ( );
  test06_np ( );
  test07_np ( );
  test08_np ( );
  test09_np ( );

  test10_np ( );
  test11_np ( );
  test12_np ( );
  test13_np ( );
  test14_np ( );
  test15_np ( );
  test16_np ( );
  test17_np ( );
  test18_np ( );
  test19_np ( );

  test20_np ( );
  test21_np ( );
  test22_np ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SANDIA_RULES_PRB\n" );
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

    TEST01 tests CHEBYSHEV1_COMPUTE against CHEBYSHEV1_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) / sqrt ( 1 - x^2 ) dx.\n" );
  printf ( "\n" );
  printf ( "  CHEBYSHEV1_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev1_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = chebyshev1_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
 
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CHEBYSHEV1_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev1_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CHEBYSHEV2_COMPUTE against CHEBYSHEV2_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) * sqrt ( 1 - x^2 ) dx.\n" );
  printf ( "\n" );
  printf ( "  CHEBYSHEV2_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev2_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = chebyshev2_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CHEBYSHEV2_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev2_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CLENSHAW_CURTIS_COMPUTE against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int n_hi;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N up to\n" );
  printf ( "    N = ORDER+1 if ORDER is odd, or\n" );
  printf ( "    N = ORDER   if ORDER is even\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    clenshaw_curtis_compute ( order, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests CLENSHAW_CURTIS_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    clenshaw_curtis_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests FEJER2_COMPUTE against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int n_hi;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  FEJER2_COMPUTE computes a Fejer Type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N up to\n" );
  printf ( "    N = ORDER+1 if ORDER is odd, or\n" );
  printf ( "    N = ORDER   if ORDER is even\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    fejer2_compute ( order, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests FEJER2_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  FEJER2_COMPUTE computes a Fejer Type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    fejer2_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests GEGENBAUER_COMPUTE against GEGENBAUER_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.\n" );
  printf ( "\n" );
  printf ( "  GEGENBAUER_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gegenbauer_compute ( order, alpha, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gegenbauer_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
        printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6f  %14.6e\n", 
          order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests GEGENBAUER_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gegenbauer_compute ( order, alpha, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests GEN_HERMITE_COMPUTE against GEN_HERMITE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n" );
  printf ( "\n" );
  printf ( "  GEN_HERMITE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_hermite_compute ( order, alpha, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gen_hermite_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
        printf ( "  %8d  %8d  %14.6f  %14.6e  %14.6e  %14.6e\n", 
          order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests GEN_HERMITE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_hermite_compute ( order, alpha, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests GEN_LAGUERRE_COMPUTE against GEN_LAGUERRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n" );
  printf ( "\n" );
  printf ( "  GEN_LAGUERRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_laguerre_compute ( order, alpha, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gen_laguerre_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
        printf ( "  %8d  %8d  %14.6f  %14.6e  %14.6e  %14.6e\n", 
          order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests GEN_LAGUERRE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_laguerre_compute ( order, alpha, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests HERMITE_COMPUTE against HERMITE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  HERMITE_COMPUTE computes a Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n" );
  printf ( "\n" );
  printf ( "  HERMITE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    hermite_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = hermite_integral ( n );
 
      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests HERMITE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  HERMITE_COMPUTE computes a Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    hermite_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests JACOBI_COMPUTE against JACOBI_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  int test1;
  int test2;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  JACOBI_COMPUTE computes a Gauss-Jacobi rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n" );
  printf ( "\n" );
  printf ( "  JACOBI_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Beta            Estimate       Exact            Error\n" );

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];

      for ( order = 1; order <= order_max; order++ )
      {
        printf ( "\n" );

        f = ( double * ) malloc ( order * sizeof ( double ) );
        w = ( double * ) malloc ( order * sizeof ( double ) );
        x = ( double * ) malloc ( order * sizeof ( double ) );

        jacobi_compute ( order, alpha, beta, x, w );

        for ( n = 0; n <= 2 * order + 2; n = n + 1 )
        {
          exact = jacobi_integral ( n, alpha, beta );
 
          if ( n == 0 )
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = 1.0;
            }
          }
          else
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = pow ( x[i], n );
            }
          }
          estimate = 0.0;
          for ( i = 0; i < order; i++ )
          {
            estimate = estimate + w[i] * f[i];
          }
 
          error = r8_abs ( exact - estimate );
  
          printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6f  %14.6f  %14.6e\n", 
            order, n, alpha, beta, estimate, exact, error );
        }
        free ( f );
        free ( w );
        free ( x );
      }
    }
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test18 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests JACOBI_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int order;
  int order_max = 10;
  int test1;
  int test2;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  JACOBI_COMPUTE computes a Gauss-Jacobi rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n" );

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];

      for ( order = 1; order <= order_max; order++ )
      {
        printf ( "\n" );
        printf ( "  Order = %d\n", order );
        printf ( "  ALPHA = %f\n", alpha );
        printf ( "  BETA =  %f\n", beta );

        w = ( double * ) malloc ( order * sizeof ( double ) );
        x = ( double * ) malloc ( order * sizeof ( double ) );

        jacobi_compute ( order, alpha, beta, x, w );

        for ( i = 0; i < order; i =i + 1 )
        {
          printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
        }
        free ( w );
        free ( x );
      }
    }
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test19 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests LAGUERRE_COMPUTE against LAGUERRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  LAGUERRE_COMPUTE computes a Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n" );
  printf ( "\n" );
  printf ( "  LAGUERRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    laguerre_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = laguerre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6e  %14.6e  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test20 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests LAGUERRE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    laguerre_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test21 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests LEGENDRE_COMPUTE against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  LEGENDRE_COMPUTE computes a Gauss-Legendre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    legendre_compute ( order, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test22 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests LEGENDRE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  LEGENDRE_COMPUTE computes a Gauss-Legendre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    legendre_compute ( order, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  return;
}
/******************************************************************************/

void test23 ( int r )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests LEVEL_GROWTH_TO_ORDER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int R, the index of the rule to be examined.
*/
{
  int dim;
  int dim_num = 11;
  int g;
  int *growth;
  int *level;
  int *order;
  int *rule;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST23\n" );
  fprintf ( stdout, "  LEVEL_GROWTH_TO_ORDER uses Level, Growth and Rule\n" );
  fprintf ( stdout, "  to determine the orders of each entry of a vector of 1D rules.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Here we examine rule %d.\n", r );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, 
    "       LEVEL:0     1     2     3     4     5     6     7     8     9    10\n" );
  fprintf ( stdout, "GROWTH\n" );

  growth = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( dim = 0; dim < dim_num; dim++ )
  {
    rule[dim] = r;
  }

  for ( g = 0; g <= 6; g++ )
  {
    if ( r == 3 || r == 11 )
    {
      if ( g == 1 || g == 2 || g == 3 )
      {
        continue;
      }
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      growth[dim] = g;
    }
    level_growth_to_order ( dim_num, level, rule, growth, order );

    fprintf ( stdout, "  %4d  ", g );
 
    for ( dim = 0; dim < dim_num; dim++ )
    {
      fprintf ( stdout, "  %4d", order[dim] );
    }
    fprintf ( stdout, "\n" );
  }

  free ( growth );
  free ( level );
  free ( order );
  free ( rule );
 
  return;
}
/******************************************************************************/

void test24 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests LEVEL_TO_ORDER_DEFAULT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  LEVEL_TO_ORDER_DEFAULT uses a default rule to\n" );
  printf ( "  determine the order of a rule from its level.\n" );
  printf ( "\n" );
  printf ( 
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n" );
  printf ( "\n" );

  level = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( r = 1; r <= 13; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    level_to_order_default ( dim_num, level, rule, order );

    printf ( "  %4d  ", r );
 
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %4d", order[dim] );
    }
    printf ( "\n" );
  }

  free ( level );
  free ( order );
  free ( rule );
 
  return;
}
/******************************************************************************/

void test25 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests LEVEL_TO_ORDER_EXPONENTIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  LEVEL_TO_ORDER_EXPONENTIAL uses an exponential rule to\n" );
  printf ( "  determine the order of a rule from its level.\n" );
  printf ( "\n" );
  printf ( 
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n" );
  printf ( "\n" );

  level = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( r = 1; r <= 10; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    level_to_order_exponential ( dim_num, level, rule, order );

    printf ( "  %4d  ", r );
 
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %4d", order[dim] );
    }
    printf ( "\n" );
  }

  free ( level );
  free ( order );
  free ( rule );
 
  return;
}
/******************************************************************************/

void test26 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests LEVEL_TO_ORDER_EXPONENTIAL_SLOW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  LEVEL_TO_ORDER_EXPONENTIAL_SLOW uses a slow exponential rule to\n" );
  printf ( "  determine the order of a rule from its level.\n" );
  printf ( "\n" );
  printf ( "  Since it is really only useful for fully nested rules,\n" );
  printf ( "  we only consider rules 11, 12 and 13.\n" );
  printf ( "\n" );
  printf ( 
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n" );
  printf ( "\n" );

  level = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( r = 11; r <= 13; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    level_to_order_exponential_slow ( dim_num, level, rule, order );

    printf ( "  %4d  ", r );
 
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %4d", order[dim] );
    }
    printf ( "\n" );
  }

  free ( level );
  free ( order );
  free ( rule );
 
  return;
}
/******************************************************************************/

void test27 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests LEVEL_TO_ORDER_LINEAR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num = 11;
  int *level;
  int *order;
  int r;
  int *rule;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  LEVEL_TO_ORDER_LINEAR uses a linear rule to\n" );
  printf ( "  determine the order of a rule from its level.\n" );
  printf ( "\n" );
  printf ( 
    "RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10\n" );
  printf ( "\n" );

  level = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( r = 1; r <= 10; r++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      rule[dim] = r;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level[dim] = dim;
    }
    level_to_order_linear ( dim_num, level, rule, order );

    printf ( "  %4d  ", r );
 
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %4d", order[dim] );
    }
    printf ( "\n" );
  }

  free ( level );
  free ( order );
  free ( rule );
 
  return;
}
/******************************************************************************/

void test28 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST28 tests PATTERSON_LOOKUP against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int level;
  int level_max = 5;
  int n;
  int order;
  int p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  PATTERSON_LOOKUP computes a Gauss-Patterson rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = (3*ORDER+1)/2\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  for ( level = 0; level <= level_max; level++ )
  {
    order = i4_power ( 2, level + 1 ) - 1;

    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    patterson_lookup ( order, x, w );

    if ( order == 1 )
    {
      p = 1;
    }
    else
    {
      p = ( 3 * order + 1 ) / 2;
    }

    for ( n = 0; n <= p + 3; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14f  %14f  %14e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void test29 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST29 tests R8COL_TOL_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0, 
    2.0,  6.0, 10.0, 
    4.0,  8.0, 12.0, 
    1.0,  5.0,  9.0, 
    3.0,  7.0, 11.0, 
    2.0,  6.0,  0.0, 
    2.0,  0.0, 10.1, 
    2.0,  0.1, 10.0, 
    3.0,  4.0, 18.0, 
    1.9,  8.0, 10.0, 
    0.0,  0.0,  0.0, 
    0.0,  6.0, 10.0, 
    2.1,  0.0, 10.0, 
    2.0,  6.0, 10.0, 
    3.0,  7.0, 11.0, 
    2.0,  0.0, 10.0, 
    2.0,  0.0, 10.0, 
    2.0,  6.0, 10.0, 
    1.0,  5.0,  9.0, 
    2.0,  0.0, 10.1, 
    1.0,  5.0,  9.1, 
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  R8COL_TOL_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the tolerably unique columns of an R8COL,\n" );
  printf ( "  and a map from the original R8COL to the (implicit)\n" );
  printf ( "  R8COL of sorted tolerably unique elements.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  tol = 0.25;

  printf ( "\n" );
  printf ( "  Using tolerance = %f\n", tol );

  n_unique = r8col_tol_unique_count ( m, n, a, tol );

  printf ( "\n" );
  printf ( "  Number of tolerably unique columns is %d\n", n_unique );

  au = ( double * ) malloc ( m * n_unique * sizeof ( double ) );
  undx = ( int * ) malloc ( n_unique * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( n ) );

  r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  printf ( "\n" );
  printf ( "  XDNU points to the representative for each item.\n" );
  printf ( "  UNDX selects the representatives.\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  UNDX\n" );
  printf ( "\n" );
  for ( i = 0; i < n_unique; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, xdnu[i], undx[i] );
  }
  for ( i = n_unique; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, xdnu[i] );
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au, 
    "  The tolerably unique R8COL (transposed):" );

  free ( au );
  free ( undx );
  free ( xdnu );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test30 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST30 tests R8VEC_SORT_HEAP_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int i;
  int *indx;
  int n = 20;
  int seed;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  R8VEC_SORT_HEAP_INDEX_A creates an ascending\n" );
  printf ( "  sort index for a R8VEC.\n" );

  seed = 123456789;

  a = r8vec_uniform_01_new ( n, &seed );

  r8vec_print ( n, a, "  The unsorted array:" );

  indx = r8vec_sort_heap_index_a ( n, a );

  i4vec_print ( n, indx, "  The index vector:" );

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[indx[i]];
  }

  r8vec_print ( n, b, "  The sorted array A(INDX(:)):" );

  free ( a );
  free ( b );
  free ( indx );
}
/******************************************************************************/

void test01_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01_NP tests CHEBYSHEV1_COMPUTE_NP against CHEBYSHEV1_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) / sqrt ( 1 - x^2 ) dx.\n" );
  printf ( "\n" );
  printf ( "  CHEBYSHEV1_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev1_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = chebyshev1_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test02_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02_NP tests CHEBYSHEV1_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev1_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test03_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03_NP tests CHEBYSHEV2_COMPUTE_NP against CHEBYSHEV2_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) * sqrt ( 1 - x^2 ) dx.\n" );
  printf ( "\n" );
  printf ( "  CHEBYSHEV2_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev2_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = chebyshev2_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test04_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04_NP tests CHEBYSHEV2_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    chebyshev2_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test05_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05_NP tests CLENSHAW_CURTIS_COMPUTE_NP against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int n_hi;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N up to\n" );
  printf ( "    N = ORDER+1 if ORDER is odd, or\n" );
  printf ( "    N = ORDER   if ORDER is even\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    clenshaw_curtis_compute_np ( order, np, p, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test06_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06_NP tests CLENSHAW_CURTIS_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    clenshaw_curtis_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test07_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07_NP tests FEJER2_COMPUTE_NP against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int n_hi;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N up to\n" );
  printf ( "    N = ORDER+1 if ORDER is odd, or\n" );
  printf ( "    N = ORDER   if ORDER is even\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    fejer2_compute_np ( order, np, p, x, w );

    if ( ( order % 2 ) == 0 )
    {
      n_hi = order + 2;
    }
    else
    {
      n_hi = order + 3;
    }

    for ( n = 0; n <= n_hi; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test08_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08_NP tests FEJER2_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    fejer2_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test09_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09_NP tests GEGENBAUER_COMPUTE_NP against GEGENBAUER_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.\n" );
  printf ( "\n" );
  printf ( "  GEGENBAUER_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gegenbauer_compute_np ( order, np, p, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gegenbauer_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6f  %14.6e\n", 
        order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test10_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10_NP tests GEGENBAUER_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gegenbauer_compute_np ( order, np, p, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }

  free ( p );
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test11_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11_NP tests GEN_HERMITE_COMPUTE_NP against GEN_HERMITE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n" );
  printf ( "\n" );
  printf ( "  GEN_HERMITE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_hermite_compute_np ( order, np, p, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gen_hermite_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
        printf ( "  %8d  %8d  ^14.6f  %14.6e  %14.6e  %14.6e\n", 
          order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test12_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12_NP tests GEN_HERMITE_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_hermite_compute_np ( order, np, p, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }
 
  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test13_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13_NP tests GEN_LAGUERRE_COMPUTE_NP against GEN_LAGUERRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n" );
  printf ( "\n" );
  printf ( "  GEN_LAGUERRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );

      f = ( double * ) malloc ( order * sizeof ( double ) );
      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_laguerre_compute_np ( order, np, p, x, w );
 
      for ( n = 0; n <= 2 * order + 2; n = n + 1 )
      {
        exact = gen_laguerre_integral ( n, alpha );

        if ( n == 0 )
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = 1.0;
          }
        }
        else
        {
          for ( i = 0; i < order; i++ )
          {
            f[i] = pow ( x[i], n );
          }
        }
        estimate = 0.0;
        for ( i = 0; i < order; i++ )
        {
          estimate = estimate + w[i] * f[i];
        }

        error = r8_abs ( exact - estimate );
  
        printf ( "  %8d  %8d  %14.6f  %14.6e  %14.6e  %14.6e\n", 
          order, n, alpha, estimate, exact, error );
      }
      free ( f );
      free ( w );
      free ( x );
    }
  }
 
  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test14_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14_NP tests GEN_LAGUERRE_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 1;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    alpha = alpha_test[test];
    p[0] = alpha;

    for ( order = 1; order <= order_max; order++ )
    {
      printf ( "\n" );
      printf ( "  Order = %d\n", order );
      printf ( "  ALPHA = %f\n", alpha );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );

      gen_laguerre_compute_np ( order, np, p, x, w );
 
      for ( i = 0; i < order; i =i + 1 )
      {
        printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
      }
      free ( w );
      free ( x );
    }
  }
 
  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test15_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15_NP tests HERMITE_COMPUTE_NP against HERMITE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n" );
  printf ( "\n" );
  printf ( "  HERMITE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    hermite_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = hermite_integral ( n );
 
      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test16_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16_NP tests HERMITE_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    hermite_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test17_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17_NP tests JACOBI_COMPUTE_NP against JACOBI_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 2;
  int order;
  int order_max = 10;
  double *p;
  int test1;
  int test2;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n" );
  printf ( "\n" );
  printf ( "  JACOBI_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Alpha           Beta            Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];
    p[0] = alpha;

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];
      p[1] = beta;

      for ( order = 1; order <= order_max; order++ )
      {
        printf ( "\n" );

        f = ( double * ) malloc ( order * sizeof ( double ) );
        w = ( double * ) malloc ( order * sizeof ( double ) );
        x = ( double * ) malloc ( order * sizeof ( double ) );

        jacobi_compute_np ( order, np, p, x, w );

        for ( n = 0; n <= 2 * order + 2; n = n + 1 )
        {
          exact = jacobi_integral ( n, alpha, beta );
 
          if ( n == 0 )
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = 1.0;
            }
          }
          else
          {
            for ( i = 0; i < order; i++ )
            {
              f[i] = pow ( x[i], n );
            }
          }
          estimate = 0.0;
          for ( i = 0; i < order; i++ )
          {
            estimate = estimate + w[i] * f[i];
          }
 
          error = r8_abs ( exact - estimate );
  
          printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6f  %14.6f  %14.6e\n", 
            order, n, alpha, beta, estimate, exact, error );
        }
        free ( f );
        free ( w );
        free ( x );
      }
    }
  }

  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test18_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST18_NP tests JACOBI_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int i;
  int np = 2;
  int order;
  int order_max = 10;
  double *p;
  int test1;
  int test2;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];
    p[0] = alpha;

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];
      p[1] = beta;

      for ( order = 1; order <= order_max; order++ )
      {
        printf ( "\n" );
        printf ( "  Order = %d\n", order );
        printf ( "  ALPHA = %f\n", alpha );
        printf ( "  BETA =  %f\n", beta );

        w = ( double * ) malloc ( order * sizeof ( double ) );
        x = ( double * ) malloc ( order * sizeof ( double ) );

        jacobi_compute_np ( order, np, p, x, w );

        for ( i = 0; i < order; i =i + 1 )
        {
          printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
        }
        free ( w );
        free ( x );
      }
    }
  }

  free ( p );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test19_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19_NP tests LAGUERRE_COMPUTE_NP against LAGUERRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  LAGUERRE_COMPUTE_NP computes a Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n" );
  printf ( "\n" );
  printf ( "  LAGUERRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    laguerre_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = laguerre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6e  %14.6e  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test20_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20_NP tests LAGUERRE_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    laguerre_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
/******************************************************************************/

void test21_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST21_NP tests LEGENDRE_COMPUTE_NP against LEGENDRE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int n;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x)  dx.\n" );
  printf ( "\n" );
  printf ( "  LEGENDRE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "  A rule of order ORDER should be exact for monomials X^N\n" );
  printf ( "  up to N = 2*ORDER-1\n" );
  printf ( "\n" );
  printf ( "  In the following table, for each order, the LAST THREE estimates\n" );
  printf ( "  are made on monomials that exceed the exactness limit for the rule.\n" );
  printf ( "\n" );
  printf ( "     Order         N       Estimate       Exact            Error\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );

    f = ( double * ) malloc ( order * sizeof ( double ) );
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    legendre_compute_np ( order, np, p, x, w );

    for ( n = 0; n <= 2 * order + 2; n = n + 1 )
    {
      exact = legendre_integral ( n );

      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f[i] = pow ( x[i], n );
        }
      }
      estimate = 0.0;
      for ( i = 0; i < order; i++ )
      {
        estimate = estimate + w[i] * f[i];
      }

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14.6f  %14.6f  %14.6e\n", 
        order, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }

  free ( p );

  return;
}
/******************************************************************************/

void test22_np ( void )

/******************************************************************************/
/*
  Purpose:

    TEST22_NP tests LEGENDRE_COMPUTE_NP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt
*/
{
  int i;
  int np = 0;
  int order;
  int order_max = 10;
  double *p;
  int test;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx.\n" );

  p = ( double * ) malloc ( np * sizeof ( double ) );

  for ( order = 1; order <= order_max; order++ )
  {
    printf ( "\n" );
    printf ( "  Order = %d\n", order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    legendre_compute_np ( order, np, p, x, w );
 
    for ( i = 0; i < order; i =i + 1 )
    {
      printf ( "  %8d  %24.16f  %24.16f\n", i, x[i], w[i] );
    }
    free ( w );
    free ( x );
  }
 
  free ( p );

  return;
}
