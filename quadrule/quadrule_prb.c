# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "quadrule.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test0725 ( void );
void test075 ( void );
void test076 ( void );
void test078 ( void );
void test079 ( int order, double alpha );
void test08 ( void );
void test085 ( void );
void test087 ( void );
void test09 ( void );
void test095 ( void );
void test10 ( void );
void test105 ( void );
void test108 ( int order );
void test11 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
void test16 ( void );
void test165 ( int order, double alpha );
void test17 ( void );
void test18 ( int n );
void test185 ( void );
void test19 ( void );
void test20 ( void );
void test21 ( void );
void test22 ( void );
void test23 ( void );
void test24 ( void );
void test25 ( void );
void test26 ( void );
void test27 ( void );
void test28 ( void );
void test29 ( void );
void test30 ( void );
void test31 ( void );
void test32 ( void );
void test33 ( void );
void test34 ( void );
void test345 ( void );
void test35 ( void );
void test36 ( void );
void test37 ( void );
void test38 ( void );
void test39 ( void );
void test40 ( void );
void test401 ( void );
void test402 ( void );
void test403 ( void );
void test404 ( void );
void test41 ( void );
double f1sd1 ( double x );
char *function_name ( int function_index );
void function_set ( char *action, int *i );
double function_value ( double x );
double fx1sd1 ( double x );
double fx2sd1 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QUADRULE_PRB.

  Discussion:

    QUADRULE_PRB calls a set of tests for the QUADRULE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  int n;

  timestamp ( );

  printf ( "\n" );
  printf ( "QUADRULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the QUADRULE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test0725 ( );
  test075 ( );
  test076 ( );
  test078 ( );

  n = 5;
  alpha = 0.5;
  test079 ( n, alpha );

  n = 10;
  alpha = - 0.5;
  test079 ( n, alpha );

  test08 ( );
  test085 ( );
  test087 ( );
  test09 ( );
  test095 ( );

  test10 ( );
  test105 ( );
  n = 10;
  test108 ( n );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );

  n = 11;
  alpha = 0.0;
  test165 ( n, alpha );

  n = 11;
  alpha = 0.5;
  test165 ( n, alpha );

  n = 11;
  alpha = 2.0;
  test165 ( n, alpha );

  test17 ( );
/*
  Compare computed and lookup versions of Gauss-Legendre rules.
*/
  n = 31;
  test18 ( n );
  n = 64;
  test18 ( n );
  n = 129;
  test18 ( n );
  n = 255;
  test18 ( n );

  test185 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test401 ( );
  test402 ( );
  test403 ( );
  test404 ( );
  test41 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUADRULE_PRB\n" );
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

    TEST01 tests BASHFORTH_SET and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int n;
  int n_max = 10;
  double *result;
  double *w;
  double *x;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  BASHFORTH_SET sets up an Adams-Bashforth rule;\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [0,1].\n" );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "  Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( n = 1; n <= n_max; n++ )
    {
      x = ( double * ) malloc ( n * sizeof ( double ) );
      w = ( double * ) malloc ( n * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        bashforth_set ( n, x, w );
 
        result[i] = summer ( function_value, n, x, w );
 
      }
      printf ( "  %2d  ", n );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8e", result[i] );
      }
      printf ( "\n" );

      free ( x );
      free ( w );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BDFC_SET and BDF_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  BDFC_SET sets up a Backward Difference Corrector rule;\n" );
  printf ( "  BDF_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [0,1].\n" );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        bdfc_set ( order, xtab, weight );
 
        result[i] = bdf_sum ( function_value, order, xtab, weight );
 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests BDFP_SET and BDF_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  BDFP_SET sets up a Backward Difference Predictor rule;\n" );
  printf ( "  BDF_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [0,1].\n" );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        bdfp_set ( order, xtab, weight );
 
        result[i] = bdf_sum ( function_value, order, xtab, weight );
 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CHEB_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CHEB_SET sets up a Chebyshev rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        cheb_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CHEBYSHEV1_COMPUTE and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 6;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -1.0;
  b =  1.0;
  nsub = 1;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CHEBYSHEV1_COMPUTE sets up a Gauss-Chebyshev type 1 rule,\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is 1 / sqrt ( 1 - X**2 )\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        chebyshev1_compute ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.6f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests CHEBYSHEV2_COMPUTE and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 4;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -1.0;
  b =  1.0;
  nsub = 1;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  CHEBYSHEV2_COMPUTE sets up a Gauss-Chebyshev type 2 rule,\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is 1 / sqrt ( 1 - X**2 )\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        chebyshev2_compute ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.6f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests CHEBYSHEV3_COMPUTE and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 6;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -1.0;
  b =  1.0;
  nsub = 1;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  CHEBYSHEV3_COMPUTE sets up a Gauss-Chebyshev type 3 rule,\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is 1 / sqrt ( 1 - X**2 )\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        chebyshev3_compute ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.6f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );

  return;
}
/******************************************************************************/

void test0725 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0725 tests CLENSHAW_CURTIS_COMPUTE

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2006

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
  printf ( "TEST0725\n" );
  printf ( "  CLENSHAW_CURTIS_COMPUTE computes\n" );
  printf ( "  a Clenshaw-Curtis quadrature rule over [-1,1]\n" );
  printf ( "  of given order.\n" );

  printf ( "\n" );
  printf ( "    Order  W             X\n" );
  printf ( "\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( order * sizeof ( double ) );

    clenshaw_curtis_compute ( order, x, w );

    printf ( "\n" );
    printf ( "  %8d\n", order );

    for ( i = 0; i < order; i++ )
    {
      printf ( "            %14f  %14f\n", w[i], x[i] );
    }
    free ( w );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test075 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST075 tests CLENSHAW_CURTIS_SET and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 16;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST075\n" );
  printf ( "  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [-1,1].\n" );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        clenshaw_curtis_set ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight );
 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test076 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST076 compares FEJER1_COMPUTE and FEJER1_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 March 2007

  Author:

    John Burkardt
*/
{
  int i;
  int order;
  int order_max = 10;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST076\n" );
  printf ( "  FEJER1_COMPUTE computes a Fejer type 1 quadrature rule;\n" );
  printf ( "  FEJER1_SET sets a Fejer type 1 quadrature rule;\n" );
  printf ( "\n" );
  printf ( "  Compare:\n" );
  printf ( "    (W1,X1) from FEJER1_SET,\n" );
  printf ( "    (W2,X2) from FEJER1_COMPUTE.\n" );
  printf ( "\n" );
  printf ( "     Order        W1              W2              X1             X2\n" );
  printf ( "\n" );

  for ( order = 1; order <= order_max; order++ )
  {
    w1 = ( double * ) malloc ( order * sizeof ( double ) );
    x1 = ( double * ) malloc ( order * sizeof ( double ) );

    fejer1_set ( order, x1, w1 );

    w2 = ( double * ) malloc ( order * sizeof ( double ) );
    x2 = ( double * ) malloc ( order * sizeof ( double ) );

    fejer1_compute ( order, x2, w2 );

    printf ( "\n" );
    printf ( "  %8d\n", order );

    for ( i = 0; i < order; i++ )
    {
      printf ( "          %14f  %14f  %14f  %14f\n", w1[i], w2[i], x1[i], x2[i]  );
    }
    free ( w1 );
    free ( w2 );
    free ( x1 );
    free ( x2 );
  }

  return;
}
/******************************************************************************/

void test078 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST078 compares FEJER2_COMPUTE and FEJER2_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 March 2007

  Author:

    John Burkardt
*/
{
# define ORDER_MAX 10

  int i;
  int order;
  double w1[ORDER_MAX];
  double w2[ORDER_MAX];
  double x1[ORDER_MAX];
  double x2[ORDER_MAX];

  printf ( "\n" );
  printf ( "TEST078\n" );
  printf ( "  FEJER2_COMPUTE computes a Fejer type 2 quadrature rule;\n" );
  printf ( "  FEJER2_SET sets a Fejer type 2 quadrature rule;\n" );
  printf ( "\n" );
  printf ( "  Compare:\n" );
  printf ( "    (W1,X1) from FEJER2_SET,\n" );
  printf ( "    (W2,X2) from FEJER2_COMPUTE.\n" );
  printf ( "\n" );
  printf ( "     Order        W1              W2              X1             X2\n" );
  printf ( "\n" );

  for ( order = 1; order <= ORDER_MAX; order++ )
  {
    fejer2_set ( order, x1, w1 );
    fejer2_compute ( order, x2, w2 );

    printf ( "\n" );
    printf ( "  %8d\n", order );

    for ( i = 0; i < order; i++ )
    {
      printf ( "          %14f  %14f  %14f  %14f\n", w1[i], w2[i], x1[i], x2[i]  );
    }
  }

  return;
# undef ORDER_MAX
}
/******************************************************************************/

void test079 ( int order, double alpha )

/******************************************************************************/
/*
  Purpose:

    TEST079 tests GEGENBAUER_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the parameter.
*/
{
  int i;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST079\n" );
  printf ( "  GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer rule;\n" );
  printf ( "\n" );
  printf ( "  The printed output of this routine can be inserted into\n" );
  printf ( "  a C++ program.\n" );

  w = ( double * ) malloc ( order * sizeof ( double ) );
  x = ( double * ) malloc ( order * sizeof ( double ) );

  gegenbauer_compute ( order, alpha, x, w );

  printf ( "\n" );
  printf ( "  Abscissas X and weights W for a Gauss Gegenbauer rule\n" );
  printf ( "  of ORDER   = %d\n", order );
  printf ( "  with ALPHA = %f\n", alpha );
  printf ( "\n" );

  for ( i = 0; i < order; i++ )
  {
    printf ( "    x[%2d] = %24.16e\n", i, x[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    x[%2d] = %24.16e\n", i, w[i] );
  }

  free ( w );
  free ( x );

  return;
# undef ORDER
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests HERMITE_SS_COMPUTE and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  HERMITE_SS_COMPUTE computes a Gauss-Hermite rule;\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is ( -Infinity, +Infinity ).\n" );
  printf ( "  The weight function is exp ( - X**2 )\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        hermite_ss_compute ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test085 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST085 tests HERMITE_SS_COMPUTE against HERMITE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2007

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double *f_vec;
  int i;
  int n;
  int order;
  int order_max = 10;
  double *weight;
  double *xtab;

  printf ( "\n" );
  printf ( "TEST085\n" );
  printf ( "  HERMITE_SS_COMPUTE computes a Gauss-Hermite rule\n" );
  printf ( "  which is appropriate for integrands of the form\n" );
  printf ( "    f(x) * exp(-x**2) from -infinity to infinity.\n" );
  printf ( "\n" );
  printf ( "  HERMITE_INTEGRAL determines the exact value of\n" );
  printf ( "  this integal when f(x) = x^n.\n" );
  printf ( "\n" );
  printf ( "         N     Order       Estimate       Exact            Error\n" );

  for ( n = 0; n <= 10; n = n + 2 )
  {
    exact = hermite_integral ( n );

    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      f_vec = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );
      xtab = ( double * ) malloc ( order * sizeof ( double ) );

      hermite_ss_compute ( order, xtab, weight );
 
      if ( n == 0 )
      {
        for ( i = 0; i < order; i++ )
        {
          f_vec[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < order; i++ )
        {
          f_vec[i] = pow ( xtab[i], n );
        }
      }
      estimate = r8vec_dot_product ( order, weight, f_vec );

      error = r8_abs ( exact - estimate );
  
      printf ( "  %8d  %8d  %14f  %14f  %14e\n", n, order, estimate, exact, error );

      free ( f_vec );
      free ( weight );
      free ( xtab );
    }
  }
 
  return;
}
/******************************************************************************/

void test087 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST087 tests HERMITE_SS_COMPUTE.

  Discussion:

    I used this test to generate tabular values of weights and
    abscissas for Gauss-Hermite quadrature.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 July 2007

  Author:

    John Burkardt
*/
{
  int i;
  int order = 31;
  double *weight;
  double *xtab;

  printf ( "\n" );
  printf ( "TEST087\n" );
  printf ( "  HERMITE_SS_COMPUTE computes a Gauss-Hermite rule;\n" );
  printf ( "\n" );
  printf ( "  Compute the data for ORDER = %d\n", order );

  weight = ( double * ) malloc ( order * sizeof ( double ) );
  xtab = ( double * ) malloc ( order * sizeof ( double ) );

  hermite_ss_compute ( order, xtab, weight );
 
  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    xtab[%2d] = %24.16f\n", i, xtab[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    weight[%2d] = %24.16f\n", i, weight[i] );
  }
  free ( weight );
  free ( xtab );

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests HERMITE_SET and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  HERMITE_SET sets up a Gauss-Hermite rule;\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is ( -Infinity, +Infinity ).\n" );
  printf ( "  The weight function is exp ( - X**2 )\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        hermite_set ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test095 ( )

/******************************************************************************/
/*
  Purpose:

    TEST095 tests HERMITE_GENZ_KEISTER_SET and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 June 2010

  Author:

    John Burkardt
*/
{
# define L_MAX 8

  double error;
  double estimate;
  double exact;
  double *f;
  int i;
  int l;
  int m;
  int n;
  int n_list[L_MAX+1] = { 1, 3, 7, 9, 17, 19, 31, 33, 35 };
  int p;
  int p_list[L_MAX+1] = { 1, 5, 7, 15, 17, 29, 31, 33, 51 };
  double *w;
  double *x;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST095\n" );
  fprintf ( stdout, "  HERMITE_GENZ_KEISTER_SET sets up a nested rule\n" );
  fprintf ( stdout, "  for the Hermite integration problem.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  The integration interval is ( -oo, +oo ).\n" );
  fprintf ( stdout, "  The weight function is exp ( - x * x )\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  HERMITE_INTEGRAL determines the exact value of\n" );
  fprintf ( stdout, "  the integal when f(x) = x^m.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "         M         N       Estimate       Exact            Error\n" );

  for ( l = 0; l <= L_MAX; l++ )
  {
    fprintf ( stdout, "\n" );
    n = n_list[l];

    f = ( double * ) malloc ( n * sizeof ( double ) );
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    p = p_list[l];

    hermite_genz_keister_set ( n, x, w );

    for ( m = 0; m <= i4_min ( p + 2, 20 ); m = m + 2 )
    {
      exact = hermite_integral ( m );

      if ( m == 0 )
      {
        for ( i = 0; i < n; i++ )
        {
          f[i] = 1.0;
        }
      }
      else
      {
        for ( i = 0; i < n; i++ )
        {
          f[i] = pow ( x[i], m );
        }
      }

      estimate = r8vec_dot_product ( n, w, f );

      error = r8_abs ( exact - estimate );
  
      fprintf ( stdout, "  %8d  %8d  %14.6e  %14.6e  %14.6e\n",
        m, n, estimate, exact, error );
    }
    free ( f );
    free ( w );
    free ( x );
  }
  return;
# undef L_MAX
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests JACOBI_COMPUTE and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2007

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int k;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  a = -1.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  JACOBI_COMPUTE computes a Gauss-Jacobi rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( k = 1; k <= 2; k++ )
  {
    result = ( double * ) malloc ( function_num * sizeof ( double ) );

    if ( k == 1 )
    {
      alpha = 0.0;
      beta = 0.0;
    }
    else
    {
      alpha = 1.0;
      beta = 0.0;
    }
    printf ( "\n" );
    printf ( "  ALPHA = %f\n", alpha );
    printf ( "  BETA =  %f\n", beta );

    for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
    {
      ihi = i4_min ( ilo + 4, function_num - 1 );

      printf ( "\n" );
      printf ( "Order  " );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "%10s    ", function_name ( i ) );
      }
      printf ( "\n" );
      printf ( "\n" );

      for ( order = 1; order <= order_max; order++ )
      {
        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        for ( i = ilo; i <= ihi; i++ )
        {
          function_set ( "SET", &i );

          jacobi_compute ( order, alpha, beta, xtab, weight );
 
          result[i] = sum_sub ( function_value, a, b, nsub, order, 
            xlo, xhi, xtab, weight ); 
        }
        printf ( "%2d  ", order );
        for ( i = ilo; i <= ihi; i++ )
        {
          printf ( "  %12.8f", result[i] );
        }
        printf ( "\n" );

        free ( xtab );
        free ( weight );
      }
    }
    free ( result );
  }
 
  return;
}
/******************************************************************************/

void test105 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST105 tests JACOBI_COMPUTE.

  Discussion:

    Compare with tabular values on page 178 of Stroud and Secrest.

     In particular,

             X              W

     1  -0.9833999115   0.4615276287E-03
     2  -0.9447138932   0.2732249104E-02
     3  -0.8849310847   0.8045830455E-02
    ..  .............   ................
    19   0.9656375637   0.7613987785E-01
    20   0.9934477866   0.3348337670E-01

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2007

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  int i;
  int order = 20;
  double *weight;
  double *xtab;

  a = -1.0;
  b =  1.0;

  printf ( "\n" );
  printf ( "TEST105\n" );
  printf ( "  JACOBI_COMPUTE computes a Gauss-Jacobi rule;\n" );
  printf ( "  Here, we simply compute a single rule and\n" );
  printf ( "  print it, to check for accuracy.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );

  alpha = 0.0;
  beta = 1.0;

  printf ( "\n" );
  printf ( "  N = %d\n", order );
  printf ( "  ALPHA = %f\n", alpha );
  printf ( "  BETA =  %f\n", beta );

  xtab = ( double * ) malloc ( order * sizeof ( double ) );
  weight = ( double * ) malloc ( order * sizeof ( double ) );
  
  jacobi_compute ( order, alpha, beta, xtab, weight );
 
  printf ( "\n" );
  printf ( "     I        X(I)            W(I)\n" );
  printf ( "\n" );

  for ( i = 0; i < order; i++ )
  {
    printf ( "  %4d  %14f  %14f\n", i+1, xtab[i], weight[i] );
  }

  free ( xtab );
  free ( weight );
 
  return;
}
/******************************************************************************/

void test108 ( int order )

/******************************************************************************/
/*
  Purpose:

    TEST108 tests JACOBI_COMPUTE.

  Discussion:

    I used this test to generate tabular values of weights and
    abscissas for Gauss-Jacobi quadrature.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2008

  Author:

    John Burkardt
*/
{
  double alpha = 0.5;
  double beta  = 2.0;
  int i;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST108\n" );
  printf ( "  JACOBI_COMPUTE computes a Gauss-Jacobi rule;\n" );
  printf ( "\n" );
  printf ( "  The printed output of this test can be inserted into\n" );
  printf ( "  a C++ program.\n" );

  w = ( double * ) malloc ( order * sizeof ( double ) );
  x = ( double * ) malloc ( order * sizeof ( double ) );

  jacobi_compute ( order, alpha, beta, x, w );

  printf ( "\n" );
  printf ( "  Abscissas X and weights W for a Gauss Jacobi rule\n" );
  printf ( "  of ORDER   = %d\n", order );
  printf ( "  with ALPHA = %f\n", alpha );
  printf ( "  and  BETA  = %f\n", beta );
  printf ( "\n" );

  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    x[%2d] = %24.16f\n", i, x[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    w[[%2d] = %24.16f\n", i, w[i] );
  }

  free ( w );
  free ( x );

  return;
# undef ORDER
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests KRONROD_SET, LEGENDRE_SET and SUMMER_GK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt
*/
{
# define ORDERG 10
# define ORDERK 2 * ORDERG + 1

  double resultg;
  double resultk;
  double weightg[ORDERG];
  double weightk[ORDERK];
  double xtabg[ORDERG];
  double xtabk[ORDERK];

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  KRONROD_SET sets up a Kronrod rule;\n" );
  printf ( "  LEGENDRE_SET sets up Gauss-Legendre rule;\n" );
  printf ( "  SUMMER_GK carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [-1, 1].\n" );
  printf ( "  Integrand is X**2 / SQRT ( 1.1 - X**2 ).\n" );
  printf ( "\n" );

  legendre_set ( ORDERG, xtabg, weightg );

  kronrod_set ( ORDERK, xtabk, weightk );

  summer_gk ( fx2sd1, ORDERG, weightg, &resultg, 
    ORDERK, xtabk, weightk, &resultk );

  printf ( "  %2d  %16.10f\n", ORDERG, resultg );
  printf ( "  %2d  %16.10f\n", ORDERK, resultk );
  printf ( "      %16.10e\n", resultg - resultk );

  return;
# undef ORDERG
# undef ORDERK
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests KRONROD_SET, LEGENDRE_SET and SUM_SUB_GK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
# define ORDERG 7
# define ORDERK 2 * ORDERG + 1

  double a;
  double b;
  double error;
  int nsub;
  double resultg;
  double resultk;
  double weightg[ORDERG];
  double weightk[ORDERK];
  double xtabg[ORDERG];
  double xtabk[ORDERK];

  a = -1.0;
  b =   1.0;
  nsub = 5;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  KRONROD_SET sets up a Kronrod rule;\n" );
  printf ( "  LEGENDRE_SET sets up Gauss-Legendre rule;\n" );
  printf ( "  SUM_SUB_GK carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  Number of subintervals is %d\n", nsub );
  printf ( "  Integrand is X**2 / SQRT ( 1.1 - X**2 ).\n" );
  printf ( "\n" );

  legendre_set ( ORDERG, xtabg, weightg );

  kronrod_set ( ORDERK, xtabk, weightk );

  sum_sub_gk ( fx2sd1, a, b, nsub, ORDERG, weightg, &resultg, 
    ORDERK, xtabk, weightk, &resultk, &error );

  printf ( "  %2d  %16.10f\n", ORDERG, resultg );
  printf ( "  %2d  %16.10f\n", ORDERK, resultk );
  printf ( "      %16.10e\n", error );

  return;
# undef ORDERG
# undef ORDERK
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests LAGUERRE_COMPUTE and LAGUERRE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 1.0;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  LAGUERRE_COMPUTE computes a Gauss-Laguerre rule;\n" );
  printf ( "  LAGUERRE_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f, +oo).\n", a );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is EXP ( - X ).\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_compute ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests LAGUERRE_COMPUTE and LAGUERRE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  LAGUERRE_COMPUTE sets up a Gauss-Laguerre rule;\n" );
  printf ( "  LAGUERRE_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f, +oo).\n", a );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is EXP ( - X ).\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_compute ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests GEN_LAGUERRE_COMPUTE and LAGUERRE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  alpha = 1.0;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule;\n" );
  printf ( "  LAGUERRE_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f, +oo).\n", a );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is EXP ( - X ) * X^%f.\n", alpha );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        gen_laguerre_compute ( order, alpha, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests GEN_LAGUERRE_COMPUTE and LAGUERRE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  alpha = 2.0;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule;\n" );
  printf ( "  LAGUERRE_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f, +Infinity).\n", a );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is EXP ( - X ) * X^%f.\n", alpha );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        gen_laguerre_compute ( order, alpha, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test165 ( int order, double alpha )

/******************************************************************************/
/*
  Purpose:

    TEST165 tests GEN_LAGUERRE_COMPUTE.

  Discussion:

    I used this test to generate tabular values of weights and
    abscissas for generalized Gauss-Laguerre quadrature.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2007

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the parameter.
*/
{
  int i;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST165\n" );
  printf ( "  LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule;\n" );
  printf ( "\n" );
  printf ( "  The printed output of this routine can be inserted into\n" );
  printf ( "  a C++ program.\n" );

  w = ( double * ) malloc ( order * sizeof ( double ) );
  x = ( double * ) malloc ( order * sizeof ( double ) );

  gen_laguerre_compute ( order, alpha, x, w );

  printf ( "\n" );
  printf ( "  Abscissas X and weights W for a Gauss Laguerre rule\n" );
  printf ( "  of ORDER   = %d\n", order );
  printf ( "  with ALPHA = %f\n", alpha );
  printf ( "\n" );

  for ( i = 0; i < order; i++ )
  {
    printf ( "    x[%2d] = %24.16f\n", i, x[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < order; i++ )
  {
    printf ( "    w[%2d] = %24.16f\n", i, w[i] );
  }

  free ( w );
  free ( x );

  return;
# undef ORDER
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests LAGUERRE_SET and LAGUERRE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  1.0;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  LAGUERRE_SET sets up a Gauss-Laguerre rule;\n" );
  printf ( "  LAGUERRE_SUM carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f, +Infinity).\n", a );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "  The weight function is EXP ( - X ).\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        laguerre_set ( order, xtab, weight );
 
        result[i] = laguerre_sum ( function_value, a, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test18 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST18 compares LEGENDRE_COMPUTE_DR and LEGENDRE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2010

  Author:

    John Burkardt
*/
{
  int i;
  int iwdifmax;
  int ixdifmax;
  double wdifmax;
  double *w1;
  double *w2;
  double xdifmax;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule;\n" );
  printf ( "  LEGENDRE_SET looks up the same data.\n" );
  printf ( "\n" );
  printf ( "  Compare the data for N = %d\n", n );

  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  legendre_compute_dr ( n, x1, w1 );

  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  w2 = ( double * ) malloc ( n * sizeof ( double ) );

  legendre_set ( n, x2, w2 );

  xdifmax = 0.0;
  ixdifmax = -1;

  wdifmax = 0.0;
  iwdifmax = -1;

  for ( i = 0; i < n; i++ )
  {
    if ( xdifmax < r8_abs ( x1[i] - x2[i] ) )
    {
      xdifmax = r8_abs ( x1[i] - x2[i] );
      ixdifmax = i;
    }

    if ( wdifmax < r8_abs ( w1[i] - w2[i] ) )
    {
      wdifmax = r8_abs ( w1[i] - w2[i] );
      iwdifmax = i;
    }

  }

  if ( -1 < ixdifmax )
  {
    printf ( "\n" );
    printf ( "  Maximum abscissa difference is %e\n", xdifmax );
    printf ( "  for index I = %d\n", ixdifmax );
    printf ( "  Computed: %f\n", x1[ixdifmax] );
    printf ( "  Stored:   %f\n", x2[ixdifmax] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  The computed and stored abscissas are identical.\n" );
  }

  if ( -1 < iwdifmax )
  {
    printf ( "\n" );
    printf ( "  Maximum weight difference is   %e\n", wdifmax );
    printf ( "  for index I = %d\n", iwdifmax );
    printf ( "  Computed: %f\n", w1[iwdifmax] );
    printf ( "  Stored:   %f\n", w2[iwdifmax] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  The computed and stored weights are identical.\n" );
  }

  free ( w1 );
  free ( w2 );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test185 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST185 tests LEGENDRE_COMPUTE_DR.

  Discussion:

    I used this test to generate tabular values of weights and
    abscissas for Gauss-Legendre quadrature.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
# define ORDER 31

  int i;
  double weight[ORDER];
  double xtab[ORDER];

  printf ( "\n" );
  printf ( "TEST185\n" );
  printf ( "  LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule;\n" );
  printf ( "\n" );
  printf ( "  Compute the data for ORDER = %d\n", ORDER );

  legendre_compute_dr ( ORDER, xtab, weight );
 
  printf ( "\n" );
  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "    xtab[%2d] = %24.16f\n", i, xtab[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "    weight[%2d] = %24.16f\n", i, weight[i] );
  }

  return;
# undef ORDER
}
/******************************************************************************/

void test19 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests LEGENDRE_COMPUTE_DR and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 May 2006

  Author:

    John Burkardt
*/
{
# define ORDER 2

  double a;
  double b;
  int iexp;
  int nsub;
  double result;
  double weight[ORDER];
  double xhi;
  double xlo;
  double xtab[ORDER];

  a = 0.0;
  b = 1.0;

  xlo = -1.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule;\n" );
  printf ( "  SUM_SUB carries it out over subintervals.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  Here, we use a fixed order ORDER = %d\n", ORDER );
  printf ( "  and use more and more subintervals.\n" );
  printf ( "\n" );
  printf ( "  NSUB     Integral\n" );
  printf ( "\n" );
 
  legendre_compute_dr ( ORDER, xtab, weight );
 
  for ( iexp = 0; iexp <= 9; iexp++ )
  {
    nsub = i4_power ( 2, iexp );

    result = sum_sub ( fx2sd1, a, b, nsub, ORDER, xlo, xhi, 
      xtab, weight );

    printf ( "  %4d  %16.10f\n", nsub, result );
  } 
  return;
# undef ORDER
}
/******************************************************************************/

void test20 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests LEGENDRE_COMPUTE_DR and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  LEGENDRE_COMPUTE_DR sets up a Gauss-Legendre rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_compute_dr ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test21 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests LEGENDRE_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  LEGENDRE_SET sets up a Gauss-Legendre rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test22 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests LEGENDRE_SET, LEGENDRE_SET_X0_01 and RULE_ADJUST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 May 2006

  Author:

    John Burkardt
*/
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b = +1.0;
  c =  0.0;
  d =  1.0;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  LEGENDRE_SET sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating F(X) over [-1,1];\n" );
  printf ( "  RULE_ADJUST adjusts a rule for a new interval.\n" );
  printf ( "  LEGENDRE_SET_X0_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating F(X) over [0,1];\n" );
  printf ( "\n" );
  printf ( "  We will use LEGENDRE_SET to get a rule for [-1,1],\n" );
  printf ( "  adjust it to [0,1] using RULE_ADJUST,\n" );
  printf ( "  and compare the results of LEGENDRE_SET_X0_01.\n" );
  printf ( "\n" );

  legendre_set ( ORDER, xtab1, weight1 );

  r8vec_copy ( ORDER, xtab1, xtab2 );
  r8vec_copy ( ORDER, weight1, weight2 );

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );

  legendre_set_x0_01 ( ORDER, xtab3, weight3 );

  printf ( "\n" );
  printf ( "  Abscissas:\n" );
  printf ( "\n" );
  printf ( "          Original          Adjusted            Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n",  i, xtab1[i], xtab2[i], xtab3[i] );
  }

  printf ( "\n" );
  printf ( "  Weights:\n" );
  printf ( "\n" );
  printf ( "          Original          Adjusted            Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n", i, weight1[i], weight2[i], weight3[i] );
  }

  return;
# undef ORDER
}
/******************************************************************************/

void test23 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests LEGENDRE_SET_COS and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double pi = 3.141592653589793;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -0.5 * pi;
  b = +0.5 * pi;

  nsub = 1;

  xlo = -0.5 * pi;
  xhi = +0.5 * pi;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  LEGENDRE_SET_COS sets up a Gauss-Legendre rule\n" );
  printf ( "    over [-PI/2,PI/2] with weight function COS(X);\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( iexp = 0; iexp <= 4; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_cos ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test24 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests LEGENDRE_SET_SQRTX_01 and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  LEGENDRE_SET_SQRTX_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    over [0,1] with weight function SQRT(X);\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( iexp = 0; iexp <= 3; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_sqrtx_01 ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test25 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests LEGENDRE_SET_SQRTX2_01 and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int order;
  int order_max = 20;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  LEGENDRE_SET_SQRTX2_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    over [0,1] with weight function 1/SQRT(X);\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( iexp = 0; iexp <= 3; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_sqrtx2_01 ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test26 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests LEGENDRE_SET_COS2 and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int iexp;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  double pi = 3.141592653589793;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -0.5 * pi;
  b = +0.5 * pi;

  nsub = 1;

  xlo = -0.5 * pi;
  xhi = +0.5 * pi;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  LEGENDRE_SET_COS2 sets up a Gauss-Legendre rule\n" );
  printf ( "    over [0,PI/2] with weight function COS(X);\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( iexp = 1; iexp <= 4; iexp++ )
    {
      order = i4_power ( 2, iexp );

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_cos2 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test27 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests LEGENDRE_SET_LOG and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt
*/
{
# define TEST_NUM 9

  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 20;
  int order_test[TEST_NUM] = { 1, 2, 3, 4, 5, 6, 7, 8, 16 };
  double *result;
  int test;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  LEGENDRE_SET_LOG sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating -LOG(X) * F(X) over [0,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( test = 0; test < TEST_NUM; test++ )
    {
      order = order_test[test];
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_log ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test28 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST28 tests LEGENDRE_SET_X0_01 and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  LEGENDRE_SET_X0_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating F(X) over [0,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f].\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x0_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test29 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST29 tests LEGENDRE_SET_X1 and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  LEGENDRE_SET_X1 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating ( 1 + X ) * F(X) over [-1,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x1 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test30 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST30 tests LEGENDRE_SET_X1, LEGENDRE_SET_X1_01 and RULE_ADJUST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b =  1.0;
  c =  0.0;
  d =  1.0;

  printf ( "\n" );
  printf ( "TEST30:\n" );
  printf ( "  LEGENDRE_SET_X1 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating ( 1 + X ) * F(X) over [-1,1];\n" );
  printf ( "  RULE_ADJUST adjusts a rule for a new interval.\n" );
  printf ( "  LEGENDRE_SET_X1_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating X * F(X) over [0,1];\n" );
  printf ( "\n" );
  printf ( "  We will use LEGENDRE_SET_X1 to get a rule for [-1,1],\n" );
  printf ( "  adjust it to [0,1] using RULE_ADJUST,\n" );
  printf ( "  make further adjustments because the weight function\n" );
  printf ( "  is not 1,\n" );
  printf ( "  and compare the results of LEGENDRE_SET_X1_01.\n" );
  printf ( "\n" );

  legendre_set_x1 ( ORDER, xtab1, weight1 );

  for ( i = 0; i < ORDER; i++ )
  {
    xtab2[i] = xtab1[i];
  }
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight1[i];
  }

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );
/*
  Because the weight function W(X) is not 1, we need to do more
  adjustments to the weight vector.
*/
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight2[i] / 2.0;
  }
  legendre_set_x1_01 ( ORDER, xtab3, weight3 );

  printf ( "\n" );
  printf ( "  Abscissas:\n" );
  printf ( "\n" );
  printf ( "  Original  Adjusted Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n", i, xtab1[i], xtab2[i], xtab3[i] );
  }

  printf ( "\n" );
  printf ( "  Weights:\n" );
  printf ( "\n" );
  printf ( "  Original  Adjusted Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n", i, weight1[i], weight2[i], weight3[i] );
  }

  return;
# undef ORDER
}
/******************************************************************************/

void test31 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests LEGENDRE_SET_X1_01 and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  LEGENDRE_SET_X1_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating X * F(X) over [0,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x1_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test32 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST32 tests LEGENDRE_SET_Xs and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = -1.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  LEGENDRE_SET_X2 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating (1 + X)**2 * F(X) over [-1,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d.\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x2 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test33 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST33 tests LEGENDRE_SET_X2, LEGENDRE_SET_X2_01 and RULE_ADJUST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
# define ORDER 5

  double a;
  double b;
  double c;
  double d;
  int i;
  double weight1[ORDER];
  double weight2[ORDER];
  double weight3[ORDER];
  double xtab1[ORDER];
  double xtab2[ORDER];
  double xtab3[ORDER];

  a = -1.0;
  b =  1.0;
  c =  0.0;
  d =  1.0;

  printf ( "\n" );
  printf ( "TEST33:\n" );
  printf ( "  LEGENDRE_SET_X2 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating ( 1 + X )^2 * F(X) over [-1,1];\n" );
  printf ( "  RULE_ADJUST adjusts a rule for a new interval.\n" );
  printf ( "  LEGENDRE_SET_X2_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating X^2 * F(X) over [0,1];\n" );
  printf ( "\n" );
  printf ( "  We will use LEGENDRE_SET_X2 to get a rule for [-1,1],\n" );
  printf ( "  adjust it to [0,1] using RULE_ADJUST,\n" );
  printf ( "  make further adjustments because the weight function\n" );
  printf ( "  is not 1,\n" );
  printf ( "  and compare the results of LEGENDRE_SET_X2_01.\n" );
  printf ( "\n" );

  legendre_set_x2 ( ORDER, xtab1, weight1 );

  for ( i = 0; i < ORDER; i++ )
  {
    xtab2[i] = xtab1[i];
  }
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight1[i];
  }

  rule_adjust ( a, b, c, d, ORDER, xtab2, weight2 );
/*
  Because the weight function W(X) is not 1, we need to do more
  adjustments to the weight vector.
*/
  for ( i = 0; i < ORDER; i++ )
  {
    weight2[i] = weight2[i] / 4.0;
  }
  legendre_set_x2_01 ( ORDER, xtab3, weight3 );

  printf ( "\n" );
  printf ( "  Abscissas:\n" );
  printf ( "\n" );
  printf ( "  Original  Adjusted Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n", i, xtab1[i], xtab2[i], xtab3[i] );
  }

  printf ( "\n" );
  printf ( "  Weights:\n" );
  printf ( "\n" );
  printf ( "  Original  Adjusted Stored\n" );
  printf ( "\n" );

  for ( i = 0; i < ORDER; i++ )
  {
    printf ( "  %2d  %16.12f  %16.12f  %16.12f\n", i, weight1[i], weight2[i], weight3[i] );
  }

  return;
# undef ORDER
}
/******************************************************************************/

void test34 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST34 tests LEGENDRE_SET_X2_01 and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 8;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = 0.0;
  b = 1.0;

  nsub = 1;

  xlo = 0.0;
  xhi = +1.0;

  printf ( "\n" );
  printf ( "TEST34\n" );
  printf ( "  LEGENDRE_SET_X2_01 sets up a Gauss-Legendre rule\n" );
  printf ( "    for integrating X*X * F(X) over [0,1];\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        legendre_set_x2_01 ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, xlo, xhi,
          xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }
  free ( result );

  return;
}
/******************************************************************************/

void test345 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST345 tests LOBATTO_COMPUTE and LOBATTO_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 February 2007

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST345\n" );
  printf ( "  LOBATTO_COMPUTE computes a Lobatto rule;\n" );
  printf ( "  LOBATTO_SET sets a rule from a table.\n" );
  printf ( "\n" );
  printf ( "         I      X1            X2            W1            W2\n" );

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w1 = ( double * ) malloc ( n * sizeof ( double ) );
    w2 = ( double * ) malloc ( n * sizeof ( double ) );
    x1 = ( double * ) malloc ( n * sizeof ( double ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );

    lobatto_compute ( n, x1, w1 );
    lobatto_set ( n, x2, w2 );

    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %12f  %12f  %12f  %12f\n", i+1, x1[i], x2[i], w1[i], w2[i] );
    }
    free ( w1 );
    free ( w2 );
    free ( x1 );
    free ( x2 );
  }
  return;
}
/******************************************************************************/

void test35 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST35 tests LOBATTO_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 15;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a = -1.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST35\n" );
  printf ( "  LOBATTO_SET sets up a Lobatto rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 1 )
      {
        continue;
      }

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        lobatto_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test36 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST36 tests MOULTON_SET and SUMMER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  int function_num;
  int i;
  int ihi;
  int ilo;
  int order;
  int order_max = 10;
  double *result;
  double *weight;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST36\n" );
  printf ( "  MOULTON_SET sets up an Adams-Moulton rule;\n" );
  printf ( "  SUMMER carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [0,1].\n" );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        moulton_set ( order, xtab, weight );
 
        result[i] = summer ( function_value, order, xtab, weight );
 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test37 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST37 tests NCC_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 21;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST37\n" );
  printf ( "  NCC_SET computes a closed Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncc_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test38 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST38 tests NCC_COMPUTE and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 21;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST38\n" );
  printf ( "  NCC_COMPUTE computes a closed Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncc_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test39 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST39 tests NCO_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST39\n" );
  printf ( "  NCO_SET sets up an open Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        nco_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test40 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST40 tests NCO_COMPUTE and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST40\n" );
  printf ( "  NCO_COMPUTE computes an open Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        nco_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test401 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST401 tests NCOH_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST401\n" );
  printf ( "  NCOH_SET sets up an open half Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncoh_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test402 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST402 tests NCOH_COMPUTE and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST402\n" );
  printf ( "  NCOH_COMPUTE computes an open half Newton-Cotes rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        ncoh_compute ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test403 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST403 tests PATTERSON_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int level;
  int level_max = 7;
  int nsub;
  int order;
  int order_max = 9;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST403\n" );
  printf ( "  PATTERSON_SET sets up a Patterson rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order   " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( level = 1; level <= level_max; level++ )
    {
      order = i4_power ( 2, level ) - 1;

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        patterson_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%3d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

void test404 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST404 tests RADAU_COMPUTE and RADAU_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2007

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST404\n" );
  printf ( "  RADAU_COMPUTE computes a Radau rule;\n" );
  printf ( "  RADAU_SET sets a rule from a table.\n" );
  printf ( "\n" );
  printf ( "         I      X1            X2            W1            W2\n" );

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w1 = ( double * ) malloc ( n * sizeof ( double ) );
    w2 = ( double * ) malloc ( n * sizeof ( double ) );
    x1 = ( double * ) malloc ( n * sizeof ( double ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );

    radau_compute ( n, x1, w1 );
    radau_set ( n, x2, w2 );

    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %12f  %12f  %12f  %12f\n", i+1, x1[i], x2[i], w1[i], w2[i] );
    }
    free ( w1 );
    free ( w2 );
    free ( x1 );
    free ( x2 );
  }
  return;
}
/******************************************************************************/

void test41 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST41 tests RADAU_SET and SUM_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int function_num;
  int i;
  int ihi;
  int ilo;
  int nsub;
  int order;
  int order_max = 15;
  double *result;
  double *weight;
  double xhi;
  double xlo;
  double *xtab;

  function_set ( "COUNT", &function_num );

  result = ( double * ) malloc ( function_num * sizeof ( double ) );

  a =  0.0;
  b =  1.0;

  nsub = 1;

  xlo = -1.0;
  xhi =  1.0;

  printf ( "\n" );
  printf ( "TEST41\n" );
  printf ( "  RADAU_SET sets up a Radau rule;\n" );
  printf ( "  SUM_SUB carries it out.\n" );
  printf ( "\n" );
  printf ( "  The integration interval is [%f,%f]\n", a, b );
  printf ( "  The number of subintervals is %d\n", nsub );
  printf ( "  Quadrature order will vary.\n" );
  printf ( "  Integrand will vary.\n" );
  printf ( "\n" );

  for ( ilo = 0; ilo < function_num; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, function_num - 1 );

    printf ( "\n" );
    printf ( "Order  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%10s    ", function_name ( i ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( order = 1; order <= order_max; order++ )
    {
      if ( order == 8 )
      {
        continue;
      }

      xtab = ( double * ) malloc ( order * sizeof ( double ) );
      weight = ( double * ) malloc ( order * sizeof ( double ) );

      for ( i = ilo; i <= ihi; i++ )
      {
        function_set ( "SET", &i );

        radau_set ( order, xtab, weight );
 
        result[i] = sum_sub ( function_value, a, b, nsub, order, 
          xlo, xhi, xtab, weight ); 
      }
      printf ( "%2d  ", order );
      for ( i = ilo; i <= ihi; i++ )
      {
        printf ( "  %12.8f", result[i] );
      }
      printf ( "\n" );

      free ( xtab );
      free ( weight );
    }
  }

  free ( result );
 
  return;
}
/******************************************************************************/

double f1sd1 ( double x )

/******************************************************************************/
/*
  Purpose:

    F1SD1 evaluates the function 1.0D+00/ sqrt ( 1.1 - x**2 ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double F1SD1, the value of the function.
*/
{
  double value;

  value = 1.0 / sqrt ( 1.1 - x * x );
 
  return value;
}
/******************************************************************************/

char *function_name ( int function_index )

/******************************************************************************/
/*
  Purpose:

    FUNCTION_NAME returns the name of the function evaluated in FUNCTION_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, integer FUNCTION_INDEX, the index of the function.

    Output, char *FUNCTION_NAME, the name of the function.
*/
{
  char *value;

  value = ( char * ) malloc ( 10 * sizeof ( char ) );

  if ( function_index == 0 )
  {
    strcpy ( value, "         1" );
  }
  else if ( function_index == 1 )
  {
    strcpy ( value, "         X" );
  }
  else if ( function_index == 2 )
  {
    strcpy ( value, "       X^2" );
  }
  else if ( function_index == 3 )
  {
    strcpy ( value, "       X^3" );
  }
  else if ( function_index == 4 )
  {
    strcpy ( value, "       X^4" );
  }
  else if ( function_index == 5 )
  {
    strcpy ( value, "       X^5" );
  }
  else if ( function_index == 6 )
  {
    strcpy ( value, "       X^6" );
  }
  else if ( function_index == 7 )
  {
    strcpy ( value, "       X^7" );
  }
  else if ( function_index == 8 )
  {
    strcpy ( value, "    SIN(X)" );
  }
  else if ( function_index == 9 )
  {
    strcpy ( value, "    EXP(X)" );
  }
  else if ( function_index == 10 )
  {
    strcpy ( value, " SQRT(|X|)" );
  }
  else
  {
    strcpy ( value, "??????????" );
  }
  return value;
}
/******************************************************************************/

void function_set ( char *action, int *i )

/******************************************************************************/
/*
  Purpose:

    FUNCTION_SET sets the function to be returned by FUNCTION_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, char *ACTION, the action to be carried out.
    "COUNT" means the call is made to count the number of functions available.
    "GET" means the call is made to find out the current function index.
    "SET" means the call is made to set the current function index.

    Input/output, int *I.
    For "COUNT", I is output as the number of functions available;
    For "GET", I is output as the currently chosen function;
    For "SET", I is input as the user's new choice for the function.
*/
{
  static int function_index = -1;

  if ( s_eqi ( action, "COUNT" ) )
  {
    *i = 11;
  }
  else if ( s_eqi ( action, "GET" ) )
  {
    *i = function_index;
  }
  else if ( s_eqi ( action, "SET" ) )
  {
    function_index = *i;
  }
  else
  {
    printf ( "\n" );
    printf ( "FUNCTION_SET - Warning!\n" );
    printf ( "  Unrecognized action = \"%s\"\n", action );
  }

  return;
}
/******************************************************************************/

double function_value ( double x )

/******************************************************************************/
/*
  Purpose:

    FUNCTION_VALUE evaluates a function of X, as chosen by the user.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double FUNCTION_VALUE, the value of the function.
*/
{
  int function_index;
  double value;

  function_set ( "GET", &function_index );

  if ( function_index == 0 )
  {
    value = 1.0;
  }
  else if ( function_index == 1 )
  {
    value = x;
  }
  else if ( function_index == 2 )
  {
    value = pow ( x, 2 );
  }
  else if ( function_index == 3 )
  {
    value = pow ( x, 3 );
  }
  else if ( function_index == 4 )
  {
    value = pow ( x, 4 );
  }
  else if ( function_index == 5 )
  {
    value = pow ( x, 5 );
  }
  else if ( function_index == 6 )
  {
    value = pow ( x, 6 );
  }
  else if ( function_index == 7 )
  {
    value = pow ( x, 7 );
  }
  else if ( function_index == 8 )
  {
    value = sin ( x );
  }
  else if ( function_index == 9 )
  {
    value = exp ( x );
  }
  else if ( function_index == 10 )
  {
    value = sqrt ( fabs ( x ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/******************************************************************************/

double fxsd1 ( double x )

/******************************************************************************/
/*
  Purpose:

    FXSD1 evaluates the function x / sqrt ( 1.1 - x**2 ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double FXSD1, the value of the function.
*/
{
  double value;

  value = x / sqrt ( 1.1 - x * x );
 
  return value;
}
/******************************************************************************/

double fx2sd1 ( double x )

/******************************************************************************/
/*
  Purpose;

    FX2SD1 evaluates the function x**2 / sqrt ( 1.1 - x**2 ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double FX2SD1, the value of the function.
*/
{
  double value;

  value = x * x / sqrt ( 1.1 - x * x );
 
  return value;
}
