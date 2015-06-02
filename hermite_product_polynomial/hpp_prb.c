# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "hpp.h"

int main ( );
void hpp_test01 ( );
void hpp_test015 ( );
void hpp_test02 ( );
void hpp_test03 ( );
void hpp_test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HPP_PRB.

  Discussion:

    HPP_PRB tests the HERMITE_PRODUCT_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HPP_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HERMITE_PRODUCT_POLYNOMIAL library.\n" );

  hpp_test01 ( );
  hpp_test015 ( );
  hpp_test02 ( );
  hpp_test03 ( );
  hpp_test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HPP_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void hpp_test01 ( )

/******************************************************************************/
/*
  Purpose:

    HPP_TEST01 tests routines for the GRLEX ordering of compositions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 September 2014

  Author:

    John Burkardt
*/
{
  int i;
  int k = 2;
  int rank;
  int rank1;
  int rank2;
  int seed;
  int test;
  int *x;
  int x_sum;
  int x_sum_old;

  x = ( int * ) malloc ( k * sizeof ( int ) );

  printf ( "\n" );
  printf ( "HPP_TEST01:\n" );
  printf ( "  COMP_NEXT_GRLEX is given a composition, and computes the \n" );
  printf ( "  next composition in grlex order.\n" );

  printf ( "\n" );
  printf ( "  Rank   Sum   Components\n" );

  for ( i = 0; i < k; i++ )
  {
    x[i] = 0;
  }
  x_sum_old = -1;
  rank = 1;

  for ( ; ; )
  {
    x_sum = i4vec_sum ( k, x );

    if ( x_sum_old < x_sum )
    {
      x_sum_old = x_sum;
      printf ( "\n" );
    }

    printf ( "  %4d  %4d:", rank, x_sum );
    for ( i = 0; i < k; i++ )
    {
      printf ( "%4d", x[i] );
    }
    printf ( "\n" );

    if ( 20 <= rank )
    {
      break;
    }

    comp_next_grlex ( k, x );
    rank = rank + 1;
  }
  free ( x );

  printf ( "\n" );
  printf ( "  COMP_UNRANK_GRLEX is given a rank and returns the\n" );
  printf ( "  corresponding set of multinomial exponents.\n" );
  printf ( "\n" );
  printf ( "  Rank   Sum   Components\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    rank = i4_uniform_ab ( 1, 20, &seed );
    x = comp_unrank_grlex ( k, rank );
    x_sum = i4vec_sum ( k, x );
    printf ( "  %4d  %4d:", rank, x_sum );
    for ( i = 0; i < k; i++ )
    {
      printf ( "%4d", x[i] );
    }
    printf ( "\n" );
    free ( x );
  }

  printf ( "\n" );
  printf ( "  COMP_RANDOM_GRLEX randomly selects a composition\n" );
  printf ( "  between given lower and upper ranks.\n" );
  printf ( "\n" );
  printf ( "  Rank   Sum   Components\n" );
  printf ( "\n" );

  seed = 123456789;
  rank1 = 5;
  rank2 = 20;

  for ( test = 1; test <= 5; test++ )
  {
    x = comp_random_grlex ( k, rank1, rank2, &seed, &rank );
    x_sum = i4vec_sum ( k, x );
    printf ( "  %4d  %4d:", rank, x_sum );
    for ( i = 0; i < k; i++ )
    {
      printf ( "%4d", x[i] );
    }
    printf ( "\n" );
    free ( x );
  }

  printf ( "\n" );
  printf ( "  COMP_RANK_GRLEX returns the rank of a given composition.\n" );
  printf ( "\n" );
  printf ( "  Rank   Sum   Components\n" );
  printf ( "\n" );

  x = ( int * ) malloc ( k * sizeof ( int ) );

  x[0] = 4;
  x[1] = 0; 
  rank = comp_rank_grlex ( k, x );
  x_sum = i4vec_sum ( k, x );
  printf ( "  %4d  %4d:", rank, x_sum );
  for ( i = 0; i < k; i++ )
  {
    printf ( "%4d", x[i] );
  }
  printf ( "\n" );

  x[0] = 11;
  x[1] = 5; 
  rank = comp_rank_grlex ( k, x );
  x_sum = i4vec_sum ( k, x );
  printf ( "  %4d  %4d:", rank, x_sum );
  for ( i = 0; i < k; i++ )
  {
    printf ( "%4d", x[i] );
  }
  printf ( "\n" );;

  free ( x );

  return;
}
/******************************************************************************/

void hpp_test015 ( )

/******************************************************************************/
/*
  Purpose:

    HPP_TEST015 tests HEP_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int *f;
  int l[1];
  int m;
  int n;
  int o;
  int o_max;
  char title[255];

  m = 1;

  printf ( "\n" );
  printf ( "HPP_TEST015:\n" );
  printf ( "  HEP_COEFFICIENTS computes the coefficients and\n" );
  printf ( "  exponents of the Hermite polynomial He(n,x).\n" );

  for ( n = 1; n <= 5; n++ )
  {
    o = ( n + 2 ) / 2;
    c = ( double * ) malloc ( o * sizeof ( double ) );
    e = ( int * ) malloc ( o * sizeof ( int ) );
    f = ( int * ) malloc ( o * sizeof ( int ) );

    hep_coefficients ( n, &o, c, f );

    l[0] = n;
    o_max = o;

    hepp_to_polynomial ( m, l, o_max, &o, c, e );

    printf ( "\n" );
    sprintf ( title, "  He(%d,x) =", n );
    polynomial_print ( m, o, c, e, title );

    free ( c );
    free ( e );
    free ( f );
  }

  return;
}
/******************************************************************************/

void hpp_test02 ( )

/******************************************************************************/
/*
  Purpose:

    HPP_TEST02 tests HEP_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2014

  Author:

    John Burkardt
*/
{
  double e;
  int n;
  int n_data;
  int o;
  double x;
  double xvec[1];
  double fx1;
  double *fx2;

  n = 1;

  printf ( "\n" );
  printf ( "HPP_TEST02:\n" );
  printf ( "  HEP_VALUES stores values of\n" );
  printf ( "  the Hermite polynomial He(o,x).\n" );
  printf ( "  HEP_VALUE evaluates a Hermite polynomial.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     O        X          He(O,X)                   He(O,X)" );
  printf ( "                   Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hep_values ( &n_data, &o, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    xvec[0] = x;

    fx2 = hep_value ( n, o, xvec );

    e = fx1 - fx2[0];

    printf ( "  %4d  %12.8f  %24.16g  %24.16g  %8.2g\n",
      o, x, fx1, fx2[0], e );

    free ( fx2 );
  }

  return;
}
/******************************************************************************/

void hpp_test03 ( )

/******************************************************************************/
/*
  Purpose:

    HPP_TEST03 tests HEPP_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int i;
  int *l;
  int m = 3;
  int n = 1;
  int o;
  int o_max;
  int rank;
  int seed;
  double *v1;
  double *v2;
  double *x;
  double xhi;
  double xlo;

  printf ( "\n" );
  printf ( "HPP_TEST03:\n" );
  printf ( "  HEPP_VALUE evaluates a Hermite product polynomial.\n" );
  printf ( "  POLYNOMIAL_VALUE evaluates a polynomial.\n" );

  xlo = -1.0;
  xhi = +1.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, xlo, xhi, &seed );

  printf ( "\n" );
  printf ( "  Evaluate at X = " );
  for ( i = 0; i < m; i++ )
  {
    printf ( "%g  ", x[i+0*m] );
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  Rank  I1  I2  I3: He(I1,X1)*He(I2,X2)*He(I3,X3)    P(X1,X2,X3)\n" );
  printf ( "\n" );

  for ( rank = 1; rank <= 20; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );
/*
  Evaluate the LPP directly.
*/
    v1 = hepp_value ( m, n, l, x );
/*
  Convert the HePP to a polynomial.
*/
    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = ( double * ) malloc ( o_max * sizeof ( double ) );
    e = ( int * ) malloc ( o_max * sizeof ( int ) );
 
    hepp_to_polynomial ( m, l, o_max, &o, c, e );
/*
  Evaluate the polynomial.
*/
    v2 = polynomial_value ( m, o, c, e, n, x );
/*
  Compare results.
*/
    printf ( "  %4d  %2d  %2d  %2d  %14.6g  %14.6g\n",
      rank, l[0], l[1], l[2], v1[0], v2[0] );
 
    free ( c );
    free ( e );
    free ( l );
    free ( v1 );
    free ( v2 );
  }

  free ( x );

  return;
}
/******************************************************************************/

void hpp_test04 ( )

/******************************************************************************/
/*
  Purpose:

    HPP_TEST04 tests HEPP_TO_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int i;
  int *l;
  char label[255];
  int m = 2;
  int o;
  int o_max;
  int rank;

  printf ( "\n" );
  printf ( "HPP_TEST04:\n" );
  printf ( "  HEPP_TO_POLYNOMIAL is given a Hermite product polynomial\n" );
  printf ( "  and determines its polynomial representation.\n" );

  printf ( "\n" );
  printf ( "  Using spatial dimension M = %d\n", m );

  for ( rank = 1; rank <= 11; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );

    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = ( double * ) malloc ( o_max * sizeof ( double ) );
    e = ( int * ) malloc ( o_max * sizeof ( int ) );

    hepp_to_polynomial ( m, l, o_max, &o, c, e );

    sprintf ( label, "  HePP #%d = He(%d,X)*He(%d,Y) = ", rank, l[0], l[1] );

    printf ( "\n" );
    polynomial_print ( m, o, c, e, label );

    free ( c );
    free ( e );
    free ( l );
  }

  return;
}

