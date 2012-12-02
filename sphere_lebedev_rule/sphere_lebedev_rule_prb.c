# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sphere_lebedev_rule.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_LEBEDEV_RULE_PRB.

  Discussion:

    SPHERE_LEBEDEV_RULE_PRB calls the SPHERE_LEBEDEV_RULE tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 September 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_LEBEDEV_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_LEBEDEV_RULE library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_LEBEDEV_RULE_PRB\n" );
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

    TEST01 tests AVAILABLE_TABLE, ORDER_TABLE, PRECISION_TABLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

   12 September 2010

  Author:

    John Burkardt
*/
{
  int available;
  int order;
  int precision;
  int rule;
  int rule_max = 65;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  List Lebedev rule properties.\n" );
  printf ( "\n" );
  printf ( "  Rule Avail Order  Prec\n" );
  printf ( "\n" );
  for ( rule = 1; rule <= rule_max; rule++ )
  {
    available = available_table ( rule );
    order = order_table ( rule );
    precision = precision_table ( rule );
    printf ( "  %4d  %4d  %4d  %4d\n", rule, available, order, precision );
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests the SPHERE_LEBEDEV_RULE functions.

  Modified:

    11 September 2010

  Author:

    Dmitri Laikov

  Reference:

    Vyacheslav Lebedev, Dmitri Laikov,
    A quadrature formula for the sphere of the 131st
    algebraic order of accuracy,
    Russian Academy of Sciences Doklady Mathematics,
    Volume 59, Number 3, 1999, pages 477-481.
*/
{
# define nmax 65
# define mmax ((nmax*2+3)*(nmax*2+3)/3)

  double alpha;
  int available;
  double beta;
  double err;
  double err_max;
  int i;
  double integral_exact;
  double integral_approx;
  int j;
  int k;
  int m;
  int n;
  int order;
  double s[nmax+2];
  double *w;
  double *x;
  double *xn;
  double *y;
  double *yn;
  double *z;
  double *zn;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Generate each available rule and test for accuracy.\n" );
  
  for ( n = 1; n <= nmax; n++ )
  {
    available = available_table ( n );

    if ( available )
    {
      order = order_table ( n );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      x = ( double * ) malloc ( order * sizeof ( double ) );
      y = ( double * ) malloc ( order * sizeof ( double ) );
      z = ( double * ) malloc ( order * sizeof ( double ) );

      xn = ( double * ) malloc ( ( order + 1 ) * ( n + 1 ) * sizeof ( double ) );
      yn = ( double * ) malloc ( ( order + 1 ) * ( n + 1 ) * sizeof ( double ) );
      zn = ( double * ) malloc ( ( order + 1 ) * ( n + 1 ) * sizeof ( double ) );

      ld_by_order ( order, x, y, z, w );

      s[0] = 1.0;
      for ( k = 1; k <= n + 1; k++ ) 
      {
        s[k] = ( 2 * k - 1 ) * s[k-1];
      }
/*
  For each abscissa X(M), compute the values 1, X(M)^2, X(M)^4, ..., X(M)^2*N.
*/
      for ( m = 0; m < order; m++ )
      {
        xn[m*(n+1)] = 1.0;
        yn[m*(n+1)] = 1.0;
        zn[m*(n+1)] = 1.0;
        for ( k = 1; k <= n; k++ ) 
        {
          xn[k+m*(n+1)] = xn[k-1+m*(n+1)] * x[m] * x[m];
          yn[k+m*(n+1)] = yn[k-1+m*(n+1)] * y[m] * y[m];
          zn[k+m*(n+1)] = zn[k-1+m*(n+1)] * z[m] * z[m];
        }
      }

      err_max = 0.0;
      for ( i = 0; i <= n; i++ ) 
      {
        for ( j = 0; j <= n - i; j++ ) 
        {
          k = n - i - j;
/*
  Apply Lebedev rule to x^2i y^2j z^2k.
*/
          integral_approx = 0.0;
          for ( m = 0; m < order; m++ ) 
          {
            integral_approx = integral_approx 
              + w[m] * xn[i+m*(n+1)] * yn[j+m*(n+1)] * zn[k+m*(n+1)];
          }
/*
  Compute exact value of integral (aside from factor of 4 pi!).
*/
          integral_exact = s[i] * s[j] * s[k] / s[1+i+j+k];
/*
  Record the maximum error for this rule.
*/
          err = fabs ( ( integral_approx - integral_exact ) / integral_exact );
          if ( err_max < err ) 
          {
            err_max = err;
          }
        }
      }
      printf ( "\n" );
      printf ( "  Order = %4d  LMAXW = %4d  max error = %8.1e\n", 
        order, precision_table ( n ), err_max );
/*
  Convert (X,Y,Z) to (Theta,Phi) and print the data.
*/
      if ( order <= 50 )
      {
        for ( m = 0; m < order; m++ ) 
        {
          xyz_to_tp ( x[m], y[m], z[m], &alpha, &beta );
          printf ( "%20.15lf  %20.15lf  %20.15lf\n", alpha, beta, w[m] );
        }
      }
      free ( w );
      free ( x );
      free ( xn );
      free ( y );
      free ( yn );
      free ( z );
      free ( zn );
    }
  }
  return;

# undef nmax
# undef mmax
}
