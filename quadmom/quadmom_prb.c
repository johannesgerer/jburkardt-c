# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "quadmom.h"
# include "toms655.h"

int main ( );
void quadmom_prb01 ( );
void quadmom_prb02 ( );
void quadmom_prb03 ( );
void quadmom_prb04 ( );
void quadmom_prb05 ( );
void quadmom_prb06 ( );
void quadmom_prb07 ( );
void quadmom_prb08 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QUADMOM_PRB.

  Discussion:

    QUADMOM_PRB tests the QUADMOM library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 October 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "QUADMOM_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the QUADMOM library.\n" );

  quadmom_prb01 ( );
  quadmom_prb02 ( );
  quadmom_prb03 ( );
  quadmom_prb04 ( );
  quadmom_prb05 ( );
  quadmom_prb06 ( );
  quadmom_prb07 ( );
  quadmom_prb08 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUADMOM_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

 return 0;
}
/******************************************************************************/

void quadmom_prb01 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB01 tests the QUADMOM procedure for the Legendre weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "QUADMOM_PRB01:\n" );
  printf ( "  Compute the points and weights of a quadrature rule\n" );
  printf ( "  for the Legendre weight, rho(x)=1, over [-1,+1],\n" );
  printf ( "  using Golub and Welsch's moment method.\n" );
  printf ( "  Compare with a standard calculation.\n" );
/*
  N is the order of the rule we want to compute.
*/
  n = 5;
/*
  Compute M = 2*N+1 moments for the Legendre weight on [-1,+1].
*/
  m = 2 * n + 1;
  a = -1.0;
  b = 1.0;

  moment = moments_legendre ( m, a, b );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Compute the points and weights the usual way.
*/
  kind = 1;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = -1.0;
  b = +1.0;
  lo = 0;
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  w2 = ( double * ) malloc ( n * sizeof ( double ) );

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
/*
  Compare the results.
*/
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( w2 );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void quadmom_prb02 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB02 tests the QUADMOM procedure for the standard Gaussian weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int i;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  const double pi = 3.141592653589793;
  double sigma;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "QUADMOM_PRB02:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  the standard Gaussian weight, rho(x)=exp(-x^2/2)/sqrt(2pi),\n" );
  printf ( "  over (-oo,+oo), using Golub and Welsch's moment method.\n" );
  printf ( "  Compare with a standard calculation.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 5;
/*
  Compute the M = 2 * N + 1 moments for the standard Gaussian weight on (-oo,+oo).
*/
  m = 2 * n + 1;
  moment = moments_normal_01 ( m );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Compute the points and weights the usual way.
*/
  kind = 6;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 0.0;
  b = +0.5;
  lo = 0;
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  w2 = ( double * ) malloc ( n * sizeof ( double ) );

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
/*
  The CGQF weights need to be normalized by sigma * sqrt ( 2 * pi )
  because they don't divide the Gaussian PDF by that factor.
*/
  sigma = 1.0;
  for ( i = 0; i < n; i++ )
  {
    w2[i] = w2[i] / sigma / sqrt ( 2.0 * pi );
  }
/*
  Compare the results.
*/
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( w2 );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void quadmom_prb03 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB03 tests the QUADMOM procedure for the general Gaussian weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int i;
  int kind;
  int lo;
  int m;
  double *moment;
  double mu;
  int n;
  const double pi = 3.141592653589793;
  double sigma;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "QUADMOM_PRB03:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  a general Gaussian weight,\n" );
  printf ( "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n" );
  printf ( "  over (-oo,+oo), using Golub and Welsch''s moment method.\n" );
  printf ( "  Compare with a standard calculation.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 5;
/*
  Compute the M = 2 * N + 1 moments for a general Gaussian weight on (-oo,+oo).
*/
  m = 2 * n + 1;
  mu = 1.0;
  sigma = 2.0;

  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
  printf ( "  SIGMA = %g\n", sigma );

  moment = moments_normal ( m, mu, sigma );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Compute the points and weights the usual way.
*/
  kind = 6;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 1.0;
  b = 0.5 / sigma / sigma;
  lo = 0;

  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  w2 = ( double * ) malloc ( n * sizeof ( double ) );

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
/*
  The CGQF weights need to be normalized by sigma * sqrt ( 2 * pi )
  because they don't divide the Gaussian PDF by that factor.
*/
  for ( i = 0; i < n; i++ )
  {
    w2[i] = w2[i] / sigma / sqrt ( 2.0 * pi );
  }
/*
  Compare the results.
*/
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( w2 );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void quadmom_prb04 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB04 tests the QUADMOM procedure for the Laguerre weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  double alpha2;
  double b;
  double beta2;
  int kind;
  int lo;
  int m;
  double *moment;
  int n;
  const double pi = 3.141592653589793;
  double *w1;
  double *w2;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "QUADMOM_PRB04:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  the Laguerre weight, rho(x)=exp(-x),\n" );
  printf ( "  over [0,+oo), using Golub and Welsch's moment method.\n" );
  printf ( "  Compare with a standard calculation.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 5;
/*
  Compute the M = 2 * N + 1 moments for the Laguerre weight on [0,+oo).
*/
  m = 2 * n + 1;
  moment = moments_laguerre ( m );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Compute the points and weights the usual way.
*/
  kind = 5;
  alpha2 = 0.0;
  beta2 = 0.0;
  a = 0.0;
  b = +1.0;
  lo = 0;

  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  w2 = ( double * ) malloc ( n * sizeof ( double ) );

  cgqf ( n, kind, alpha2, beta2, a, b, lo, x2, w2 );
/*
  Compare the results.
*/
  r8vec2_print ( n, x1, x2, 
    "  Points from GW moment and orthogonal polynomial methods:" );

  r8vec2_print ( n, w1, w2, 
    "  Weights from GW moment and orthogonal polynomial methods:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( w2 );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void quadmom_prb05 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB05 tests the QUADMOM procedure for the truncated normal weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  double b;
  int m;
  double *moment;
  double mu;
  int n;
  double sigma;
  double *w1;
  double *x1;

  printf ( "\n" );
  printf ( "QUADMOM_PRB05:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  a truncated normal weight,\n" );
  printf ( "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n" );
  printf ( "  over [a,b], using Golub and Welsch's moment method.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 5;
/*
  Compute the M = 2 * N + 1 moments.
*/
  m = 2 * n + 1;
  mu = 100.0;
  sigma = 25.0;
  a = 50.0;
  b = 150.0;

  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
  printf ( "  SIGMA = %g\n", sigma );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );

  moment = moments_truncated_normal_ab ( m, mu, sigma, a, b );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Print the results.
*/
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( x1 );

  return;
}
/******************************************************************************/

void quadmom_prb06 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB06 tests QUADMOM for the lower truncated normal weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  int m;
  double *moment;
  double mu;
  int n;
  double sigma;
  double *w1;
  double *x1;

  printf ( "\n" );
  printf ( "QUADMOM_PRB06:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  a lower truncated normal weight,\n" );
  printf ( "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n" );
  printf ( "  over [a,+oo), using Golub and Welsch's moment method.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 9;
/*
  Compute the M = 2 * N + 1 moments.
*/
  m = 2 * n + 1;
  mu = 2.0;
  sigma = 0.5;
  a = 0.0;

  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
  printf ( "  SIGMA = %g\n", sigma );
  printf ( "  A = %g\n", a );

  moment = moments_truncated_normal_a ( m, mu, sigma, a );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Print the results.
*/
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( x1 );

  return;
}
/******************************************************************************/

void quadmom_prb07 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB07 tests QUADMOM procedure for the upper truncated normal weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double b;
  int m;
  double *moment;
  double mu;
  int n;
  double sigma;
  double *w1;
  double *x1;

  printf ( "\n" );
  printf ( "QUADMOM_PRB07:\n" );
  printf ( "  Compute the points and weights of a quadrature rule for\n" );
  printf ( "  a upper truncated normal weight,\n" );
  printf ( "  rho(mu,s;x)=exp(-((x-mu)/sigma)^2/2)/sigma^2/sqrt(2pi),\n" );
  printf ( "  over (-oo,b], using Golub and Welsch's moment method.\n" );
/*
  N is the order of the quadrature rule.
*/
  n = 9;
/*
  Compute the M = 2 * N + 1 moments.
*/
  m = 2 * n + 1;
  mu = 2.0;
  sigma = 0.5;
  b = 3.0;

  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
  printf ( "  SIGMA = %g\n", sigma );
  printf ( "  B = %g\n", b );

  moment = moments_truncated_normal_b ( m, mu, sigma, b );
/*
  Compute the points and weights by the moment method.
*/
  x1 = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  moment_method ( n, moment, x1, w1 );
/*
  Print the results.
*/
  r8vec_print ( n, x1, "  Points from GW moment method:" );

  r8vec_print ( n, w1, "  Weights from GW moment method:" );
/*
  Free memory.
*/
  free ( moment );
  free ( w1 );
  free ( x1 );

  return;
}
/******************************************************************************/

void quadmom_prb08 ( )

/******************************************************************************/
/*
  Purpose:

    QUADMOM_PRB08 integrates sin(x) against a lower truncated normal weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
*/
{
  double a;
  int i;
  int m;
  double *moment;
  double mu;
  int n;
  double q;
  double sigma;
  double *w1;
  double *x1;

  printf ( "\n" );
  printf ( "QUADMOM_PRB08:\n" );
  printf ( "  Integrate sin(x) against a lower truncated normal weight.\n" );

  mu = 0.0;
  sigma = 1.0;
  a = -3.0;

  printf ( "\n" );
  printf ( "  MU = %g\n", mu );
  printf ( "  SIGMA = %g\n", sigma );
  printf ( "  A = %g\n", a );
  printf ( "\n" );
  printf ( "   N   Estimate\n" );
  printf ( "\n" );
/*
  N is the order of the quadrature rule.
*/
  for ( n = 1; n <= 9; n++ )
  {
/*
  Compute the M = 2 * N + 1 moments.
*/
    m = 2 * n + 1;

    moment = moments_truncated_normal_a ( m, mu, sigma, a );
/*
  Compute the points and weights by the moment method.
*/
    x1 = ( double * ) malloc ( n * sizeof ( double ) );
    w1 = ( double * ) malloc ( n * sizeof ( double ) );

    moment_method ( n, moment, x1, w1 );

    q = 0.0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w1[i] * sin ( x1[i] );
    }
    printf ( "  %2d  %14.6g\n", n, q );
/*
  Free memory.
*/
    free ( moment );
    free ( w1 );
    free ( x1 );
  }

  return;
}


