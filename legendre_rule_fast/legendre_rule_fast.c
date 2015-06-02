# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
void legendre_compute_glr ( int n, double x[], double w[] );
void legendre_compute_glr0 ( int n, double *p, double *pp );
void legendre_compute_glr1 ( int n, double *roots, double *ders );
void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
void legendre_handle ( int n, double a, double b );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
double rk2_leg ( double t, double tn, double x, int n );
void timestamp ( void );
double ts_mult ( double *u, double h, int n );
double wtime ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LEGENDRE_RULE_FAST.

  Discussion:

    This program computes a Gauss-Legendre quadrature rule
    and writes it to a file.

  Usage:

    legendre_rule_fast n a b

    where

    * n is the number of points in the rule;
    * a is the left endpoint;
    * b is the right endpoint.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int n;

  timestamp ( );
  printf ( "\n" );
  printf ( "LEGENDRE_RULE\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Compute a Gauss-Legendre rule for approximating\n" );
  printf ( "\n" );
  printf ( "    Integral ( a <= x <= b ) f(x) dx\n" );
  printf ( "\n" );
  printf ( "  of order N.\n" );
  printf ( "\n" );
  printf ( "  The user specifies N, A, and B.\n" );
  printf ( "\n" );
  printf ( "  The rule is stored in 3 files:\n" );
  printf ( "\n" );
  printf ( "  * leg_oN_w.txt - the weight file\n" );
  printf ( "  * leg_oN_x.txt - the abscissa file.\n" );
  printf ( "  * leg_oN_r.txt - the region file.\n" );
/*
  Get N.
*/
  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter N:\n" );
    scanf ( "%d", &n );
  }
  printf ( "\n" );
  printf ( "  N = %d\n", n );
/*
  Get A:
*/
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter A.\n" );
    scanf ( "%lf", &a );
  }
  printf ( "\n" );
  printf ( "  A = %g\n", a );
/*
  Get B:
*/
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter B.\n" );
    scanf ( "%lf", &b );
  }
  printf ( "\n" );
  printf ( "  B = %g\n", b );
/*
  Construct the rule and output it.
*/
  legendre_handle ( n, a, b );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LEGENDRE_RULE_FAST:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void legendre_compute_glr ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 October 2009

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Reference:

    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
    A fast algorithm for the calculation of the roots of special functions, 
    SIAM Journal on Scientific Computing,
    Volume 29, Number 4, pages 1420-1438, 2007.

  Parameters:

    Input, int N, the order.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  int i;
  double p;
  double pp;
  double w_sum;
/*
  Get the value and derivative of the N-th Legendre polynomial at 0.
*/
  legendre_compute_glr0 ( n, &p, &pp );
/*
  Either zero is a root, or we have to call a function to find the first root.
*/  
  if ( n % 2 == 1 )
  {
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }
  else
  {
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
/*
  Get the complete set of roots and derivatives.
*/
  legendre_compute_glr1 ( n, x, w );
/*
  Compute the weights.
*/
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr0 ( int n, double *p, double *pp )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 October 2009

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Reference:

    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
    A fast algorithm for the calculation of the roots of special functions, 
    SIAM Journal on Scientific Computing,
    Volume 29, Number 4, pages 1420-1438, 2007.

  Parameters:

    Input, int N, the order of the Legendre polynomial.

    Output, double *P, *PP, the value of the N-th Legendre polynomial
    and its derivative at 0.
*/
{
  double dk;
  int k;
  double pm1;
  double pm2;
  double ppm1;
  double ppm2;

  pm2 = 0.0;
  pm1 = 1.0;
  ppm2 = 0.0;
  ppm1 = 0.0;

  for ( k = 0; k < n; k++ )
  {
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr1 ( int n, double *x, double *ders )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.

  Discussion:

    This routine requires that a starting estimate be provided for one
    root and its derivative.  This information will be stored in entry
    (N+1)/2 if N is odd, or N/2 if N is even, of ROOTS and DERS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 October 2009

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Reference:

    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
    A fast algorithm for the calculation of the roots of special functions, 
    SIAM Journal on Scientific Computing,
    Volume 29, Number 4, pages 1420-1438, 2007.

  Parameters:

    Input, int N, the order of the Legendre polynomial.

    Input/output, double X[N].  On input, a starting value
    has been set in one entry.  On output, the roots of the Legendre 
    polynomial.

    Input/output, double DERS[N].  On input, a starting value
    has been set in one entry.  On output, the derivatives of the Legendre 
    polynomial at the zeros.

  Local Parameters:

    Local, int M, the number of terms in the Taylor expansion.
*/
{
  double dk;
  double dn;
  double h;
  int j;
  int k;
  int l;
  int m = 30;
  int n2;
  const double pi = 3.141592653589793;
  int s;
  double *u;
  double *up;
  double xp;

  if ( n % 2 == 1 )
  {
    n2 = ( n - 1 ) / 2;
    s = 1;
  }
  else
  {
    n2 = n / 2;
    s = 0;
  }

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;

  for ( j = n2; j < n - 1; j++ )
  {
    xp = x[j];

    h = rk2_leg ( pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = ders[j];

    up[0] = 0.0;
    up[1] = u[2];

    for ( k = 0; k <= m - 2; k++ )
    {
      dk = ( double ) k;

      u[k+3] = 
      ( 
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 )
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }

    for ( l = 0; l < 5; l++ )
    { 
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    ders[j+1] = ts_mult ( up, h, m-1 );
  }

  free ( u );
  free ( up );

  for ( k = 0; k < n2 + s; k++ )
  {
    x[k] = - x[n-k-1];
    ders[k] = ders[n-k-1];
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr2 ( double pn0, int n, double *x1,  double *d1 )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_GLR2 finds the first real root.

  Discussion:

    This routine is only called if N is even.

    Thanks to Morten Welinder, for pointing out a typographical error
    in indexing, 17 May 2013.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2013

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Reference:

    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
    A fast algorithm for the calculation of the roots of special functions, 
    SIAM Journal on Scientific Computing,
    Volume 29, Number 4, pages 1420-1438, 2007.

  Parameters:

    Input, double PN0, the value of the N-th Legendre polynomial at 0.

    Input, int N, the order of the Legendre polynomial.

    Output, double *X1, the first real root.

    Output, double *D1, the derivative at X1.

  Local Parameters:

    Local, int M, the number of terms in the Taylor expansion.
*/
{
  double dk;
  double dn;
  int k;
  int l;
  int m = 30;
  const double pi = 3.141592653589793;
  double t;
  double *u;
  double *up;

  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;
/*
  U[0] and UP[0] are never used.
  U[M+1] is set, but not used, and UP[M] is set and not used.
  What gives?
*/
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;
 
  for ( k = 0; k <= m - 2; k = k + 2 )
  {
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / ( dk + 1.0 ) / ( dk + 2.0 );
 
    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }
  
  for ( l = 0; l < 5; l++ )
  {
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  free ( u );
  free ( up) ;

  return;
}
/******************************************************************************/

void legendre_handle ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 October 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the rule.

    Input, double A, B, the left and right endpoints.
*/ 
{
  int i;
  char output_r[255];
  char output_w[255];
  char output_x[255];
  double *r;
  double t;
  double *w;
  double *x;

  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );

  r[0] = a;
  r[1] = b;
/*
  Compute the rule.
*/
  t = wtime ( );
  legendre_compute_glr ( n, x, w );
  t = wtime ( ) - t;

  printf ( "\n" );
  printf ( "  Elapsed time during computation was %g seconds.\n", t );
/*
  Rescale the rule to [A,B].
*/
  rescale ( a, b, n, x, w );
/*
  Write the rule to 3 files.
*/
  sprintf ( output_w, "leg_o%d_w.txt", n );
  sprintf ( output_x, "leg_o%d_x.txt", n );
  sprintf ( output_r, "leg_o%d_r.txt", n );

  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", output_w );
  printf ( "  Abscissa file will be \"%s\".\n", output_x );
  printf ( "  Region file will be   \"%s\".\n", output_r );
            
  r8mat_write ( output_w, 1, n, w );
  r8mat_write ( output_x, 1, n, x );
  r8mat_write ( output_r, 1, 2, r );

  free ( r );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the table data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void rescale ( double a, double b, int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 October 2009

  Author:

    Original MATLAB version by Nick Hale.
    C version by John Burkardt.

  Reference:

    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
    A fast algorithm for the calculation of the roots of special functions, 
    SIAM Journal on Scientific Computing,
    Volume 29, Number 4, pages 1420-1438, 2007.

  Parameters:

    Input, double A, B, the endpoints of the new interval.

    Input, int N, the order.

    Input/output, double X[N], on input, the abscissas for [-1,+1].
    On output, the abscissas for [A,B].

    Input/output, double W[N], on input, the weights for [-1,+1].
    On output, the weights for [A,B].
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
/******************************************************************************/

double rk2_leg ( double t1, double t2, double x, int n )

/******************************************************************************/
/*
  Purpose:

    RK2_LEG advances the value of X(T) using a Runge-Kutta method.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 October 2009

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Parameters:

    Input, double T1, T2, the range of the integration interval.

    Input, double X, the value of X at T1.

    Input, int N, the number of steps to take.

    Output, double RK2_LEG, the value of X at T2.
*/
{
  double f;
  double h;
  int j;
  double k1;
  double k2;
  int m = 10;
  double snn1;
  double t;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );

  t = t1;

  for ( j = 0; j < m; j++ )
  {
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;

    t = t + h;

    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );   
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
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
/******************************************************************************/

double ts_mult ( double *u, double h, int n )

/******************************************************************************/
/*
  Purpose:

    TS_MULT evaluates a polynomial.

  Discussion:

    TS_MULT = U[1] + U[2] * H + ... + U[N] * H^(N-1).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2013

  Author:

    Original C version by Nick Hale.
    This C version by John Burkardt.

  Parameters:

    Input, double U[N+1], the polynomial coefficients.
    U[0] is ignored.

    Input, double H, the polynomial argument.

    Input, int N, the number of terms to compute.

    Output, double TS_MULT, the value of the polynomial.
*/
{
  double hk;
  int k;
  double ts;
  
  ts = 0.0;
  hk = 1.0;
  for ( k = 1; k<= n; k++ )
  {
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}
/******************************************************************************/

double wtime ( void )

/******************************************************************************/
/*
  Purpose:

    WTIME estimates the elapsed wall clock time.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 October 2009

  Author:

    John Burkardt

  Parameters:

    Output, double WTIME, the current elapsed wall clock time.
*/
{
  double now;

  now = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC; 

  return now;
}  
