# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] );
void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] );
double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] );
void imtqlx ( int n, double d[], double e[], double z[] );
void parchk ( int kind, int m, double alpha, double beta );
double r8_epsilon ( );
double r8_gamma ( double x );
double r8_sign ( double x );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rule_write ( int order, char *filename, double x[], double w[], 
  double r[] );
void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b );
void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LEGENDRE_RULE.

  Discussion:

    This program computes a standard Gauss-Legendre quadrature rule
    and writes it to a file.

    The user specifies:
    * the ORDER (number of points) in the rule
    * A, the left endpoint;
    * B, the right endpoint;
    * FILENAME, the root name of the output files.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double alpha;
  double b;
  double beta;
  char filename[255];
  int kind;
  int order;
  double *r;
  double *w;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "LEGENDRE_RULE\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Compute a Gauss-Legendre rule for approximating\n" );
  printf ( "    Integral ( A <= x <= B ) f(x) dx\n" );
  printf ( "  of order ORDER.\n" );
  printf ( "\n" );
  printf ( "  The user specifies ORDER, A, B and FILENAME.\n" );
  printf ( "\n" );
  printf ( "  ORDER is the number of points.\n" );
  printf ( "  A is the left endpoint.\n" );
  printf ( "  B is the right endpoint.\n" );
  printf ( "  FILENAME is used to generate 3 files:\n" );
  printf ( "    filename_w.txt - the weight file\n" );
  printf ( "    filename_x.txt - the abscissa file.\n" );
  printf ( "    filename_r.txt - the region file.\n" );
/*
  Initialize parameters;
*/
  alpha = 0.0;
  beta = 0.0;
/*
  Get ORDER.
*/
  if ( 1 < argc )
  {
    order = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of ORDER (1 or greater)\n" );
    scanf ( "%d", &order );
  }
/*
  Get A.
*/
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of A:\n" );
    scanf ( "%lf", &a );
  }
/*
  Get B.
*/
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of B:\n" );
    scanf ( "%lf", &b );
  }
/*
  Get FILENAME:
*/
  if ( 4 < argc )
  {
    strcpy ( filename, argv[4] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter FILENAME, the \"root name\" of the quadrature files).\n" );
    scanf ( "%s", filename );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  ORDER = %d\n", order );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  FILENAME = \"%s\".\n", filename );
/*
  Construct the rule.
*/
  w = ( double * ) malloc ( order * sizeof ( double ) );
  x = ( double * ) malloc ( order * sizeof ( double ) );
  
  kind = 1;
  cgqf ( order, kind, alpha, beta, a, b, x, w );
/*
  Write the rule.
*/
  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  r[0] = a;
  r[1] = b;

  rule_write ( order, filename, x, w, r );
/*
  Free memory.
*/
  free ( r );
  free ( w );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LEGENDRE_RULE:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] )

/******************************************************************************/
/*
  Purpose:

    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.

  Discussion:

    This routine computes all the knots and weights of a Gauss quadrature
    formula with a classical weight function with default values for A and B,
    and only simple knots.

    There are no moments checks and no printing is done.

    Use routine EIQFS to evaluate a quadrature computed by CGQFS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  double *aj;
  double *bj;
  double zemu;

  parchk ( kind, 2 * nt, alpha, beta );
/*
  Get the Jacobi matrix and zero-th moment.
*/
  aj = ( double * ) malloc ( nt * sizeof ( double ) );
  bj = ( double * ) malloc ( nt * sizeof ( double ) );

  zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );
/*
  Compute the knots and weights.
*/
  sgqf ( nt, aj, bj, zemu, t, wts );

  free ( aj );
  free ( bj );

  return;
}
/******************************************************************************/

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] )

/******************************************************************************/
/*
  Purpose:

    CGQF computes knots and weights of a Gauss quadrature formula.

  Discussion:

    The user may specify the interval (A,B).

    Only simple knots are produced.

    Use routine EIQFS to evaluate this quadrature formula.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 September 2013

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Input, double A, B, the interval endpoints.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  int i;
  int *mlt;
  int *ndx;
/*
  Compute the Gauss quadrature formula for default values of A and B.
*/
  cdgqf ( nt, kind, alpha, beta, t, wts );
/*
  Prepare to scale the quadrature formula to other weight function with 
  valid A and B.
*/
  mlt = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 1;
  }
  ndx = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    ndx[i] = i + 1;
  }
  scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b );

  free ( mlt );
  free ( ndx );

  return;
}
/******************************************************************************/

double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] )

/******************************************************************************/
/*
  Purpose:

    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.

  Discussion:

    This routine computes the diagonal AJ and sub-diagonal BJ
    elements of the order M tridiagonal symmetric Jacobi matrix
    associated with the polynomials orthogonal with respect to
    the weight function specified by KIND.

    For weight functions 1-7, M elements are defined in BJ even
    though only M-1 are needed.  For weight function 8, BJ(M) is
    set to zero.

    The zero-th moment of the weight function is returned in ZEMU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, int M, the order of the Jacobi matrix.

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
    of the Jacobi matrix.

    Output, double CLASS_MATRIX, the zero-th moment.
*/
{
  double a2b2;
  double ab;
  double aba;
  double abi;
  double abj;
  double abti;
  double apone;
  int i;
  const double pi = 3.14159265358979323846264338327950;
  double temp;
  double temp2;
  double zemu;

  temp = r8_epsilon ( );

  parchk ( kind, 2 * m - 1, alpha, beta );

  temp2 = 0.5;

  if ( 500.0 * temp < fabs ( pow ( r8_gamma ( temp2 ), 2 ) - pi ) )
  {
    printf ( "\n" );
    printf ( "CLASS_MATRIX - Fatal error!\n" );
    printf ( "  Gamma function does not match machine parameters.\n" );
    exit ( 1 );
  }

  if ( kind == 1 )
  {
    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 2 )
  {
    zemu = pi;

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] =  sqrt ( 0.5 );
    for ( i = 1; i < m; i++ )
    {
      bj[i] = 0.5;
    }
  }
  else if ( kind == 3 )
  {
    ab = alpha * 2.0;
    zemu = pow ( 2.0, ab + 1.0 ) * pow ( r8_gamma ( alpha + 1.0 ), 2 )
      / r8_gamma ( ab + 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
    for ( i = 2; i <= m; i++ )
    {
      bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * pow ( i + alpha, 2 ) - 1.0 ) );
    }
  }
  else if ( kind == 4 )
  {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow ( 2.0, ab + 1.0 ) * r8_gamma ( alpha + 1.0 ) 
      * r8_gamma ( beta + 1.0 ) / r8_gamma ( abi );
    aj[0] = ( beta - alpha ) / abi;
    bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
      / ( ( abi + 1.0 ) * abi * abi ) );
    a2b2 = beta * beta - alpha * alpha;

    for ( i = 2; i <= m; i++ )
    {
      abi = 2.0 * i + ab;
      aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi * abi;
      bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
        / ( ( abi - 1.0 ) * abi ) );
    }
  }
  else if ( kind == 5 )
  {
    zemu = r8_gamma ( alpha + 1.0 );

    for ( i = 1; i <= m; i++ )
    {
      aj[i-1] = 2.0 * i - 1.0 + alpha;
      bj[i-1] = sqrt ( i * ( i + alpha ) );
    }
  }
  else if ( kind == 6 )
  {
    zemu = r8_gamma ( ( alpha + 1.0 ) / 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
    }
  }
  else if ( kind == 7 )
  {
    ab = alpha;
    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 8 )
  {
    ab = alpha + beta;
    zemu = r8_gamma ( alpha + 1.0 ) * r8_gamma ( - ( ab + 1.0 ) ) 
      / r8_gamma ( - beta );
    apone = alpha + 1.0;
    aba = ab * apone;
    aj[0] = - apone / ( ab + 2.0 );
    bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for ( i = 2; i <= m; i++ )
    {
      abti = ab + 2.0 * i;
      aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
    }

    for ( i = 2; i <= m - 1; i++ )
    {
      abti = ab + 2.0 * i;
      bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
        / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
    }
    bj[m-1] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      bj[i] =  sqrt ( bj[i] );
    }
  }
  else
  {
    printf ( "\n" );
    printf ( "CLASS_MATRIX - Fatal error!\n" );
    printf ( "  Illegal value of KIND = %d.\n", kind );
    exit ( 1 );
  }

  return zemu;
}
/******************************************************************************/

void imtqlx ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to 
    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 

    The authors thank the authors of EISPACK for permission to use this
    routine. 

    It has been modified to produce the product Q' * Z, where Z is an input 
    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
    The changes consist (essentially) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the 
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        printf ( "\n" );
        printf ( "IMTQLX - Fatal error!\n" );
        printf ( "  Iteration limit exceeded\n" );
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
/*
  Sorting.
*/
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
/******************************************************************************/

void parchk ( int kind, int m, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    PARCHK checks parameters ALPHA and BETA for classical weight functions. 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, int M, the order of the highest moment to
    be calculated.  This value is only needed when KIND = 8.

    Input, double ALPHA, BETA, the parameters, if required
    by the value of KIND.
*/
{
  double tmp;

  if ( kind <= 0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  KIND <= 0.\n" );
    exit ( 1 );
  }
/*
  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
*/
  if ( 3 <= kind && alpha <= -1.0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  3 <= KIND and ALPHA <= -1.\n" );
    exit ( 1 );
  }
/*
  Check BETA for Jacobi.
*/
  if ( kind == 4 && beta <= -1.0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  KIND == 4 and BETA <= -1.0.\n" );
    exit ( 1 );
  }
/*
  Check ALPHA and BETA for rational.
*/
  if ( kind == 8 )
  {
    tmp = alpha + beta + m + 1.0;
    if ( 0.0 <= tmp || tmp <= beta )
    {
      printf ( "\n" );
      printf ( "PARCHK - Fatal error!\n" );
      printf ( "  KIND == 8 but condition on ALPHA and BETA fails.\n" );
      exit ( 1 );
    }
  }
  return;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  const double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

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

    Input, double TABLE[M*N], the data.
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
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
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

void rule_write ( int order, char *filename, double x[], double w[],
  double r[] )

/******************************************************************************/
/*
  Purpose:

    RULE_WRITE writes a quadrature rule to three files.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double A, the left endpoint.

    Input, double B, the right endpoint.

    Input, char *FILENAME, specifies the output filenames.
    "filename_w.txt", "filename_x.txt", "filename_r.txt"
    defining weights, abscissas, and region.
*/
{
  char filename_r[80];
  char filename_w[80];
  char filename_x[80];

  strcpy ( filename_r, filename );
  strcat ( filename_r, "_r.txt" );
  strcpy ( filename_w, filename );
  strcat ( filename_w, "_w.txt" );
  strcpy ( filename_x, filename );
  strcat ( filename_x, "_x.txt" );

  printf ( "\n" );
  printf ( "  Creating quadrature files.\n" );
  printf ( "\n" );
  printf ( "  Root file name is     \"%s\".\n", filename );
  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", filename_w );
  printf ( "  Abscissa file will be \"%s\".\n", filename_x );
  printf ( "  Region file will be   \"%s\".\n", filename_r );

  r8mat_write ( filename_w, 1, order, w );
  r8mat_write ( filename_x, 1, order, x );
  r8mat_write ( filename_r, 1, 2,     r );

  return;
}
/******************************************************************************/

void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b )

/******************************************************************************/
/*
  Purpose:

    SCQF scales a quadrature formula to a nonstandard interval.

  Discussion:

    The arrays WTS and SWTS may coincide.

    The arrays T and ST may coincide.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, double T[NT], the original knots.

    Input, int MLT[NT], the multiplicity of the knots.

    Input, double WTS[NWTS], the weights.

    Input, int NWTS, the number of weights.

    Input, int NDX[NT], used to index the array WTS.  
    For more details see the comments in CAWIQ.

    Output, double SWTS[NWTS], the scaled weights.

    Output, double ST[NT], the scaled knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Input, double A, B, the interval endpoints.
*/
{
  double al;
  double be;
  int i;
  int k;
  int l;
  double p;
  double shft;
  double slp;
  double temp;
  double tmp;

  temp = r8_epsilon ( );

  parchk ( kind, 1, alpha, beta );

  if ( kind == 1 )
  {
    al = 0.0;
    be = 0.0;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 2 )
  {
    al = -0.5;
    be = -0.5;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 3 )
  {
    al = alpha;
    be = alpha;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 4 )
  {
    al = alpha;
    be = beta;

    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 5 )
  {
    if ( b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B <= 0\n" );
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 6 )
  {
    if ( b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B <= 0.\n" );
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / sqrt ( b );
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 7 )
  {
    al = alpha;
    be = 0.0;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 8 )
  {
    if ( a + b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  A + B <= 0.\n" );
      exit ( 1 );
    }
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  }

  p = pow ( slp, al + be + 1.0 );

  for ( k = 0; k < nt; k++ )
  {
    st[k] = shft + slp * t[k];
    l = abs ( ndx[k] );

    if ( l != 0 )
    {
      tmp = p;
      for ( i = l - 1; i <= l - 1 + mlt[k] - 1; i++ )
      {
        swts[i] = wts[i] * tmp;
        tmp = tmp * slp;
      }
    }
  }
  return;
}
/******************************************************************************/

void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] )

/******************************************************************************/
/*
  Purpose:

    SGQF computes knots and weights of a Gauss Quadrature formula.

  Discussion:

    This routine computes all the knots and weights of a Gauss quadrature
    formula with simple knots from the Jacobi matrix and the zero-th
    moment of the weight function, using the Golub-Welsch technique.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, double AJ[NT], the diagonal of the Jacobi matrix.

    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.

    Input, double ZEMU, the zero-th moment of the weight function.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  int i;
/*
  Exit if the zero-th moment is not positive.
*/
  if ( zemu <= 0.0 )
  {
    printf ( "\n" );
    printf ( "SGQF - Fatal error!\n" );
    printf ( "  ZEMU <= 0.\n" );
    exit ( 1 );
  }
/*
  Set up vectors for IMTQLX.
*/
  for ( i = 0; i < nt; i++ )
  {
    t[i] = aj[i];
  }
  wts[0] = sqrt ( zemu );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
/*
  Diagonalize the Jacobi matrix.
*/
  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i];
  }

  return;
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
