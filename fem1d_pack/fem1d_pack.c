# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "fem1d_pack.h"

/******************************************************************************/

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.

  Discussion:

    The quantity computed here is the "geometric" bandwidth determined
    by the finite element mesh alone.

    If a single finite element variable is associated with each node
    of the mesh, and if the nodes and variables are numbered in the
    same way, then the geometric bandwidth is the same as the bandwidth
    of a typical finite element matrix.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where 

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    Because the finite element node adjacency relationship is symmetric,
    we are guaranteed that ML = MU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

void legendre_com ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose: 

    LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.

  Integration interval:

    [ -1, 1 ]

  Weight function:

    1.

  Integral to approximate:

    Integral ( -1 <= X <= 1 ) F(X) dX.

  Approximate integral:

    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2005

  Author:

    John Burkardt

  Parameters:

    Input, int NORDER, the order of the rule.
    NORDER must be greater than 0.

    Output, double XTAB[NORDER], the abscissas of the rule.

    Output, double WEIGHT[NORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
# define PI 3.141592653589793

  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_COM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of NORDER = %d\n", order );
    exit ( 1 );
  }
 
  e1 = ( double ) ( order * ( order + 1 ) );
 
  m = ( order + 1 ) / 2;
 
  for ( i = 1; i <= ( order + 1 ) / 2; i++ )
  {
    mp1mi = m + 1 - i;
    t = PI * ( double ) ( 4 * i - 1 ) / ( double ) ( 4 * order + 2 );
    x0 = cos(t) * ( 1.0 - ( 1.0 - 1.0 / 
      ( double ) ( order ) ) / ( double ) ( 8 * order * order ) );
 
    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }
 
    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
/*
  Initial approximation H:
*/
    h = - u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn 
      / ( 3.0 * dpn ) ) ) );
/*
  Refine H using one step of Newton's method:
*/
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;
 
    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;
 
    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx ); 
  }
 
  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
/*
  Shift the data up.
*/
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
/*
  Reflect values for the negative abscissas.
*/
  for ( i = 0; i < order - nmove; i++ )
  {
    xtab[i] = - xtab[order-1-i];
    weight[i] = weight[order-1-i];
  }
 
  return;

# undef PI
}
/******************************************************************************/

double *local_basis_1d ( int order, double node_x[], double x )

/******************************************************************************/
/*
  Purpose:

    LOCAL_BASIS_1D evaluates the basis functions in an element.

  Discussion:

    PHI(I)(X) = product ( J ~= I ) ( X         - NODE_X(I) ) 
                                 / ( NODE_X(J) - NODE_X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the element.
    0 <= ORDER.  ORDER = 1 means piecewise linear.

    Input, double NODE_X[ORDER], the element nodes.  
    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
    and 0 when X is equal to any other node.

    Input, double X, the point at which the basis functions are to 
    be evaluated.

    Output, double LOCAL_BASIS_1D[ORDER], the basis functions.
*/
{
  int i;
  int j;
  double *phi;

  phi = ( double * ) malloc ( order * sizeof ( double ) );

  for ( j = 0; j < order; j++ )
  {
    phi[j] = 1.0;
    for ( i = 0; i < order; i++ )
    {
      if ( j != i )
      {
        phi[j] = ( phi[j] * ( x - node_x[i] ) ) / ( node_x[j] - node_x[i] );
      }
    }
  }

  return phi;
}
/******************************************************************************/

double *local_basis_prime_1d ( int order, double node_x[], double x )

/******************************************************************************/
/*
  Purpose:

    LOCAL_BASIS_PRIME_1D evaluates the basis function derivatives in an element.

  Discussion:

    PHI(I)(X) = product ( J ~= I ) ( X - NODE_X(I) ) 
                                 / ( NODE_X(J) - NODE_X(I) )

    dPHIdx(I)(X) = sum ( J ~= I ) ( 1 / ( NODE_X(J) - NODE_X(I) ) *
      product ( K ~= ( J, I ) ) ( X - NODE_X(I) ) / ( NODE_X(J) - NODE_X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the element.
    0 <= ORDER.  ORDER = 1 means piecewise linear.

    Input, double NODE_X[ORDER], the element nodes.  
    These must be distinct.  Basis function I is 1 when X = NODE_X(I) 
    and 0 when X is equal to any other node.

    Input, double X, the point at which the basis functions are to 
    be evaluated.

    Output, double LOCAL_BASIS_PRIME_1D[ORDER], the basis functions.
*/
{
  double *dphidx;
  int i;
  int j;
  int k;
  double term;

  dphidx = ( double * ) malloc ( order * sizeof ( double ) );

  for ( i = 0; i < order; i++ )
  {
    dphidx[i] = 0.0;
    for ( j = 0; j < order; j++ )
    {
      if ( j != i )
      {
        term = 1.0 / ( node_x[j] - node_x[i] );
        for ( k = 0; k < order; k++ )
        {
          if ( k != i && k != j )
          {
            term = term * ( x - node_x[i] ) / ( node_x[k] - node_x[i] );
          }
        }
        dphidx[i] = dphidx[i] + term;
      }
    }
  }
  return dphidx;
}
/******************************************************************************/

double *local_fem_1d ( int order, double node_x[], double node_v[], 
  int sample_num, double sample_x[] )

/******************************************************************************/
/*
  Purpose:

    LOCAL_FEM_1D evaluates a local finite element function.

  Discussion:

    A local finite element function is a finite element function
    defined over a single element.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the element.
    0 <= ORDER.  ORDER = 1 means piecewise linear.

    Input, double NODE_X[ORDER], the element nodes.  
    These must be distinct.  Basis function I is 1 when X = NODE_X(I) and 0 
    when X is equal to any other node.

    Input, double NODE_V[ORDER], the value of the finite element 
    function at each node.

    Input, int SAMPLE_NUM, the number of sample points.

    Input, double SAMPLE_X[SAMPLE_NUM], the sample points at which 
    the local finite element function is to be evaluated.

    Output, double LOCAL_FEM_1D[SAMPLE_NUM], the values of the local 
    finite element basis functions.
*/
{
  double *phi;
  int sample;
  double *sample_v;
  double x;

  sample_v = ( double * ) malloc ( sample_num * sizeof ( double ) );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    x = sample_x[sample];
    phi = local_basis_1d ( order, node_x, x );
    sample_v[sample] = r8vec_dot_product ( order, node_v, phi );
    free ( phi );
  }

  return sample_v;
}
/******************************************************************************/

double r8_uniform_ab ( double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_AB returns a pseudorandom R8 scaled to [A,B].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double R8_UNIFORM_AB, a number strictly between A and B.
*/
{
  int i4_huge = 2147483647;
  int k;
  double r;
  double value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  value = a + ( b - a ) * r;

  return value;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double r8vec_sum ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SUM returns the sum of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_SUM, the sum of the vector.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

