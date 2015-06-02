# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "sphere_stereograph.h"

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

void plane_normal_basis_3d ( double pp[3], double pn[3], double pq[3], 
  double pr[3] )

/******************************************************************************/
/*
  Purpose:

    PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.

  Discussion:

    The normal form of a plane in 3D is:

      PP is a point on the plane,
      N is a normal vector to the plane.

    The two vectors to be computed, PQ and PR, can be regarded as
    the basis of a Cartesian coordinate system for points in the plane.
    Any point in the plane can be described in terms of the "origin" 
    point PP plus a weighted sum of the two vectors PQ and PR:

      P = PP + a * PQ + b * PR.

    The vectors PQ and PR have unit length, and are perpendicular to N
    and to each other.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double PP[3], a point on the plane.

    Input, double PN[3], a normal vector to the plane.  The
    vector must not have zero length, but it is not necessary for PN
    to have unit length.

    Output, double PQ[3], a vector of unit length, perpendicular
    to the vector PN and the vector PR.

    Output, double PR[3], a vector of unit length, perpendicular
    to the vector PN and the vector PQ.
*/
{
# define DIM_NUM 3

  int i;
  double normal_norm;
  double pr_norm;
  double *temp;
/*
  Compute the length of NORMAL.
*/
  normal_norm = r8vec_norm ( DIM_NUM, pn );

  if ( normal_norm == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PLANE_NORMAL_BASIS_3D - Fatal error!\n" );
    fprintf ( stderr, "  The normal vector is 0.\n" );
    exit ( 1 );
  }
/*
  Find a vector PQ that is normal to PN and has unit length.
*/
  temp = r8vec_any_normal ( DIM_NUM, pn );
  r8vec_copy ( DIM_NUM, temp, pq );
  free ( temp );
/*
  Now just take the cross product PN x PQ to get the PR vector.
*/
  temp = r8vec_cross_product_3d ( pn, pq );

  pr_norm = r8vec_norm ( DIM_NUM, temp );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    pr[i] = temp[i] / pr_norm;
  }
  free ( temp );

  return;
# undef DIM_NUM
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
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

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14g", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

double *r8vec_any_normal ( int dim_num, double v1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ANY_NORMAL returns some normal vector to V1.

  Discussion:

    An R8VEC is a vector of R8's.

    If DIM_NUM < 2, then no normal vector can be returned.

    If V1 is the zero vector, then any unit vector will do.

    No doubt, there are better, more robust algorithms.  But I will take
    just about ANY reasonable unit vector that is normal to V1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double V1[DIM_NUM], the vector.

    Output, double R8VEC_ANY_NORMAL[DIM_NUM], a vector that is
    normal to V2, and has unit Euclidean length.
*/
{
  int i;
  int j;
  int k;
  double *v2;
  double vj;
  double vk;

  if ( dim_num < 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_ANY_NORMAL - Fatal error!\n" );
    fprintf ( stderr, "  Called with DIM_NUM < 2.\n" );
    exit ( 1 );
  }

  v2 = ( double * ) malloc ( dim_num * sizeof ( double ) );

  if ( r8vec_norm ( dim_num, v1 ) == 0.0 )
  {
    r8vec_zero ( dim_num, v2 );
    v2[0] = 1.0;
    return v2;
  }
/*
  Seek the largest entry in V1, VJ = V1(J), and the
  second largest, VK = V1(K).

  Since V1 does not have zero norm, we are guaranteed that
  VJ, at least, is not zero.
*/
  j = -1;
  vj = 0.0;

  k = -1;
  vk = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( r8_abs ( vk ) < r8_abs ( v1[i] ) || k == -1 )
    {
      if ( r8_abs ( vj ) < r8_abs ( v1[i] ) || j == -1 )
      {
        k = j;
        vk = vj;
        j = i;
        vj = v1[i];
      }
      else
      {
        k = i;
        vk = v1[i];
      }
    }
  }
/*
  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
  will just about do the trick.
*/
  r8vec_zero ( dim_num, v2 );

  v2[j] = -vk / sqrt ( vk * vk + vj * vj );
  v2[k] =  vj / sqrt ( vk * vk + vj * vj );

  return v2;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

double *r8vec_cross_product_3d ( double v1[3], double v2[3] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V1[3], V2[3], the coordinates of the vectors.

    Output, double R8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
*/
{
  double *v3;

  v3 = ( double * ) malloc ( 3 * sizeof ( double ) );

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
/******************************************************************************/

double r8vec_norm ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM returns the L2 norm of an R8VEC.

  Discussion:

    The vector L2 norm is defined as:

      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_NORM, the L2 norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
/******************************************************************************/

double r8vec_norm_affine ( int n, double v0[], double v1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.

  Discussion:

    The affine vector L2 norm is defined as:

      R8VEC_NORM_AFFINE(V0,V1)
        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double V0[N], the base vector.

    Input, double V1[N], the vector whose affine L2 norm is desired.

    Output, double R8VEC_NORM_AFFINE, the affine L2 norm of V1.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

void r8vec_normal_01 ( int n, int *seed, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    This routine can generate a vector of values on one call.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values desired.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double X[N], a sample of the standard normal PDF.

  Local parameters:

    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.
*/
{
  int i;
  int m;
  const double pi = 3.141592653589793;
  double *r;
  int x_hi;
  int x_lo;
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
  }

  return;
}
/******************************************************************************/

void r8vec_transpose_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".

  Discussion:

    An R8VEC is a vector of R8's.

  Example:

    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
    TITLE = 'My vector:  '

    My vector:
        1.0    2.1    3.2    4.3    5.4
        6.5    7.6    8.7    9.8   10.9
       11.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;
  int ihi;
  int ilo;

  printf ( "\n" );
  printf ( "%s\n", title );

  if ( n <= 0 )
  {
    printf ( "  (Empty)\n" );
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      printf ( "  %12g", a[i] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void r8vec_uniform_01 ( int n, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r8vec_zero ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ZERO zeroes an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double A[N], a vector of zeroes.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
/******************************************************************************/

double *sphere_stereograph ( int m, int n, double p[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.

  Discussion:

    We start with a sphere of radius 1 and center (0,0,0).

    The north pole N = (0,0,1) is the point of tangency to the sphere
    of a plane, and the south pole S = (0,0,-1) is the focus for the
    stereographic projection.

    For any point P on the sphere, the stereographic projection Q of the
    point is defined by drawing the line from S through P, and computing
    Q as the intersection of this line with the plane.

    Actually, we allow the spatial dimension M to be arbitrary.  Values
    of M make sense starting with 2.  The north and south poles are
    selected as the points (0,0,...,+1) and (0,0,...,-1).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    John Burkardt

  Reference:

    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double P[M*N], a set of points on the unit sphere.

    Output, double SPHERE_STEREOGRAPH[M*N], the coordinates of the
    image points.
*/
{
  int i;
  int j;
  double *q;

  q = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      q[i+j*m] = 2.0 * p[i+j*m] / ( 1.0 + p[m-1+j*m] );
    }
    q[m-1+j*m] = 1.0;
  }

  return q;
}
/******************************************************************************/

double *sphere_stereograph_inverse ( int m, int n, double q[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.

  Discussion:

    We start with a sphere of radius 1 and center (0,0,0).

    The north pole N = (0,0,1) is the point of tangency to the sphere
    of a plane, and the south pole S = (0,0,-1) is the focus for the
    stereographic projection.

    For any point Q on the plane, the stereographic inverse projection
    P of the point is defined by drawing the line from S through Q, and
    computing P as the intersection of this line with the sphere.

    Actually, we allow the spatial dimension M to be arbitrary.  Values
    of M make sense starting with 2.  The north and south poles are
    selected as the points (0,0,...,+1) and (0,0,...,-1).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    John Burkardt

  Reference:

    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double Q[M*N], the points, which are presumed to lie
    on the plane Z = 1.

    Output, double SPHERE_STEREOGRAPH_INVERSE[M*N], the stereographic
    inverse projections of the points.
*/
{
  int i;
  int j;
  double *p;
  double qn;

  p = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    qn = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      qn = qn + pow ( q[i+j*m], 2 );
    }
    for ( i = 0; i < m - 1; i++ )
    {
      p[i+j*m] = 4.0 * q[i+j*m] / ( 4.0 + qn );
    }
    p[m-1+j*m] = ( 4.0 - qn ) / ( 4.0 + qn );
  }

  return p;
}
/******************************************************************************/

double *sphere_stereograph2 ( int m, int n, double p[], double focus[],
  double center[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.

  Discussion:

    We start with a sphere of center C.

    F is a point on the sphere which is the focus of the mapping,
    and the antipodal point 2*C-F is the point of tangency
    to the sphere of a plane.

    For any point P on the sphere, the stereographic projection Q of the
    point is defined by drawing the line from F through P, and computing
    Q as the intersection of this line with the plane.

    The spatial dimension M is arbitrary, but should be at least 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    John Burkardt

  Reference:

    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double P[M*N], a set of points on the unit sphere.

    Input, double FOCUS[M], the coordinates of the focus point.

    Input, double CENTER[M], the coordinates of the center of the sphere.

    Output, double SPHERE_STEREOGRAPH2[M*N], the coordinates of the
    image points,
*/
{
  double cf_dot_pf;
  double cf_normsq;
  int i;
  int j;
  double *q;
  double s;

  q = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    cf_normsq = 0.0;
    cf_dot_pf = 0.0;
    for ( i = 0; i < m; i++ )
    {
      cf_normsq = cf_normsq + pow ( center[i] - focus[i], 2 );
      cf_dot_pf = cf_dot_pf + ( center[i] - focus[i] ) * ( p[i+j*m] - focus[i] );
    }
    s = 2.0 * cf_normsq / cf_dot_pf;
    for ( i = 0; i < m; i++ )
    {
      q[i+j*m] = s * p[i+j*m] + ( 1.0 - s ) * focus[i];
    }
  }
  return q;
}
/******************************************************************************/

double *sphere_stereograph2_inverse ( int m, int n, double q[], double focus[],
  double center[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.

  Discussion:

    We start with a sphere of center C.

    F is a point on the sphere which is the focus of the mapping,
    and the antipodal point 2*C-F is the point of tangency
    to the sphere of a plane.

    For any point Q on the plane, the stereographic inverse projection
    P of the point is defined by drawing the line from F through Q, and
    computing P as the intersection of this line with the sphere.

    The spatial dimension M is arbitrary, but should be at least 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    John Burkardt

  Reference:

    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double Q[M*N], the points, which are presumed to lie
    on the plane.

    Input, double FOCUS[M], the coordinates of the focus point.

    Input, double CENTER[M], the coordinates of the center of the sphere.

    Output, double SPHERE_STEREOGRAPH2_INVERSE[M*N], the stereographic
    inverse projections of the points.
*/
{
  double cf_dot_qf;
  int i;
  int j;
  double *p;
  double qf_normsq;
  double s;

  p = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    cf_dot_qf = 0.0;
    qf_normsq = 0.0;
    for ( i = 0; i < m; i++ )
    {
      cf_dot_qf = cf_dot_qf + ( center[i] - focus[i] ) * ( q[i+j*m] - focus[i] );
      qf_normsq = qf_normsq + pow ( q[i+j*m] - focus[i], 2 );
    }
    s = 2.0 * cf_dot_qf / qf_normsq;
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = s * q[i+j*m] + ( 1.0 - s ) * focus[i];
    }
  }
  return p;
}
/******************************************************************************/

void timestamp ( )

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
/******************************************************************************/

double *uniform_on_sphere01_map ( int dim_num, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.

  Discussion:

    The sphere has center 0 and radius 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.

    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.

  Parameters:

    Input, int DIM_NUM, the dimension of the space.

    Input, int N, the number of points.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double UNIFORM_ON_SPHERE01_MAP[DIM_NUM*N], the points.
*/
{
  int i;
  int j;
  double norm;
  double *u;
  double *x;

  u = ( double * ) malloc ( dim_num * sizeof ( double ) );
  x = ( double * ) malloc ( dim_num * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
/*
  Fill a vector with normally distributed values.
*/
    r8vec_normal_01 ( dim_num, seed, u );
/*
  Compute the length of the vector.
*/
    norm = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      norm = norm + u[i] * u[i];
    }
    norm = sqrt ( norm );
/*
  Normalize the vector.
*/
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = u[i] / norm;
    }

  }

  free ( u );

  return x;
}
