# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "dream_user.h"
# include "pdflib.h"
# include "rnglib.h"
# include "problem1_covariance.h"

int main ( );
void test01 ( int par_num, int sample_num );
double r8_huge ( void );
double *r8mat_covariance ( int m, int n, double x[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PROBLEM1_MAIN.

  Discussion:

    The coding of PROBLEM1 is tricky enough that I want to be able to
    try it out independently of the DREAM code.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int chain_num;
  int cr_num;
  int gen_num;
  int pair_num;
  int par_num;
  int sample_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "PROBLEM1_MAIN\n" );
  printf ( "  C version\n" );
/*
  Initialize the random number generator library.
*/
  initialize ( );
/*
  By calling PROBLEM_SIZE, we implicitly set up the covariance as well.
*/
  problem_size ( &chain_num, &cr_num, &gen_num, &pair_num, &par_num );

  sample_num = 10000;
  test01 ( par_num, sample_num );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PROBLEM1_MAIN\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int par_num, int sample_num )

/******************************************************************************/
/*
  Purpose:

    TEST01 calls the sampling function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  double *cov_sample;
  int i;
  int j;
  double z;
  double *zp;
  double *zp_array;
  double *zp_ave;
  double *zp_max;
  double *zp_min;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Call PRIOR_SAMPLE many times.\n" );
  printf ( "  Compare statistics to PDF parameters.\n" );
  printf ( "  Note that the covariance estimate can be very bad\n" );
  printf ( "  unless the matrix is strongly diagonal.\n" );
  printf ( "\n" );
  printf ( "  Parameter dimension is %d\n", par_num );
  printf ( "  Number of samples is %d\n", sample_num );
/*
  Compute N multinormal samples.
*/
  zp_array = ( double * ) malloc ( par_num * sample_num * sizeof ( double ) );
  for ( j = 0; j < sample_num; j++ )
  {
    zp = prior_sample ( par_num );
    for ( i = 0; i < par_num; i++ )
    {
      zp_array[i+j*par_num] = zp[i];
    }
    free ( zp );
  }

  zp_ave = ( double * ) malloc ( par_num * sizeof ( double ) );
  zp_max = ( double * ) malloc ( par_num * sizeof ( double ) );
  zp_min = ( double * ) malloc ( par_num * sizeof ( double ) );

  for ( i = 0; i < par_num; i++ )
  {
    zp_min[i] =   r8_huge ( );
    zp_max[i] = - r8_huge ( );
    zp_ave[i] = 0.0;
    for ( j = 0; j < sample_num; j++ )
    {
      z = zp_array[i+j*par_num];
      if ( z < zp_min[i] )
      {
        zp_min[i] = z;
      }
      if ( zp_max[i] < z )
      {
        zp_max[i] = z;
      }
      zp_ave[i] = zp_ave[i] + z;
    }
    zp_ave[i] = zp_ave[i] / ( double ) ( sample_num );
  }
  printf ( "\n" );
  printf ( " Index       Min            Ave              Max             MU\n" );
  printf ( "\n" );
  for ( i = 0; i < par_num; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g  %14g\n",
      i, zp_min[i], zp_ave[i], zp_max[i], zp_mean[i] );
  }

  cov_sample = r8mat_covariance ( par_num, sample_num, zp_array );

  r8mat_print ( par_num, par_num, cov_sample, "  Sample covariance:" );

  r8mat_print ( par_num, par_num, c, "  PDF covariance:" );

  free ( cov_sample );
  free ( zp_array );
  free ( zp_ave );
  free ( zp_max );
  free ( zp_min );

  return;
}
/******************************************************************************/

double r8_huge ( void )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

double *r8mat_covariance ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COVARIANCE computes the sample covariance of a set of vector data.

  Discussion:

    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt.

  Parameters:

    Input, int M, the size of a single data vectors.

    Input, int N, the number of data vectors.
    N should be greater than 1.

    Input, double X[M*N], an array of N data vectors, each
    of length M.

    Output, double C[M*M], the covariance matrix for the data.
*/
{
  double *c;
  int i;
  int j;
  int k;
  double *x_mean;

  c = ( double * ) malloc ( m * m * sizeof ( double ) );
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = 0.0;
    }
  }
/*
  Special case of N = 1.
*/
  if ( n == 1 )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+i*m] = 1.0;
    }
    return c;
  }
/*
  Determine the sample means.
*/
  x_mean = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    x_mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      x_mean[i] = x_mean[i] + x[i+j*m];
    }
    x_mean[i] = x_mean[i] / ( double ) ( n );
  }
/*
  Determine the sample covariance.
*/
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      for ( k = 0; k < n; k++ )
      {
        c[i+j*m] = c[i+j*m] 
          + ( x[i+k*m] - x_mean[i] ) * ( x[j+k*m] - x_mean[j] );
      }
    }
  }

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = c[i+j*m] / ( double ) ( n - 1 );
    }
  }

  free ( x_mean );

  return c;
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
