# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cell.h"

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

int i4vec_max ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MAX returns the value of the maximum element in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int IVEC_MAX, the value of the maximum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

double r8cvv_iget ( int mn, double a[], int m, int roff[], int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8CVV_IGET gets item J from row I in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row of the item.
    0 <= I < M.

    Input, int J, the column of the item.
    0 <= J < NR[I].

    Output, double R8CVV_IGET, the value of item A(I,J).
*/
{
  double aij;
  int k;

  k = roff[i] + j;
  aij = a[k];

  return aij;
}
/******************************************************************************/

void r8cvv_iinc ( int mn, double a[], int m, int roff[], int i, int j, 
  double daij )

/******************************************************************************/
/*
  Purpose:

    R8CVV_IINC increments item J from row I in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A(MN), the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row of the item.
    0 <= I < M.

    Input, int J, the column of the item.
    0 <= J < NR(I).

    Input, double DAIJ, the increment to the value of item A(I,J).
*/
{
  int k;

  k = roff[i] + j;
  a[k] = a[k] + daij;

  return;
}
/******************************************************************************/

void r8cvv_iset ( int mn, double a[], int m, int roff[], int i, int j, 
  double aij )

/******************************************************************************/
/*
  Purpose:

    R8CVV_ISET sets item J from row I in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A(MN), the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row of the item.
    0 <= I < M.

    Input, int J, the column of the item.
    0 <= J < NR[I].

    Input, double AIJ, the new value of item A(I,J).
*/
{
  int k;

  k = roff[i] + j;
  a[k] = aij;

  return;
}
/******************************************************************************/

double *r8cvv_nget_new ( int mn, double a[], int m, int roff[], int nn, 
  int in[], int jn[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_NGET_NEW gets N items JN(*) from row IN(*) in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input, double A(MN), the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int NN, the number of items.

    Input, int IN[NN], the rows of the items.
    0 <= IN(*) < M.

    Input, int JN[NN], the columns of the items.
    0 <= JN(*) < NR(IN(*)).

    Output, double R8CVV_NGET[NN], the value of items A(IN(*),JN(*)).
*/
{
  int i;
  int k;
  double *vn;

  vn = ( double * ) malloc ( nn * sizeof ( double ) );

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    vn[i] = a[k];
  }
  return vn;
}
/******************************************************************************/

void r8cvv_ninc ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double dvn[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_NINC increments items JN(*) from row IN(*) in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int NN, the number of items.

    Input, int IN[NN], the rows of the items.
    0 <= IN(*) < M.

    Input, int JN[NN], the columns of the items.
    0 <= JN(*) < NR(IN(*)).

    Input, double DVN[NN], the increments of items A(IN(*),JN(*)).
*/
{
  int i;
  int k;

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    a[k] = a[k] + dvn[i];
  }
  return;
}
/******************************************************************************/

void r8cvv_nset ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double vn[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_NSET sets items JN(*) from row IN(*) in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int NN, the number of items.

    Input, int IN[NN], the rows of the items.
    0 <= IN(*) < M.

    Input, int JN[NN], the columns of the items.
    0 <= JN(*) < NR(IN(*)).

    Input, double VN[NN], the new value of items A(IN(*),JN(*)).
*/
{
  int i;
  int k;

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    a[k] = vn[i];
  }
  return;
}
/******************************************************************************/

int *r8cvv_offset ( int m, int nr[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_OFFSET determines the row offsets of an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the array.

    Input, int NR(M), the row sizes.

    Output, int R8CVV_OFFSET[M+1], the row offsets.
*/
{
  int i;
  int *roff;

  roff = ( int * ) malloc ( ( m + 1 ) * sizeof ( int ) );

  roff[0] = 0;
  for ( i = 0; i < m; i++ )
  {
    roff[i+1] = roff[i] + nr[i];
  }

  return roff;
}
/******************************************************************************/

void r8cvv_print ( int mn, double a[], int m, int roff[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8CVV_PRINT prints an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int k;
  int k1;
  int k2;
  int khi;
  int klo;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    k1 = roff[i];
    k2 = roff[i+1];

    for ( klo = k1; klo < k2; klo = klo + 5 )
    {
      khi = i4_min ( klo + 5, k2 );
      if ( klo == k1 )
      {
        printf ( "%5d", i );
      }
      else
      {
        printf ( "     " );
      }
      printf ( "  " );
      for ( k = klo; k < khi; k++ )
      {
        printf ( "%14g", a[k] );
      }
      printf ( "\n" );
    }
  }
  return;
}
/******************************************************************************/

double *r8cvv_rget_new ( int mn, double a[], int m, int roff[], int i )

/******************************************************************************/
/*
  Purpose:

    R8CVV_RGET_NEW gets row I from an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row.
    0 <= I < M.

    Output, double R8CVV_GET[NR[I]], the value of A(I,*).
*/
{
  double *ai;
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  ai = ( double * ) malloc ( nv * sizeof ( double ) );
  for ( j = 0; j < nv; j++ )
  {
    ai[j] = a[k1+j];
  }

  return ai;
}
/******************************************************************************/

void r8cvv_rinc ( int mn, double a[], int m, int roff[], int i, double dai[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_RINC increments row I in an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row.
    0 <= I < M.

    Input, double DAI[NR[I]], the increment for A(I,*).
*/
{
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  for ( j = 0; j < nv; j++ )
  {
    a[k1+j] = a[k1+j] + dai[j];
  }

  return;
}
/******************************************************************************/

void r8cvv_rset ( int mn, double a[], int m, int roff[], int i, double ai[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_RSET sets row I from an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int MN, the size of the cell array.

    Input/output, double A[MN], the cell array.

    Input, int M, the number of rows in the array.

    Input, int ROFF[M+1], the row offsets.

    Input, int I, the row.
    0 <= I < M.

    Input, double AI[NR[I]], the new value of A(I,*).
*/
{
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  for ( j = 0; j < nv; j++ )
  {
    a[k1+j] = ai[j];
  }

  return;
}
/******************************************************************************/

int r8cvv_size ( int m, int nr[] )

/******************************************************************************/
/*
  Purpose:

    R8CVV_SIZE determines the size of an R8CVV.

  Discussion:

    An R8CVV is a "vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the array.

    Input, int NR[M], the size of each row.

    Output, int R8CVV_SIZE, the size of the cell array.
*/
{
  int i;
  int mn;

  mn = 0;
  for ( i = 0; i < m; i++ )
  {
    mn = mn + nr[i];
  }
  return mn;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
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
  printf ( "\n" );

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
      printf ( "  %12f", a[i] );
    }
    printf ( "\n" );
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

