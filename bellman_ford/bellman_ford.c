# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "bellman_ford.h"

/******************************************************************************/

void bellman_ford ( int v_num, int e_num, int source, int e[],
  double e_weight[], double v_weight[], int predecessor[] )

/******************************************************************************/
/*
  Purpose:

    BELLMAN_FORD finds shortest paths from a given vertex of a weighted directed graph.

  Discussion:

    The Bellman-Ford algorithm is used.

    Each edge of the graph has a weight, which may be negative.  However,
    it should not be the case that there is any negative loop, that is,
    a circuit whose total weight is negative.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2014

  Author:

    John Burkardt

  Parameters:

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of edges.

    Input, int SOURCE, the vertex from which distances will be calculated.

    Input, int E[2*E_NUM], the edges, given as pairs of vertex indices.

    Input, double E_WEIGHT[E_NUM], the weight of each edge.

    Output, double V_WEIGHT[V_NUM], the weight of each node, that is,
    its minimum distance from SOURCE.

    Output, int PREDECESSOR[V_NUM], a list of predecessors, which can be
    used to recover the shortest path from any node back to SOURCE.
*/    
{
  int i;
  int j;
  const double r8_big = 1.0E+30;
  double t;
  int u;
  int v;
  double w;
/*
  Step 1: initialize the graph.
*/
  for ( i = 0; i < v_num; i++ )
  {
    if ( i == source )
    {
      v_weight[i] = 0.0;
    }
    else
    {
      v_weight[i] = r8_big;
    }
  }

  for ( i = 0; i < v_num; i++ )
  {
    predecessor[i] = -1;
  }
/*
  Step 2: Relax edges repeatedly.
*/
  for ( i = 1; i < v_num; i++ )
  {
    for ( j = 0; j < e_num; j++ )
    {
      u = e[1+j*2];
      v = e[0+j*2];
      t = v_weight[u] + e_weight[j];
      if ( t < v_weight[v] )
      {
        v_weight[v] = t;
        predecessor[v] = u;
      }
    }
  }
/*
  Step 3: check for negative-weight cycles
*/
  for ( j = 0; j < e_num; j++ )
  {
    u = e[1+j*2];
    v = e[0+j*2];
    if ( v_weight[u] + e_weight[j] < v_weight[v] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "BELLMAN_FORD - Fatal error!\n" );
      fprintf ( stderr, "  Graph contains a cycle with negative weight.\n" );
      exit ( 1 );
    }
  }

  return;
}
/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

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
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
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

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < jhi )
    {
      j2hi = n;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
  }

  return;
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

