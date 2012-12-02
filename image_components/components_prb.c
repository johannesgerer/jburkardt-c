# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "components.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for COMPONENTS_PRB.

  Discussion:

    COMPONENTS_PRB tests the COMPONENTS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "COMPONENTS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the COMPONENTS library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "COMPONENTS_PRB\n" );
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

    TEST01 tests I4VEC_COMPONENTS on a simple case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
# define N 28

  int a[N] = {
    0, 0, 1, 2, 4, 0, 0, 4, 0, 0,
    0, 8, 9, 9, 1, 2, 3, 0, 0, 5,
    0, 1, 6, 0, 0, 0, 4, 0 };
  int c[N];
  int component_num;
  int j;
  int n = N;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  I4VEC_COMPONENTS finds and labels connected\n" );
  printf ( "  components in a 1D integer vector.\n" );

  printf ( "\n" );
  printf ( "  A:\n" );
  printf ( "\n" );
  printf ( "    " );
  for ( j = 0; j < n; j++ )
  {
    printf ( "%d", a[j] );
  }
  printf ( "\n" );

  component_num = i4vec_components ( n, a, c );

  printf ( "\n" );
  printf ( "  Number of components = %d\n", component_num );
  printf ( "\n" );
  printf ( "  C:\n" );
  printf ( "\n" );
  printf ( "    " );
  for ( j = 0; j < n; j++ )
  {
    printf ( "%d", c[j] );
  }
  printf ( "\n" );

  return;
# undef N
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests I4MAT_COMPONENTS on a simple case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
# define M 9
# define N 17

  int a[M*N] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 1, 0, 1, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 1, 1, 1, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 0, 1, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 1, 1, 0,
    0, 1, 0, 1, 1, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 1, 0, 1, 1, 0,
    0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int *c;
  int component_num;
  int i;
  int j;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  I4MAT_COMPONENTS finds and labels connected\n" );
  printf ( "  components in a 2D integer array.\n" );

  printf ( "\n" );
  printf ( "  A:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "    " );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%d", a[i+j*m] );
    }
    printf ( "\n" );
  }

  c = ( int * ) malloc ( m * n * sizeof ( int ) );

  component_num = i4mat_components ( m, n, a, c );

  printf ( "\n" );
  printf ( "  Number of components = %d\n", component_num );
  printf ( "\n" );
  printf ( "  C:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "    " );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%d", c[i+j*m] );
    }
    printf ( "\n" );
  }

  free ( c );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests I4BLOCK_COMPONENTS on a simple case.

  Discussion:

    This calculation is also done by a program called REGION.
    The two programs differ in the number of components discovered
    because REGION uses the full 3x3 block of pixels, resulting
    in 26 potential neighbors, whereas I4BLOCK_COMPONENTS uses only
    the north/south, east/west, up/down directions for 8 neighbors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
# define L 64
# define M 64
# define N 26

  int a[L*M*N];
  int c[L*M*N];
  int component_num;
  char filename[30] = "indices.txt";
  int i;
  int i1;
  int *indices;
  int j;
  int j1;
  int k;
  int l = L;
  int m = M;
  int m1;
  int n = N;
  int n1;
  int *s;
  int s_total;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  I4BLOCK_COMPONENTS finds and labels connected\n" );
  printf ( "  components in a 3D integer block.\n" );

  printf ( "\n" );
  printf ( "  A is a 3D block of order %d * %d * %d\n", l, m, n );

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        a[i+j*l+k*l*m] = 0;
      }
    }
  }
/*
  Retrieve the indices of nonzero data in A by reading a file.
*/
  i4mat_header_read ( filename, &m1, &n1 );

  indices = i4mat_data_read ( filename, m1, n1 );

  for ( j1 = 0; j1 < n1; j1++ )
  {
    i = indices[0+j1*3] - 1;
    j = indices[1+j1*3] - 1;
    k = indices[2+j1*3] - 1;
    a[i+j*l+k*l*m] = 1;
  }

  free ( indices );

  printf ( "\n" );
  printf ( "  Number of nonzero A values is %d\n", n1 );
/*
  Determine the components.
*/
  component_num = i4block_components ( l, m, n, a, c );

  printf ( "\n" );
  printf ( "  Number of components = %d\n", component_num );

  s = ( int * ) malloc ( component_num * sizeof ( int ) );

  for ( i = 0; i < component_num; i++ )
  {
    s[i] = 0;
  }

  printf ( "\n" );
  printf ( "  Number of components = %d\n", component_num );

  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        if ( c[i+j*l+k*l*m] != 0 )
        {
          s[c[i+j*l+k*l*m]-1] = s[c[i+j*l+k*l*m]-1] + 1;
        }
      }
    }
  }

  printf ( "\n" );
  printf ( "  Component  Size\n" );
  printf ( "\n" );
  s_total = 0;
  for ( i = 0; i < component_num; i++ )
  {
    printf ( "  %4d  %8d\n", i + 1, s[i] );
    s_total = s_total + s[i];
  }
  printf ( "------  --------\n" );
  printf ( " Total  %8d\n", s_total );

  free ( s );

  return;
# undef L
# undef M
# undef N
}
