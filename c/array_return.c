# include <stdlib.h>
# include <stdio.h>

int main ( );
void make_arrays ( int m, int **a, int n, int **b );
void i4vec_print ( int n, int a[], char *title );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ARRAY_RETURN.

  Discussion:

    The correct form of this program was worked out with the somewhat
    bemused assistance of Miro Stoyanov.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 February 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int *b;
  int m;
  int n;

  printf ( "\n" );
  printf ( "ARRAY_RETURN:\n" );
  printf ( "  C version\n" );
  printf ( "  Create two arrays in a function, \n" );
  printf ( "  return them in the argument list.\n" );
/*
  Specify the size of the arrays to be created.
*/
  m = 10;
  n = 5;

  make_arrays ( m, &a, n, &b );
/*
  Verify that the arrays were created and transferred properly.
*/
  i4vec_print ( m, a, "  A as received by main:" );
  i4vec_print ( n, b, "  B as received by main:" );
/*
  Free memory.
*/
  free ( a );
  free ( b );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ARRAY_RETURN:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void make_arrays ( int m, int **a, int n, int **b )

/******************************************************************************/
/*
  Purpose:

    MAKE_ARRAYS creates, sets, and returns two arrays using the argument list.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 February 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the desired size of the first array.

    Output, int **A, the first array.

    Input, int N, the desired size of the second array.

    Output, int **B, the second array.
*/
{
  int i;
/*
  We create PA and PB simply as a convenience.
  They make the code a little more readable.
*/
  int *pa;
  int *pb;
  
  ( *a ) = ( int * ) malloc ( m * sizeof ( int ) );
  pa = *a;

  for ( i = 0; i < m; i++ )
  {
    pa[i] = 10 + i;
  }
  i4vec_print ( m, pa, "  A as defined in MAKE_ARRAYS:" );

  ( *b ) = ( int * ) malloc ( n * sizeof ( int ) );
  pb = *b;

  for ( i = 0; i < n; i++ )
  {
    pb[i] = 100 + 2 * i;
  }
  i4vec_print ( n, pb, "  B as defined in MAKE_ARRAYS:" );

  return;
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
