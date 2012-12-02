# include <stdlib.h>
# include <stdio.h>

int main ( );
void i4pp_delete ( int **a, int m, int n );
int **i4pp_new ( int m, int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DYNAMIC_ARRAY_2D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2011

  Author:

    John Burkardt
*/
{
  int **a;
  int **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  printf ( "\n" );
  printf ( "DYNAMIC_ARRAY_2D:\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrate the allocation and use of dynamic arrays.\n" );
/*
  Request the array.
*/
  m = 4;
  n = 5;

  a = i4pp_new ( m, n  );
/*
  Assign values.
*/
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 10 * i + j;
    }
  }
  printf ( "\n" );
  printf ( "  A = a dynamically allocated matrix:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %2d", a[i][j] );
    }
    printf ( "\n" );
  }
/*
  Multiply A' * A.
*/
  b = i4pp_new ( n, n );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
  printf ( "\n" );
  printf ( "  B = A' * A\n" );
  printf ( "  (also dynamically allocated)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %2d", b[i][j] );
    }
    printf ( "\n" );
  }
/*
  Free memory.
*/
  i4pp_delete ( a, m, n );
  i4pp_delete ( b, n, n );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DYNAMIC_ARRAY_2D:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void i4pp_delete ( int **a, int m, int n )

/******************************************************************************/
/*
  Purpose:

    I4PP_DELETE frees the memory set aside by I4PP_NEW.

  Discussion:

    This function releases the memory associated with an array that was 
    created by a command like
      int **a;
      a = i4pp_new ( m, n );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int **A, a pointer to the pointers to the array.

    Input, int M, N, the number of rows and columns.
*/
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    free ( a[i] );
  }

  free ( a );

  return;
}
/******************************************************************************/

int **i4pp_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    I4PP_NEW sets up an integer matrix of the desired dimensions.

  Discussion:

    A declaration of the form
      int **a;
    is necesary.  Then an assignment of the form:
      a = i4pp_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17;
      y = a[1][0];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, int **I4PP_NEW, a pointer to the pointers to the data.
*/
{
  int **a;
  int i;

  a = ( int ** ) malloc ( m * n * sizeof ( int * ) );

  for ( i = 0; i < m; i++ )
  {
    a[i] = ( int * ) malloc ( n * sizeof ( int ) );
  }
  return a;
}

