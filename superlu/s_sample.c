# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "slu_sdefs.h"

int main ( );
float *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], float acc[], 
  float x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], float acc[], 
  char *title );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    S_SAMPLE tests the SUPERLU solver with a 5x5 single precision real matrix.

  Discussion:

    The general (GE) representation of the matrix is:

      [ 19  0 21 21  0
        12 21  0  0  0
         0 12 16  0  0 
         0  0  0  5 21
        12 12  0  0 18 ]

    The (0-based) compressed column (CC) representation of this matrix is:

      I  CC   A
     --  --  --
      0   0  19
      1      12
      4      12

      1   3  21
      2      12
      4      12

      0   6  21
      2      16

      0   8  21
      3       5

      3  10  21
      4      18

      *  12   *

    The right hand side B and solution X are

      #   B     X
     --  --  ----------
      0   1  -0.03125
      1   1   0.0654762
      2   1   0.0133929
      3   1   0.0625
      4   1   0.0327381 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2014

  Author:

    John Burkardt

  Reference:

    James Demmel, John Gilbert, Xiaoye Li,
    SuperLU Users's Guide.
*/
{
# define N 5
# define NCC 12

  SuperMatrix A;
  float acc[NCC] = { 
    19.0, 12.0, 12.0, 
    21.0, 12.0, 12.0,
    21.0, 16.0, 
    21.0,  5.0, 
    21.0, 18.0 };
  float *b;
  float *b2;
  SuperMatrix B;
  int ccc[N+1] = { 0, 3, 6, 8, 10, 12 };
  int i;
  int icc[NCC] = { 
    0, 1, 4,
    1, 2, 4,
    0, 2, 
    0, 3, 
    3, 4 };
  int info;
  int j;
  SuperMatrix L;
  int m = N;
  int n = N;
  int nrhs = 1;
  int ncc = NCC;
  superlu_options_t options;
  int *perm_c;
  int permc_spec;
  int *perm_r;
  SuperLUStat_t stat;
  SuperMatrix U;

  timestamp ( );
  printf ( "\n" );
  printf ( "S_SAMPLE:\n" );
  printf ( "  C version\n" );
  printf ( "  SUPERLU solves a single precision real linear system.\n" );
/*
  Print the matrix.
*/
  cc_print ( m, n, ncc, icc, ccc, acc, "  0-based CC matrix SAMPLE:" );
/*
  Convert the compressed column (CC) matrix into a SuperMatrix A. 
*/
  sCreate_CompCol_Matrix ( &A, m, n, ncc, acc, icc, ccc, SLU_NC, SLU_S, SLU_GE );    
/*
  Create the right-hand side matrix.
*/
  b = ( float * ) malloc ( m * sizeof ( float ) );
  for ( i = 0; i < m; i++ )
  {
    b[i] = 1.0;
  }
  printf ( "\n" );
  printf ( "  Right hand side:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "%g\n", b[i] );
  }
/*
  Create Super Right Hand Side.
*/
  sCreate_Dense_Matrix ( &B, m, nrhs, b, m, SLU_DN, SLU_S, SLU_GE );
/*
  Set space for the permutations.
*/
  perm_r = ( int * ) malloc ( m * sizeof ( int ) );
  perm_c = ( int * ) malloc ( n * sizeof ( int ) );
/*
  Set the input options. 
*/
  set_default_options ( &options );
  options.ColPerm = NATURAL;
/*
  Initialize the statistics variables. 
*/
  StatInit ( &stat );
/*
  Solve the linear system. 
*/
  sgssv ( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
    
  sPrint_CompCol_Matrix ( ( char * ) "A", &A );
  sPrint_CompCol_Matrix ( ( char * ) "U", &U );
  sPrint_SuperNode_Matrix ( ( char * ) "L", &L );
  print_int_vec ( ( char * ) "\nperm_r", m, perm_r );
/*
  By some miracle involving addresses, 
  the solution has been put into the B vector.
*/
  printf ( "\n" );
  printf ( "  Computed solution:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %g\n", b[i] );
  }
/*
  B now contains the solution X.
  Set B2 = A * X.
*/
  b2 = cc_mv ( m, n, ncc, icc, ccc, acc, b );
  printf ( "\n" );
  printf ( "  Product A*X:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %g\n", b2[i] );
  }
/*
  Free memory.
*/
  free ( b );
  free ( b2 );
  free ( perm_c );
  free ( perm_r );

  Destroy_SuperMatrix_Store ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperNode_Matrix ( &L );
  Destroy_CompCol_Matrix ( &U );
  StatFree ( &stat );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "S_SAMPLE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;

# undef N
# undef NCC
}
/******************************************************************************/

float *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], float acc[], 
  float x[] )

/******************************************************************************/
/*
  Purpose:

    CC_MV multiplies a CC matrix by a vector

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Reference:

    Iain Duff, Roger Grimes, John Lewis,
    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
    October 1992

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NCC, the number of CC values.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns

    Input, float ACC[NCC], the CC values.

    Input, float X[N], the vector to be multiplied.

    Output, float CC_MV[M], the product A*X.
*/
{
  float *b;
  int i;
  int j;
  int k;

  b = ( float * ) malloc ( m * sizeof ( float ) );

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = ccc[j]; k < ccc[j+1]; k++ )
    {
      i = icc[k];
      b[i] = b[i] + acc[k] * x[j];
    }
  }

  return b;
}
/******************************************************************************/

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], float acc[], 
  char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT prints a sparse matrix in CC format.

  Discussion:

    The index data in ICC and CCC is assumed to be 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the matrix.

    Input, int N, the number of columns in the matrix.

    Input, int NCC, the number of CC elements.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, float ACC[NCC], the CC values.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jnext;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J         A\n" );
  printf ( "  ----  ----  ----  ----------------\n" );
  printf ( "\n" );

  j = 0;
  jnext = ccc[1];

  for ( k = 0; k < ncc; k++ )
  {
    i = icc[k];
    while ( jnext <= k )
    {
      j = j + 1;
      jnext = ccc[j+1];
    }
 
    printf ( "  %4d  %4d  %4d  %16.8g\n", k, i, j, acc[k] );
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

