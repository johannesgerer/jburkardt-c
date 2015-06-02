# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "slu_zdefs.h"

int main ( );
doublecomplex *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], 
  doublecomplex acc[], doublecomplex x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], doublecomplex acc[], 
  char *title );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    Z_SAMPLE_ST tests the SUPERLU solver with a 5x5 double precision complex matrix.

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

    17 July 2014

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
  doublecomplex *acc;
  SuperMatrix B;
  doublecomplex *b;
  doublecomplex *b2;
  int *ccc;
  int i;
  int *icc;
  int info;
  int j;
  SuperMatrix L;
  int m = 5;
  int n = N;
  int ncc = NCC;
  int nrhs = 1;
  superlu_options_t options;
  int *perm_c;
  int permc_spec;
  int *perm_r;
  SuperLUStat_t stat;
  SuperMatrix U;
  doublecomplex *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "Z_SAMPLE_ST:\n" );
  printf ( "  C version\n" );
  printf ( "  SUPERLU solves a double precision complex linear system.\n" );
  printf ( "  The matrix is read from a Sparse-Triplet (ST) file.\n" );
/*
  Read the matrix from a file associated with standard input,
  in Sparse Triplet format, into compressed column (CC) format.
*/
  zreadtriple ( &m, &n, &ncc, &acc, &icc, &ccc );
/*
  Print the matrix.
*/
  cc_print ( m, n, ncc, icc, ccc, acc, "  CC matrix:" );
/*
  Convert the compressed column (CC) matrix into a SuperMatrix A. 
*/
  zCreate_CompCol_Matrix ( &A, m, n, ncc, acc, icc, ccc, SLU_NC, SLU_Z, SLU_GE );    
/*
  Create right-hand side matrix.
*/
  b = ( doublecomplex * ) malloc ( m * sizeof ( doublecomplex ) );
  for ( i = 0; i < m; i++ )
  {
    b[i] = ( doublecomplex ) { 1.0, 0.0 };
  }
/*
  Create right hand side matrix B in the format expected by SuperLU. 
*/
  zCreate_Dense_Matrix ( &B, m, nrhs, b, m, SLU_DN, SLU_Z, SLU_GE );
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
  zgssv ( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
    
  zPrint_CompCol_Matrix ( ( char * ) "A", &A );
  zPrint_CompCol_Matrix ( ( char * ) "U", &U );
  zPrint_SuperNode_Matrix ( ( char * ) "L", &L );
  print_int_vec ( ( char * ) "\nperm_r", m, perm_r );
/*
  Get a pointer to the solution.
*/
  x = ( doublecomplex * ) ( ( DNformat * ) B.Store )->nzval; 
/*
  Print the soluton.
*/
  printf ( "\n" );
  printf ( "  Computed solution:\n" );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    printf ( "%g  %g\n", x[i].r, x[i].i );
  }
/*
  Demonstrate that X is really the solution.
*/
  b2 = cc_mv ( m, n, ncc, icc, ccc, acc, x );
  printf ( "\n" );
  printf ( "  Product A*X:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "%g  %g\n", b2[i].r, b2[i].i );
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
  printf ( "Z_SAMPLE_ST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;

# undef N
# undef NCC
}
/******************************************************************************/

doublecomplex *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], 
  doublecomplex acc[], doublecomplex x[] )

/******************************************************************************/
/*
  Purpose:

    CC_MV multiplies a CC matrix by a vector

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2014

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

    Input, doublecomplex ACC[NCC], the CC values.

    Input, doublecomplex X[N], the vector to be multiplied.

    Output, doublecomplex CC_MV[M], the product A*X.
*/
{
  doublecomplex *b;
  int i;
  int j;
  int k;

  b = ( doublecomplex * ) malloc ( m * sizeof ( doublecomplex ) );

  for ( i = 0; i < m; i++ )
  {
    b[i].r = 0.0;
    b[i].i = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = ccc[j]; k < ccc[j+1]; k++ )
    {
      i = icc[k];
      b[i].r = b[i].r + acc[k].r * x[j].r - acc[k].i * x[j].i;
      b[i].i = b[i].i + acc[k].r * x[j].i + acc[k].i * x[j].r;
    }
  }

  return b;
}
/******************************************************************************/

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], doublecomplex acc[], 
  char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT prints a sparse matrix in CC format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the matrix.

    Input, int N, the number of columns in the matrix.

    Input, int NCC, the number of CC elements.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, doublecomplex ACC[NCC], the CC values.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jnext;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J         A.r               A.i\n" );
  printf ( "  ----  ----  ----  ----------------  ----------------\n" );
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
    printf ( "  %4d  %4d  %4d  %16g  %16g\n", k, i, j, acc[k].r, acc[k].i );
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

