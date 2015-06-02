# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "linpack_d.h"
# include "blas1_d.h"

int main ( int argc, char *argv[] );
int get_seed ( void );
double *pseudo_inverse ( int m, int n, double u[], double s[], 
  double v[] );
void pseudo_linear_solve_test ( int m, int n, double a[], 
  double a_pseudo[], int *seed );
void pseudo_product_test ( int m, int n, double a[], double a_pseudo[] );
int r8_nint ( double x );
double r8mat_dif_fro ( int m, int n, double a[], double b[] );
double r8mat_norm_fro ( int m, int n, double a[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
void r8mat_svd_linpack ( int m, int n, double a[], double u[], double s[], 
  double v[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double r8vec_norm_l2 ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
void rank_one_print_test ( int m, int n, double a[], double u[], 
  double s[], double v[] );
void rank_one_test ( int m, int n, double a[], double u[], double s[], 
  double v[] );
int s_len_trim ( char *s );
void svd_product_test ( int m, int n, double a[], double u[], 
  double s[], double v[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SVD_DEMO.

  Discussion:

    SVD_DEMO demonstrates the singular value decomposition.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Usage:

    svd_demo m n

  Command Parameters:

    Command parameter, integer M, N, the number of rows and columns
    of the matrix.

  Local Parameters:

    Local, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Local, double S[M*N], the diagonal factor
    in the singular value decomposition of A.

    Local, int SEED, a seed used to define the random number generator.

    Output, double U[M*M], the first orthogonal factor
    in the singular value decomposition of A.

    Output, double V[N*N], the second orthogonal factor
    in the singular value decomposition of A.
*/
{
  double *a;
  double *a_pseudo;
  int i;
  int j;
  int m;
  int n;
  double *s;
  int seed;
  char string[80];
  double *u;
  double *v;

  timestamp ( );

  printf ( "\n" );
  printf ( "SVD_DEMO:\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Demonstrate the singular value decomposition (SVD)\n" );
  printf ( "\n" );
  printf ( "  A real MxN matrix A can be factored as:\n" );
  printf ( "\n" );
  printf ( "    A = U * S * V'\n" );
  printf ( "\n" );
  printf ( "  where\n" );
  printf ( "\n" );
  printf ( "    U = MxM orthogonal,\n" );
  printf ( "    S = MxN zero except for diagonal,\n" );
  printf ( "    V = NxN orthogonal.\n" );
  printf ( "\n" );
  printf ( "  The diagonal of S contains only nonnegative numbers\n" );
  printf ( "  and these are arranged in descending order.\n" );
/*
  If M was not on the command line, get it now.
*/
  if ( argc < 2 ) 
  {
    printf ( "\n" );
    printf ( "SVD_DEMO:\n" );
    printf ( "  Please enter the value of M:\n" );
    printf ( "  (Number of rows in matrix A).\n" );
    printf ( "  (We prefer M <= 10!).\n" );
    scanf ( "d", &m );
  }
  else 
  {
    strcpy ( string, argv[1] );
    m = atoi ( string );
  }
  printf ( "\n" );
  printf ( "  Matrix row order    M = %d\n", m );
/*
  If N was not on the command line, get it now.
*/
  if ( argc < 3 ) 
  {
    printf ( "\n" );
    printf ( "SVD_DEMO:\n" );
    printf ( "  Please enter the value of N:\n" );
    printf ( "  (Number of columns in matrix A).\n" );
    printf ( "  (We prefer N <= 10!).\n" );
    scanf ( "%d", &n );
  }
  else 
  {
    strcpy ( string, argv[2] );
    n = atoi ( string );
  }
  printf ( "  Matrix column order N = %d\n", n );
/*
  If SEED was not on the command line, use GET_SEED.
*/
  if ( argc < 4 ) 
  {
    seed = get_seed ( );
    printf ( "  Random number SEED    = %d\n", seed );
    printf ( "  (Chosen by the program.)\n" );
  }
  else 
  {
    strcpy ( string, argv[3] );
    seed = atoi ( string );
    printf ( "  Random number SEED    = %d\n", seed );
    printf ( "  (Chosen by the user.)\n" );
  }
/*
  Set aside space for the arrays.
*/
  u = ( double * ) malloc ( m * m * sizeof ( double ) );
  s = ( double * ) malloc ( m * n * sizeof ( double ) );
  v = ( double * ) malloc ( n * n * sizeof ( double ) );
/*
  Generate the matrix.
*/
  printf ( "\n" );
  printf ( "  We choose a \"random\" matrix A, with integral\n" );
  printf ( "  values between 0 and 10.\n" );

  a = r8mat_uniform_01_new ( m, n, &seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+m*j] = r8_nint ( 10.0 * a[i+m*j] );
    }
  }
  r8mat_print ( m, n, a, "  The matrix A:\n" );
/*
  Get the SVD from LINPACK.
*/
  r8mat_svd_linpack ( m, n, a, u, s, v );
/*
  Print the SVD.
*/
  r8mat_print ( m, m, u, "  The orthogonal factor U:" );

  r8mat_print ( m, n, s, "  The diagonal factor S:" );

  r8mat_print ( n, n, v, "  The orthogonal factor V:" );
/*
  Check that A = U * S * V'.
*/
  svd_product_test ( m, n, a, u, s, v );
/*
  Compute the norm of the difference between A and the successive
  sums of rank one approximants.
*/
  rank_one_test ( m, n, a, u, s, v );
/*
  Actually print the sums of rank one approximants.
*/
  rank_one_print_test ( m, n, a, u, s, v );
/*
  Compute the pseudoinverse.
*/
  a_pseudo = pseudo_inverse ( m, n, u, s, v );

  r8mat_print ( n, m, a_pseudo, "  The pseudoinverse of A:" );
/*
  Test A*A+ = I+, A+*A = I+
*/
  pseudo_product_test ( m, n, a, a_pseudo );
/*
  Demonstrate the use of the pseudoinverse for linear systems.
*/
  pseudo_linear_solve_test ( m, n, a, a_pseudo, &seed );
/*
  Free memory.
*/
  free ( a );
  free ( a_pseudo );
  free ( s );
  free ( u );
  free ( v );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SVD_DEMO:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int get_seed ( void )

/******************************************************************************/
/*
  Purpose:

    GET_SEED returns a random seed for the random number generator.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2004

  Author:

    John Burkardt

  Parameters:

    Output, int GET_SEED, a random seed value.
*/
{
  time_t clock;
  int i;
  int i4_huge = 2147483647;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
/*
  If the internal seed is 0, generate a value based on the time.
*/
  clock = time ( &tloc );
  lt = localtime ( &clock );
/*
  Hours is 1, 2, ..., 12.
*/
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
/*
  Move Hours to 0, 1, ..., 11
*/
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
/*
  We want values in [1,43200], not [0,43199].
*/
  seed = seed + 1;
/*
  Remap SEED from [1,43200] to [1,HUGE].
*/
  seed = ( int ) 
    ( ( ( double ) seed )
    * ( ( double ) i4_huge ) / ( 60.0 * 60.0 * 12.0 ) );
/*
  Never use a seed of 0.
*/
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
}
/******************************************************************************/

double *pseudo_inverse ( int m, int n, double u[], double s[], 
  double v[] )

/******************************************************************************/
/*
  Purpose:

    PSEUDO_INVERSE computes the pseudoinverse.

  Discussion:

    Given the singular value decomposition of a real MxN matrix A:

      A = U * S * V'

    where 

      U is MxM orthogonal,
      S is MxN, and entirely zero except for the diagonal;
      V is NxN orthogonal.

    the pseudo inverse is the NxM matrix A+ with the form

      A+ = V * S+ * U'

    where 

      S+ is the NxM matrix whose nonzero diagonal elements are
      the inverses of the corresponding diagonal elements of S.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double U[M*M], S[M*N], V[N*N], the factors
    that form the singular value decomposition of A.

    Output, double PSEUDO_INVERSE[N*M], the pseudo_inverse of A.
*/
{
  double *a_pseudo;
  int i;
  int j;
  int k;
  double *sp;
  double *sput;

  sp = ( double * ) malloc ( n * m * sizeof ( double ) );
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j && s[i+i*m] != 0.0 )
      {
        sp[i+j*n] = 1.0 / s[i+i*m];
      }
      else
      {
        sp[i+j*n] = 0.0;
      }
    }
  }

  sput = ( double * ) malloc ( n * m * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      sput[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        sput[i+j*n] = sput[i+j*n] + sp[i+k*n] * u[j+k*m];
      }
    }
  }

  a_pseudo = ( double * ) malloc ( n * m * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      a_pseudo[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        a_pseudo[i+j*n] = a_pseudo[i+j*n] + v[i+k*n] * sput[k+j*n];
      }
    }
  }

  free ( sp );

  return a_pseudo;
}
/******************************************************************************/

void pseudo_linear_solve_test ( int m, int n, double a[], 
  double a_pseudo[], int *seed )

/******************************************************************************/
/*
  Purpose:

    PSEUDO_LINEAR_SOLVE_TEST uses the pseudoinverse for linear systems.

  Discussion:

    Given an MxN matrix A, and its pseudoinverse A+:

      "Solve" A  * x = b by x = A+  * b.

      "Solve" A' * x = b by x = A+' * b.

    When the system is overdetermined, the solution minimizes the
    L2 norm of the residual.  

    When the system is underdetermined, the solution
    is the solution of minimum L2 norm.     

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double A_PSEUDO[N*M], the pseudo_inverse of A.

    Input/output, int *SEED, a seed for the random number generator.
*/
{
  double *bm;
  double *bn;
  int i;
  int j;
  double *rm;
  double *rn;
  double *xm1;
  double *xm2;
  double *xn1;
  double *xn2;

  printf ( "\n" );
  printf ( "PSEUDO_LINEAR_SOLVE_TEST\n" );
/*
  A * x = b, b in range of A.
*/
  xn1 = r8vec_uniform_01_new ( n, seed );
  for ( i = 0; i < n; i++ )
  {
    xn1[i] = r8_nint ( 10.0 * xn1[i] );
  }

  bm = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    bm[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      bm[i] = bm[i] + a[i+j*m] * xn1[j];
    }
  }

  xn2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    xn2[i] = 0.0;
    for ( j = 0; j < m; j++ )
    {
      xn2[i] = xn2[i] + a_pseudo[i+j*n] * bm[j];
    }
  }

  rm = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    rm[i] = bm[i];
    for ( j = 0; j < n; j++ )
    {
      rm[i] = rm[i] - a[i+j*m] * xn2[j];
    }
  }

  printf ( "\n" );
  printf ( "  Given:\n" );
  printf ( "    b = A * x1\n" );
  printf ( "  so that b is in the range of A, solve\n" );
  printf ( "    A * x = b\n" );
  printf ( "  using the pseudoinverse:\n" );
  printf ( "    x2 = A+ * b.\n" );
  printf ( "\n" );
  printf ( "  Norm of x1 = %g\n", r8vec_norm_l2 ( n, xn1 ) );
  printf ( "  Norm of x2 = %g\n", r8vec_norm_l2 ( n, xn2 ) );
  printf ( "  Norm of residual = %g\n", r8vec_norm_l2 ( m, rm ) );

  free ( bm );
  free ( rm );
  free ( xn1 );
  free ( xn2 );
/*
  A * x = b, b not in range of A.
*/
  if ( n < m )
  {
    printf ( "\n" );
    printf ( "  For N < M, most systems A*x=b will not be\n" );
    printf ( "  exactly and uniquely solvable, except in the\n" );
    printf ( "  least squares sense.\n" );
    printf ( "\n" );
    printf ( "  Here is an example:\n" );

    bm = r8vec_uniform_01_new ( m, seed );

    xn2 = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      xn2[i] = 0.0;
      for ( j = 0; j < m; j++ )
      {
        xn2[i] = xn2[i] + a_pseudo[i+j*n] * bm[j];
      }
    }

    rm = ( double * ) malloc ( m * sizeof ( double ) );
    for ( i = 0; i < m; i++ )
    {
      rm[i] = bm[i];
      for ( j = 0; j < n; j++ )
      {
        rm[i] = rm[i] - a[i+j*m] * xn2[j];
      }
    }

    printf ( "\n" );
    printf ( "  Given b is NOT in the range of A, solve\n" );
    printf ( "    A * x = b\n" );
    printf ( "  using the pseudoinverse:\n" );
    printf ( "    x2 = A+ * b.\n" );
    printf ( "\n" );
    printf ( "  Norm of x2 = %g\n", r8vec_norm_l2 ( n, xn2 ) );
    printf ( "  Norm of residual = %g\n", r8vec_norm_l2 ( m, rm ) );

    free ( bm );
    free ( rm );
    free ( xn2 );
  }
/*
  A' * x = b, b is in the range of A'.
*/
  xm1 = r8vec_uniform_01_new ( m, seed );
  for ( i = 0; i < m; i++ )
  {
    xm1[i] = r8_nint ( 10.0 * xm1[i] );
  }

  bn = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    bn[i] = 0.0;
    for ( j = 0; j < m; j++ )
    {
      bn[i] = bn[i] + a[j+i*m] * xm1[j];
    }
  }

  xm2 = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    xm2[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      xm2[i] = xm2[i] + a_pseudo[j+i*n] * bn[j];
    }
  }

  rn = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    rn[i] = bn[i];
    for ( j = 0; j < m; j++ )
    {
      rn[i] = rn[i] - a[j+i*m] * xm2[j];
    }
  }
  printf ( "\n" );
  printf ( "  Given:\n" );
  printf ( "    b = A' * x1\n" );
  printf ( "  so that b is in the range of A', solve\n" );
  printf ( "    A' * x = b\n" );
  printf ( "  using the pseudoinverse:\n" );
  printf ( "    x2 = A+' * b.\n" );
  printf ( "\n" );
  printf ( "  Norm of x1 = %g\n", r8vec_norm_l2 ( m, xm1 ) );
  printf ( "  Norm of x2 = %g\n", r8vec_norm_l2 ( m, xm2 ) );
  printf ( "  Norm of residual = %g\n", r8vec_norm_l2 ( n, rn ) );

  free ( bn );
  free ( rn );
  free ( xm1 );
  free ( xm2 );
/*
  A' * x = b, b is not in the range of A'.
*/
  if ( m < n )
  {
    printf ( "\n" );
    printf ( "  For M < N, most systems A'*x=b will not be\n" );
    printf ( "  exactly and uniquely solvable, except in the\n" );
    printf ( "  least squares sense.\n" );
    printf ( "\n" );
    printf ( "  Here is an example:\n" );

    bn = r8vec_uniform_01_new ( n, seed );

    xm2 = ( double * ) malloc ( m * sizeof ( double ) );
    for ( i = 0; i < m; i++ )
    {
      xm2[i] = 0.0;
      for ( j = 0; j < n; j++ )
      {
        xm2[i] = xm2[i] + a_pseudo[j+i*n] * bn[j];
      }
    }

    rn = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      rn[i] = bn[i];
      for ( j = 0; j < m; j++ )
      {
        rn[i] = rn[i] - a[j+i*m] * xm2[j];
      }
    }

    printf ( "\n" );
    printf ( "  Given b is NOT in the range of A', solve\n" );
    printf ( "    A' * x = b\n" );
    printf ( "  using the pseudoinverse:\n" );
    printf ( "    x2 = A+ * b.\n" );
    printf ( "\n" );
    printf ( "  Norm of x2 = %g\n", r8vec_norm_l2 ( m, xm2 ) );
    printf ( "  Norm of residual = %g\n", r8vec_norm_l2 ( n, rn ) );

    free ( bn );
    free ( rn );
    free ( xm2 );
  }

  return;
}
/******************************************************************************/

void pseudo_product_test ( int m, int n, double a[], double a_pseudo[] )

/******************************************************************************/
/*
  Purpose:

    PSEUDO_PRODUCT_TEST examines pseudoinverse products.

  Discussion:

    Given an MxN matrix A, and its pseudoinverse A+, we must have

      A+ * A * A+ = A+
      A * A+ * A = A
      ( A * A+ )' = A * A+ (MxM symmetry)
      ( A+ * A )' = A+ * A (NxN symmetry)

    If M <= N, A * A+ may be "interesting" (equal to or "like" the identity),
    if N <= M, A+ * A may be "interesting" (equal to or "like" the identity).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double A_PSEUDO[N*M], the pseudo_inverse of A.
*/
{
  double *bmm;
  double *bmn;
  double *bnm;
  double *bnn;
  double dif1;
  double dif2;
  double dif3;
  double dif4;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "PSEUDO_PRODUCT_TEST\n" );
  printf ( "\n" );
  printf ( "  The following relations MUST hold:\n" );
  printf ( "\n" );
  printf ( "   A  * A+ * A  = A\n" );
  printf ( "   A+ * A  * A+ = A+\n" );
  printf ( " ( A  * A+ ) is MxM symmetric;\n" );
  printf ( " ( A+ * A  ) is NxN symmetric\n" );
/*
  Compute A * A+ * A.
*/
  bnn = ( double * ) malloc ( n * n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      bnn[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        bnn[i+j*n] = bnn[i+j*n] + a_pseudo[i+k*n] * a[k+j*m];
      }
    }
  }
  bmn = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      bmn[i+j*m] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        bmn[i+j*m] = bmn[i+j*m] + a[i+k*m] * bnn[k+j*n];
      }
    }
  }
  dif1 = r8mat_dif_fro ( m, n, a, bmn );

  free ( bmn );
  free ( bnn );
/*
  Compute A+ * A * A+.
*/
  bmm = ( double * ) malloc ( m * m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      bmm[i+j*m] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        bmm[i+j*m] = bmm[i+j*m] + a[i+k*m] * a_pseudo[k+j*n];
      }
    }
  }

  bnm = ( double * ) malloc ( n * m * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      bnm[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        bnm[i+j*n] = bnm[i+j*n] + a_pseudo[i+k*n] * bmm[k+j*m];
      }
    }
  }

  dif2 = r8mat_dif_fro ( n, m, a_pseudo, bnm );

  free ( bnm );
  free ( bmm );
/*
  Compute norm of A * A+ - (A * A+)'.
*/
  bmm = ( double * ) malloc ( m * m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      bmm[i+j*m] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        bmm[i+j*m] = bmm[i+j*m] + a[i+k*m] * a_pseudo[k+j*n];
      }
    }
  }
  dif3 = 0.0;
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif3 = dif3 + pow ( bmm[i+j*m] - bmm[j+i*m], 2 );
    }
  }
  dif3 = sqrt ( dif3 );

  free ( bmm );
/*
  Compute norm of A+ * A - (A+ * A)'
*/
  bnn = ( double * ) malloc ( n * n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      bnn[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        bnn[i+j*n] = bnn[i+j*n] + a_pseudo[i+k*n] * a[k+j*m];
      }
    }
  }
  dif4 = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      dif4 = dif4 + pow ( bnn[i+j*n] - bnn[j+i*n], 2 );
    }
  }
  dif4 = sqrt ( dif4 );

  free ( bnn );
/*
  Report.
*/
  printf ( "\n" );
  printf ( "  Here are the Frobenius norms of the errors\n" );
  printf ( "  in these relationships:\n" );
  printf ( "\n" );
  printf ( "   A  * A+ * A  = A            %g\n", dif1 );
  printf ( "   A+ * A  * A+ = A+           %g\n", dif2 );
  printf ( " ( A  * A+ ) is MxM symmetric; %g\n", dif3 );
  printf ( " ( A+ * A  ) is NxN symmetric; %g\n", dif4 );

  printf ( "\n" );
  printf ( "  In some cases, the matrix A * A+\n" );
  printf ( "  may be interesting (if M <= N, then\n" );
  printf ( "  it MIGHT look like the identity.)\n" );
  printf ( "\n" );
  bmm = ( double * ) malloc ( m * m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      bmm[i+j*m] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        bmm[i+j*m] = bmm[i+j*m] + a[i+k*m] * a_pseudo[k+j*n];
      }
    }
  }
  r8mat_print ( m, m, bmm, "  A * A+:" );

  free ( bmm );

  printf ( "\n" );
  printf ( "  In some cases, the matrix A+ * A\n" );
  printf ( "  may be interesting (if N <= M, then\n" );
  printf ( "  it MIGHT look like the identity.)\n" );
  printf ( "\n" );

  bnn = ( double * ) malloc ( n * n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      bnn[i+j*n] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        bnn[i+j*n] = bnn[i+j*n] + a_pseudo[i+k*n] * a[k+j*m];
      }
    }
  }

  r8mat_print ( n, n, bnn, "  A+ * A" );

  free ( bnn );

  return;
}
/******************************************************************************/

int r8_nint ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_NINT returns the nearest integer to an R8.

  Example:

        X         R8_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value.

    Output, int R8_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

double r8mat_dif_fro ( int m, int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIF_FRO returns the Frobenius norm of the difference of R8MAT's.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of double precision values, which
    may be stored as a vector in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )

    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], double B[M*N], the matrices for which we
    want the Frobenius norm of the difference.

    Output, double R8MAT_DIF_FRO, the Frobenius norm of ( A - B ).
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - b[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double r8mat_norm_fro ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose Frobenius
    norm is desired.

    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

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

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

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
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

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
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

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

void r8mat_svd_linpack ( int m, int n, double a[], double u[], double s[], 
  double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.

  Discussion:

    The singular value decomposition of a real MxN matrix A has the form:

      A = U * S * V'

    where 

      U is MxM orthogonal,
      S is MxN, and entirely zero except for the diagonal;
      V is NxN orthogonal.

    Moreover, the nonzero entries of S are positive, and appear
    in order, from largest magnitude to smallest.

    This routine calls the LINPACK routine DSVDC to compute the
    factorization.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Output, double U[M*M], S[M*N], V[N*N], the factors
    that form the singular value decomposition of A.
*/
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int lwork;
  double *sdiag;
  double *work;
/*
  The correct size of E and SDIAG is min ( m+1, n).
*/
  a_copy = ( double * ) malloc ( m * n * sizeof ( double ) );
  e = ( double * ) malloc ( ( m + n ) * sizeof ( double ) );
  sdiag = ( double * ) malloc ( ( m + n )  * sizeof ( double ) );
  work = ( double * ) malloc ( m * sizeof ( double ) );
/*
  Compute the eigenvalues and eigenvectors.
*/
  job = 11;
  lda = m;
  ldu = m;
  ldv = n;
/*
  The input matrix is destroyed by the routine.  Since we need to keep
  it around, we only pass a copy to the routine.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    { 
      a_copy[i+j*m] = a[i+j*m];
    }
  }
  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );
 
  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "R8MAT_SVD_LINPACK - Failure!\n" );
    printf ( "  The SVD could not be calculated.\n" );
    printf ( "  LINPACK routine DSVDC returned a nonzero\n" );
    printf ( "  value of the error flag, INFO = %d\n", info );
    return;
  }
/*
  Make the MxN matrix S from the diagonal values in SDIAG.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        s[i+j*m] = sdiag[i];
      }
      else
      {
        s[i+j*m] = 0.0;
      }
    }
  }
/*
  Note that we do NOT need to transpose the V that comes out of LINPACK!
*/
  free ( a_copy );
  free ( e );
  free ( sdiag );
  free ( work );

  return;
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with unit pseudorandom values.

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

double r8vec_norm_l2 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.

  Discussion:

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_NORM_L2, the L2 norm of A.
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

void rank_one_print_test ( int m, int n, double a[], double u[], 
  double s[], double v[] )

/******************************************************************************/
/*
  Purpose:

    RANK_ONE_PRINT_TEST prints the sums of rank one matrices.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double U[M*M], S[M*N], V[N*N], the factors
    that form the singular value decomposition of A.
*/
{
  double a_norm;
  double dif_norm;
  int i;
  int j;
  int k;
  int r;
  double *svt;
  char title[100];
  double *usvt;

  a_norm = r8mat_norm_fro ( m, n, a );

  printf ( "\n" );
  printf ( "RANK_ONE_PRINT_TEST:\n" );
  printf ( "  Print the sums of R rank one matrices.\n" );

  for ( r = 0; r <= i4_min ( m, n ); r++ )
  {
    svt = ( double * ) malloc ( r * n * sizeof ( double ) );
    for ( i = 0; i < r; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        svt[i+j*r] = 0.0;
        for ( k = 0; k < r; k++ )
        {
          svt[i+j*r] = svt[i+j*r] + s[i+k*m] * v[j+k*n];
        }
      }
    }
    usvt = ( double * ) malloc ( m * n * sizeof ( double ) );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        usvt[i+j*m] = 0.0;
        for ( k = 0; k < r; k++ )
        { 
          usvt[i+j*m] = usvt[i+j*m] + u[i+k*m] * svt[k+j*r];
        }
      }
    }

    sprintf ( title, "  Rank R = %d", r );
    r8mat_print ( m, n, usvt, title );


    free ( svt );
    free ( usvt );
  }

  r8mat_print ( m, n, a, "  Original matrix A:" );

  return;
}
/******************************************************************************/

void rank_one_test ( int m, int n, double a[], double u[], double s[], 
  double v[] )

/******************************************************************************/
/*
  Purpose:

    RANK_ONE_TEST compares A to the sum of rank one matrices.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double U[M*M], S[M*N], V[N*N], the factors
    that form the singular value decomposition of A.
*/
{
  double a_norm;
  double dif_norm;
  int i;
  int j;
  int k;
  int r;
  double *svt;
  double *usvt;

  a_norm = r8mat_norm_fro ( m, n, a );

  printf ( "\n" );
  printf ( "RANK_ONE_TEST:\n" );
  printf ( "  Compare A to the sum of R rank one matrices.\n" );
  printf ( "\n" );
  printf ( "         R    Absolute      Relative\n" );
  printf ( "              Error         Error\n" );
  printf ( "\n" );

  for ( r = 0; r <= i4_min ( m, n ); r++ ) 
  {
    svt = ( double * ) malloc ( r * n * sizeof ( double ) );
    for ( i = 0; i < r; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        svt[i+j*r] = 0.0;
        for ( k = 0; k < r; k++ )
        {
          svt[i+j*r] = svt[i+j*r] + s[i+k*m] * v[j+k*n];
        }
      }
    }
    usvt = ( double * ) malloc ( m * n * sizeof ( double ) );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        usvt[i+j*m] = 0.0;
        for ( k = 0; k < r; k++ )
        { 
          usvt[i+j*m] = usvt[i+j*m] + u[i+k*m] * svt[k+j*r];
        }
      }
    }
    dif_norm = r8mat_dif_fro ( m, n, a, usvt );

    printf ( "  %8d  %14g  %14g\n", r, dif_norm, dif_norm / a_norm );

    free ( svt );
    free ( usvt );
  }
  return;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

void svd_product_test ( int m, int n, double a[], double u[], 
  double s[], double v[] )

/******************************************************************************/
/*
  Purpose:

    SVD_PRODUCT_TEST tests that A = U * S * V'.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix A.

    Input, double A[M*N], the matrix whose singular value
    decomposition we are investigating.

    Input, double U[M*M], S[M*N], V[N*N], the factors
    that form the singular value decomposition of A.
*/
{
  double a_norm;
  double dif_norm;
  int i;
  int j;
  int k;
  double *svt;
  double *usvt;

  a_norm = r8mat_norm_fro ( m, n, a );

  svt = ( double * ) malloc ( m * n * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      svt[i+j*m] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        svt[i+j*m] = svt[i+j*m] + s[i+k*m] * v[j+k*n];
      }
    }
  }
  usvt = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      usvt[i+j*m] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        usvt[i+j*m] = usvt[i+j*m] + u[i+k*m] * svt[k+j*m];
      }
    }
  }

  r8mat_print ( m, n, usvt, "  The product U * S * V':" );

  dif_norm = r8mat_dif_fro ( m, n, a, usvt );

  printf ( "\n" );
  printf ( "  Frobenius Norm of A, A_NORM = %g\n", a_norm );
  printf ( "\n" );
  printf ( "  ABSOLUTE ERROR for A = U*S*V'\n" );
  printf ( "  Frobenius norm of difference A-U*S*V' = %g\n", dif_norm );
  printf ( "\n" );
  printf ( "  RELATIVE ERROR for A = U*S*V':\n" );
  printf ( "  Ratio of DIF_NORM / A_NORM = %g\n", dif_norm / a_norm );

  free ( svt );
  free ( usvt );

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
