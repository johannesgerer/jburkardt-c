# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "vandermonde.h"

int main ( );
void bivand1_test ( );
void bivand2_test ( );
void dvand_test ( );
void dvandprg_test ( );
void pvand_test ( );
void pvandprg_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for VANDERMONDE_PRB.

  Discussion:

    VANDERMONDE_TEST tests the VANDERMONDE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "VANDERMONDE_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the VANDERMONDE library.\n" );

  bivand1_test ( );
  bivand2_test ( );
  dvand_test ( );
  dvandprg_test ( );
  pvand_test ( );
  pvandprg_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "VANDERMONDE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void bivand1_test ( )

/******************************************************************************/
/*
  Purpose:

    BIVAND1_TEST tests BIVAND1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2014

  Author:

    John Burkardt
*/
{
# define N 3

  double *a;
  double alpha[N] = { 1.0, 2.0, 3.0 };
  double beta[N] = { 10.0, 20.0, 30.0 };
  int n = N;
  int n2;

  printf ( "\n" );
  printf ( "BIVAND1_TEST:\n" );
  printf ( "  Compute a bidimensional Vandermonde matrix\n" );
  printf ( "  associated with polynomials of total degree\n" );
  printf ( "  less than N.\n" );

  r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );
  r8vec_print ( n, beta, "  Vandermonde vector BETA:" );

  a = bivand1 ( n, alpha, beta );

  n2 = ( n * ( n + 1 ) ) / 2;
  r8mat_print ( n2, n2, a, "  Bidimensional Vandermonde matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void bivand2_test ( )

/******************************************************************************/
/*
  Purpose:

    BIVAND2_TEST tests BIVAND2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 May 2014

  Author:

    John Burkardt
*/
{
# define N 3

  double *a;
  double alpha[N] = { 1.0, 2.0, 3.0 };
  double beta[N] = { 10.0, 20.0, 30.0 };
  int n = N;
  int n2;

  printf ( "\n" );
  printf ( "BIVAND2_TEST:\n" );
  printf ( "  Compute a bidimensional Vandermonde matrix\n" );
  printf ( "  associated with polynomials of maximum degree\n" );
  printf ( "  less than N.\n" );

  r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );
  r8vec_print ( n, beta, "  Vandermonde vector BETA:" );

  a = bivand2 ( n, alpha, beta );

  n2 = n * n;
  r8mat_print ( n2, n2, a, "  Bidimensional Vandermonde matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void dvand_test ( )

/******************************************************************************/
/*
  Purpose:

    DVAND_TEST tests DVAND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2014

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  int n = N;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  printf ( "\n" );
  printf ( "DVAND_TEST:\n" );
  printf ( "  Solve a Vandermonde linear system A'*x=b\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 ) 
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, &seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mtv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    free ( x );

    x = dvand ( n, alpha, b );
    r8vec_print ( n, x, "  Solution X:" );

    free ( a );
    free ( alpha );
    free ( b );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void dvandprg_test ( )

/******************************************************************************/
/*
  Purpose:

    DVANDPRG_TEST tests DVANDPRG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2014

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  double *c;
  double *m;
  int n = N;
  int nsub;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  printf ( "\n" );
  printf ( "DVANDPRG_TEST:\n" );
  printf ( "  Solve a Vandermonde linear system A'*x=b\n" );
  printf ( "  progressively.\n" );
  printf ( "  First we use ALPHA = 0, 1, 2, 3, 4.\n" );
  printf ( "  Then we choose ALPHA at random.\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 ) 
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, &seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mtv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    free ( x );

    x = ( double * ) malloc ( n * sizeof ( double ) );
    c = ( double * ) malloc ( n * sizeof ( double ) );
    m = ( double * ) malloc ( n * sizeof ( double ) );

    for ( nsub = 1; nsub <= n; nsub++ )
    {
      dvandprg ( nsub, alpha, b, x, c, m );
      r8vec_print ( nsub, x, "  Solution X:" );
    }

    free ( a );
    free ( alpha );
    free ( b );
    free ( c );
    free ( m );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void pvand_test ( )

/******************************************************************************/
/*
  Purpose:

    PVAND_TEST tests PVAND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2014

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  int n = N;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  printf ( "\n" );
  printf ( "PVAND_TEST:\n" );
  printf ( "  Solve a Vandermonde linear system A*x=b\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, &seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    free ( x );

    x = pvand ( n, alpha, b );
    r8vec_print ( n, x, "  Solution X:" );

    free ( a );
    free ( alpha );
    free ( b );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void pvandprg_test ( )

/******************************************************************************/
/*
  Purpose:

    PVANDPRG_TEST tests PVANDPRG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 April 2014

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  double *c;
  double *m;
  int n = N;
  int nsub;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  printf ( "\n" );
  printf ( "PVANDPRG_TEST:\n" );
  printf ( "  Solve a Vandermonde linear system A*x=b\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, &seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    free ( x );

    x = ( double * ) malloc ( n * sizeof ( double ) );
    c = ( double * ) malloc ( n * sizeof ( double ) );
    m = ( double * ) malloc ( n * sizeof ( double ) );

    for ( nsub = 1; nsub <= n; nsub++ )
    {
      pvandprg ( nsub, alpha, b, x, c, m );
      r8vec_print ( nsub, x, "  Solution X:" );
    }

    free ( a );
    free ( alpha );
    free ( b );
    free ( c );
    free ( m );
    free ( x );
  }

  return;
# undef N
}

