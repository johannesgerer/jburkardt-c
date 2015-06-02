# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cg_rc.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CG_RC_PRB.

  Discussion:

    CG_RC_PRB tests the CG_RC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CG_RC_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CG_RC library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CG_RC_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 uses CG_RC for the simple 1, -2, 1 matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2013

  Author:

    John Burkardt
*/
{
  double angle;
  double *b;
  double bnrm2;
  double err;
  int i;
  int it;
  int it_max;
  int job;
  int n = 21;
  double *p;
  double pi = 3.141592653589793;
  double *q;
  double *r;
  double rnrm2;
  double t;
  double tol;
  double *x;
  double *x_exact;
  double *z;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use CG_RC on the 1, -2, 1 matrix.\n" );
/*
  In order to specify the right hand side, pick an exact solution,
  and multiply by the matrix.
*/
  x_exact = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    angle = 2.0 * pi * ( double ) ( i ) / ( double ) ( n - 1 );
    x_exact[i] = sin ( angle );
  }

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = - 2.0 * x_exact[i];
  }
  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = b[i] + x_exact[i+1];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + x_exact[i-1];
  }
/*
  Here is the initial guess for the solution.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
/*
  Parameters for the stopping test.
*/
  it = 0;
  it_max = 30;
  tol = 1.0E-05;
  bnrm2 = 0.0;
  for ( i = 0; i < n; i++ )
  {
    bnrm2 = bnrm2 + b[i] * b[i];
  }
  bnrm2 = sqrt ( bnrm2 );
/*
  Set parameters for the CG_RC code.
*/
  r = ( double * ) malloc ( n * sizeof ( double ) );
  z = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );
  q = ( double * ) malloc ( n * sizeof ( double ) );

  job = 1;
/*
  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
*/
  for ( ; ; )
  {
    job = cg_rc ( n, b, x, r, z, p, q, job );
/*
  Compute q = A * p.
*/
    if ( job == 1 )
    {
      for ( i = 0; i < n; i++ )
      {
        q[i] = - 2.0 * p[i];
      }
      for ( i = 0; i < n - 1; i++ )
      {
        q[i] = q[i] + p[i+1];
      }
      for ( i = 1; i < n; i++ )
      {
        q[i] = q[i] + p[i-1];
      }
    }
/*
  Solve M * z = r.
*/
    else if ( job == 2 )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = r[i] / ( - 2.0 );
      }
    }
/*
  Compute r = r - A * x.
*/
    else if ( job == 3 )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i] = r[i] + 2.0 * x[i];
      }
      for ( i = 0; i < n - 1; i++ )
      {
        r[i] = r[i] - x[i+1];
      }
      for ( i = 1; i < n; i++ )
      {
        r[i] = r[i] - x[i-1];
      }
    }
/*
  Stopping test on R.
*/
    else if ( job == 4 )
    {
      rnrm2 = 0.0;
      for ( i = 0; i < n; i++ )
      {
        rnrm2 = rnrm2 + r[i] * r[i];
      }
      rnrm2 = sqrt ( rnrm2 );

      if ( bnrm2 == 0.0 )
      {
        if ( rnrm2 <= tol )
        {
          break;
        }
      }
      else
      {
        if ( rnrm2 <= tol * bnrm2 )
        {
          break;
        }
      }

      it = it + 1;

      if ( it_max <= it )
      {
        printf ( "\n" );
        printf ( "  Iteration limit exceeded.\n" );
        printf ( "  Terminating early.\n" );
        break;
      }
    }
    job = 2;
  }
  
  printf ( "\n" );
  printf ( "  Number of iterations was %d\n", it );
  printf ( "  Estimated error is %g\n", rnrm2 );
  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = fabs ( x_exact[i] - x[i] );
    if ( err < t )
    {
      err = t;
    }
  }
  printf ( "  Loo error is %g\n", err );

  printf ( "\n" );
  printf ( "     I      X(I)         X_EXACT(I)        B(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n",
      i, x[i], x_exact[i], b[i] );
  }

  free ( b );
  free ( p );
  free ( q );
  free ( r );
  free ( x );
  free ( x_exact );
  free ( z );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CG_RC with the Wathen matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *ax;
  double *b;
  double bnrm2;
  double err;
  int i;
  int ii;
  int it;
  int it_max;
  int job;
  int n;
  int nx;
  int ny;
  double *p;
  double *q;
  double *r;
  double rnrm2;
  int seed;
  double t;
  double tol;
  double *x;
  double *x_exact;
  double *z;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use CG_RC to solve a linear system\n" );
  printf ( "  involving the Wathen matrix.\n" );

  nx = 5;
  ny = 4;
 
  printf ( "\n" );
  printf ( "  NX = %d\n", nx );
  printf ( "  NY = %d\n", ny );
  printf ( "  N  = %d\n", n );

  n = wathen_order ( nx, ny );

  a = wathen ( nx, ny, n );

  seed = 123456789;
  x_exact = r8vec_uniform_01_new ( n, &seed );

  b = ( double * ) malloc ( n * sizeof ( double ) );
  r8mat_mv ( n, n, a, x_exact, b );
/*
  Here is the initial guess for the solution.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  ax = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Parameters for the stopping test.
*/
  it = 0;
  it_max = 30;
  tol = 1.0E-05;
  bnrm2 = 0.0;
  for ( i = 0; i < n; i++ )
  {
    bnrm2 = bnrm2 + b[i] * b[i];
  }
  bnrm2 = sqrt ( bnrm2 );
/*
  Set parameters for the CG_RC code.
*/
  r = ( double * ) malloc ( n * sizeof ( double ) );
  z = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );
  q = ( double * ) malloc ( n * sizeof ( double ) );
  job = 1;
/*
  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
*/
  for ( ; ; )
  {
    job = cg_rc ( n, b, x, r, z, p, q, job );
/*
  Compute q = A * p.
*/
    if ( job == 1 )
    {
      r8mat_mv ( n, n, a, p, q );
    }
/*
  Solve M * z = r.
*/
    else if ( job == 2 )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = r[i] / a[i+i*n];
      }
    }
/*
  Compute r = r - A * x.
*/
    else if ( job == 3 )
    {
      r8mat_mv ( n, n, a, x, ax );
      for ( i = 0; i < n; i++ )
      {
        r[i] = r[i] - ax[i];
      }
    }
/*
  Stopping test.
*/
    else if ( job == 4 )
    {
      rnrm2 = 0.0;
      for ( i = 0; i < n; i++ )
      {
        rnrm2 = rnrm2 + r[i] * r[i];
      }
      rnrm2 = sqrt ( rnrm2 );

      if ( bnrm2 == 0.0 )
      {
        if ( rnrm2 <= tol )
        {
          break;
        }
      }
      else
      {
        if ( rnrm2 <= tol * bnrm2 )
        {
          break;
        }
      }

      it = it + 1;

      if ( it_max <= it )
      {
        printf ( "\n" );
        printf ( "  Iteration limit exceeded.\n" );
        printf ( "  Terminating early.\n" );
        break;
      }
    }
    job = 2;
  }
  
  printf ( "\n" );
  printf ( "  Number of iterations was %d\n", it );
  printf ( "  Estimated error is %g\n", rnrm2 );
  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = fabs ( x_exact[i] - x[i] );
    if ( err < t )
    {
      err = t;
    }
  }
  printf ( "  Loo error is %g\n", err );

  printf ( "\n" );
  printf ( "     I      X(I)         X_EXACT(I)        B(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n",
      i, x[i], x_exact[i], b[i] );
  }

  free ( a );
  free ( ax );
  free ( b );
  free ( p );
  free ( q );
  free ( r );
  free ( x );
  free ( x_exact );
  free ( z );

  return;
}

