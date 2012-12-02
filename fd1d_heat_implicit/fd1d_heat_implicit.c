# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

int main ( void );
void f ( double a, double b, double t0, double t, int n, double x[], 
  double value[] );
int r83_np_fa ( int n, double a[] );
double *r83_np_sl ( int n, double a_lu[], double b[], int job );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void timestamp ( void );
void u0 ( double a, double b, double t0, int n, double x[], double value[] );
double ua ( double a, double b, double t0, double t );
double ub ( double a, double b, double t0, double t );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FD1D_HEAT_IMPLICIT.

  Discussion:

    FD1D_HEAT_IMPLICIT solves the 1D heat equation with an implicit method.

    This program solves

      dUdT - k * d2UdX2 = F(X,T)

    over the interval [A,B] with boundary conditions

      U(A,T) = UA(T),
      U(B,T) = UB(T),

    over the time interval [T0,T1] with initial conditions

      U(X,T0) = U0(X)

    The code uses the finite difference method to approximate the
    second derivative in space, and an implicit backward Euler approximation
    to the first derivative in time.

    The finite difference form can be written as

      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
      ------------------  = F(X,T) + k *  ------------------------------------
               dt                                   dx * dx

    so that we have the following linear system for the values of U at time T+dt:

            -     k * dt / dx / dx   * U(X-dt,T+dt)
      + ( 1 + 2 * k * dt / dx / dx ) * U(X,   T+dt)
            -     k * dt / dx / dx   * U(X+dt,T+dt)
      =               dt             * F(X,   T+dt)
      +                                U(X,   T)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 June 2009

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *fvec;
  int i;
  int info;
  int j;
  int job;
  double k;
  double *t;
  double t_delt;
  char *t_file = "t.txt";
  double t_max;
  double t_min;
  int t_num;
  double *u;
  char *u_file = "u.txt";
  double w;
  double *x;
  double x_delt;
  char *x_file = "x.txt";
  double x_max;
  double x_min;
  int x_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "FD1D_HEAT_IMPLICIT\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Finite difference solution of\n" );
  printf ( "  the time dependent 1D heat equation\n" );
  printf ( "\n" );
  printf ( "    Ut - k * Uxx = F(x,t)\n" );
  printf ( "\n" );
  printf ( "  for space interval A <= X <= B with boundary conditions\n" );
  printf ( "\n" );
  printf ( "    U(A,t) = UA(t)\n" );
  printf ( "    U(B,t) = UB(t)\n" );
  printf ( "\n" );
  printf ( "  and time interval T0 <= T <= T1 with initial condition\n" );
  printf ( "\n" );
  printf ( "    U(X,T0) = U0(X).\n" );
  printf ( "\n" );
  printf ( "  A second order difference approximation is used for Uxx.\n" );
  printf ( "\n" );
  printf ( "  A first order backward Euler difference approximation\n" );
  printf ( "  is used for Ut.\n" );

  k = 5.0E-07;
/*
  Set X values.
*/
  x_min = 0.0;
  x_max = 0.3;
  x_num = 11;
  x_delt = ( x_max - x_min ) / ( double ) ( x_num - 1 );

  x = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( i = 0; i < x_num; i++ )
  {
    x[i] = ( ( double ) ( x_num - i - 1 ) * x_min   
           + ( double ) (         i     ) * x_max ) 
           / ( double ) ( x_num     - 1 );
  }
/*
  Set T values.
*/
  t_min = 0.0;
  t_max = 22000.0;
  t_num = 51;
  t_delt = ( t_max - t_min ) / ( double ) ( t_num - 1 );

  t = ( double * ) malloc ( t_num * sizeof ( double ) );

  for ( j = 0; j < t_num; j++ )
  {
    t[j] = ( ( double ) ( t_num - j - 1 ) * t_min   
           + ( double ) (         j     ) * t_max ) 
           / ( double ) ( t_num     - 1 );
  }
/*
  Set the initial data, for time T_MIN.
*/
  u = ( double * ) malloc ( x_num * t_num * sizeof ( double ) );

  u0 ( x_min, x_max, t_min, x_num, x, u );
/*
  The matrix A does not change with time.  We can set it once,
  factor it once, and solve repeatedly.
*/
  w = k * t_delt / x_delt / x_delt;

  a = ( double * ) malloc ( 3 * x_num * sizeof ( double ) );

  a[0+0*3] = 0.0;

  a[1+0*3] = 1.0;
  a[0+1*3] = 0.0;

  for ( i = 1; i < x_num - 1; i++ )
  {
    a[2+(i-1)*3] =           - w;
    a[1+ i   *3] = 1.0 + 2.0 * w;
    a[0+(i+1)*3] =           - w;
  }

  a[2+(x_num-2)*3] = 0.0;
  a[1+(x_num-1)*3] = 1.0;

  a[2+(x_num-1)*3] = 0.0;
/*
  Factor the matrix.
*/
  info = r83_np_fa ( x_num, a );

  b = ( double * ) malloc ( x_num * sizeof ( double ) );
  fvec = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( j = 1; j < t_num; j++ )
  {
/*
  Set the right hand side B.
*/
    b[0] = ua ( x_min, x_max, t_min, t[j] );

    f ( x_min, x_max, t_min, t[j-1], x_num, x, fvec );

    for ( i = 1; i < x_num - 1; i++ )
    {
      b[i] = u[i+(j-1)*x_num] + t_delt * fvec[i];
    }

    b[x_num-1] = ub ( x_min, x_max, t_min, t[j] );

    free ( fvec );

    job = 0;
    fvec = r83_np_sl ( x_num, a, b, job );

    for ( i = 0; i < x_num; i++ )
    {
      u[i+j*x_num] = fvec[i];
    }
  }

  r8mat_write ( x_file, 1, x_num, x );
  r8mat_write ( t_file, 1, t_num, t );
  r8mat_write ( u_file, x_num, t_num, u );

  printf ( "\n" );
  printf ( "  X data written to \"%s\".\n", x_file );
  printf ( "  T data written to \"%s\".\n", t_file );
  printf ( "  U data written to \"%s\".\n", u_file );

  free ( a );
  free ( b );
  free ( fvec );
  free ( t );
  free ( u );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FD1D_HEAT_IMPLICIT\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void f ( double a, double b, double t0, double t, int n, double x[], 
  double value[] )

/******************************************************************************/
/*
  Purpose:

    F returns the right hand side of the heat equation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the left and right endpoints.

    Input, double T0, the initial time.

    Input, double T, the current time.

    Input, int N, the number of points.

    Input, double X[N], the current spatial positions.

    Output, double VALUE[N], the prescribed value of U(X(:),T0).
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    value[i] = 0.0;
  }
  return;
}
/******************************************************************************/

int r83_np_fa ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83_NP_FA factors a R83 system without pivoting.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    Because this routine does not use pivoting, it can fail even when
    the matrix is not singular, and it is liable to make larger
    errors.

    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
    in one step, and does not save the factorization.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input/output, double A[3*N].
    On input, the tridiagonal matrix.  On output, factorization information.

    Output, int R83_NP_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
    if ( a[1+(i-1)*3] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R83_NP_FA - Fatal error!\n" );
      printf ( "  Zero pivot on step %d\n", i );
      return i;
    }
/*
  Store the multiplier in L.
*/
    a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];
/*
  Modify the diagonal entry in the next column.
*/
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
  }

  if ( a[1+(n-1)*3] == 0.0 )
  {
    printf ( "\n" );
    printf ( "R83_NP_FA - Fatal error!\n" );
    printf ( "  Zero pivot on step %d\n", n );
    return n;
  }

  return 0;
}
/******************************************************************************/

double *r83_np_sl ( int n, double a_lu[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R83_NP_SL solves a R83 system factored by R83_NP_FA.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A_LU[3*N], the LU factors from R83_NP_FA.

    Input, double B[N], the right hand side of the linear system.
    On output, B contains the solution of the linear system.

    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.

    Output, double R83_NP_SL[N], the solution of the linear system.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
/*
  Solve L * Y = B.
*/
    for ( i = 1; i < n; i++ )
    {
      x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];
    }
/*
  Solve U * X = Y.
*/
    for ( i = n; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( 1 < i )
      {
        x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
      }
    }
  }
  else
  {
/*
  Solve U' * Y = B
*/
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( i < n )
      {
        x[i] = x[i] - a_lu[0+i*3] * x[i-1];
      }
    }
/*
  Solve L' * X = Y.
*/
    for ( i = n-1; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] - a_lu[2+(i-1)*3] * x[i];
    }
  }

  return x;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the table data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16e", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void u0 ( double a, double b, double t0, int n, double x[], double value[] )

/******************************************************************************/
/*
  Purpose:

    U0 returns the initial condition at the starting time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the left and right endpoints

    Input, double T0, the initial time.

    Input, double T, the current time.

    Input, int N, the number of points where initial data is needed.

    Input, double X[N], the positions where initial data is needed.

    Output, double VALUE[N], the prescribed value of U(X,T0).
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    value[i] = 100.0;
  }
  return;
}
/******************************************************************************/

double ua ( double a, double b, double t0, double t )

/******************************************************************************/
/*
  Purpose:

    UA returns the Dirichlet boundary condition at the left endpoint.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the left and right endpoints

    Input, double T0, the initial time.

    Input, double T, the current time.

    Output, double UA, the prescribed value of U(A,T).
*/
{
  double value;

  value = 20.0;

  return value;
}
/******************************************************************************/

double ub ( double a, double b, double t0, double t )

/******************************************************************************/
/*
  Purpose:

    UB returns the Dirichlet boundary condition at the right endpoint.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the left and right endpoints

    Input, double T0, the initial time.

    Input, double T, the current time.

    Output, double UB, the prescribed value of U(B,T).
*/
{
  double value;

  value = 20.0;

  return value;
}
