# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "test_opt_con.h"

/******************************************************************************/

void p00_ab ( int problem, int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P00_AB returns bounds for a problem specified by index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  if ( problem == 1 )
  {
    p01_ab ( m, a, b );
  }
  else if ( problem == 2 )
  {
    p02_ab ( m, a, b );
  }
  else if ( problem == 3 )
  {
    p03_ab ( m, a, b );
  }
  else if ( problem == 4 )
  {
    p04_ab ( m, a, b );
  }
  else if ( problem == 5 )
  {
    p05_ab ( m, a, b );
  }
  else if ( problem == 6 )
  {
    p06_ab ( m, a, b );
  }
  else if ( problem == 7 )
  {
    p07_ab ( m, a, b );
  }
  else if ( problem == 8 )
  {
    p08_ab ( m, a, b );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "P00_AB - Fatal error!\n" );
    fprintf ( stderr, "  Problem index out of bounds.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

double *p00_f ( int problem, int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P00_F returns the objective function value for a problem specified by index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P00_F[N], the function values.
//
{
  double *f;

  if ( problem == 1 )
  {
    f = p01_f ( m, n, x );
  }
  else if ( problem == 2 )
  {
    f = p02_f ( m, n, x );
  }
  else if ( problem == 3 )
  {
    f = p03_f ( m, n, x );
  }
  else if ( problem == 4 )
  {
    f = p04_f ( m, n, x );
  }
  else if ( problem == 5 )
  {
    f = p05_f ( m, n, x );
  }
  else if ( problem == 6 )
  {
    f = p06_f ( m, n, x );
  }
  else if ( problem == 7 )
  {
    f = p07_f ( m, n, x );
  }
  else if ( problem == 8 )
  {
    f = p08_f ( m, n, x );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "P00_F - Fatal error!\n" );
    fprintf ( stderr, "  Problem index out of bounds.\n" );
    exit ( 1 );
  }
  return f;
}
/******************************************************************************/

int p00_m ( int problem )

/******************************************************************************/
//
//  Purpose:
//
//    P00_M returns the spatial dimension for a problem specified by index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, int P00_M, the spatial dimension.
//
{
  int m;

  if ( problem == 1 )
  {
    m = p01_m ( );
  }
  else if ( problem == 2 )
  {
    m = p02_m ( );
  }
  else if ( problem == 3 )
  {
    m = p03_m ( );
  }
  else if ( problem == 4 )
  {
    m = p04_m ( );
  }
  else if ( problem == 5 )
  {
    m = p05_m ( );
  }
  else if ( problem == 6 )
  {
    m = p06_m ( );
  }
  else if ( problem == 7 )
  {
    m = p07_m ( );
  }
  else if ( problem == 8 )
  {
    m = p08_m ( );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "P00_M - Fatal error!\n" );
    fprintf ( stderr, "  Problem index out of bounds.\n" );
    exit ( 1 );
  }
  return m;
}
/******************************************************************************/

int p00_problem_num ( )

/******************************************************************************/
//
//  Purpose:
//
//    P00_PROBLEM_NUM returns the number of problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROBLEM_NUM, the number of defined problems.
//
{
  int problem_num;

  problem_num = 8;

  return problem_num;
}
/******************************************************************************/

double *p00_sol ( int problem, int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P00_SOL returns known solutions for a problem specified by index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Input, int M, the order of the problem.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P00_SOL[M], the solution.
//
{
  double *x;

  if ( problem == 1 )
  {
    x = p01_sol ( m, know );
  }
  else if ( problem == 2 )
  {
    x = p02_sol ( m, know );
  }
  else if ( problem == 3 )
  {
    x = p03_sol ( m, know );
  }
  else if ( problem == 4 )
  {
    x = p04_sol ( m, know );
  }
  else if ( problem == 5 )
  {
    x = p05_sol ( m, know );
  }
  else if ( problem == 6 )
  {
    x = p06_sol ( m, know );
  }
  else if ( problem == 7 )
  {
    x = p07_sol ( m, know );
  }
  else if ( problem == 8 )
  {
    x = p08_sol ( m, know );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "P00_SOL - Fatal error!\n" );
    fprintf ( stderr, "  Problem index out of bounds.\n" );
    exit ( 1 );
  }
  return x;
}
/******************************************************************************/

void p00_title ( int problem, char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P00_TITLE returns a title for a problem specified by index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROBLEM, the problem index.
//
//    Output, char *TITLE, a title for the problem.
//
{
  if ( problem == 1 )
  {
    p01_title ( title );
  }
  else if ( problem == 2 )
  {
    p02_title ( title );
  }
  else if ( problem == 3 )
  {
    p03_title ( title );
  }
  else if ( problem == 4 )
  {
    p04_title ( title );
  }
  else if ( problem == 5 )
  {
    p05_title ( title );
  }
  else if ( problem == 6 )
  {
    p06_title ( title );
  }
  else if ( problem == 7 )
  {
    p07_title ( title );
  }
  else if ( problem == 8 )
  {
    p08_title ( title );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "P00_TITLE - Fatal error!\n" );
    fprintf ( stderr, "  Problem number out of bounds.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void p01_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P01_AB returns bounds for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p01_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P01_F returns the objective function value for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P01_F[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double p;
  double s;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    p = 1.0;
    s = 0.0;
    for ( i = 0; i < m; i++ )
    {
      p = p * x[i+j*m];
      s = s + x[i+j*m];
    }
    f[j] = - exp ( p ) * sin ( s );
  }
  return f;
}
/******************************************************************************/

int p01_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P01_M returns the spatial dimension for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P01_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p01_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P01_SOL returns known solutions for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P01_SOL[M], the solution.
//
{
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    x[0] = 0.409887209247642;
    x[1] = 0.409887209247642;
    x[2] = 0.409887209247642;
    x[3] = 0.409887209247642;
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p01_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P01_TITLE returns a title for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = - exp(prod(x)) * sin(sum(x))." );

  return;
}
/******************************************************************************/

void p02_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P02_AB returns bounds for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p02_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P02_F returns the objective function value for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P02_F[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double p;
  double s;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    p = x[0+j*m] * pow ( x[1+j*m], 2 ) * pow ( x[2+j*m], 3 ) * pow ( x[3+j*m], 4 );
    s = 0.0;
    for ( i = 0; i < m; i++ )
    {
      s = s + x[i+j*m];
    }
    f[j] = - exp ( p ) * sin ( s );
  }
  return f;
}
/******************************************************************************/

int p02_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P02_M returns the spatial dimension for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P02_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p02_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P02_SOL returns known solutions for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P02_SOL[M], the solution.
//
{
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    x[0] = 0.390500591228663;
    x[1] = 0.392051909813608;
    x[2] = 0.393601661544812;
    x[3] = 0.395149843840982;
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p02_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P02_TITLE returns a title for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = - exp(x(1)*x(2)^2*x(3)^3*x(4)^4) * sin(sum(x))." );

  return;
}
/******************************************************************************/

void p03_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P03_AB returns bounds for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p03_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P03_F returns the objective function value for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P03_F[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double p;
  double s;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    s = - x[0+j*m] - 2.0 * x[1+j*m] - 3.0 * x[2+j*m] - 4.0 * x[3+j*m];
    p = 1.0;
    for ( i = 0; i < m; i++ )
    {
      p = p * x[i+j*m];
    }
    f[j] = - 10000.0 * p * exp ( s );
  }
  return f;
}
/******************************************************************************/

int p03_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P03_M returns the spatial dimension for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P03_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p03_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P03_SOL returns known solutions for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P03_SOL[M], the solution.
//
{
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    x[0] = 0.999980569087140;
    x[1] = 0.500000721280566;
    x[2] = 0.333341891834645;
    x[3] = 0.249997266604697;
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p03_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P03_TITLE returns a title for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = -1000 * product(x) * exp(-x(1)-2x(2)-3x(3)-4x(4))." );

  return;
}
/******************************************************************************/

void p04_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P04_AB returns bounds for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p04_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P04_F returns the objective function value for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P04_F[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double p;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    p = 1.0;
    for ( i = 0; i < m; i++ )
    {
      p = p * x[i+j*m];
    }
    f[j] = - 100.0 * p * exp ( - x[3+j*m] ) 
      / pow ( 1.0 + x[0+j*m] * x[1+j*m] * x[2+j*m], 2 );
  }
  return f;
}
/******************************************************************************/

int p04_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P04_M returns the spatial dimension for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P04_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p04_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P04_SOL returns known solutions for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P04_SOL[M], the solution.
//
{
  int i;
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    for ( i = 0; i < m; i++ )
    {
      x[i] = 1.0;
    }
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p04_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P04_TITLE returns a title for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = -100 * product(x) * exp(-x(4)) / (1+x(1)+x(2)+x(3))." );

  return;
}
/******************************************************************************/

void p05_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P05_AB returns bounds for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p05_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P05_F returns the objective function value for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P05_F[N], the function values.
//
{
  double *f;
  int j;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    f[j] = pow ( x[0+j*m] -  3.0 / 11.0, 2 )
         + pow ( x[1+j*m] -  6.0 / 13.0, 2 )
         + pow ( x[2+j*m] - 12.0 / 23.0, 4 )
         + pow ( x[3+j*m] -  8.0 / 37.0, 6 );
  }
  return f;
}
/******************************************************************************/

int p05_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P05_M returns the spatial dimension for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P05_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p05_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P05_SOL returns known solutions for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P05_SOL[M], the solution.
//
{
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    x[0] =  3.0 / 11.0;
    x[1] =  6.0 / 13.0;
    x[2] = 12.0 / 23.0;
    x[3] =  8.0 / 37.0;
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p05_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P05_TITLE returns a title for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = (x(1)-3/11)^2+(x(2)-6/13)^2+(x(3)-12/23)^4+(x(4)-8/37)^6" );

  return;
}
/******************************************************************************/

void p06_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P06_AB returns bounds for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
/******************************************************************************/

double *p06_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P06_F returns the objective function value for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P06_F[N], the function values.
//
{
  double arg;
  double *f;
  int j;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    arg = 
        1.0 / x[0+j*m]
      + 1.0 / x[1+j*m]
      + 1.0 / x[2+j*m]
      + 1.0 / x[3+j*m];
    f[j] = - sin ( arg );
  }
  return f;
}
/******************************************************************************/

int p06_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P06_M returns the spatial dimension for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, int P06_M, the spatial dimension.
//
{
  int m;

  m = 4;

  return m;
}
/******************************************************************************/

double *p06_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P06_SOL returns known solutions for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int *KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P06_SOL[M], the solution.
//
{
  double *x;

  if ( *know == 0 )
  {
    *know = 1;
    x = ( double * ) malloc ( m * sizeof ( double ) );
    x[0] = 0.509282516910744;
    x[1] = 0.509282516910744;
    x[2] = 0.509282516910746;
    x[3] = 0.509282516910744;
  }
  else
  {
    *know = 0;
    x = NULL;
  }
  return x;
}
/******************************************************************************/

void p06_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P06_TITLE returns a title for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Harald Niederreiter, Kevin McCurley,
//    Optimization of functions by quasi-random search methods,
//    Computing,
//    Volume 22, Number 2, 1979, pages 119-123.
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = - sin(1/x(1)+1/x(2)+1/x(3)+1/x(4))" );

  return;
}
/******************************************************************************/

void p07_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P07_AB returns bounds for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 10.0;
  }
  return;
}
/******************************************************************************/

double *p07_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P07_F returns the objective function value for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Langerman10 reference?
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P07_F[N], the function values.
//
{
  double a[2*5] = {
    3.0, 5.0, 
    5.0, 2.0, 
    2.0, 1.0, 
    1.0, 4.0, 
    7.0, 9.0 };
  double arg;
  double c[5] = { 1.0, 2.0, 5.0, 2.0, 3.0 };
  double *f;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( k = 0; k < 5; k++ )
    {
      arg = 0.0;
      for ( i = 0; i < m; i++ )
      {
        arg = arg + pow ( x[i+j*m] - a[i+k*m], 2 );
      }
      f[j] = f[j] - c[k] * exp ( - arg / pi ) * cos ( pi * arg );
    }
  }
  return f;
}
/******************************************************************************/

int p07_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P07_M returns the spatial dimension for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P07_M, the spatial dimension.
//
{
  int m;

  m = 2;

  return m;
}
/******************************************************************************/

double *p07_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P07_SOL returns known solutions for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P07_SOL[M], the solution.
//
{
  double *x;

  *know = 0;
  x = NULL;

  return x;
}
/******************************************************************************/

void p07_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P07_TITLE returns a title for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = Langerman2 function" );

  return;
}
/******************************************************************************/

void p08_ab ( int m, double a[], double b[] )

/******************************************************************************/
//
//  Purpose:
//
//    P08_AB returns bounds for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 10.0;
  }
  return;
}
/******************************************************************************/

double *p08_f ( int m, int n, double x[] )

/******************************************************************************/
//
//  Purpose:
//
//    P08_F returns the objective function value for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Langerman10 reference?
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P08_F[N], the function values.
//
{
  double a[10*30] = {
    9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020, 
    9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374, 
    8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982, 
    2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426, 
    8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567, 
    7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208, 
    1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448, 
    8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762, 
    0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637, 
    7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247, 
    0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016, 
    2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789, 
    8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109, 
    2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564, 
    4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670, 
    8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826, 
    8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591, 
    4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740, 
    2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675, 
    6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258, 
    0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070, 
    5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234, 
    3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027, 
    8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064, 
    1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224, 
    0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644, 
    0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229, 
    4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506, 
    9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732, 
    4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500 };
  double arg;
  double c[30] = {
    0.806, 0.517, 1.500, 0.908, 0.965, 
    0.669, 0.524, 0.902, 0.531, 0.876, 
    0.462, 0.491, 0.463, 0.714, 0.352, 
    0.869, 0.813, 0.811, 0.828, 0.964, 
    0.789, 0.360, 0.369, 0.992, 0.332, 
    0.817, 0.632, 0.883, 0.608, 0.326 };
  double *f;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( k = 0; k < 30; k++ )
    {
      arg = 0.0;
      for ( i = 0; i < m; i++ )
      {
        arg = arg + pow ( x[i+j*m] - a[i+k*m], 2 );
      }
      f[j] = f[j] - c[k] * exp ( - arg / pi ) * cos ( pi * arg );
    }
  }
  return f;
}
/******************************************************************************/

int p08_m ( )

/******************************************************************************/
//
//  Purpose:
//
//    P08_M returns the spatial dimension for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P08_M, the spatial dimension.
//
{
  int m;

  m = 10;

  return m;
}
/******************************************************************************/

double *p08_sol ( int m, int *know )

/******************************************************************************/
//
//  Purpose:
//
//    P08_SOL returns known solutions for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P08_SOL[M], the solution.
//
{
  double *x;

  *know = 0;
  x = NULL;

  return x;
}
/******************************************************************************/

void p08_title ( char *title )

/******************************************************************************/
//
//  Purpose:
//
//    P08_TITLE returns a title for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P08_TITLE, a title for the problem.
//
{
  strcpy ( title, "f(x) = Langerman10 function" );

  return;
}
/******************************************************************************/

double *r8col_uniform_new ( int m, int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIFORM_NEW fills an R8COL with scaled pseudorandom numbers.

  Discussion:

    An R8COL is an array of R8 values, regarded as a set of column vectors.

    The user specifies a minimum and maximum value for each row.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 December 2011

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

    Input, double A[M], B[M], the upper and lower limits.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8COL_UNIFORM_NEW[M*N], a matrix of pseudorandom values.
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
      r[i+j*m] = a[i] 
        + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  double *r8vec_pointer;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
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

