# include <stdlib.h>
# include <stdio.h>
# include <math.h>

int main ( );
void jacobi_test01 ( );

# include "gnuplot_i.h"
# include "jacobi.h"

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for JACOBI_PRB.

  Discussion:

    JACOBI_PRB tests the JACOBI library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "JACOBI_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the JACOBI library.\n" );

  jacobi_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "JACOBI_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void jacobi_test01 ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_TEST01 tests JACOBI1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int i;
  int it;
  int it_num;
  int j;
  double *m_plot;
  double *ml_plot;
  int n;
  double *r_plot;
  double *rl_plot;
  double *s_plot;
  double t;
  gnuplot_ctrl *window1;
  gnuplot_ctrl *window2;
  double *x;
  double *x_exact;
  double *x_new;
  double *x_plot;

  double w = 0.5;

  printf ( "\n" );
  printf ( "JACOBI_TEST01:\n" );

  it_num = 2000;
  n = 33;
/*
  Set the matrix A.
*/
  a = dif2 ( n, n );
/*
  Determine the right hand side vector B.
*/
  x_exact = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    t = ( double ) i / ( double ) ( n - 1 );
    x_exact[i] = exp ( t ) * ( t - 1 ) * t;
//  x_exact[i] = ( double ) ( i + 1 );
  }
  b = r8mat_mv_new ( n, n, a, x_exact );
//
//  Set the initial estimate for the solution.
//
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
/*
  Allocate plot arrays.
*/
  m_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );
  r_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );
  s_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );
  x_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );
/*
  Initialize plot arrays.
*/
  r_plot[0] = r8mat_residual_norm ( n, n, a, x, b );
  m_plot[0] = 1.0;
  for ( i = 0; i < n; i++ )
  {
    x_plot[i+0*n] = x[i];
  }
  for ( j = 0; j <= it_num; j++ )
  {
    s_plot[j] = ( double ) j;
  }
/*
  Carry out the iteration.
*/
  for ( it = 1; it <= it_num; it++ )
  {
    x_new = jacobi1 ( n, a, b, x );

    r_plot[it] = r8mat_residual_norm ( n, n, a, x_new, b );
/*
  Compute the average point motion.
*/
    m_plot[it] = r8vec_diff_norm_squared ( n, x, x_new ) / ( double ) n;
/*
  Update the solution
*/
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( 1.0 - w ) * x[i] + w * x_new[i];
    }
/*
    r8vec_copy ( n, x_new, x );
*/
    for ( i = 0; i < n; i++ )
    {
      x_plot[i+0*n] = x[i];
    }

    free ( x_new );
  }
    r8vec_print ( n, x, "Solution" );
/*
  Plot the residual.
*/
  rl_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );

  for ( j = 0; j <= it_num; j++ )
  {
    rl_plot[j] = log ( r_plot[j] );
  }
  window1 = gnuplot_init ( );
  gnuplot_setstyle ( window1, "lines" );
  gnuplot_cmd ( window1, "set grid" );
  gnuplot_set_xlabel ( window1, "Step" );
  gnuplot_set_ylabel ( window1, "Log ( Residual )" );
  gnuplot_setstyle ( window1, "lines" );
  gnuplot_plot1d_var2v ( window1, s_plot, rl_plot, it_num + 1, "Log(Residual)" );
  sleep ( 1 );
//
//  Plot the average point motion.
//
  ml_plot = ( double * ) malloc ( ( it_num + 1 ) * sizeof ( double ) );
  for ( j = 0; j <= it_num; j++ )
  {
    ml_plot[j] = log ( m_plot[j] );
//  printf ( " j << "  " << ml_plot[j] << "\n" );
  }
  window2 = gnuplot_init ( );
  gnuplot_setstyle ( window2, "lines" );
  gnuplot_cmd ( window2, "set grid" );
  gnuplot_set_xlabel ( window2, "Step" );
  gnuplot_set_ylabel ( window2, "Log ( Motion )" );
  gnuplot_setstyle ( window2, "lines" );
  gnuplot_plot1d_var2v ( window2, s_plot, ml_plot, it_num + 1, "Log(Motion)" );
  sleep ( 10 );
//
//  Plot the evolution of the locations of the generators.
//
//figure ( 3 )

//y = ( 0 : it_num );
//for k = 1 : n
//  plot ( x_plot(k,1:it_num+1), y )
//  hold on;
//end
//grid on
//hold off;

//title ( "Generator evolution." );
//xlabel ( "Generator positions" );
//ylabel ( "Iterations" ); 

  gnuplot_close ( window1 );
  gnuplot_close ( window2 );

  free ( a );
  free ( b );
  free ( m_plot );
  free ( ml_plot );
  free ( r_plot );
  free ( rl_plot );
  free ( s_plot );
  free ( x );
  free ( x_exact );
  free ( x_plot );

  return;
}
