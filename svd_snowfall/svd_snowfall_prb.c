# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "svd_snowfall.h"

int main ( );
void svd_snowfall_test02 ( int m, int n, double x[] );
void svd_snowfall_test03 ( int m, int n, double x[] );
void svd_snowfall_test04 ( int m, int n, double x[] );
void svd_snowfall_test05 ( int m, int n, double x[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SVD_SNOWFALL_PRB.

  Discussion:

    SVD_SNOWFALL_PRB tests the SVD_SNOWFALL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2013

  Author:

    John Burkardt
*/
{
  char *filename = "snowfall.txt";
  int m;
  int n;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "SVD_SNOWFALL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SVD_SNOWFALL library.\n" );
/*
  Retrieve the data.
  It's really easier to do this in the main program.
*/
  printf ( "\n" );
  printf ( "SVD_SNOWFALL_TEST01\n" );
  printf ( "  Read, process, and return snowfall data in \"%s\".\n", filename );
/*
  Determine the size of the data.
*/
  r8mat_header_read ( filename, &m, &n );

  printf ( "\n" );
  printf ( "  Number of data rows    M = %d\n", m );
  printf ( "  Number of data columns N = %d\n", n );

  x = r8mat_data_read ( filename, m, n );

  printf ( "\n" );
  printf ( "  Data has been read from the file.\n" );

  svd_snowfall_test02 ( m, n, x );
  svd_snowfall_test03 ( m, n, x );
  svd_snowfall_test04 ( m, n, x );
  svd_snowfall_test05 ( m, n, x );
/*
  Free memory.
*/
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SVD_SNOWFALL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void svd_snowfall_test02 ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    SVD_SNOWFALL_TEST02 looks at the singular values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double X[M*N], the snowfall data.
*/
{
  char command_filename[] = "singular_values_commands.txt";
  FILE *command;
  char data_filename[] = "singular_values_data.txt";
  FILE *data;
  double *e;
  double *e_cum;
  double e_sum;
  int i;
  int mn;
  double *s;
  double *s_diag;
  double *u;
  double *v;

  printf ( "\n" );
  printf ( "SVD_SNOWFALL_TEST02\n" );
  printf ( "  Look at the singular values.\n" );
  printf ( "  If the singular values are close, then the data is\n" );
  printf ( "  well spread out.  If the singular values decay rapidly,\n" );
  printf ( "  then the data exhibits patterns, or is constrained to\n" );
  printf ( "  a lower-dimensional subspace.\n" );
/*
  Compute the SVD.
*/
  u = ( double * ) malloc ( m * m * sizeof ( double ) );
  s = ( double * ) malloc ( m * n * sizeof ( double ) );
  v = ( double * ) malloc ( n * n * sizeof ( double ) );

  r8mat_svd_linpack ( m, n, x, u, s, v );
/*
  Extract the diagonal of S.
*/
  mn = i4_min ( m, n );
  s_diag = ( double * ) malloc ( mn * sizeof ( double ) );

  for ( i = 0; i < mn; i++ )
  {
    s_diag[i] = s[i+i*m];
  }
/*
  Print the singular values.
*/
  r8vec_print ( mn, s_diag, "  The singular values:" );
/*
  Plot the singular values.
*/
  data = fopen ( data_filename, "wt" );
  for ( i = 0; i < mn; i++ )
  {
    fprintf ( data, "  %4d  %14.6g\n", i, s_diag[i] );
  }
  fclose ( data );
  printf ( "\n" );
  printf ( "  Created data file \"%s\".\n", data_filename );

  command = fopen ( command_filename, "wt" );
  fprintf ( command, "# %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'singular_values.png'\n" );
  fprintf ( command, "set xlabel 'Index I'\n" );
  fprintf ( command, "set ylabel 'S(I)'\n" );
  fprintf ( command, "set title 'Snowfall Singular Values'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue'\n", 
    data_filename );
  fprintf ( command, "quit\n" );
  fclose ( command );
  printf ( "  Created command file \"%s\".\n", command_filename );
/*
  Print the cumulative "energy" of the singular values.
*/
  e = ( double * ) malloc ( mn * sizeof ( double ) );

  for ( i = 0; i < mn; i++ )
  {
    e[i] = pow ( s_diag[i], 2 );
  }
  e_sum = r8vec_sum ( mn, e );
  for ( i = 0; i < mn; i++ )
  {
    e[i] = e[i] / e_sum;
  }
  e_cum = r8vec_cum0_new ( mn, e );

  r8vec_print ( mn + 1, e_cum, "  The cumulative energy:" );

  free ( e );
  free ( e_cum );
  free ( s );
  free ( s_diag );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void svd_snowfall_test03 ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    SVD_SNOWFALL_TEST03 computes low rank approximations to the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double X[M*N], the snowfall data.
*/
{
  double *a1;
  double *a2;
  double *a3;
  double *a4;
  double *a5;
  char command_filename[] = "approx_commands.txt";
  FILE *command;
  char data_filename[] = "approx_data.txt";
  FILE *data;
  int i;
  double *s;
  double *u;
  double *v;

  printf ( "\n" );
  printf ( "SVD_SNOWFALL_TEST03\n" );
  printf ( "  Compute the rank 1 through rank 5 approximations to the data.\n" );
  printf ( "  Compare each of these to the 2012 snowfall data.\n" );
/*
  Compute the SVD.
*/
  u = ( double * ) malloc ( m * m * sizeof ( double ) );
  s = ( double * ) malloc ( m * n * sizeof ( double ) );
  v = ( double * ) malloc ( n * n * sizeof ( double ) );

  r8mat_svd_linpack ( m, n, x, u, s, v );
/*
  Form the rank 1, 2, 3, 4, 5 approximants to A.
*/
  a1 = r8mat_svd_low_rank ( m, n, 1, u, s, v );
  a2 = r8mat_svd_low_rank ( m, n, 2, u, s, v );
  a3 = r8mat_svd_low_rank ( m, n, 3, u, s, v );
  a4 = r8mat_svd_low_rank ( m, n, 4, u, s, v );
  a5 = r8mat_svd_low_rank ( m, n, 5, u, s, v );
/*
  Column 1 of X is the 2012 snowfall.
  Column 1 of A1 is the rank 1 approximant to 2012 snowfall.
*/
  data = fopen ( data_filename, "wt" );
  for ( i = 0; i < m; i++ )
  {
    fprintf ( data, "  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g\n",
      i, x[i+0*m], a1[i+0*m], a2[i+0*m], a3[i+0*m], a4[i+0*m], a5[i+0*m] );
  }
  fclose ( data );
  printf ( "  Created data file \"%s\".\n", data_filename );

  command = fopen ( command_filename, "wt" );
  fprintf ( command, "# %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'approx0.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title '2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue'\n", data_filename );

  fprintf ( command, "set output 'approx1.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Rank 1 Approx to 2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:3 lw 3 linecolor rgb 'red'\n", data_filename );

  fprintf ( command, "set output 'approx2.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Rank 2 Approx to 2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:3 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:4 lw 3 linecolor rgb 'red'\n", data_filename );

  fprintf ( command, "set output 'approx3.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Rank 3 Approx to 2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:3 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:4 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:5 lw 3 linecolor rgb 'red'\n", data_filename );

  fprintf ( command, "set output 'approx4.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Rank 4 Approx to 2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:3 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:4 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:5 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:6 lw 3 linecolor rgb 'red'\n", data_filename );

  fprintf ( command, "set output 'approx5.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Rank 5 Approx to 2012 Snowfall'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:3 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:4 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:5 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:6 lw 3 linecolor rgb 'gray',\\\n", data_filename );
  fprintf ( command, "     '%s' using 1:7 lw 3 linecolor rgb 'red'\n", data_filename );

  fprintf ( command, "quit\n" );
  fclose ( command );
  printf ( "  Created command file '%s'.\n", command_filename );

  free ( a1 );
  free ( a2 );
  free ( a3 );
  free ( a4 );
  free ( a5 );
  free ( s );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void svd_snowfall_test04 ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    SVD_SNOWFALL_TEST04 looks at the first 6 modes in the U matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double X[M*N], the snowfall data.
*/
{
  char command_filename[] = "umode_commands.txt";
  FILE *command;
  char data_filename[] = "umode_data.txt";
  FILE *data;
  int i;
  double *s;
  double *u;
  double *v;

  printf ( "\n" );
  printf ( "SVD_SNOWFALL_TEST04\n" );
  printf ( "  Look at the first 6 modes in the U matrix.\n" );
  printf ( "  Each of these represents a pattern for snowfall over a year.\n" );
  printf ( "  The first mode is the pattern that is strongest in the data.\n" );
/*
  Compute the SVD.
*/
  u = ( double * ) malloc ( m * m * sizeof ( double ) );
  s = ( double * ) malloc ( m * n * sizeof ( double ) );
  v = ( double * ) malloc ( n * n * sizeof ( double ) );

  r8mat_svd_linpack ( m, n, x, u, s, v );
/*
  Normalize the patterns so that each column has maximum entry 1.
*/
  r8col_normalize_li ( m, m, u );
/*
  Plot the U modes.
*/
  data = fopen ( data_filename, "wt" );
  for ( i = 0; i < m; i++ )
  {
    fprintf ( data, "  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g\n",
    i, u[i+0*m], u[i+1*m], u[i+2*m], u[i+3*m], u[i+4*m], u[i+5*m] );
  }
  fclose ( data );
  printf ( "  Created data file \"%s\".\n", data_filename );

  command = fopen ( command_filename, "wt" );
  fprintf ( command, "# %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'umode1.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 1'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:2 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "set output 'umode2.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 2'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:3 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "set output 'umode3.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 3'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:4 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "set output 'umode4.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 4'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:5 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "set output 'umode5.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 5'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:6 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "set output 'umode6.png'\n" );
  fprintf ( command, "set xlabel 'Month'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Monthly Snowfall Mode 6'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "plot '%s' using 1:7 lw 3 linecolor rgb 'blue'\n", 
    data_filename );

  fprintf ( command, "quit\n" );
  fclose ( command );
  printf ( "  Created command file '%s'.\n", command_filename );

  free ( s );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void svd_snowfall_test05 ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    SVD_SNOWFALL_TEST05 looks at the first 6 modes in the V matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double X[M*N], the snowfall data.
*/
{
  char command_filename[] = "vmode_commands.txt";
  FILE *command;
  char data_filename[] = "vmode_data.txt";
  FILE *data;
  int i;
  double *s;
  double *u;
  double *v;

  printf ( "\n" );
  printf ( "SVD_SNOWFALL_TEST05\n" );
  printf ( "  Look at the first 6 modes in the V matrix.\n" );
  printf ( "  Each of these represents a pattern shared by all the months,\n" );
  printf ( "  and extending across the 123 sampling years.\n" );
/*
  Compute the SVD.
*/
  u = ( double * ) malloc ( m * m * sizeof ( double ) );
  s = ( double * ) malloc ( m * n * sizeof ( double ) );
  v = ( double * ) malloc ( n * n * sizeof ( double ) );

  r8mat_svd_linpack ( m, n, x, u, s, v );
/*
  Normalize the patterns so that each column has maximum entry 1.
*/
  r8col_normalize_li ( n, n, v );
/*
  Reverse the row ordering.
*/
  r8row_reverse ( n, n, v );
/*
  Plot the V modes.
*/
  data = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( data, "  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g\n",
    i, v[i+0*n], v[i+1*n], v[i+2*n], v[i+3*n], v[i+4*n], v[i+5*n] );
  }
  fclose ( data );
  printf ( "  Created data file \"%s\".\n", data_filename );

  command = fopen ( command_filename, "wt" );
  fprintf ( command, "# %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'vmode1.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 1'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:2 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "set output 'vmode2.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 2'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:3 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "set output 'vmode3.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 3'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:4 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "set output 'vmode4.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 4'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:5 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "set output 'vmode5.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 5'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:6 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "set output 'vmode6.png'\n" );
  fprintf ( command, "set xlabel 'Year'\n" );
  fprintf ( command, "set ylabel 'Snowfall'\n" );
  fprintf ( command, "set title 'Yearly Snowfall Mode 6'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style data points\n" );
  fprintf ( command, "plot '%s' using 1:7 with points lt 3 pt 3\n",
    data_filename );

  fprintf ( command, "quit\n" );
  fclose ( command );
  printf ( "  Created command file '%s'.\n", command_filename );

  free ( s );
  free ( u );
  free ( v );

  return;
}
