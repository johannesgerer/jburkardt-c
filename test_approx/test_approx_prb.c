# include <stdlib.h>
# include <stdio.h>

# include "test_approx.h"
# include "spline.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( void );
void test11 ( void );
void test12 ( void );
void test13 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_APPROX_PRB.

  Discussion:

    TEST_APPROX_PRB calls the TEST_APPROX tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  timestamp (  );

  printf ( "\n" );
  printf ( "TEST_APPROX_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_APPROX library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_APPROX_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 shows how P00_TITLE can be called.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;
  char title[100];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Demonstrate some of the bookkeeping routines.\n" );
  printf ( "  P00_PROB_NUM returns the number of problems.\n" );
  printf ( "  P00_TITLE returns the problem title.\n" );
  printf ( "  P00_LIMIT returns the problem limits.\n" );

  prob_num = p00_prob_num ( );

  printf ( "\n" );
  printf ( "  Number of problems = %d\n", prob_num );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );
    printf ( "  %2d  \"%s\"\n", prob, title );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 shows how P00_STORY can be called.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  P00_STORY prints the problem \"story\".\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    p00_story ( prob );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses polynomial interpolation on data vector problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double *diftab;
  int i;
  int j;
  int jhi;
  char mark;
  int max_tab = 12;
  int ntab;
  int data_num;
  int prob;
  int prob_num;
  char title[100];
  double x;
  double *xdata;
  double yapprox;
  double *ydata;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Polynomial interpolation to a vector of data.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    printf ( "  %s\n", title );

    data_num = p00_data_num ( prob );

    printf ( "  DATA_NUM = %d\n", data_num );

    if ( max_tab < data_num )
    {
      printf ( "\n" );
      printf ( "  Skipped problem %d\n", prob );
      printf ( "  Too big.\n" );
    }
    else
    {
      xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
      ydata = ( double * ) malloc ( data_num * sizeof ( double ) );
      diftab = ( double * ) malloc ( data_num * sizeof ( double ) );

      p00_dat ( prob, data_num, xdata, ydata );

      ntab = data_num;

      printf ( "\n" );
      printf ( "  Interpolating polynomial order = %d\n", ntab );
      printf ( "\n" );
/*
  Construct the interpolating polynomial via finite differences.
*/
      data_to_dif ( ntab, xdata, ydata, diftab );
/*
  Print out the approximation, including midpoints of the intervals.
*/
      for ( i = 1; i <= ntab; i++ )
      {
        if ( i < ntab )
        {
          jhi = 2;
        }
        else
        {
          jhi = 1;
        }

        for ( j = 1; j <= jhi; j++ )
        {
          if ( i < ntab )
          {
            x = ( ( double ) ( jhi - j + 1 ) * xdata[i-1]   
                + ( double ) (       j - 1 ) * xdata[i] ) 
                / ( double ) ( jhi         );
          }
          else
          {
            x = xdata[ntab-1];
          }

          if ( j == 1 )
          {
            mark = '*';
          }
          else
          {
            mark = ' ';
          }

          yapprox = dif_val ( ntab, xdata, diftab, x );

          printf ( "  %c  %14g  %14g\n", mark, x, yapprox );
        }
      }
      free ( diftab );
      free ( xdata );
      free ( ydata );
    }
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 uses linear spline interpolation on all problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int imax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  char title[100];
  double *xdata;
  double xval;
  double *ydata;
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Linear spline interpolation.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );

    data_num = p00_data_num ( prob );
    xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
    ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    printf ( "  %s\n", title );
    printf ( "\n" );
    printf ( "       X          Y          Y'\n" );
    printf ( "\n" );
/*
  Evaluate the interpolation function.
*/
    imax = 2 * data_num - 1;

    for ( i = 1; i <= imax; i++ )
    {
      xval = ( ( double ) ( imax - i     ) * a   
             + ( double ) (        i - 1 ) * b ) 
             / ( double ) ( imax     - 1 );

      spline_linear_val ( data_num, xdata, ydata, xval, &yval, &ypval );

      if ( ( i % 2 ) == 1 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %14g  %14g  %14g\n", mark, xval, yval, ypval );
    }
    free ( xdata );
    free ( ydata );
  }

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 uses Overhauser spline interpolation on all problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int num_dim = 1;
  int prob;
  int prob_num;
  char title[100];
  double *xdata;
  double xval;
  double *ydata;
  double yval;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Overhauser spline interpolation.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );

    data_num = p00_data_num ( prob );

    xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
    ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    printf ( "  %s\n", title );
    printf ( "\n" );
    printf ( "  X   Y\n" );
    printf ( "\n" );
/*
  Evaluate the interpolation function.
*/
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]     
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, &yval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        printf ( "  %c  %14g  %14g\n", mark, xval, yval );
      }

    }
    free ( xdata );
    free ( ydata );
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 uses cubic spline interpolation on all problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  char title[100];
  double *xdata;
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Cubic spline interpolation.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );

    data_num = p00_data_num ( prob );

    xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
    ydata = ( double * ) malloc ( data_num * sizeof ( double ) );
    ypp = ( double * ) malloc ( data_num * sizeof ( double ) );

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];
/*
  Set up the interpolation function.
*/
    ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, 
      ibcend, ybcend );

    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    printf ( "  %s\n", title );
    printf ( "\n" );
    printf ( "    X   Y\n" );
    printf ( "\n" );
/*
  Evaluate the interpolation function.
*/
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, &yppval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        printf ( "  %c  %14g  %14g", mark, xval, yval );
      }
    }
    free ( xdata );
    free ( ydata );
    free ( ypp );
  }

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 plots an Overhauser spline interpolant for problem 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  char *approx_filename = "test07_approx.txt";
  char *data_filename = "test07_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int num_dim = 1;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Plot an Overhauser spline interpolant for problem 7.\n" );
/*
  Get the problem data.
*/
  prob = 7;

  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Evaluate the approximating function.
*/
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = ( double * ) malloc ( nplot * sizeof ( double ) );
  yplot = ( double * ) malloc ( nplot * sizeof ( double ) );

  plot = 0;

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1] 
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, &yval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( xplot );
  free ( ydata );
  free ( yplot );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 plots a cubic spline interpolant for problem 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  char *approx_filename = "test08_approx.txt";
  char *data_filename = "test08_data.txt";
  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double *yplot;
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Plot a cubic spline interpolant for problem 7.\n" );

  prob = 7;
/*
  Get the data.
*/
  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ypp = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Set up the interpolation function.
*/
  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, ybcend );
/*
  Evaluate the interpolation function.
*/
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = ( double * ) malloc ( nplot * sizeof ( double ) );
  yplot = ( double * ) malloc ( nplot * sizeof ( double ) );

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }
    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, 
        &yppval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( xplot );
  free ( ydata );
  free ( yplot );
  free ( ypp );

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 uses B spline approximation on all problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  int jhi;
  int jmax;
  char mark;
  int data_num;
  int prob;
  int prob_num;
  char title[100];
  double *xdata;
  double xval;
  double *ydata;
  double yval;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  B spline approximation.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    p00_title ( prob, title );

    data_num = p00_data_num ( prob );

    xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
    ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

    p00_dat ( prob, data_num, xdata, ydata );

    a = xdata[0];
    b = xdata[data_num-1];

    printf ( "\n" );
    printf ( "  Problem %d\n", prob );
    printf ( "  %s\n", title );
    printf ( "\n" );
    printf ( "       X        Y\n" );
    printf ( "\n" );
/*
  Evaluate the interpolation function.
*/
    for ( i = 1; i < data_num; i++ )
    {
      jmax = 3;

      if ( i == data_num - 1 )
      {
        jhi = jmax;
      }
      else
      {
        jhi = jmax - 1;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        xval = ( ( double ) ( jmax - j     ) * xdata[i-1]
               + ( double ) (        j - 1 ) * xdata[i] ) 
               / ( double ) ( jmax     - 1 );

        yval = spline_b_val ( data_num, xdata, ydata, xval );

        if ( j == 1 || j == 3 )
        {
          mark = '*';
        }
        else
        {
          mark = ' ';
        }
        printf ( "  %c  %14g  %14g\n", mark, xval, yval );
      }
    }
    free ( xdata );
    free ( ydata );
  }

  return;
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 plots a B spline approximant for problem 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  char *approx_filename = "test10_approx.txt";
  char *data_filename = "test10_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  char title[100];
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  Plot a B spline approximant for problem 7\n" );

  prob = 7;

  p00_title ( prob, title );
/*
  Get the data.
*/
  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Evaluate the approximation function.
*/
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = ( double * ) malloc ( nplot * sizeof ( double ) );
  yplot = ( double * ) malloc ( nplot * sizeof ( double ) );

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }
    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1]  
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_b_val ( data_num, xdata, ydata, xval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( xplot );
  free ( ydata );
  free ( yplot );

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 plots a beta spline approximant for problem 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  char *approx_filename = "test11_approx.txt";
  double beta1;
  double beta2;
  char *data_filename = "test11_data.txt";
  int i;
  int j;
  int jhi;
  int jmax = 7;
  int data_num;
  int nplot;
  int plot;
  int prob;
  char title[100];
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  beta1 = 100.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  Plot a beta spline approximant for problem 7\n" );
  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );

  prob = 7;

  p00_title ( prob, title );
/*
  Get the data.
*/
  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Evaluate the interpolation function.
*/
  plot = 0;
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1;

  xplot = ( double * ) malloc ( nplot * sizeof ( double ) );
  yplot = ( double * ) malloc ( nplot * sizeof ( double ) );

  for ( i = 1; i < data_num; i++ )
  {
    if ( i == data_num - 1 )
    {
      jhi = jmax;
    }
    else
    {
      jhi = jmax - 1;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      xval = ( ( double ) ( jmax - j     ) * xdata[i-1] 
             + ( double ) (        j - 1 ) * xdata[i] ) 
             / ( double ) ( jmax     - 1 );

      yval = spline_beta_val ( beta1, beta2, data_num, xdata, ydata, xval );

      xplot[plot] = xval;
      yplot[plot] = yval;
      plot = plot + 1;
    }
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( xplot );
  free ( ydata );
  free ( yplot );

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 plots a Bernstein spline approximant for problem 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
  double a;
  char *approx_filename = "test12_approx.txt";
  double b;
  char *data_filename = "test12_data.txt";
  int i;
  int data_num;
  int nplot = 101;
  int plot;
  int prob;
  double *xdata;
  double *xplot;
  double xval;
  double *ydata;
  double *yplot;
  double yval;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  Plot a Bernstein approximant for problem 5.\n" );
  printf ( "  Note that the Bernstein approximant requires equally\n" );
  printf ( "  spaced data!\n" );

  prob = 5;
/*
  Get the data.
*/
  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Evaluate the approximant function.
*/
  xplot = ( double * ) malloc ( nplot * sizeof ( double ) );
  yplot = ( double * ) malloc ( nplot * sizeof ( double ) );

  a = xdata[0];
  b = xdata[data_num-1];

  for ( plot = 1; plot <= nplot; plot++ )
  {
    xval = ( ( double ) ( nplot - plot     ) * a     
           + ( double ) (         plot - 1 ) * b ) 
           / ( double ) ( nplot        - 1 );

    yval = bpab_approx ( data_num - 1, a, b, ydata, xval );

    xplot[plot-1] = xval;
    yplot[plot-1] = yval;
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( xplot );
  free ( ydata );
  free ( yplot );

  return;
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 plots a cubic spline interpolant for problem 5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt
*/
{
# define NPLOT 101

  char *approx_filename = "test13_approx.txt";
  char *data_filename = "test13_data.txt";
  int ibcbeg;
  int ibcend;
  int j;
  int data_num;
  int nplot = NPLOT;
  int prob;
  double *xdata;
  double xplot[NPLOT];
  double xval;
  double ybcbeg;
  double ybcend;
  double *ydata;
  double yplot[NPLOT];
  double *ypp;
  double yppval;
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  Plot a cubic spline interpolant for problem 5\n" );

  prob = 5;

  data_num = p00_data_num ( prob );

  xdata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ydata = ( double * ) malloc ( data_num * sizeof ( double ) );
  ypp = ( double * ) malloc ( data_num * sizeof ( double ) );

  p00_dat ( prob, data_num, xdata, ydata );

  r8vec2_write ( data_filename, data_num, xdata, ydata );

  printf ( "\n" );
  printf ( "  Data values stored in \"%s\".\n", data_filename );
/*
  Set up the interpolation function.
*/
  ibcbeg = 0;
  ibcend = 0;
  ybcbeg = 0.0;
  ybcend = 0.0;

  ypp = spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, ybcend );
/*
  Evaluate the interpolation function.
*/
  for ( j = 1; j <= nplot; j++ )
  {
    xval = ( ( double ) ( nplot - j     ) * xdata[0]
           + ( double ) (         j - 1 ) * xdata[data_num-1] ) 
           / ( double ) ( nplot     - 1 );

    yval = spline_cubic_val ( data_num, xdata, ydata, ypp, xval, &ypval, 
      &yppval );

    xplot[j-1] = xval;
    yplot[j-1] = yval;
  }

  r8vec2_write ( approx_filename, nplot, xplot, yplot );

  printf ( "  Approximant values stored in \"%s\".\n", approx_filename ); 

  free ( xdata );
  free ( ydata );
  free ( ypp );

  return;
# undef NPLOT
}
