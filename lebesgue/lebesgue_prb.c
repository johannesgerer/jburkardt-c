# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "lebesgue.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LEBESGUE_PRB.

  Discussion:

    LEBESGUE_PRB tests the LEBESGUE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LEBESGUE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LEBESGUE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LEBESGUE_PRB\n" );
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

    LEBESGUE_TEST01 looks at Chebyshev1 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "chebyshev1";
  double *l;
  char label[] = "Chebyshev1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST01:\n" );
  printf ( "  Analyze Chebyshev1 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l,
    "  Chebyshev1 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = chebyshev1 ( n );
  r8vec_print ( n, x, "  Chebyshev1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST02 looks at Chebyshev2 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "chebyshev2";
  double *l;
  char label[] = "Chebyshev2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST02:\n" );
  printf ( "  Analyze Chebyshev2 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l,
    "  Chebyshev2 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = chebyshev2 ( n );
  r8vec_print ( n, x, "  Chebyshev2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST03 looks at Chebyshev3 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2014

  Author:

    John Burkardt
*/
{
  char filename[] = "chebyshev3";
  double *l;
  char label[] = "Chebyshev3 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST03:\n" );
  printf ( "  Analyze Chebyshev3 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev3 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l, 
    "  Chebyshev3 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = chebyshev3 ( n );
  r8vec_print ( n, x, "  Chebyshev3 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST04 looks at Chebyshev4 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "chebyshev4";
  double *l;
  char label[] = "Chebyshev4 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST04:\n" );
  printf ( "  Analyze Chebyshev4 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev4 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l, 
    "  Chebyshev4 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = chebyshev4 ( n );
  r8vec_print ( n, x, "  Chebyshev4 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST05 looks at Equidistant1 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "equidistant1";
  double *l;
  char label[] = "Equidistant1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST05:\n" );
  printf ( "  Analyze Equidistant1 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l, 
    "  Equidistant1 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = equidistant1 ( n );
  r8vec_print ( n, x, "  Equidistant1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST06 looks at Equidistant2 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "equidistant2";
  double *l;
  char label[] = "Equidistant2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST06:\n" );
  printf ( "  Analyze Equidistant2 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l, 
    "  Equidistant2 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = equidistant2 ( n );
  r8vec_print ( n, x, "  Equidistant2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST07 looks at Equidistant3 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{

  char filename[] = "equidistant3";
  double *l;
  char label[] = "Equidistant3 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST07:\n" );
  printf ( "  Analyze Equidistant3 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant3 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l,
    "  Equidistant3 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = equidistant3 ( n );
  r8vec_print ( n, x, "  Equidistant3 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST08 looks at Fejer 1 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "fejer1";
  double *l;
  char label[] = "Fejer1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST08:\n" );
  printf ( "  Analyze Fejer1 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = fejer1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l,
    "  Fejer1 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = fejer1 ( n );
  r8vec_print ( n, x, "  Fejer1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_TEST09 looks at Fejer2 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  char filename[] = "fejer2";
  double *l;
  char label[] = "Fejer2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  printf ( "\n" );
  printf ( "LEBESGUE_TEST09:\n" );
  printf ( "  Analyze Fejer2 points.\n" );

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = ( double * ) malloc ( n_max * sizeof ( double ) );

  for ( n = 1; n <= n_max; n++ )
  {
    x = fejer2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    free ( x );
  }

  r8vec_print ( n_max, l,
    "  Fejer2 Lebesgue constants for N = 1 to 11:" );
/*
  Examine one case more closely.
*/
  n = 11;
  x = fejer2 ( n );
  r8vec_print ( n, x, "  Fejer2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  free ( l );
  free ( x );
  free ( xfun );

  return;
}

