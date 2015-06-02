# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "tetrahedron_arbq_rule.h"

int main ( );
void test01 ( int degree, int n );
void test02 ( int degree, int n, char *header );
void test03 ( int degree, int n, char *header );
void test04 ( int degree, int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TETRAHEDRON_ARBQ_RULE_PRB.

  Discussion:

    TETRAHEDRON_ARBQ_RULE_PRB tests the TETRAHEDRON_ARBQ_RULE library.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    10 July 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.

  Reference:

    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature
    rules in two and higher dimensions,
    Computers and Mathematics with Applications,
    Volume 59, 2010, pages 663-676.
*/
{
  int degree;
  char header[255];
  int n;

  timestamp ( );
  printf ( "\n" );
  printf ( "TETRAHEDRON_ARBQ_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TETRAHEDRON_ARBQ_RULE library.\n" );

  degree = 8;
  n = tetrahedron_arbq_size ( degree );
  strcpy ( header, "tetrahedron08" );

  test01 ( degree, n );

  test02 ( degree, n, header );

  test03 ( degree, n, header );

  test04 ( degree, n );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TETRAHEDRON_ARBQ_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int degree, int n )

/******************************************************************************/
/*
  Purpose:

    TEST01 calls TETRAHEDRON_ARBQ for a quadrature rule of given order.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    08 July 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.

  Reference:

    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature
    rules in two and higher dimensions,
    Computers and Mathematics with Applications,
    Volume 59, 2010, pages 663-676.

  Parameters:

    Input, int DEGREE, the desired total polynomial degree exactness
    of the quadrature rule.  1 <= DEGREE <= 15.

    Input, int N, the number of nodes.
*/
{
  double d;
  int j;
  double volume;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Quadrature rule for a tetrahedron.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );

  volume = sqrt ( 8.0 ) / 3.0;
/*
  Retrieve and print a symmetric quadrature rule.
*/
  x = ( double * ) malloc ( 3 * n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  tetrahedron_arbq ( degree, n, x, w );

  printf ( "\n" );
  printf ( "  Number of nodes N = %d\n", n );

  printf ( "\n" );
  printf ( "     J  W       X       Y      Z\n" );
  printf ( "\n" );
  for ( j = 0; j < n; j++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g  %14.6g\n",
      j, w[j], x[0+j*3], x[1+j*3], x[2+j*3] );
  }

  d = r8vec_sum ( n, w );

  printf ( "   Sum    %g\n", d );
  printf ( "  Volume  %g\n", volume );

  return;
}
/******************************************************************************/

void test02 ( int degree, int n, char *header )

/******************************************************************************/
/*
  Purpose:

    TEST02 gets a rule and writes it to a file.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    10 July 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.

  Reference:

    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature
    rules in two and higher dimensions,
    Computers and Mathematics with Applications,
    Volume 59, 2010, pages 663-676.

  Parameters:

    Input, int DEGREE, the desired total polynomial degree exactness
    of the quadrature rule.  0 <= DEGREE <= 15.

    Input, int N, the number of nodes to be used by the rule.

    Input, char *HEADER, an identifier for the filenames.
*/
{
  int i;
  FILE *rule_unit;
  char rule_filename[255];
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Get a quadrature rule for the tetrahedron.\n" );
  printf ( "  Then write it to a file.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );
/*
  Retrieve a quadrature rule.
*/
  x = ( double * ) malloc ( 3 * n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  tetrahedron_arbq ( degree, n, x, w );
/*
  Write the points and weights to a file.
*/
  strcpy ( rule_filename, header );
  strcat ( rule_filename, ".txt" );

  rule_unit = fopen ( rule_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( rule_unit, "%g  %g  %g  %g\n", x[0+i*3], x[1+i*3], x[2+i*3], w[i] );
  }
  fclose ( rule_unit );
  printf ( "\n" );
  printf ( "  Quadrature rule written to file '%s'\n", rule_filename );

  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( int degree, int n, char *header )

/******************************************************************************/
/*
  Purpose:

    TEST03 gets a rule and creates GNUPLOT input files.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    08 July 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.

  Reference:

    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature
    rules in two and higher dimensions,
    Computers and Mathematics with Applications,
    Volume 59, 2010, pages 663-676.

  Parameters:

    Input, int DEGREE, the desired total polynomial degree exactness
    of the quadrature rule.  0 <= DEGREE <= 15.

    Input, int N, the number of nodes to be used by the rule.

    Input, char *HEADER, an identifier for the filenames.
*/
{
  int i;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Get a quadrature rule for the tetrahedron.\n" );
  printf ( "  Set up GNUPLOT graphics input.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );
/*
  Retrieve a symmetric quadrature rule.
*/
  x = ( double * ) malloc ( 3 * n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  tetrahedron_arbq ( degree, n, x, w );
/*
  Create files for input to GNUPLOT.
*/
  tetrahedron_arbq_gnuplot ( n, x, header );

  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test04 ( int degree, int n )

/******************************************************************************/
/*
  Purpose:

    TEST04 gets a rule and tests its accuracy.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    08 July 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.

  Reference:

    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature
    rules in two and higher dimensions,
    Computers and Mathematics with Applications,
    Volume 59, 2010, pages 663-676.

  Parameters:

    Input, int DEGREE, the desired total polynomial degree exactness
    of the quadrature rule.  0 <= DEGREE <= 15.

    Input, int N, the number of nodes to be used by the rule.
*/
{
  double d;
  int i;
  int j;
  int npols;
  double *pols;
  double *rints;
  double volume;
  double *w;
  double *x;
  double z[3];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Get a quadrature rule for the tetrahedron.\n" );
  printf ( "  Test its accuracy.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );
/*
  Retrieve a quadrature rule.
*/
  x = ( double * ) malloc ( 3 * n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  tetrahedron_arbq ( degree, n, x, w );

  npols = ( ( degree + 1 ) * ( degree + 2 ) * ( degree + 3 ) ) / 6;
  rints = ( double * ) malloc ( npols * sizeof ( double ) );

  for ( j = 0; j < npols; j++ )
  {
    rints[j] = 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    z[0] = x[0+i*3];
    z[1] = x[1+i*3];
    z[2] = x[2+i*3];

    pols = ortho3eva ( degree, z );

    for ( j = 0; j < npols; j++ )
    {
      rints[j] = rints[j] + w[i] * pols[j];
     }
    free ( pols );
  }

  volume = sqrt ( 8.0 ) / 3.0;

  d = 0.0;
  d = pow ( rints[0] - sqrt ( volume ), 2 );
  for ( i = 1; i < npols; i++ )
  {
    d = d + pow ( rints[i], 2 );
  }
  d = sqrt ( d ) / ( double ) ( npols );

  printf ( "\n" );
  printf ( "  RMS error = %g\n", d );

  free ( rints );
  free ( w );
  free ( x );

  return;
}
