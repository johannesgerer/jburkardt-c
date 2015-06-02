# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "triangle_symq_rule.h"

int main ( );
void test01 ( );
void test02 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] );
void test03 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], char *header );
void test04 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], char *header );
void test05 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_SYMQ_RULE_PRB.

  Discussion:

    TRIANGLE_SYMQ_RULE_PRB tests the TRIANGLE_SYMQ_RULE library.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    28 June 2014

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
  int itype;
  int numnodes;
  double vert1[2];
  double vert2[2];
  double vert3[2];

  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_SYMQ_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_SYMQ_RULE library.\n" );

  test01 ( );

  for ( itype = 0; itype <= 2; itype++ )
  {
    if ( itype == 0 )
    {
      printf ( "\n" );
      printf ( "  Region is user-defined triangle.\n" );
      vert1[0] = 1.0;
      vert1[1] = 0.0;
      vert2[0] = 4.0;
      vert2[1] = 4.0;
      vert3[0] = 0.0;
      vert3[1] = 3.0;
      strcpy ( header, "user08" );
      degree = 8;
    }
    else if ( itype == 1 )
    {
      printf ( "\n" );
      printf ( "  Region is standard equilateral triangle.\n" );
      vert1[0] = -1.0;
      vert1[1] = -1.0 / sqrt ( 3.0 );
      vert2[0] = +1.0;
      vert2[1] = -1.0 / sqrt ( 3.0 );
      vert3[0] =  0.0;
      vert3[1] =  2.0 / sqrt ( 3.0 );
      strcpy ( header, "equi08" );
      degree = 8;
    }
    else if ( itype == 2 )
    {
      printf ( "\n" );
      printf ( "  Region is the simplex (0,0),(1,0),(0,1).\n" );
      vert1[0] = 0.0;
      vert1[1] = 0.0;
      vert2[0] = 1.0;
      vert2[1] = 0.0;
      vert3[0] = 0.0;
      vert3[1] = 1.0;
      strcpy ( header, "simp08" );
      degree = 8;
    }

    printf ( "\n" );
    printf ( "  Triangle:\n" );
    printf ( "\n" );
    printf ( "  %14.6g  %14.6g\n", vert1[0], vert1[1] );
    printf ( "  %14.6g  %14.6g\n", vert2[0], vert2[1] );
    printf ( "  %14.6g  %14.6g\n", vert3[0], vert3[1] );
/*
  Determine the size of the rule.
*/
    numnodes = rule_full_size ( degree );
/*
  Retrieve a rule and print it.
*/
    test02 ( degree, numnodes, vert1, vert2, vert3 );
/*
  Get a rule, and write data files that gnuplot can use to plot the points.
*/
    test03 ( degree, numnodes, vert1, vert2, vert3, header );

    test04 ( degree, numnodes, vert1, vert2, vert3, header );

    test05 ( degree, numnodes, vert1, vert2, vert3 );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_SYMQ_RULE_PRB\n" );
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

    TEST01 tests TRIANGLE_TO_SIMPLEX, TRIANGLE_TO_REF, REF_TO_TRIANGLE, SIMPLEX_TO_TRIANGLE.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    30 June 2014

  Author:

    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
    This C version by John Burkardt.
*/
{
  int i;
  double *rp1;
  double rv1[2];
  double rv2[2];
  double rv3[2];
  int seed;
  double *sp1;
  double *sp2;
  double sv1[2];
  double sv2[2];
  double sv3[2];
  double *tp1;
  double *tp2;
  double tv1[2];
  double tv2[2];
  double tv3[2];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Map points from one triangle to another.\n" );
  printf ( "\n" );
  printf ( "  R = reference triangle\n" );
  printf ( "  S = simplex\n" );
  printf ( "  T = user-defined triangle.\n" );
  printf ( "  REF_TO_TRIANGLE:     R => T\n" );
  printf ( "  SIMPLEX_TO_TRIANGLE: S => T\n" );
  printf ( "  TRIANGLE_TO_REF:     T => R\n" );
  printf ( "  TRIANGLE_TO_SIMPLEX: T => S\n" );
/*
  Reference triangle
*/
  rv1[0] = -1.0;
  rv1[1] = -1.0 / sqrt ( 3.0 );
  rv2[0] = +1.0;
  rv2[1] = -1.0 / sqrt ( 3.0 );
  rv3[0] =  0.0;
  rv3[1] =  2.0 / sqrt ( 3.0 );
/*
  Simplex
*/
  sv1[0] = 0.0;
  sv1[1] = 0.0;
  sv2[0] = 1.0;
  sv2[1] = 0.0;
  sv3[0] = 0.0;
  sv3[1] = 1.0;
/*
  User triangle.
*/
  tv1[0] = 1.0;
  tv1[1] = 0.0;
  tv2[0] = 4.0;
  tv2[1] = 4.0;
  tv3[0] = 0.0;
  tv3[1] = 3.0;

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    sp1 = r8vec_uniform_01_new ( 2, &seed );

    if ( 1.0 < sp1[0] + sp1[1] )
    {
      sp1[0] = 1.0 - sp1[0];
      sp1[1] = 1.0 - sp1[1];
    }

    tp1 = simplex_to_triangle ( tv1, tv2, tv3, sp1 );
    rp1 = triangle_to_ref ( tv1, tv2, tv3, tp1 );
    tp2 = ref_to_triangle ( tv1, tv2, tv3, rp1 );
    sp2 = triangle_to_simplex ( tv1, tv2, tv3, tp2 );

    printf ( "\n" );
    printf ( "  SP1: %14.6g  %14.6g\n", sp1[0], sp1[1] );
    printf ( "  TP1: %14.6g  %14.6g\n", tp1[0], tp1[1] );
    printf ( "  RP1: %14.6g  %14.6g\n", rp1[0], rp1[1] );
    printf ( "  TP2: %14.6g  %14.6g\n", tp2[0], tp2[1] );
    printf ( "  SP2: %14.6g  %14.6g\n", sp2[0], sp2[1] );

    free ( rp1 );
    free ( sp1 );
    free ( sp2 );
    free ( tp1 );
    free ( tp2 );
  }

  return;
}
/******************************************************************************/

void test02 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] )

/******************************************************************************/
/*
  Purpose:

    TEST02 calls TRIASYMQ for a quadrature rule of given order and region.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    30 June 2014

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
    of the quadrature rule.  0 <= DEGREE <= 50.

    Input, int NUMNODES, the number of nodes to be used by the rule.

    Input, double VERT1[2], VERT2[2], VERT3[2], the
    vertices of the triangle.
*/
{
  double area;
  double d;
  int j;
  double *rnodes;
  double *weights;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Symmetric quadrature rule for a triangle.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );

  area = triangle_area ( vert1, vert2, vert3 );
/*
  Retrieve and print a symmetric quadrature rule.
*/
  rnodes = ( double * ) malloc ( 2 * numnodes * sizeof ( double ) );
  weights = ( double * ) malloc ( numnodes * sizeof ( double ) );

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );

  printf ( "\n" );
  printf ( "  NUMNODES = %d\n", numnodes );

  printf ( "\n" );
  printf ( "     J      W               X               Y\n" );
  printf ( "\n" );
  for ( j = 0; j < numnodes; j++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n",
      j, weights[j], rnodes[0+j*2], rnodes[1+j*2] );
  }

  d = r8vec_sum ( numnodes, weights );

  printf ( "   Sum  %14.6g\n", d );
  printf ( "  Area  %14.6g\n", area );

  free ( rnodes );
  free ( weights );

  return;
}
/******************************************************************************/

void test03 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], char *header )

/******************************************************************************/
/*
  Purpose:

    TEST03 calls TRIASYMQ_GNUPLOT to generate graphics files.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    28 June 2014

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
    of the quadrature rule.  0 <= DEGREE <= 50.

    Input, int NUMNODES, the number of nodes to be used by the rule.

    Input, double VERT1[2], VERT2[2], VERT3[2], the
    vertices of the triangle.

    Input, char *HEADER, an identifier for the graphics filenames.
*/
{
  double *rnodes;
  double *weights;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  TRIASYMQ_GNUPLOT creates gnuplot graphics files.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );

  rnodes = ( double * ) malloc ( 2 * numnodes * sizeof ( double ) );
  weights = ( double * ) malloc ( numnodes * sizeof ( double ) );

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );

  printf ( "  Number of nodes = %d\n", numnodes );

  triasymq_gnuplot ( vert1, vert2, vert3, numnodes, rnodes, header );

  free ( rnodes );
  free ( weights );

  return;
}
/******************************************************************************/

void test04 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[], char *header )

/******************************************************************************/
/*
  Purpose:

    TEST04 gets a rule and writes it to a file.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    20 June 2014

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
    of the quadrature rule.  0 <= DEGREE <= 50.

    Input, int NUMNODES, the number of nodes to be used by the rule.

    Input, double VERT1[2], VERT2[2], VERT3[2], the
    vertices of the triangle.

    Input, char *HEADER, an identifier for the filenames.
*/
{
  int j;
  double *rnodes;
  FILE *rule_unit;
  char rule_filename[255];
  double *weights;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Get a quadrature rule for a triangle.\n" );
  printf ( "  Then write it to a file.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );
/*
  Retrieve a symmetric quadrature rule.
*/
  rnodes = ( double * ) malloc ( 2 * numnodes * sizeof ( double ) );
  weights = ( double * ) malloc ( numnodes * sizeof ( double ) );

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );
/*
  Write the points and weights to a file.
*/
  strcpy ( rule_filename, header );
  strcat ( rule_filename, ".txt" );
 
  rule_unit = fopen ( rule_filename, "wt" );
  for ( j = 0; j < numnodes; j++ )
  {
    fprintf ( rule_unit, "21.15e  21.15e  21.15e\n",
      rnodes[0+j*2], rnodes[1+j*2], weights[j] );
  }
  fclose ( rule_unit );
  printf ( "\n" );
  printf ( "  Quadrature rule written to file '%s'\n", rule_filename );

  free ( rnodes );
  free ( weights );

  return;
}
/******************************************************************************/

void test05 ( int degree, int numnodes, double vert1[], double vert2[], 
  double vert3[] )

/******************************************************************************/
/*
  Purpose:

    TEST05 calls TRIASYMQ for a quadrature rule of given order and region.

  Licensing:

    This code is distributed under the GNU GPL license.

  Modified:

    28 June 2014

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

    Input, int DEGREE, the desired total polynomial degree 
    exactness of the quadrature rule.  0 <= DEGREE <= 50.

    Input, int NUMNODES, the number of nodes to be used by the rule.

    Input, double VERT1[2], VERT2[2], VERT3[2], the
    vertices of the triangle.
*/
{
  double area;
  double d;
  int i;
  int j;
  int npols;
  double *pols;
  double *r;
  double *rints;
  double *rnodes;
  double scale;
  double *weights;
  double z[2];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Compute a quadrature rule for a triangle.\n" );
  printf ( "  Check it by integrating orthonormal polynomials.\n" );
  printf ( "  Polynomial exactness degree DEGREE = %d\n", degree );

  area = triangle_area ( vert1, vert2, vert3 );
/*
  Retrieve a symmetric quadrature rule.
*/
  rnodes = ( double * ) malloc ( 2 * numnodes * sizeof ( double ) );
  weights = ( double * ) malloc ( numnodes * sizeof ( double ) );

  triasymq ( degree, vert1, vert2, vert3, rnodes, weights, numnodes );
/*
  Construct the matrix of values of the orthogonal polynomials
  at the user-provided nodes
*/
  npols = ( degree + 1 ) * ( degree + 2 ) / 2;
  rints = ( double * ) malloc ( npols * sizeof ( double ) );

  for ( j = 0; j < npols; j++ )
  {
    rints[j] = 0.0;
  }

  for ( i = 0; i < numnodes; i++ )
  {
    z[0] = rnodes[0+i*2];
    z[1] = rnodes[1+i*2];
    r = triangle_to_ref ( vert1, vert2, vert3, z );
    pols = ortho2eva ( degree, r );
    for ( j = 0; j < npols; j++ )
    {
      rints[j] = rints[j] + weights[i] * pols[j];
    }
    free ( pols );
    free ( r );
  }

  scale = sqrt ( sqrt ( 3.0 ) ) / sqrt ( area );

  for ( j = 0; j < npols; j++ )
  {
    rints[j] = rints[j] * scale;
  }

  d = pow ( rints[0] - sqrt ( area ), 2 );
  for ( j = 1; j < npols; j++ )
  {
    d = d + rints[j] * rints[j];
  }
  d = sqrt ( d ) / ( double ) ( npols );

  printf ( "\n" );
  printf ( "  RMS integration error = %g\n", d );

  free ( rints );
  free ( rnodes );
  free ( weights );

  return;
}
