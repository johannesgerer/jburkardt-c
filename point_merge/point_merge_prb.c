# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "point_merge.h"

int main ( void );
void test01 ( int m, int n, int n_unique, int seed );
void test02 ( int m, int n, int n_unique, double tol, int seed );
void test03 ( int m, int n, int n_unique, double tol, int seed );
void test04 ( int m, int n, int n_unique, double tol, int seed );
void test05 ( int m, int n, int n_unique, double tol, int seed );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POINT_MERGE_PRB.

  Discussion:

    POINT_MERGE_PRB tests the POINT_MERGE library.

    Compare correctness of the codes.

    Compare speed of the codes.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int n_unique;
  int seed;
  double tol;

  timestamp ( );
  printf (  " \n" );
  printf (  "POINT_MERGE_PRB\n" );
  printf (  "  C version\n" );
  printf (  "  Test the POINT_MERGE library.\n" );
/*
  TEST01 gives me some confidence that, at least for zero-tolerance,
  the radial approach is accurate, as compared to the "Unique count"
  (which cannot be extended to a tolerance version in multidimensions)
  and the "Tol Unique Count", which is an O(N^2) algorithm.
*/
  m = 3;
  n = 10;
  n_unique = 7;
  seed = 123456789;
  test01 ( m, n, n_unique, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  seed = 987654321;
  test01 ( m, n, n_unique, seed );
/*
  In TEST02, I want to compute the same data, but with "blurred"
  duplicates, and a tolerance version of the radial approach,
  compared to "Tol Unique Count".
*/
  m = 3;
  n = 10;
  n_unique = 7;
  tol = 0.00001;
  seed = 123456789;
  test02 ( m, n, n_unique, tol, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  tol = 0.00001;
  seed = 987654321;
  test02 ( m, n, n_unique, tol, seed );
/*
  In TEST03, I want to measure the time required for a sequence
  of increasingly hard problems.
*/
  m = 3;
  n = 100;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 1000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 10000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  if ( 0 )
  {
    m = 3;
    n = 100000;
    n_unique = n / 2;
    tol = 0.00001;
    seed = 123456789;
    test03 ( m, n, n_unique, tol, seed );
  }
/*
  In TEST04, repeat TEST02, but now compute the index vector.
*/
  m = 3;
  n = 10;
  n_unique = 7;
  tol = 0.00001;
  seed = 123456789;
  test04 ( m, n, n_unique, tol, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  tol = 0.00001;
  seed = 987654321;
  test04 ( m, n, n_unique, tol, seed );
/*
  In TEST05, I want to measure the time required for a sequence
  of increasingly hard problems.
*/
  m = 3;
  n = 100;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 1000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 10000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  if ( 0 )
  {
    m = 3;
    n = 100000;
    n_unique = n / 2;
    tol = 0.00001;
    seed = 123456789;
    test05 ( m, n, n_unique, tol, seed );
  }
/*
  Terminate.
*/
  printf (  " \n" );
  printf (  "POINT_MERGE_PRB\n" );
  printf (  "  Normal end of execution.\n" );
  printf (  " \n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int m, int n, int n_unique, int seed )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests uniqueness counting with no tolerance. 

  Discussion:

    POINT_UNIQUE_COUNT uses an O(N) algorithm.
    POINT_RADIAL_UNIQUE_COUNT uses an algorithm that should be,
      in general, O(N);
    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.

    For this test, we just want to make sure the algorithms agree
    in the counting.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double tol;
  int unique_num;

  printf (  " \n" );
  printf (  "TEST01\n" );
  printf (  "  To count the unique columns in an R8COL, we call\n" );
  printf (  "  POINT_UNIQUE_COUNT,\n" );
  printf (  "  POINT_RADIAL_UNIQUE_COUNT, (with random center)\n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n" );
  printf (  " \n" );
  printf (  "  M =     %d\n", m );
  printf (  "  N =     %d\n", n );
  printf (  "  SEED =  %d\n",seed );

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );

  printf (  " \n" );
  printf (  "  N_UNIQUE =                  %d\n", n_unique );

  unique_num = point_unique_count ( m, n, a );
  printf (  "  POINT_UNIQUE_COUNT =        %d\n", unique_num );

  unique_num = point_radial_unique_count ( m, n, a, &seed );
  printf (  "  POINT_RADIAL_UNIQUE_COUNT = %d\n", unique_num );

  tol = 0.0;
  unique_num = point_tol_unique_count ( m, n, a, tol );
  printf (  "  POINT_TOL_UNIQUE_COUNT =    %d\n", unique_num );

  free ( a );

  return;
}
/******************************************************************************/

void test02 ( int m, int n, int n_unique, double tol, int seed )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests uniqueness counting with a tolerance. 

  Discussion:

    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
      in general, O(N);
    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.

    For this test, we just want to make sure the algorithms agree
    in the counting.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int j;
  double *r;
  double r_norm;
  int unique_num;

  printf (  " \n" );
  printf (  "TEST02\n" );
  printf (  "  To count the unique columns in an R8COL, we call\n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n" );
  printf (  " \n" );
  printf (  "  M =     %d\n", m );
  printf (  "  N =     %d\n", n );
  printf (  "  TOL =  %e\n", tol );
  printf (  "  SEED =  %d\n", seed );

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );
/*
  The form of the tolerance test means that if two vectors are initially
  equal, they remain "tolerably equal" after the addition of random
  perturbation vectors whose 2-norm is no greater than TOL/2.
*/
  r = ( double * ) malloc ( m * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  free ( r );

  r8mat_transpose_print ( m, n, a, "  Blurred matrix:" );

  printf (  " \n" );
  printf (  "  N_UNIQUE =                      %d\n", n_unique );

  unique_num = point_radial_tol_unique_count ( m, n, a, tol, &seed );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT = %d\n", unique_num );

  unique_num = point_tol_unique_count ( m, n, a, tol );
  printf (  "  POINT_TOL_UNIQUE_COUNT =        %d\n", unique_num );

  free ( a );

  return;
}
/******************************************************************************/

void test03 ( int m, int n, int n_unique, double tol, int seed )

/******************************************************************************/
/*
  Purpose:

    TEST03 compares timings for two uniqueness counters.

  Discussion:

    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
      in general, O(N);
    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double ctime;
  int i;
  int j;
  double *r;
  double r_norm;
  int unique_num;

  printf (  " \n" );
  printf (  "TEST03\n" );
  printf (  "  To count the unique columns in an R8COL, we call\n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n" );
  printf (  " \n" );
  printf (  "  M =     %d\n", m );
  printf (  "  N =     %d\n", n );
  printf (  "  TOL =  %e\n", tol );
  printf (  "  SEED =  %d\n", seed );

  a = r8col_duplicates ( m, n, n_unique, &seed );
/*
  The form of the tolerance test means that if two vectors are initially
  equal, they remain "tolerably equal" after the addition of random
  perturbation vectors whose 2-norm is no greater than TOL/2.
*/
  r = ( double * ) malloc ( m * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  free ( r );

  printf (  " \n" );
  printf (  "  N_UNIQUE =                      %d\n", n_unique );

  ctime = cpu_time ( );
  unique_num = point_radial_tol_unique_count ( m, n, a, tol, &seed );
  ctime = cpu_time ( ) - ctime;
  printf (  " \n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT = %d\n", unique_num );
  printf (  "  CPU_TIME = %f\n", ctime );

  ctime = cpu_time ( );
  unique_num = point_tol_unique_count ( m, n, a, tol );
  ctime = cpu_time ( ) - ctime;
  printf (  " \n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT =        %d\n", unique_num );
  printf (  "  CPU_TIME = %f\n", ctime );

  free ( a );

  return;
}
/******************************************************************************/

void test04 ( int m, int n, int n_unique, double tol, int seed )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests uniqueness indexing with a tolerance. 

  Discussion:

    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
      in general, O(N);
    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.

    For this test, we just want to make sure the algorithms agree
    in the counting.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double dist;
  int i;
  int j;
  int k;
  double *r;
  double r_norm;
  int *undx;
  int unique_num;
  int *xdnu;

  printf (  " \n" );
  printf (  "TEST04\n" );
  printf (  "  To index the unique columns in an R8COL, we call\n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n" );
  printf (  " \n" );
  printf (  "  M =     %d\n", m );
  printf (  "  N =     %d\n", n );
  printf (  "  TOL =  %e\n", tol );
  printf (  "  SEED =  %d\n", seed );

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );
/*
  The form of the tolerance test means that if two vectors are initially
  equal, they remain "tolerably equal" after the addition of random
  perturbation vectors whose 2-norm is no greater than TOL/2.
*/
  r = ( double * ) malloc ( m * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  free ( r );

  r8mat_transpose_print ( m, n, a, "  Blurred matrix:" );

  printf (  " \n" );
  printf (  "  N_UNIQUE =                      %d\n", n_unique );

  undx = ( int * ) malloc ( n * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( int ) );

  unique_num = point_radial_tol_unique_index ( m, n, a, tol, &seed, undx, 
    xdnu );

  printf (  " \n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_INDEX\n" );
  printf (  "  Unique_num = %d\n", unique_num );

  i4vec_print ( unique_num, undx, "  UNDX:" );

  i4vec_print ( n, xdnu, "  XDNU:" );

  printf (  " \n" );
  printf (  "  List of nonunique points P(J), represented by\n" );
  printf (  "  point with index I(J).\n" );
  printf (  " \n" );
  printf (  "  J, P(J)\n" );
  printf (  "  I(J), P(I(J))\n" );
  printf (  "  || P(J) - P(I(J)) || (should be <= TOL)\n" );
  printf (  " \n" );
  for ( j = 0; j < n; j++ )
  {
    k = undx[xdnu[j]];
    if ( j != k )
    {
      printf (  " \n" );
      printf (  "  %4d", j );
      for ( i = 0; i < m; i++ )
      {
        printf (  "  %10f", a[i+j*m] );
      }
      printf (  "\n" );
      printf (  "  %4d", k );
      for ( i = 0; i < m; i++ )
      {
        printf (  "  %10f", a[i+k*m] );
      }
      printf (  "\n" );
      dist = 0.0;
      for ( i = 0; i < m; i++ )
      {
        dist = dist + pow ( a[i+j*m] - a[i+k*m], 2 );
      }
      dist = sqrt ( dist );
      printf (  "          %10e\n", dist );
    }
  }
/*
  The interpretation of XDNU is simpler for POINT_TOL_UNIQUE_INDEX.
*/
  unique_num = point_tol_unique_index ( m, n, a, tol, xdnu );

  printf (  " \n" );
  printf (  "  POINT_TOL_UNIQUE_INDEX\n" );
  printf (  "  Unique_num = %d\n", unique_num );
  printf (  " \n" );
  printf (  "  List of nonunique points P(J), represented by\n" );
  printf (  "  point with index I(J).\n" );
  printf (  " \n" );
  printf (  "  J, P(J)\n" );
  printf (  "  I(J), P(I(J))\n" );
  printf (  "  || P(J) - P(I(J)) || (should be <= TOL)\n" );
  printf (  " \n" );
  for ( j = 0; j < n; j++ )
  {
    k = xdnu[j];
    if ( j != k )
    {
      printf (  " \n" );
      printf (  "  %4d", j );
      for ( i = 0; i < m; i++ )
      {
        printf (  "  %10f", a[i+j*m] );
      }
      printf (  "\n" );
      printf (  "  %4d", k );
      for ( i = 0; i < m; i++ )
      {
        printf (  "  %10f", a[i+k*m] );
      }
      printf (  "\n" );
      dist = 0.0;
      for ( i = 0; i < m; i++ )
      {
        dist = dist + pow ( a[i+j*m] - a[i+k*m], 2 );
      }
      dist = sqrt ( dist );
      printf (  "          %10e\n", dist );
    }
  }

  free ( a );
  free ( undx );
  free ( xdnu );

  return;
}
/******************************************************************************/

void test05 ( int m, int n, int n_unique, double tol, int seed )

/******************************************************************************/
/*
  Purpose:

    TEST05 times uniqueness indexing with a tolerance. 

  Discussion:

    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
      in general, O(N);
    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.

    For this test, we just want to make sure the algorithms agree
    in the counting.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 July 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double ctime;
  double dist;
  int i;
  int j;
  double *r;
  double r_norm;
  int *undx;
  int unique_num;
  int *xdnu;

  printf (  " \n" );
  printf (  "TEST05\n" );
  printf (  "  We time the computations in TEST04, calling\n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n" );
  printf (  "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n" );
  printf (  " \n" );
  printf (  "  M =     %d\n", m );
  printf (  "  N =     %d\n", n );
  printf (  "  TOL =  %e\n", tol );
  printf (  "  SEED =  %d\n", seed );

  a = r8col_duplicates ( m, n, n_unique, &seed );
/*
  The form of the tolerance test means that if two vectors are initially
  equal, they remain "tolerably equal" after the addition of random
  perturbation vectors whose 2-norm is no greater than TOL/2.
*/
  r = ( double * ) malloc ( m * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  free ( r );

  printf (  " \n" );
  printf (  "  N_UNIQUE =                      %d\n", n_unique );

  undx = ( int * ) malloc ( n * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( int ) );

  ctime = cpu_time ( );
  unique_num = point_radial_tol_unique_index ( m, n, a, tol, &seed, undx,
    xdnu );
  ctime = cpu_time ( ) - ctime;

  printf (  " \n" );
  printf (  "  POINT_RADIAL_TOL_UNIQUE_INDEX\n" );
  printf (  "  Unique_num = %d\n", unique_num );
  printf (  "  Time = %f\n", ctime );

  ctime = cpu_time ( );
  unique_num = point_tol_unique_index ( m, n, a, tol, xdnu );
  ctime = cpu_time ( ) - ctime;

  printf (  " \n" );
  printf (  "  POINT_TOL_UNIQUE_INDEX\n" );
  printf (  "  Unique_num = %d\n", unique_num );
  printf (  "  Time = %f\n", ctime );

  free ( a );
  free ( undx );
  free ( xdnu );

  return;
}
