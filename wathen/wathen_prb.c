# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "wathen.h"

int main ( );
void test01 ( );
void test02 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test10 ( );
void test11 ( );
void test115 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WATHEN_PRB.

  Discussion:

    WATHEN_PRB tests the WATHEN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WATHEN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WATHEN library.\n" );

  test01 ( );
  test02 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test10 ( );
  test11 ( );
  test115 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WATHEN_PRB\n" );
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

    TEST01 assembles, factor and solve using WATHEN_GE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int info;
  int *ipvt;
  int job;
  int n;
  int nx;
  int ny;
  int seed;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Assemble, factor and solve a Wathen system\n" );
  printf ( "  defined by WATHEN_GE.\n" );
  printf ( "\n" );

  nx = 4;
  ny = 4;
  printf ( "  Elements in X direction NX = %d\n", nx );
  printf ( "  Elements in Y direction NY = %d\n", ny );
  printf ( "  Number of elements = %d\n", nx * ny );
/*
  Compute the number of unknowns.
*/
  n = wathen_order ( nx, ny );
  printf ( "  Number of nodes N = %d\n", n );
/*
  Set up a random solution X.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
  seed = 123456789;
  a = wathen_ge ( nx, ny, n, &seed );
/*
  Compute the corresponding right hand side B.
*/
  b = mv_ge ( n, n, a, x1 );
/*
  Solve the linear system.
*/
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );
  info = dgefa ( a, n, n, ipvt );

  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = b[i];
  }
  job = 0;
  dgesl ( a, n, n, ipvt, x2, job );
/*
  Compute the maximum solution error.
*/
  e = r8vec_diff_norm_li ( n, x1, x2 );
  printf ( "  Maximum solution error is %g\n", e );
/*
  Free memory.
*/
  free ( a );
  free ( b );
  free ( ipvt );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 assembles, factors and solves using WATHEN_GB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int info;
  int *ipvt;
  int j;
  int jhi;
  int jlo;
  int job;
  int lda;
  int md;
  int ml;
  int mu;
  int n;
  int nx;
  int ny;
  int seed;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Assemble, factor and solve a Wathen system\n" );
  printf ( "  using WATHEN_GB.\n" );
  printf ( "\n" );

  nx = 4;
  ny = 4;
  printf ( "  Elements in X direction NX = %d\n", nx );
  printf ( "  Elements in Y direction NY = %d\n", ny );
  printf ( "  Number of elements = %d\n", nx * ny );
/*
  Compute the number of unknowns.
*/
  n = wathen_order ( nx, ny );
  printf ( "  Number of nodes N = %d\n", n );
/*
  Compute the bandwidth.
*/
  wathen_bandwidth ( nx, ny, &ml, &md, &mu );
  printf ( "  Lower bandwidth ML = %d\n", ml );
  printf ( "  Upper bandwidth MU = %d\n", mu );
/*
  Set up a random solution X1.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
  seed = 123456789;
  a = wathen_gb ( nx, ny, n, &seed );
/*
  Compute the corresponding right hand side B.
*/
  b = mv_gb ( n, n, ml, mu, a, x1 );
/*
  Solve the linear system.
*/
  lda = 2 * ml + mu + 1;
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );
  info = dgbfa ( a, lda, n, ml, mu, ipvt );

  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = b[i];
  }
  job = 0;
  dgbsl ( a, lda, n, ml, mu, ipvt, x2, job );
/*
  Compute the maximum solution error.
*/
  e = r8vec_diff_norm_li ( n, x1, x2 );
  printf ( "  Maximum solution error is %g\n", e );
/*
  Free memory.
*/
  free ( a );
  free ( b );
  free ( ipvt );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 measures the storage needed for the Wathen system.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int bd1;
  int bd2;
  int bl1;
  int bl2;
  int bu1;
  int bu2;
  int bw1;
  int bw2;
  int n;
  int nx;
  int ny;
  int seed;
  int storage_gb;
  int storage_ge;
  int storage_sparse;
  int test;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For various problem sizes and storage schemes,\n" );
  printf ( "  measure the storage used for the Wathen system.\n" );
  printf ( "\n" );
  printf ( "                                   Predicted  Observed\n" );
  printf ( "                              GE        Band      Band      " );
  printf ( "Band    Sparse\n" );
  printf ( "    NX  Elements   Nodes   storage     width     width   " );
  printf ( "  storage   storage\n" );
  printf ( "\n" );

  nx = 1;
  ny = 1;

  for ( test = 1; test <= 6; test++ )
  {
/*
  Compute the number of unknowns.
*/
    n = wathen_order ( nx, ny );
/*
  Predict the bandwidth.
*/
    wathen_bandwidth ( nx, ny, &bl1, &bd1, &bu1 );
    bw1 = bl1 + bd1 + bu1;
/*
  Compute the matrix.
*/
    seed = 123456789;
    a = wathen_ge ( nx, ny, n, &seed );

    storage_ge = n * n;

    bandwidth ( n, n, a, &bw2, &bl2, &bd2, &bu2 );
    storage_gb = ( 2 * bl2 + 1 + bu2 ) * n;

    storage_sparse = nonzeros ( n, n, a );
/*
  Report.
*/
    printf ( "  %4d      %4d  %6d  %8d  %8d  %8d  %8d  %8d\n",
      nx, nx * ny, n, storage_ge, bw1, bw2, storage_gb, storage_sparse );
/*
  Ready for next iteration.
*/
    nx = nx * 2;
    ny = ny * 2;

    free ( a );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 times WATHEN_GE assembly and solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int info;
  int *ipvt;
  int job;
  int n;
  int nx;
  int ny;
  int seed;
  int storage_ge;
  double t0;
  double t1;
  double t2;
  int test;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For various problem sizes,\n" );
  printf ( "  time the assembly and factorization of a Wathen system\n" );
  printf ( "  using the WATHEN_GE function.\n" );
  printf ( "\n" );
  printf ( 
    "    NX  Elements   Nodes   Storage    Assembly      Factor      Error\n" );
  printf ( "\n" );

  nx = 1;
  ny = 1;

  for ( test = 1; test <= 6; test++ )
  {
/*
  Compute the number of unknowns.
*/
    n = wathen_order ( nx, ny );
    storage_ge = n * n;
/*
  Set up a random solution X1.
*/
    seed = 123456789;
    x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix, and measure the storage required.
*/
    seed = 123456789;

    t0 = cpu_time ( );
    a = wathen_ge ( nx, ny, n, &seed );
    t1 = cpu_time ( );
    t1 = t1 - t0;
/*
  Compute the corresponding right hand side B.
*/
    b = mv_ge ( n, n, a, x1 );
/*
  Solve the system.
*/
    ipvt = ( int * ) malloc ( n * sizeof ( int ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x2[i] = b[i];
    }
    job = 0;

    t0 = cpu_time ( );
    info = dgefa ( a, n, n, ipvt );
    dgesl ( a, n, n, ipvt, x2, job );
    t2 = cpu_time ( );
    t2 = t2 - t0;
/*
  Compute the maximum solution error.
*/
    e = r8vec_diff_norm_li ( n, x1, x2 );
/*
  Report.
*/
    printf ( "  %4d      %4d  %6d  %8d  %10.2e  %10.2e  %10.2e\n",
      nx, nx * ny, n, storage_ge, t1, t2, e );
/*
  Ready for next iteration.
*/
    nx = nx * 2;
    ny = ny * 2;
/*
  Free memory.
*/
    free ( a );
    free ( b );
    free ( ipvt );
    free ( x1 );
    free ( x2 );
  }

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 times WATHEN_GB assembly and solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int info;
  int *ipvt;
  int j;
  int jhi;
  int jlo;
  int job;
  int lda;
  int md;
  int ml;
  int mu;
  int n;
  int nx;
  int ny;
  int seed;
  int storage_gb;
  double t0;
  double t1;
  double t2;
  int test;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  For various problem sizes,\n" );
  printf ( "  time the assembly and factorization of a Wathen system\n" );
  printf ( "  using the WATHEN_GB function.\n" );
  printf ( "\n" );
  printf ( "    NX  Elements   Nodes   Storage    Assembly      Factor      Error\n" );
  printf ( "\n" );

  nx = 1;
  ny = 1;

  for ( test = 1; test <= 6; test++ )
  {
/*
  Compute the number of unknowns.
*/
    n = wathen_order ( nx, ny );
/*
  Compute the bandwidth.
*/
    wathen_bandwidth ( nx, ny, &ml, &md, &mu );
    storage_gb = ( 2 * ml + mu + 1 ) * n;
/*
  Set up a random solution X1.
*/
    seed = 123456789;
    x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
    seed = 123456789;
    t0 = cpu_time ( );
    a = wathen_gb ( nx, ny, n, &seed );
    t1 = cpu_time ( );
    t1 = t1 - t0;
/*
  Compute the corresponding right hand side B.
*/
    b = mv_gb ( n, n, ml, mu, a, x1 );
/*
  Solve the system.
*/
    lda = 2 * ml + mu + 1;
    ipvt = ( int * ) malloc ( n * sizeof ( int ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x2[i] = b[i];
    }
    job = 0;

    t0 = cpu_time ( );
    info = dgbfa ( a, lda, n, ml, mu, ipvt );
    dgbsl ( a, lda, n, ml, mu, ipvt, x2, job );
    t2 = cpu_time ( );
    t2 = t2 - t0;
/*
  Compute the maximum solution error.
*/
    e = r8vec_diff_norm_li ( n, x1, x2 );
/*
  Report.
*/
    printf ( "  %4d      %4d  %6d  %8d  %10.2e  %10.2e  %10.2e\n",
      nx, nx * ny, n, storage_gb, t1, t2, e );
/*
  Ready for next iteration.
*/
    nx = nx * 2;
    ny = ny * 2;
/*
  Free memory.
*/
    free ( a );
    free ( b );
    free ( ipvt );
    free ( x1 );
    free ( x2 );
  }

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 times WATHEN_GE/WATHEN_GB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int info;
  int *ipvt;
  int j;
  int jhi;
  int jlo;
  int job;
  int lda;
  int md;
  int ml;
  int mu;
  int n;
  int nx;
  int ny;
  int seed;
  int storage_gb;
  int storage_ge;
  double t0;
  double t1;
  double t2;
  int test;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  For various problem sizes,\n" );
  printf ( "  time the assembly and factorization of a Wathen system\n" );
  printf ( "  WATHEN_GE/WATHEN_GB\n" );
  printf ( "\n" );
  printf ( "                   NX  Elements   Nodes   Storage    " );
  printf ( "  Assembly      Factor      Error\n" );

  nx = 1;
  ny = 1;

  for ( test = 1; test <= 6; test++ )
  {
/*
  Compute the number of unknowns.
*/
    n = wathen_order ( nx, ny );
    storage_ge = n * n;
/*
  Set up a random solution X1.
*/
    seed = 123456789;
    x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
    seed = 123456789;
    t0 = cpu_time ( );
    a = wathen_ge ( nx, ny, n, &seed );
    t1 = cpu_time ( );
    t1 = t1 - t0;
/*
  Compute the corresponding right hand side B.
*/
    b = mv_ge ( n, n, a, x1 );
/*
  Solve the system.
*/
    ipvt = ( int * ) malloc ( n * sizeof ( int ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x2[i] = b[i];
    }
    job = 0;

    t0 = cpu_time ( );
    info = dgefa ( a, n, n, ipvt );
    dgesl ( a, n, n, ipvt, x2, job );
    t2 = cpu_time ( );
    t2 = t2 - t0;
/*
  Compute the maximum solution error.
*/
    e = r8vec_diff_norm_li ( n, x1, x2 );
/*
  Report.
*/
    printf ( "\n" );
    printf ( 
      "  WATHEN_GE      %4d      %4d  %6d  %8d  %10.2e  %10.2e  %10.2e\n",
      nx, nx * ny, n, storage_ge, t1, t2, e );
/*
  Free memory.
*/
    free ( a );
    free ( b );
    free ( ipvt );
    free ( x1 );
    free ( x2 );
/*
  Compute the bandwidth.
*/
    wathen_bandwidth ( nx, ny, &ml, &md, &mu );
    storage_gb = ( 2 * ml + mu + 1 ) * n;
/*
  Set up a random solution X1.
*/
    seed = 123456789;
    x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
    seed = 123456789;
    t0 = cpu_time ( );
    a = wathen_gb ( nx, ny, n, &seed );
    t1 = cpu_time ( );
    t1 = t1 - t0;
/*
  Compute the corresponding right hand side B.
*/
    b = mv_gb ( n, n, ml, mu, a, x1 );
/*
  Solve the system.
*/
    lda = 2 * ml + mu + 1;
    ipvt = ( int * ) malloc ( n * sizeof ( int ) );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      x2[i] = b[i];
    }
    job = 0;

    t0 = cpu_time ( );
    info = dgbfa ( a, lda, n, ml, mu, ipvt );
    dgbsl ( a, lda, n, ml, mu, ipvt, x2, job );
    t2 = cpu_time ( );
    t2 = t2 - t0;
/*
  Compute the maximum solution error.
*/
    e = r8vec_diff_norm_li ( n, x1, x2 );
/*
  Report.
*/
    printf ( 
      "  WATHEN_GB      %4d      %4d  %6d  %8d  %10.2e  %10.2e  %10.2e\n",
      nx, nx * ny, n, storage_gb, t1, t2, e );
/*
  Free memory.
*/
    free ( a );
    free ( b );
    free ( ipvt );
    free ( x1 );
    free ( x2 );
/*
  Ready for next iteration.
*/
    nx = nx * 2;
    ny = ny * 2;
  }

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 assembles, factor and solve using WATHEN_GE and CG_GE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int n;
  int nx;
  int ny;
  int seed;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  Assemble, factor and solve a Wathen system\n" );
  printf ( "  defined by WATHEN_GE and CG_GE.\n" );
  printf ( "\n" );

  nx = 1;
  ny = 1;
  printf ( "  Elements in X direction NX = %d\n", nx );
  printf ( "  Elements in Y direction NY = %d\n", ny );
  printf ( "  Number of elements = %d\n", nx * ny );
/*
  Compute the number of unknowns.
*/
  n = wathen_order ( nx, ny );
  printf ( "  Number of nodes N = %d\n", n );
/*
  Set up a random solution X.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
  seed = 123456789;
  a = wathen_ge ( nx, ny, n, &seed );
/*
  Compute the corresponding right hand side B.
*/
  b = mv_ge ( n, n, a, x1 );
/*
  Solve the linear system.
*/
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  cg_ge ( n, a, b, x2 );
/*
  Compute the maximum solution error.
*/
  e = r8vec_diff_norm_li ( n, x1, x2 );
  printf ( "  Maximum solution error is %g\n", e ); 
/*
  Free memory.
*/
  free ( a );
  free ( b );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 assemble, factor and solve using WATHEN_ST + CG_ST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int *col;
  double e;
  int i;
  int n;
  int nx;
  int ny;
  int nz_num;
  int *row;
  int seed;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  Assemble, factor and solve a Wathen system\n" );
  printf ( "  defined by WATHEN_ST and CG_ST.\n" );
  printf ( "\n" );

  nx = 1;
  ny = 1;
  printf ( "  Elements in X direction NX = %d\n", nx );
  printf ( "  Elements in Y direction NY = %d\n", ny );
  printf ( "  Number of elements = %d\n", nx * ny );
/*
  Compute the number of unknowns.
*/
  n = wathen_order ( nx, ny );
  printf ( "  Number of nodes N = %d\n", n );
/*
  Set up a random solution X1.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix size.
*/
  nz_num = wathen_st_size ( nx, ny );
  printf ( "  Number of nonzeros NZ_NUM = %d\n", nz_num );
/*
  Compute the matrix.
*/
  seed = 123456789;
  row = ( int * ) malloc ( nz_num * sizeof ( int ) );
  col = ( int * ) malloc ( nz_num * sizeof ( int ) );
  a = wathen_st ( nx, ny, nz_num, &seed, row, col );
/*
  Compute the corresponding right hand side B.
*/
  b = mv_st ( n, n, nz_num, row, col, a, x1 );
/*
  Solve the linear system.
*/
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  cg_st ( n, nz_num, row, col, a, b, x2 );
/*
  Compute the maximum solution error.
*/
  e = r8vec_diff_norm_li ( n, x1, x2 );
  printf ( "  Maximum solution error is %g\n", e );
/*
  Free memory.
*/
  free ( a );
  free ( b );
  free ( col );
  free ( row );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test115 ( )

/******************************************************************************/
/*
  Purpose:

    TEST115 assembles, factors and solves using WATHEN_GB and CG_GB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double e;
  int i;
  int j;
  int jhi;
  int jlo;
  int lda;
  int md;
  int ml;
  int mu;
  int n;
  int nx;
  int ny;
  int seed;
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST115\n" );
  printf ( "  Assemble, factor and solve a Wathen system\n" );
  printf ( "  using WATHEN_GB and CG_GB.\n" );
  printf ( "\n" );

  nx = 4;
  ny = 4;
  printf ( "  Elements in X direction NX = %d\n", nx );
  printf ( "  Elements in Y direction NY = %d\n", ny );
  printf ( "  Number of elements = %d\n", nx * ny );
/*
  Compute the number of unknowns.
*/
  n = wathen_order ( nx, ny );
  printf ( "  Number of nodes N = %d\n", n );
/*
  Compute the bandwidth.
*/
  wathen_bandwidth ( nx, ny, &ml, &md, &mu );
  printf ( "  Lower bandwidth ML = %d\n", ml );
  printf ( "  Upper bandwidth MU = %d\n", mu );
/*
  Set up a random solution X1.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Compute the matrix.
*/
  seed = 123456789;
  a = wathen_gb ( nx, ny, n, &seed );
/*
  Compute the corresponding right hand side B.
*/
  b = mv_gb ( n, n, ml, mu, a, x1 );
/*
  Solve the linear system.
*/
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  cg_gb ( n, ml, mu, a, b, x2 );
/*
  Compute the maximum solution error.
*/
  e = r8vec_diff_norm_li ( n, x1, x2 );
  printf ( "  Maximum solution error is %g\n", e );
/*
  Free memory.
*/
  free ( a );
  free ( b );
  free ( x1 );
  free ( x2 );

  return;
}
