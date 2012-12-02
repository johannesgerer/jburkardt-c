# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "condition.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CONDITION_PRB.
//
//  Discussion:
//
//    CONDITION_PRB tests the CONDITION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  printf ( "\n" );
  printf ( "CONDITION_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CONDITION library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "CONDITION_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CONDITION_LINPACK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int n;
  char name[80];
  int *pivot;
  int seed;
  double *z;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  CONDITION_LINPACK estimates the L1 condition number.\n" );
  printf ( "\n" );
  printf ( "  Matrix               Order   Condition         Linpack\n" );
  printf ( "\n" );
//
//  Combinatorial matrix.
//
  strcpy ( name, "Combinatorial" );
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX1
//
  strcpy ( name, "CONEX1" );
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX2
//
  strcpy ( name, "CONEX2" );
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX3
//
  strcpy ( name, "CONEX3" );
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX4
//
  strcpy ( name, "CONEX4" );
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  KAHAN
//
  strcpy ( name, "KAHAN" );
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_linpack ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  Random
//
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    strcpy ( name, "RANDOM" );
    n = 4;
    a = r8mat_uniform_01_new ( n, n, &seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = ( int * ) malloc ( n * sizeof ( int ) );
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    cond = condition_linpack ( n, a );
    printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
    free ( a );
    free ( a_inverse );
    free ( a_lu );
    free ( pivot );
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CONDITION_SAMPLE1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int j;
  int m;
  int m_test[3] = { 10, 1000, 100000 };
  int n;
  char name[80];
  int *pivot;
  int seed;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  CONDITION_SAMPLE1 estimates the L1 condition number using sampling.\n" );
  printf ( "\n" );
  printf ( "  Matrix                 Samples Order   Condition        Estimate\n" );
//
//  Combinatorial matrix.
//
  strcpy ( name, "Combinatorial" );
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  CONEX1
//
  strcpy ( name, "CONEX1" );
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  CONEX2
//
  strcpy ( name, "CONEX2" );
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  CONEX3
//
  strcpy ( name, "CONEX3" );
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  CONEX4
//
  strcpy ( name, "CONEX4" );
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  KAHAN
//
  strcpy ( name, "KAHAN" );
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    m = m_test[i];
    cond = condition_sample1 ( n, a, m );
    printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
  }
  free ( a );
  free ( a_inverse );
//
//  Random
//
  seed = 123456789;

  for ( j = 1; j <= 5; j++ )
  {
    strcpy ( name, "RANDOM" );
    n = 4;
    a = r8mat_uniform_01_new ( n, n, &seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = ( int * ) malloc ( n * sizeof ( int ) );
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    printf ( "\n" );
    for ( i = 0; i < 3; i++ )
    {
      m = m_test[i];
      cond = condition_sample1 ( n, a, m );
      printf ( "  %20s  %8d  %4d  %14g  %14g\n", name, m, n, cond_l1, cond );
    }
    free ( a );
    free ( a_inverse );
    free ( a_lu );
    free ( pivot );
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CONDITION_HAGER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *a_inverse;
  double a_inverse_norm_l1;
  double *a_lu;
  double a_norm_l1;
  double alpha;
  double beta;
  double cond;
  double cond_l1;
  int i;
  int info;
  int n;
  char name[80];
  int *pivot;
  int seed;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  CONDITION_HAGER estimates the L1 condition number.\n" );
  printf ( "\n" );
  printf ( "  Matrix               Order   Condition         Hager\n" );
  printf ( "\n" );
//
//  Combinatorial matrix.
//
  strcpy ( name, "Combinatorial" );
  n = 4;
  alpha = 2.0;
  beta = 3.0;
  a = combin ( alpha, beta, n );
  a_inverse = combin_inverse ( alpha, beta, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX1
//
  strcpy ( name, "CONEX1" );
  n = 4;
  alpha = 100.0;
  a = conex1 ( alpha );
  a_inverse = conex1_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX2
//
  strcpy ( name, "CONEX2" );
  n = 3;
  alpha = 100.0;
  a = conex2 ( alpha );
  a_inverse = conex2_inverse ( alpha );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX3
//
  strcpy ( name, "CONEX3" );
  n = 5;
  a = conex3 ( n );
  a_inverse = conex3_inverse ( n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  CONEX4
//
  strcpy ( name, "CONEX4" );
  n = 4;
  a = conex4 ( );
  a_inverse = conex4_inverse ( );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  KAHAN
//
  strcpy ( name, "KAHAN" );
  n = 4;
  alpha = 0.25;
  a = kahan ( alpha, n, n );
  a_inverse = kahan_inverse ( alpha, n );
  a_norm_l1 = r8mat_norm_l1 ( n, n, a );
  a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
  cond_l1 = a_norm_l1 * a_inverse_norm_l1;
  cond = condition_hager ( n, a );
  printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
  free ( a );
  free ( a_inverse );
//
//  Random
//
  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    strcpy ( name, "RANDOM" );
    n = 4;
    a = r8mat_uniform_01_new ( n, n, &seed );
    a_lu = r8mat_copy_new ( n, n, a );
    pivot = ( int * ) malloc ( n * sizeof ( int ) );
    info = r8ge_fa ( n, a_lu, pivot );
    a_inverse = r8ge_inverse ( n, a_lu, pivot );
    a_norm_l1 = r8mat_norm_l1 ( n, n, a );
    a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse );
    cond_l1 = a_norm_l1 * a_inverse_norm_l1;
    cond = condition_hager ( n, a );
    printf ( "  %20s  %4d  %14g  %14g\n", name, n, cond_l1, cond );
    free ( a );
    free ( a_inverse );
    free ( a_lu );
    free ( pivot );
  }
  return;
}
