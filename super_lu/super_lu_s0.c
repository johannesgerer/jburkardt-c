# include <stdlib.h>
# include <stdio.h>

# include "ssp_defs.h"

int main ( int arc, char *argv[] );

/**********************************************************************/

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
  Purpose:

    SUPER_LU_S0 runs a small 5 by 5 example of the use of SUPER_LU.

  Modified:

    27 April 2004

  Reference:

    James Demmel, John Gilbert, Xiaoye Li,
    SuperLU Users's Guide,
    Sections 1 and 2.
*/
{
  float *a;
  SuperMatrix A;
  int *asub;
  SuperMatrix B;
  int i;
  int info;
  SuperMatrix L;
  int m;
  int n;
  int nnz;
  int nrhs;
  superlu_options_t options;
  int *perm_c;
  int *perm_r;
  int permc_spec;
  float *rhs;
  float sol[5];
  SuperLUStat_t stat;
  SuperMatrix U;
  int *xa;
/*
  Say hello.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_S0:\n" );
  printf ( "  Simple 5 by 5 example of SUPER_LU solver.\n" );
/* 
  Initialize parameters. 
*/
  m = 5;
  n = 5;
  nnz = 12;
/* 
  Set aside space for the arrays. 
*/
  a = floatMalloc ( nnz );
  if ( !a ) 
  {
    ABORT ( "Malloc fails for a[]." );
  }

  asub = intMalloc ( nnz );
  if ( !asub ) 
  {
    ABORT ( "Malloc fails for asub[]." );
  }

  xa = intMalloc ( n+1 );
  if ( !xa ) 
  { 
    ABORT ( "Malloc fails for xa[]." );
  }
/* 
  Initialize matrix A. 
*/
  a[0] = 19.0; 
  a[1] = 12.0; 
  a[2] = 12.0; 
  a[3] = 21.0; 
  a[4] = 12.0; 
  a[5] = 12.0;
  a[6] = 21.0; 
  a[7] = 16.0; 
  a[8] = 21.0; 
  a[9] =  5.0; 
  a[10]= 21.0; 
  a[11]= 18.0;

  asub[0] = 0; 
  asub[1] = 1; 
  asub[2] = 4; 
  asub[3] = 1;
  asub[4] = 2; 
  asub[5] = 4; 
  asub[6] = 0; 
  asub[7] = 2;
  asub[8] = 0; 
  asub[9] = 3; 
  asub[10]= 3; 
  asub[11]= 4;

  xa[0] = 0; 
  xa[1] = 3; 
  xa[2] = 6; 
  xa[3] = 8; 
  xa[4] = 10; 
  xa[5] = 12;

  sol[0] = -0.031250000;
  sol[1] =  0.065476190;
  sol[2] =  0.013392857;
  sol[3] =  0.062500000;
  sol[4] =  0.032738095;
/* 
  Create matrix A in the format expected by SuperLU. 
*/
  sCreate_CompCol_Matrix ( &A, m, n, nnz, a, asub, xa, SLU_NC, SLU_S, SLU_GE );
/* 
  Create the right-hand side matrix B. 
*/
  nrhs = 1;
  rhs = floatMalloc ( m * nrhs );
  if ( !rhs ) 
  {
    ABORT ( "Malloc fails for rhs[].");
  }

  for ( i = 0; i < m; i++ ) 
  {
    rhs[i] = 1.0;
  }

  sCreate_Dense_Matrix ( &B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE );
/* 
  Set up the arrays for the permutations. 
*/
  perm_r = intMalloc ( m );
  if ( !perm_r ) 
  {
    ABORT ( "Malloc fails for perm_r[]." );
  }

  perm_c = intMalloc ( n );
  if ( !perm_c ) 
  {
    ABORT ( "Malloc fails for perm_c[]." );
  }
/* 
  Set the default input options, and then adjust some of them.
*/
  set_default_options ( &options );
  options.ColPerm = NATURAL;
/* 
  Initialize the statistics variables. 
*/
  StatInit ( &stat );
/*
  Factor the matrix and solve the system.
*/
  sgssv ( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
/*
  Print some of the results.
*/
  sPrint_CompCol_Matrix ( "Matrix A", &A );
  sPrint_SuperNode_Matrix ( "Factor L", &L );
  sPrint_CompCol_Matrix ( "Factor U", &U );
  sPrint_Dense_Matrix ( "Solution X", &B );

  printf ( "\n" );
  printf ( "  The exact solution:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%d  %f\n", i, sol[i] );
  }

  printf ( "\n" );
  print_int_vec ( "perm_r", m, perm_r );
/* 
  De-allocate storage.
*/
  SUPERLU_FREE ( rhs );
  SUPERLU_FREE ( perm_r );
  SUPERLU_FREE ( perm_c );
  Destroy_CompCol_Matrix ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperNode_Matrix ( &L );
  Destroy_CompCol_Matrix ( &U );
  StatFree ( &stat );

  printf ( "\n" );
  printf ( "SUPER_LU_S0:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
