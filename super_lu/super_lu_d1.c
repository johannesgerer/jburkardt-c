# include <stdlib.h>
# include <stdio.h>

# include "dsp_defs.h"

int main ( int argc, char *argv[] );

/**********************************************************************/

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
  Purpose:

    SUPER_LU_D1 solves a general sparse system read from a file.

  Discussion:

    The sparse matrix is stored in a file using the Harwell-Boeing
    sparse matrix format.  The file should be assigned to the standard
    input of this program.  For instance, if the matrix is stored
    in the file "g10_rua.txt", the execution command might be:

      super_lu_d1 < g10_rua.txt

  Modified:

    27 April 2004

  Reference:

    James Demmel, John Gilbert, Xiaoye Li,
    SuperLU Users's Guide,
    Sections 1 and 2.

  Local parameters:

    SuperMatrix L, the computed L factor.

    int *perm_c, the column permutation vector.

    int *perm_r, the row permutations from partial pivoting.

    SuperMatrix U, the computed U factor.
*/
{
  SuperMatrix A;
  NCformat *Astore;
  double *a;
  int *asub;
  SuperMatrix B;
  int info;
  SuperMatrix L;
  int ldx;
  SCformat *Lstore;
  int m;
  mem_usage_t mem_usage;
  int n;
  int nnz;
  int nrhs;
  superlu_options_t options;
  int *perm_c;
  int *perm_r;
  double *rhs;
  double *sol;
  SuperLUStat_t stat;
  SuperMatrix U;
  NCformat *Ustore;
  int *xa;
  double *xact;
/*
  Say hello.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_D1:\n" );
  printf ( "  Read a general sparse matrix A from standard input,\n");
  printf ( "  stored in Harwell-Boeing Sparse Matrix format.\n" );
  printf ( "\n" );
  printf ( "  Solve a linear system A * X = B.\n" );
/* 
  Set the default input options:

  options.Fact = DOFACT;
  options.Equil = YES;
  options.ColPerm = COLAMD;
  options.DiagPivotThresh = 1.0;
  options.Trans = NOTRANS;
  options.IterRefine = NOREFINE;
  options.SymmetricMode = NO;
  options.PivotGrowth = NO;
  options.ConditionNumber = NO;
  options.PrintStat = YES;
*/
  set_default_options ( &options );
/*  
  Read the matrix from standard input, in Harwell-Boeing format. 
*/
  dreadhb ( &m, &n, &nnz, &a, &asub, &xa );
/*
  Create storage for a compressed column matrix.
*/
  dCreate_CompCol_Matrix ( &A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE );
  Astore = A.Store;

  printf ( "\n" );
  printf ( "  Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz );
/*
  Set up the right hand side.
*/    
  nrhs = 1;
  rhs = doubleMalloc ( m * nrhs );

  if ( !rhs ) 
  {
    ABORT ( "Malloc fails for rhs[]." );
  }

  dCreate_Dense_Matrix ( &B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE );

  xact = doubleMalloc ( n * nrhs );
  if ( !xact ) 
  {
    ABORT ( "Malloc fails for xact[]." );
  }

  ldx = n;
  dGenXtrue ( n, nrhs, xact, ldx );
  dFillRHS ( options.Trans, nrhs, xact, ldx, &A, &B );

  perm_c = intMalloc ( n );
  if ( !perm_c ) 
  {
    ABORT ( "Malloc fails for perm_c[]." );
  }

  perm_r = intMalloc ( m );
  if ( !perm_r )
  {
    ABORT ( "Malloc fails for perm_r[]." );
  }
/* 
  Initialize the statistics variables. 
*/
  StatInit ( &stat );
/*
  Call DGSSV to factor the matrix and solve the linear system.
*/
  dgssv ( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
 
  if ( info == 0 ) 
  {
/* 
  To conveniently access the solution matrix, you need to get a pointer to it. 
*/
    sol = (double*) ((DNformat*) B.Store)->nzval; 
/*
  Compute the infinity norm of the error.
*/
    dinf_norm_error ( nrhs, &B, xact );

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;

    printf ( "\n" );
    printf ( "  Number of nonzeros in factor L = %d\n", Lstore->nnz );
    printf ( "  Number of nonzeros in factor U = %d\n", Ustore->nnz );
    printf ( "  Number of nonzeros in L+U = %d\n", 
      Lstore->nnz + Ustore->nnz - n );
  
    dQuerySpace ( &L, &U, &mem_usage );

    printf ( "\n" );
    printf ( "  L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
      mem_usage.expansions );
  
  } 
  else
  {
    printf ( "\n" );
    printf ( "  DGSSV error returns INFO = %d\n", info );
/*
  Was the factorization completed?
*/
    if ( info <= n ) 
    { 
      dQuerySpace ( &L, &U, &mem_usage );
      printf ( "  L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
       mem_usage.expansions );
    }
  }

  if ( options.PrintStat ) 
  {
    StatPrint ( &stat );
  }

  StatFree ( &stat );
/*
  Free the memory.
*/
  SUPERLU_FREE ( rhs );
  SUPERLU_FREE ( xact );
  SUPERLU_FREE ( perm_r );
  SUPERLU_FREE ( perm_c );
  Destroy_CompCol_Matrix ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperNode_Matrix ( &L );
  Destroy_CompCol_Matrix ( &U );
/*
  Say goodbye.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_D1:\n" );
  printf ( "  Normal end of execution.\n");

  return 0;
}
