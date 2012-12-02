# include <stdlib.h>
# include <stdio.h>

# include "dsp_defs.h"

int main ( int argc, char *argv[] );

static void parse_command_line ( int argc, char *argv[], int *lwork,
  double *u, yes_no_t *equil, trans_t *trans );

/**********************************************************************/

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
  Purpose:

    SUPER_LU_D5 solves a sparse system read from a file using DGSSVX.

  Discussion:

    In this example, we show how to take advantage of the fact that
    several systems to be solved have the same sparsity pattern
    (though not, presumably, the same values.)  

    DGSSVX can reuse the sparsity pattern information from one
    system, to speed up the solution of the next, presuming, of
    course, that the systems do share a sparsity pattern!

    In particular, PERM_C and ETREE will be computed once and
    reused in subsequent calls.

    The sparse matrix is stored in a file using the Harwell-Boeing
    sparse matrix format.  The file should be assigned to the standard
    input of this program.  For instance, if the matrix is stored
    in the file "g10_rua.txt", the execution command might be:

      super_lu_d5 < g10_rua.txt

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
  double *a;
  SuperMatrix A1;
  double *a1;
  NCformat *Astore;
  int *asub;
  int *asub1;
  SuperMatrix B;
  SuperMatrix B1;
  double *berr;
  double *C;
  char equed[1];
  yes_no_t equil;
  int *etree;
  double *ferr;
  int i;
  int info;
  int j;
  SuperMatrix L;
  int ldx;
  SCformat *Lstore;
  int lwork;
  int m;
  mem_usage_t mem_usage;
  int n;
  int nnz;
  int nrhs;
  superlu_options_t options;
  int *perm_c;
  int *perm_r;
  double *R;
  double rcond;
  double *rhsb;
  double *rhsb1;
  double *rhsx;
  double rpg;
  double *sol;
  SuperLUStat_t stat;
  trans_t trans;
  SuperMatrix U;
  double u;
  NCformat *Ustore;
  void *work;
  SuperMatrix X;
  int *xa;
  int *xa1;
  double *xact;
/*
  Say hello.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_D5:\n" );
  printf ( "  Solve multiple linear systems which share a sparsity,\n");
  printf ( "  pattern.  Save time by letting DGSSVX know it can reuse.\n" );
  printf ( "  the sparsity pattern, rather than recomputing it.\n" );
/* 
  Defaults 
*/
  lwork = 0;
  nrhs = 1;
  equil = YES;  
  u = 1.0;
  trans = NOTRANS;
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
  Can use command line input to modify the defaults. 
*/
  parse_command_line ( argc, argv, &lwork, &u, &equil, &trans );
  options.Equil = equil;
  options.DiagPivotThresh = u;
  options.Trans = trans;

  printf ( "\n" );
  printf ( "  Length of work array LWORK = %d\n", lwork );
  printf ( "  Equilibration option EQUIL = %d\n", equil );
  printf ( "  Diagonal pivot threshhold value U = %f\n", u );
  printf ( "  Tranpose option TRANS = %d\n", trans );

  if ( 0 < lwork ) 
  {
    work = SUPERLU_MALLOC ( lwork );
    if ( !work )
    {
      ABORT ( "SUPERLU_MALLOC cannot allocate work[]" );
    }
  }
/* 
  Read matrix A from a file in Harwell-Boeing format.
*/
  dreadhb ( &m, &n, &nnz, &a, &asub, &xa );

  a1 = doubleMalloc ( nnz );
  if ( !a1 ) 
  {
    ABORT("Malloc fails for a1[].");
  }

  asub1 = intMalloc ( nnz );
  if ( !asub1 ) 
  {
    ABORT("Malloc fails for asub1[].");
  }

  xa1 = intMalloc ( n+1 );
  if ( !xa1 ) 
  {
    ABORT("Malloc fails for xa1[].");
  }

  for ( i = 0; i < nnz; i++ ) 
  {
    a1[i] = a[i];
    asub1[i] = asub[i];
  }
  for ( i = 0; i < n+1; i++ ) 
  {
    xa1[i] = xa[i];
  }
/*
  Create storage for a column compressed matrix.
*/
  dCreate_CompCol_Matrix ( &A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE );
  Astore = A.Store;
  printf ( "\n" );
  printf ( "  Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz );
/*
  Set up space for the right hand side and solution.
*/  
  rhsb = doubleMalloc ( m * nrhs );
  if ( !rhsb ) 
  {
    ABORT ( "Malloc fails for rhsb[]." );
  }

  rhsb1 = doubleMalloc ( m * nrhs );
  if ( !rhsb1 ) 
  {
    ABORT ( "Malloc fails for rhsb1[]." );
  }

  rhsx = doubleMalloc ( m * nrhs );
  if ( !rhsx ) 
  {
    ABORT ( "Malloc fails for rhsx[]." );
  }

  dCreate_Dense_Matrix ( &B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE );

  dCreate_Dense_Matrix ( &X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE );

  xact = doubleMalloc ( n * nrhs );
  if ( !xact ) 
  {
    ABORT ( "Malloc fails for xact." );
  }
  ldx = n;
  dGenXtrue ( n, nrhs, xact, ldx );
  dFillRHS ( trans, nrhs, xact, ldx, &A, &B );
  for ( j = 0; j < nrhs; j++ )
  {
    for ( i = 0; i < m; i++ ) 
    {
      rhsb1[i+j*m] = rhsb[i+j*m];
    }
  }
    
  perm_c = intMalloc ( n );
  if ( !perm_c ) 
  {
    ABORT("Malloc fails for perm_c[].");
  }

  perm_r = intMalloc ( m );
  if ( !perm_r ) 
  {
    ABORT("Malloc fails for perm_r[].");
  }

  etree = intMalloc(n);
  if ( !etree ) 
  {
    ABORT("Malloc fails for etree[].");
  }

  R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double));
  if ( !R )
  { 
    ABORT("SUPERLU_MALLOC fails for R[].");
  }

  C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double));
  if ( !C )
  {
    ABORT("SUPERLU_MALLOC fails for C[].");
  }

  ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double));
  if ( !ferr )
  {
    ABORT("SUPERLU_MALLOC fails for ferr[].");
  }

  berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double));
  if ( !berr ) 
  {
    ABORT("SUPERLU_MALLOC fails for berr[].");
  }
/* 
  Initialize the statistics variables. 
*/
  StatInit ( &stat );
/*
  Solve the first linear system, A * X = B.
*/
  dgssvx ( &options, &A, perm_c, perm_r, etree, equed, R, C,
    &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
    &mem_usage, &stat, &info );

  printf ( "\n" );
  printf ( "  First system: DGSSVX returns INFO = %d\n", info );

  if ( info == 0 || info == n+1 )
  {
    sol = (double*) ((DNformat*) X.Store)->nzval; 

    if ( options.PivotGrowth ) 
    {
      printf ( "\n" );
      printf ( "  Reciprocal pivot growth = %e\n", rpg );
    }

    if ( options.ConditionNumber )
    {
      printf ( "\n" );
      printf ( "  Reciprocal condition number = %e\n", rcond );
    }

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;

    printf ( "  Number of nonzeros in factor L = %d\n", Lstore->nnz );
    printf ( "  Number of nonzeros in factor U = %d\n", Ustore->nnz );
    printf ( "  Number of nonzeros in L+U = %d\n", 
      Lstore->nnz + Ustore->nnz - n );

    printf ( "\n" );
    printf ( "  L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
      mem_usage.expansions );

    if ( options.IterRefine )
    {
      printf ( "\n" );
      printf ( "  Iterative Refinement:\n");
      printf ( "%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
      for ( i = 0; i < nrhs; i++ )
      {
        printf ( "%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i] );
      }
    }
    fflush ( stdout );

  }
  else if ( 0 < info && lwork == -1 )
  {
    printf ( "\n" );
    printf ( "**Estimated memory: %d bytes\n", info - n );
  }

  if ( options.PrintStat ) 
  {
    StatPrint ( &stat );
  }
  StatFree ( &stat );

  Destroy_CompCol_Matrix ( &A );
  Destroy_Dense_Matrix ( &B );
/* 
  Deallocate storage associated with L and U. 
*/
  if ( 0 <= lwork )
  { 
    Destroy_SuperNode_Matrix ( &L );
    Destroy_CompCol_Matrix ( &U );
  }
/*
  Now solve a new system A1 * X = B1, noting, however,
  that the sparsity pattern of A1 is the same as that of A.
*/
  options.Fact = SamePattern;
/*
  Reinitialize the statistics.
*/
  StatInit ( &stat );

  dCreate_CompCol_Matrix ( &A1, m, n, nnz, a1, asub1, xa1,
    SLU_NC, SLU_D, SLU_GE );

  dCreate_Dense_Matrix ( &B1, m, nrhs, rhsb1, m, SLU_DN, SLU_D, SLU_GE );

  dgssvx ( &options, &A1, perm_c, perm_r, etree, equed, R, C,
    &L, &U, work, lwork, &B1, &X, &rpg, &rcond, ferr, berr,
    &mem_usage, &stat, &info );

  printf ( "\n" );
  printf ( "  Second system: DGSSVX returns INFO = %d\n", info );

  if ( info == 0 || info == n+1 )
  {
    sol = (double*) ((DNformat*) X.Store)->nzval; 

    if ( options.PivotGrowth ) 
    {
      printf ( "\n" );
      printf ( "  Reciprocal pivot growth = %e\n", rpg );
    }

    if ( options.ConditionNumber )
    {
      printf ( "\n" );
      printf ( "  Reciprocal condition number = %e\n", rcond );
    }

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;

    printf ( "  Number of nonzeros in factor L = %d\n", Lstore->nnz );
    printf ( "  Number of nonzeros in factor U = %d\n", Ustore->nnz );
    printf ( "  Number of nonzeros in L+U = %d\n", 
      Lstore->nnz + Ustore->nnz - n );

    printf ( "\n" );
    printf ( "  L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
      mem_usage.expansions );

    if ( options.IterRefine )
    {
      printf ( "\n" );
      printf ( "  Iterative Refinement:\n" );
      printf ( "%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR" );
      for ( i = 0; i < nrhs; i++ )
      {
        printf ( "%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i] );
      }
    }
    fflush ( stdout );
  } 
  else if ( 0 < info && lwork == -1 )
  {
    printf ( "\n" );
    printf ( "**Estimated memory: %d bytes\n", info - n );
  }

  if ( options.PrintStat ) 
  {
    StatPrint ( &stat );
  }
  StatFree ( &stat );

  SUPERLU_FREE ( xact );
  SUPERLU_FREE ( etree );
  SUPERLU_FREE ( perm_r );
  SUPERLU_FREE ( perm_c );
  SUPERLU_FREE ( R );
  SUPERLU_FREE ( C );
  SUPERLU_FREE ( ferr );
  SUPERLU_FREE ( berr );
  Destroy_CompCol_Matrix ( &A1 );
  Destroy_Dense_Matrix ( &B1 );
  Destroy_Dense_Matrix ( &X );

  if ( 0 <= lwork )
  {
    Destroy_SuperNode_Matrix ( &L );
    Destroy_CompCol_Matrix ( &U );
  }
/*
  Say goodbye.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_D4:\n" );
  printf ( "  Normal end of execution.\n");

  return 0;
}
/**********************************************************************/

static void parse_command_line ( int argc, char *argv[], int *lwork,
  double *u, yes_no_t *equil, trans_t *trans )

/**********************************************************************/
/*
  Purpose:

    PARSE_COMMAND_LINE parses the command line.

  Discussion:

    The user can set command line options, including:
      -l <int>,    length of work array;
      -u <double>,  diagonal pivot threshhold value;
      -e <0 or 1>, equilibration option;
      -t <0 or 1>, tranpose option.

  Modified:

    26 April 2004

  Reference:

    James Demmel, John Gilbert, Xiaoye Li,
    SuperLU Users's Guide,
    Sections 1 and 2.
*/
{
  int c;
  extern char *optarg;

  while ( (c = getopt(argc, argv, "hl:u:e:t:")) != EOF ) 
  {
    switch (c) 
    {
      case 'h':
        printf ( "\n" );
        printf ( "Options:\n");
        printf ( "  -l <int> - length of work[*] array\n");
        printf ( "  -u <double> - pivoting threshold\n");
        printf ( "  -e <0 or 1> - equilibrate or not\n");
        printf ( "  -t <0 or 1> - solve transposed system or not\n");
        exit ( 1 );
        break;
      case 'l': 
        *lwork = atoi ( optarg );
        break;
      case 'u': 
        *u = ( double ) atof ( optarg ); 
        break;
      case 'e': 
        *equil = atoi ( optarg ); 
        break;
      case 't': 
        *trans = atoi ( optarg );
        break;
    }
  }
  return;
}
