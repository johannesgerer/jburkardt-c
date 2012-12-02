# include <stdlib.h>
# include <stdio.h>

# include "ssp_defs.h"

int main ( int argc, char *argv[] );

static void parse_command_line ( int argc, char *argv[], int *lwork,
  float *u, yes_no_t *equil, trans_t *trans );

/**********************************************************************/

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
  Purpose:

    SUPER_LU_S4 solves a sparse system read from a file using SGSSVX.

  Discussion:

    In this example, we show how the LU factors computed by one call to
    SGSSVX can be reused to solve a second linear system with the same
    matrix by calling SGSSVX with the appropriate options set.

    The sparse matrix is stored in a file using the Harwell-Boeing
    sparse matrix format.  The file should be assigned to the standard
    input of this program.  For instance, if the matrix is stored
    in the file "g10_rua.txt", the execution command might be:

      super_lu_s4 < g10_rua.txt

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
  float *a;
  NCformat *Astore;
  int *asub;
  SuperMatrix B;
  float *berr;
  float *C;
  char equed[1];
  yes_no_t equil;
  int *etree;
  float *ferr;
  int i;
  int info;
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
  float *R;
  float rcond;
  float *rhsb;
  float *rhsx;
  float rpg;
  float *sol;
  SuperLUStat_t stat;
  trans_t trans;
  SuperMatrix U;
  float u;
  NCformat *Ustore;
  void *work;
  SuperMatrix X;
  int *xa;
  float *xact;
/*
  Say hello.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_S4:\n" );
  printf ( "  Read a sparse matrix A from standard input,\n");
  printf ( "  stored in Harwell-Boeing Sparse Matrix format.\n" );
  printf ( "\n" );
  printf ( "  Using SGSSVX:\n" );
  printf ( "    Factor the matrix A;\n" );
  printf ( "    LATER, solve a linear system A * X = B.\n" );
/* 
  Defaults 
*/
  lwork = 0;
  nrhs = 1;
  equil = YES;  
  u = 1.0;
  trans = NOTRANS;
/* 
  Set the default values for options argument:
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
  sreadhb ( &m, &n, &nnz, &a, &asub, &xa );
/*
  Create storage for a column compressed matrix.
*/
  sCreate_CompCol_Matrix ( &A, m, n, nnz, a, asub, xa, SLU_NC, SLU_S, SLU_GE );
  Astore = A.Store;

  printf ( "\n" );
  printf ( "  Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz );
/*
  Set up space for the right hand side and solution.
*/
  rhsb = floatMalloc ( m * nrhs );
  if ( !rhsb )
  {
    ABORT ( "Malloc fails for rhsb[]." );
  }

  rhsx = floatMalloc ( m * nrhs );
  if ( !rhsx ) 
  {
    ABORT ( "Malloc fails for rhsx[]." );
  }

  sCreate_Dense_Matrix ( &B, m, nrhs, rhsb, m, SLU_DN, SLU_S, SLU_GE );
  sCreate_Dense_Matrix ( &X, m, nrhs, rhsx, m, SLU_DN, SLU_S, SLU_GE );
/*
  Set up the exact solution.
*/
  xact = floatMalloc ( n * nrhs );
  if ( !xact ) 
  {
    ABORT ( "SUPERLU_MALLOC cannot allocate xact[]" );
  }
  ldx = n;
  sGenXtrue ( n, nrhs, xact, ldx );
  sFillRHS ( trans, nrhs, xact, ldx, &A, &B );
    
  etree = intMalloc ( n );
  if ( !etree )
  {
    ABORT ( "Malloc fails for etree[]." );
  }

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

  R = (float *) SUPERLU_MALLOC ( A.nrow * sizeof(float) );
  if ( !R ) 
  {
    ABORT ( "SUPERLU_MALLOC fails for R[]." );
  }

  C = (float *) SUPERLU_MALLOC ( A.ncol * sizeof(float) );
  if ( !C )
  {
    ABORT ( "SUPERLU_MALLOC fails for C[]." );
  }

  ferr = (float *) SUPERLU_MALLOC ( nrhs * sizeof(float) );
  if ( !ferr )
  {
    ABORT ( "SUPERLU_MALLOC fails for ferr[]." );
  }

  berr = (float *) SUPERLU_MALLOC ( nrhs * sizeof(float) );
  if ( !berr ) 
  {
    ABORT ( "SUPERLU_MALLOC fails for berr[]." );
  }
/* 
  Initialize the statistics variables. 
*/
  StatInit ( &stat );
/* 
  Only perform the LU decomposition.

  Setting B.ncol = 0 indicates not to solve the system now.
*/
  sgssvx ( &options, &A, perm_c, perm_r, etree, equed, R, C,
    &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
    &mem_usage, &stat, &info );

  printf ( "\n" );
  printf ( "  LU factorization: SGSSVX returns INFO = %d\n", info );

  if ( info == 0 || info == n+1 )
  {

    if ( options.PivotGrowth ) 
    {
      printf ( "\n" );
      printf ( "  Reciprocal pivot growth = %e\n", rpg);
    }

    if ( options.ConditionNumber )
    {
      printf ( "\n" );
      printf ( "  Reciprocal condition number = %e\n", rcond );
    }

    Lstore = (SCformat *) L.Store;
    Ustore = (NCformat *) U.Store;
    printf ( "\n" );
    printf ( "  Number of nonzeros in factor L = %d\n", Lstore->nnz );
    printf ( "  Number of nonzeros in factor U = %d\n", Ustore->nnz );
    printf ( "  Number of nonzeros in L+U = %d\n", 
      Lstore->nnz + Ustore->nnz - n );
    printf ( "  L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
      mem_usage.expansions );
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
/*
  Now we solve the linear system using the computed LU factors.

  options.Fact = FACTORED indicates that the LU factors are already computed.
*/
  options.Fact = FACTORED;
  B.ncol = nrhs;
/* 
  Initialize the statistics variables. 
*/
  StatInit ( &stat );

  sgssvx ( &options, &A, perm_c, perm_r, etree, equed, R, C,
    &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
    &mem_usage, &stat, &info );

  printf ( "\n" );
  printf ( "  Triangular solve: SGSSVX returns INFO = %d\n", info );

  if ( info == 0 || info == n+1 )
  {
    sol = (float*) ((DNformat*) X.Store)->nzval; 

  if ( options.IterRefine )
  {
    printf ( "\n" );
    printf ( "  Iterative Refinement:\n");
    printf ( "%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
    for (i = 0; i < nrhs; i++ )
    {
        printf ( "%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[i], berr[i] );
    }
  }
  fflush ( stdout );
    
  } 
  else if ( info > 0 && lwork == -1 ) 
  {
    printf ( "\n" );
    printf ( "**Estimated memory: %d bytes\n", info - n );
  }

  if ( options.PrintStat ) 
  {
    StatPrint ( &stat );
  }
  StatFree ( &stat );

  SUPERLU_FREE ( rhsb );
  SUPERLU_FREE ( rhsx );
  SUPERLU_FREE ( xact );
  SUPERLU_FREE ( etree );
  SUPERLU_FREE ( perm_r );
  SUPERLU_FREE ( perm_c );
  SUPERLU_FREE ( R );
  SUPERLU_FREE ( C );
  SUPERLU_FREE ( ferr );
  SUPERLU_FREE ( berr );
  Destroy_CompCol_Matrix ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperMatrix_Store ( &X );

  if ( 0 <= lwork )
  {
    Destroy_SuperNode_Matrix ( &L );
    Destroy_CompCol_Matrix ( &U );
  }
/*
  Say goodbye.
*/
  printf ( "\n" );
  printf ( "SUPER_LU_S4:\n" );
  printf ( "  Normal end of execution.\n");

  return 0;
}
/**********************************************************************/

static void parse_command_line ( int argc, char *argv[], int *lwork,
  float *u, yes_no_t *equil, trans_t *trans )

/**********************************************************************/
/*
  Purpose:

    PARSE_COMMAND_LINE parses the command line.

  Discussion:

    The user can set command line options, including:
      -l <int>,    length of work array;
      -u <float>,  diagonal pivot threshhold value;
      -e <0 or 1>, equilibration option;
      -t <0 or 1>, tranpose option.

  Modified:

    25 April 2004

  Reference:

    James Demmel, John Gilbert, Xiaoye Li,
    SuperLU Users's Guide,
    Sections 1 and 2.
*/
{
  int c;
  extern char *optarg;

  while ( ( c = getopt ( argc, argv, "hl:u:e:t:" ) ) != EOF )
  {
    switch ( c )
    {
      case 'h':
        printf ( "\n" );
        printf ( "Options:\n");
        printf ( "  -l <int> - length of work[*] array\n");
        printf ( "  -u <float> - pivoting threshold\n");
        printf ( "  -e <0 or 1> - equilibrate or not\n");
        printf ( "  -t <0 or 1> - solve transposed system or not\n");
        exit ( 1 );
        break;
      case 'l': 
        *lwork = atoi ( optarg );
        break;
      case 'u': 
        *u = atof ( optarg ); 
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
