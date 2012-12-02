# include <stdio.h>
# include <stddef.h>
# include <stdlib.h>
# include <time.h>

# include "mgmres.h"

/******************************************************************************/

void atx_cr ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    ATX_CR computes A'*x for a matrix stored in sparse compressed row form.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    For this version of MGMRES, the row and column indices are assumed
    to use the C/C++ convention, in which indexing begins at 0.

    If your index vectors IA and JA are set up so that indexing is based 
    at 1, each use of those vectors should be shifted down by 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A'.

    Output, double W[N], the value of A'*X.
*/
{
  int i;
  int k;
  int k1;
  int k2;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
    k1 = ia[i];
    k2 = ia[i+1];
    for ( k = k1; k < k2; k++ )
    {
      w[ja[k]] = w[ja[k]] + a[k] * x[i];
    }
  }
  return;
}
/******************************************************************************/

void atx_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    ATX_ST computes A'*x for a matrix stored in sparse triplet form.

  Discussion:

    The matrix A is assumed to be stored in sparse triplet format.  Only
    the nonzero entries of A are stored.  For instance, the K-th nonzero
    entry in the matrix is stored by:

      A(K) = value of entry,
      IA(K) = row of entry,
      JA(K) = column of entry.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
    of the matrix values.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A'.

    Output, double W[N], the value of A'*X.
*/
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = ia[k];
    j = ja[k];
    w[j] = w[j] + a[k] * x[i];
  }

  return;
}
/******************************************************************************/

void ax_cr ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    AX_CR computes A*x for a matrix stored in sparse compressed row form.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    For this version of MGMRES, the row and column indices are assumed
    to use the C/C++ convention, in which indexing begins at 0.

    If your index vectors IA and JA are set up so that indexing is based 
    at 1, each use of those vectors should be shifted down by 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A.

    Output, double W[N], the value of A*X.
*/
{
  int i;
  int k;
  int k1;
  int k2;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
    k1 = ia[i];
    k2 = ia[i+1];
    for ( k = k1; k < k2; k++ )
    {
      w[i] = w[i] + a[k] * x[ja[k]];
    }
  }
  return;
}
/******************************************************************************/

void ax_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    AX_ST computes A*x for a matrix stored in sparse triplet form.

  Discussion:

    The matrix A is assumed to be stored in sparse triplet format.  Only
    the nonzero entries of A are stored.  For instance, the K-th nonzero
    entry in the matrix is stored by:

      A(K) = value of entry,
      IA(K) = row of entry,
      JA(K) = column of entry.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
    of the matrix values.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A.

    Output, double W[N], the value of A*X.
*/
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = ia[k];
    j = ja[k];
    w[i] = w[i] + a[k] * x[j];
  }

  return;
}
/******************************************************************************/

void diagonal_pointer_cr ( int n, int nz_num, int ia[], int ja[], int ua[] )

/******************************************************************************/
/*
  Purpose:

    DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    The array UA can be used to locate the diagonal elements of the matrix.

    It is assumed that every row of the matrix includes a diagonal element,
    and that the elements of each row have been ascending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.  On output,
    the order of the entries of JA may have changed because of the sorting.

    Output, int UA[N], the index of the diagonal element of each row.
*/
{
  int i;
  int j;
  int j1;
  int j2;
  int k;

  for ( i = 0; i < n; i++ )
  {
    ua[i] = -1;
    j1 = ia[i];
    j2 = ia[i+1];

    for ( j = j1; j < j2; j++ )
    {
      if ( ja[j] == i ) 
      {
        ua[i] = j;
      }
    }

  }
  return;
}
/******************************************************************************/

double **dmatrix ( int nrl, int nrh, int ncl, int nch )

/******************************************************************************/
/*
  Purpose:

    DMATRIX allocates a double matrix.

  Discussion:

    The matrix will have a subscript range m[nrl...nrh][ncl...nch] .

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Output, double **DMATRIX, a doubly-dimensioned array with
    the requested row and column ranges.
*/
{
  int i;
  double **m;
  int nrow = nrh - nrl + 1;
  int ncol = nch - ncl + 1;
/* 
  Allocate pointers to the rows.
*/
  m = ( double **) malloc ( (size_t) ( ( nrow + 1 ) * sizeof ( double* ) ) );

  if ( !m ) 
  {
    printf ( "\n" );
    printf ( "DMATRIX - Fatal error!\n" );
    printf ( "  Failure allocating pointers to rows.\n");
    exit ( 1 );
  }
  m = m + 1;
  m = m - nrl;
/* 
  Allocate each row and set pointers to them.
*/
  m[nrl] = ( double * ) malloc ( (size_t) ( ( nrow * ncol + 1 ) * sizeof ( double ) ) );

  if ( !m[nrl] ) 
  {
    printf ( "\n" );
    printf ( "DMATRIX - Fatal error!\n" );
    printf ( "  Failure allocating rows.\n");
    exit ( 1 );
  }
  m[nrl] = m[nrl] + 1;
  m[nrl] = m[nrl] - ncl;

  for ( i = nrl + 1; i <= nrh; i++ ) 
  { 
    m[i] = m[i-1] + ncol;
  }
/* 
  Return the pointer to the array of pointers to the rows;
*/
  return m;
}
/******************************************************************************/

void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch )

/******************************************************************************/
/*
  Purpose:

    FREE_DMATRIX frees a double matrix allocated by DMATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Input, double **M, the pointer to the doubly-dimensioned array,
    previously created by a call to DMATRIX.
*/
{
  free ( ( char* ) ( m[nrl] + ncl - 1 ) );
  free ( ( char* ) ( m + nrl - 1 ) );

  return;
}
/******************************************************************************/

void ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], int ua[],
  double l[] )

/******************************************************************************/
/*
  Purpose:

    ILU_CR computes the incomplete LU factorization of a matrix.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input, int UA[N], the index of the diagonal element of each row.

    Output, double L[NZ_NUM], the ILU factorization of A.
*/
{
  int *iw;
  int i;
  int j;
  int jj;
  int jrow;
  int jw;
  int k;
  double tl;

  iw = ( int * ) malloc ( n * sizeof ( int ) );
/*
  Copy A.
*/
  for ( k = 0; k < nz_num; k++ ) 
  {
    l[k] = a[k];
  }

  for ( i = 0; i < n; i++ ) 
  {
/*
  IW points to the nonzero entries in row I.
*/
    for ( j = 0; j < n; j++ )
    {
      iw[j] = -1;
    }

    for ( k = ia[i]; k <= ia[i+1] - 1; k++ ) 
    {
      iw[ja[k]] = k;
    }

    j = ia[i];
    do 
    {
      jrow = ja[j];
      if ( i <= jrow )
      {
        break;
      }
      tl = l[j] * l[ua[jrow]];
      l[j] = tl;
      for ( jj = ua[jrow] + 1; jj <= ia[jrow+1] - 1; jj++ ) 
      {
        jw = iw[ja[jj]];
        if ( jw != -1 ) 
        {
          l[jw] = l[jw] - tl * l[jj];
        }
      }
      j = j + 1;
    } while ( j <= ia[i+1] - 1 );

    ua[i] = j;

    if ( jrow != i ) 
    {
      printf ( "\n" );
      printf ( "ILU_CR - Fatal error!\n" );
      printf ( "  JROW != I\n" );
      printf ( "  JROW = %d\n", jrow );
      printf ( "  I    = %d\n", i );
      exit ( 1 );
    }

    if ( l[j] == 0.0 ) 
    {
      printf ( "\n" );
      printf ( "ILU_CR - Fatal error!\n" );
      printf ( "  Zero pivot on step I = \n", i );
      printf ( "  L[%d] = 0.0\n", j );
      exit ( 1 );
    }

    l[j] = 1.0 / l[j];
  }

  for ( k = 0; k < n; k++ ) 
  {
    l[ua[k]] = 1.0 / l[ua[k]];
  }

  free ( iw );

  return;
}
/******************************************************************************/

void lus_cr ( int n, int nz_num, int ia[], int ja[], double l[], int ua[], 
  double r[], double z[] )

/******************************************************************************/
/*
  Purpose:

    LUS_CR applies the incomplete LU preconditioner.

  Discussion:

    The linear system M * Z = R is solved for Z.  M is the incomplete
    LU preconditioner matrix, and R is a vector supplied by the user.
    So essentially, we're solving L * U * Z = R.

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double L[NZ_NUM], the matrix values.

    Input, int UA[N], the index of the diagonal element of each row.

    Input, double R[N], the right hand side.

    Output, double Z[N], the solution of the system M * Z = R.
*/
{
  int i;
  int j;
  double *w;

  w = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Copy R in.
*/
  for ( i = 0; i < n; i++ )
  {
    w[i] = r[i];
  }
/*
  Solve L * w = w where L is unit lower triangular.
*/
  for ( i = 1; i < n; i++ )
  {
    for ( j = ia[i]; j < ua[i]; j++ )
    {
      w[i] = w[i] - l[j] * w[ja[j]];
    }
  }
/*
  Solve U * w = w, where U is upper triangular.
*/
  for ( i = n - 1; 0 <= i; i-- ) 
  {
    for ( j = ua[i] + 1; j < ia[i+1]; j++ ) 
    {
      w[i] = w[i] - l[j] * w[ja[j]];
    }
    w[i] = w[i] / l[ua[i]];
  }
/*
  Copy Z out.
*/
  for ( i = 0; i < n; i++ )
  {
    z[i] = w[i];
  }

  free ( w );

  return;
}
/******************************************************************************/

void mgmres_st ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel )

/******************************************************************************/
/*
  Purpose:

    MGMRES_ST applies restarted GMRES to a matrix in sparse triplet form.

  Discussion:

    The matrix A is assumed to be stored in sparse triplet form.  Only 
    the nonzero entries of A are stored.  For instance, the 
    K-th nonzero entry in the matrix is stored by:

      A(K) = value of entry,
      IA(K) = row of entry,
      JA(K) = column of entry.

    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
    corrections to the code on 31 May 2007.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994.
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 2003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the linear system.

    Input, int NZ_NUM, the number of nonzero matrix values.

    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices of 
    the matrix values.

    Input, double A[], the matrix values.

    Input/output, double X[N]; on input, an approximation to
    the solution.  On output, an improved approximation.

    Input, double RHS[N], the right hand side of the linear system.

    Input, int ITR_MAX, the maximum number of (outer) iterations to take.

    Input, int MR, the maximum number of (inner) iterations to take.
    MR must be less than N.

    Input, double TOL_ABS, an absolute tolerance applied to the
    current residual.

    Input, double TOL_REL, a relative tolerance comparing the
    current residual to the initial residual.
*/
{
  double av;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double **h;
  double htmp;
  int i;
  int itr;
  int itr_used;
  int j;
  int k;
  int k_copy;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  double *s;
  double **v;
  int verbose = 1;
  double *y;

  if ( n < mr )
  {
    printf ( "\n" );
    printf ( "MGMRES_ST - Fatal error!\n" );
    printf ( "  N < MR.\n" );
    printf ( "  N = %d\n", n );
    printf ( "  MR = %d\n", mr );
    exit ( 1 );
  }

  itr_used = 0;

  c = ( double * ) malloc ( mr * sizeof ( double ) );
  g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
  h = dmatrix(0,mr,0,mr-1);
  r = ( double * ) malloc ( n * sizeof ( double ) );
  s = ( double * ) malloc ( mr * sizeof ( double ) );
  v = dmatrix(0,mr,0,n-1);
  y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );

  for ( itr = 0; itr < itr_max; itr++ ) 
  {
    ax_st ( n, nz_num, ia, ja, a, x, r );

    for ( i = 0; i < n; i++ )
    {
      r[i] = rhs[i] - r[i];
    }

    rho = sqrt ( r8vec_dot ( n, r, r ) );

    if ( verbose )
    {
      printf ( "  ITR = %8d  Residual = %e\n", itr, rho );
    }

    if ( itr == 0 )
    {
      rho_tol = rho * tol_rel;
    }

    for ( i = 0; i < n; i++ )
    {
      v[0][i] = r[i] / rho;
    }

    g[0] = rho;
    for ( i = 1; i < mr + 1; i++ )
    {
      g[i] = 0.0;
    }

    for ( i = 0; i < mr + 1; i++ )
    {
      for ( j = 0; j < mr; j++ ) 
      {
        h[i][j] = 0.0;
      }
    }

    for ( k = 0; k < mr; k++ )
    {
      k_copy = k;

      ax_st ( n, nz_num, ia, ja, a, v[k], v[k+1] );

      av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      for ( j = 0; j < k+1; j++ )
      {
        h[j][k] = r8vec_dot ( n, v[k+1], v[j] );
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
        }
      }

      h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      if ( ( av + delta * h[k+1][k] ) == av )
      {
        for ( j = 0; j < k+1; j++ )
        {
          htmp = r8vec_dot ( n, v[k+1], v[j] );
          h[j][k] = h[j][k] + htmp;
          for ( i = 0; i < n; i++ ) 
          {
            v[k+1][i] = v[k+1][i] - htmp * v[j][i];
          }
        }
        h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
      }

      if ( h[k+1][k] != 0.0 )
      {
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] / h[k+1][k];
        }
      }

      if ( 0 < k )
      {
        for ( i = 0; i < k + 2; i++ )
        {
          y[i] = h[i][k];
        }
        for ( j = 0; j < k; j++ ) 
        {
          mult_givens ( c[j], s[j], j, y );
        }
        for ( i = 0; i < k + 2; i++ ) 
        {
          h[i][k] = y[i];
        }
      }

      mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );
      c[k] = h[k][k] / mu;
      s[k] = -h[k+1][k] / mu;
      h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
      h[k+1][k] = 0.0;
      mult_givens ( c[k], s[k], k, g );

      rho = fabs ( g[k+1] );

      itr_used = itr_used + 1;

      if ( verbose )
      {
        printf ( "  K =   %8d  Residual = %e\n", k, rho );
      }

      if ( rho <= rho_tol && rho <= tol_abs )
      {
        break;
      }
    }

    k = k_copy;

    y[k] = g[k] / h[k][k];
    for ( i = k - 1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = i+1; j < k + 1; j++ ) 
      {
        y[i] = y[i] - h[i][j] * y[j];
      }
      y[i] = y[i] / h[i][i];
    }

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k + 1; j++ )
      {
        x[i] = x[i] + v[j][i] * y[j];
      }
    }

    if ( rho <= rho_tol && rho <= tol_abs ) 
    {
      break;
    }
  }
  if ( verbose )
  {
    printf ( "\n" );
    printf ( "MGMRES_ST:\n" );
    printf ( "  Iterations = %d\n", itr_used );
    printf ( "  Final residual = %e\n", rho );
  }

  free ( c );
  free ( g );
  free_dmatrix(h,0,mr,0,mr-1);
  free ( r );
  free ( s );
  free_dmatrix(v,0,mr,0,n-1);
  free ( y );

  return;
}
/******************************************************************************/

void mult_givens ( double c, double s, int k, double *g )

/******************************************************************************/
/*
  Purpose:

    MULT_GIVENS applies a Givens rotation to two vector elements.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, double C, S, the cosine and sine of a Givens
    rotation.

    Input, int K, indicates the location of the first vector entry.

    Input/output, double G[K+2], the vector to be modified.  On output,
    the Givens rotation has been applied to entries G(K) and G(K+1).
*/
{
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k+1];
  g2 = s * g[k] + c * g[k+1];

  g[k]   = g1;
  g[k+1] = g2;

  return;
}
/******************************************************************************/

void pmgmres_ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel )

/******************************************************************************/
/*
  Purpose:

    PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.

  Discussion:

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

    This routine uses the incomplete LU decomposition for the
    preconditioning.  This preconditioner requires that the sparse
    matrix data structure supplies a storage position for each diagonal
    element of the matrix A, and that each diagonal element of the
    matrix A is not zero.

    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
    corrections to the code on 31 May 2007.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994.
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 2003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the linear system.

    Input, int NZ_NUM, the number of nonzero matrix values.

    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
    of the matrix values.  The row vector has been compressed.

    Input, double A[NZ_NUM], the matrix values.

    Input/output, double X[N]; on input, an approximation to
    the solution.  On output, an improved approximation.

    Input, double RHS[N], the right hand side of the linear system.

    Input, int ITR_MAX, the maximum number of (outer) iterations to take.

    Input, int MR, the maximum number of (inner) iterations to take.
    MR must be less than N.

    Input, double TOL_ABS, an absolute tolerance applied to the
    current residual.

    Input, double TOL_REL, a relative tolerance comparing the
    current residual to the initial residual.
*/
{
  double av;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double **h;
  double htmp;
  int i;
  int itr;
  int itr_used;
  int j;
  int k;
  int k_copy;
  double *l;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  double *s;
  int *ua;
  double **v;
  int verbose = 1;
  double *y;

  itr_used = 0;

  c = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
  g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
  h = dmatrix ( 0, mr, 0, mr-1 );
  l = ( double * ) malloc ( ( ia[n] + 1 ) * sizeof ( double ) );
  r = ( double * ) malloc ( n * sizeof ( double ) );
  s = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
  ua = ( int * ) malloc ( n * sizeof ( int ) );
  v = dmatrix ( 0, mr, 0, n-1 );
  y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );

  rearrange_cr ( n, nz_num, ia, ja, a );

  diagonal_pointer_cr ( n, nz_num, ia, ja, ua );

  ilu_cr ( n, nz_num, ia, ja, a, ua, l );

  if ( verbose )
  {
    printf ( "\n" );
    printf ( "PMGMRES_ILU_CR\n" );
    printf ( "  Number of unknowns = %d\n", n );
  }

  for ( itr = 0; itr < itr_max; itr++ ) 
  {
    ax_cr ( n, nz_num, ia, ja, a, x, r );

    for ( i = 0; i < n; i++ ) 
    {
      r[i] = rhs[i] - r[i];
    }

    lus_cr ( n, nz_num, ia, ja, l, ua, r, r );

    rho = sqrt ( r8vec_dot ( n, r, r ) );

    if ( verbose )
    {
      printf ( "  ITR = %d  Residual = %e\n", itr, rho );
    }

    if ( itr == 0 )
    {
      rho_tol = rho * tol_rel;
    }

    for ( i = 0; i < n; i++ ) 
    {
      v[0][i] = r[i] / rho;
    }

    g[0] = rho;
    for ( i = 1; i < mr + 1; i++ ) 
    {
      g[i] = 0.0;
    }

    for ( i = 0; i < mr + 1; i++ ) 
    {
      for ( j = 0; j < mr; j++ ) 
      {
        h[i][j] = 0.0;
      }
    }

    for ( k = 0; k < mr; k++ )
    {
      k_copy = k;

      ax_cr ( n, nz_num, ia, ja, a, v[k], v[k+1] ); 

      lus_cr ( n, nz_num, ia, ja, l, ua, v[k+1], v[k+1] );

      av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      for ( j = 0; j <= k; j++ ) 
      {
        h[j][k] = r8vec_dot ( n, v[k+1], v[j] );
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
        }
      }
      h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      if ( ( av + delta * h[k+1][k]) == av ) 
      {
        for ( j = 0; j < k + 1; j++ )
        {
          htmp = r8vec_dot ( n, v[k+1], v[j] );
          h[j][k] = h[j][k] + htmp;
          for ( i = 0; i < n; i++ ) 
          {
            v[k+1][i] = v[k+1][i] - htmp * v[j][i];
          }
        }
        h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
      }

      if ( h[k+1][k] != 0.0 )
      {
        for ( i = 0; i < n; i++ )
        {
          v[k+1][i] = v[k+1][i] / h[k+1][k];
        } 
      }

      if ( 0 < k )  
      {
        for ( i = 0; i < k + 2; i++ ) 
        {
          y[i] = h[i][k];
        }
        for ( j = 0; j < k; j++ ) 
        {
          mult_givens ( c[j], s[j], j, y );
        }
        for ( i = 0; i < k + 2; i++ )
        {
          h[i][k] = y[i];
        }
      }
      mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );

      c[k] = h[k][k] / mu;
      s[k] = -h[k+1][k] / mu;
      h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
      h[k+1][k] = 0.0;
      mult_givens ( c[k], s[k], k, g );

      rho = fabs ( g[k+1] );

      itr_used = itr_used + 1;

      if ( verbose )
      {
        printf ( "  K   = %d  Residual = %e\n", k, rho );
      }

      if ( rho <= rho_tol && rho <= tol_abs )
      {
        break;
      }
    }

    k = k_copy;

    y[k] = g[k] / h[k][k];
    for ( i = k - 1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = i + 1; j < k + 1; j++ ) 
      {
        y[i] = y[i] - h[i][j] * y[j];
      }
      y[i] = y[i] / h[i][i];
    }
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k + 1; j++ )
      {
        x[i] = x[i] + v[j][i] * y[j];
      }
    }

    if ( rho <= rho_tol && rho <= tol_abs )
    {
      break;
    }
  }

  if ( verbose )
  {
    printf ( "\n" );
    printf ( "PMGMRES_ILU_CR:\n" );
    printf ( "  Iterations = %d\n", itr_used );
    printf ( "  Final residual = %e\n", rho );
  }

  free ( c );
  free ( g );
  free_dmatrix ( h, 0, mr, 0, mr-1 );
  free ( l );
  free ( r );
  free ( s );
  free ( ua );
  free_dmatrix ( v, 0, mr, 0, n-1 );
  free ( y );

  return;
}
/******************************************************************************/

double r8vec_dot ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double *r8vec_uniform_01 ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01 fills a double precision vector with unit pseudorandom values.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      unif = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
*/
{
  int i;
  int k;
  double *r;

  r = malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void rearrange_cr ( int n, int nz_num, int ia[], int ja[], double a[] )

/******************************************************************************/
/*
  Purpose:

    REARRANGE_CR sorts a sparse compressed row matrix.

  Discussion:

    This routine guarantees that the entries in the CR matrix
    are properly sorted.

    After the sorting, the entries of the matrix are rearranged in such
    a way that the entries of each column are listed in ascending order
    of their column values.

    The matrix A is assumed to be stored in compressed row format.  Only
    the nonzero entries of A are stored.  The vector JA stores the
    column index of the nonzero value.  The nonzero values are sorted
    by row, and the compressed row vector IA then has the property that
    the entries in A and JA that correspond to row I occur in indices
    IA[I] through IA[I+1]-1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[N+1], the compressed row index vector.

    Input/output, int JA[NZ_NUM], the column indices of the matrix values.
    On output, the order of the entries of JA may have changed because of
    the sorting.

    Input/output, double A[NZ_NUM], the matrix values.  On output, the
    order of the entries may have changed because of the sorting.
*/
{
  double dtemp;
  int i;
  int is;
  int itemp;
  int j;
  int j1;
  int j2;
  int k;

  for ( i = 0; i < n; i++ )
  {
    j1 = ia[i];
    j2 = ia[i+1];
    is = j2 - j1;

    for ( k = 1; k < is; k++ ) 
    {
      for ( j = j1; j < j2 - k; j++ ) 
      {
        if ( ja[j+1] < ja[j] ) 
        {
          itemp = ja[j+1];
          ja[j+1] =  ja[j];
          ja[j] =  itemp;

          dtemp = a[j+1];
          a[j+1] =  a[j];
          a[j] = dtemp;
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
