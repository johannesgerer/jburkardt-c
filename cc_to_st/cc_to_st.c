# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cc_to_st.h"

/******************************************************************************/

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT prints a sparse matrix in CC format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the matrix.

    Input, int N, the number of columns in the matrix.

    Input, int NCC, the number of CC elements.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, double ACC[NCC], the CC values.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J           A\n" );
  printf ( "  ----  ----  ----  ----------------\n" );
  printf ( "\n" );

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k )
      {
        j = j + 1;
      }
      printf ( "  %4d  %4d  %4d  %16.8g\n", k, i, j, acc[k] );
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j] <= k + 1 )
      {
        j = j + 1;
      }
      printf ( "  %4d  %4d  %4d  %16.8g\n", k + 1, i, j, acc[k] );
    }
  }

  return;
}
/******************************************************************************/

void cc_to_st ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  int *nst, int ist[], int jst[], double ast[] )

/******************************************************************************/
/*
  Purpose:

    CC_TO_ST converts sparse matrix information from CC to ST format.

  Discussion:

    Only JST actually needs to be computed.  The other three output 
    quantities are simply copies.  
    
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NCC, the number of CC elements.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the CC compressed columns.

    Input, double ACC[NCC], the CC values.

    Output, int NST, the number of ST elements.

    Output, int IST[NST], JST[NST], the ST rows and columns.

    Output, double AST[NST], the ST values.
*/
{
  int j;
  int jhi;
  int jlo;
  int k;
  int khi;
  int klo;

  *nst = 0;

  if ( ccc[0] == 0 )
  {
    jlo = 0;
    jhi = n - 1;
  
    for ( j = jlo; j <= jhi; j++ )
    {
      klo = ccc[j];
      khi = ccc[j+1] - 1;

      for ( k = klo; k <= khi; k++ )
      {
        ist[*nst] = icc[k];
        jst[*nst] = j;
        ast[*nst] = acc[k];
        *nst = *nst + 1;
      }
    }
  }
  else
  {
    jlo = 1;
    jhi = n;
  
    for ( j = jlo; j <= jhi; j++ )
    {
      klo = ccc[j-1];
      khi = ccc[j] - 1;

      for ( k = klo; k <= khi; k++ )
      {
        ist[*nst] = icc[k-1];
        jst[*nst] = j;
        ast[*nst] = acc[k-1];
        *nst = *nst + 1;
      }
    }
  }

  return;
}
/******************************************************************************/

void st_print ( int m, int n, int nst, int ist[], int jst[], double ast[], 
  char *title )

/******************************************************************************/
/*
  Purpose:

    ST_PRINT prints a sparse matrix in ST format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NST, the number of ST elements.

    Input, int IST[NST], JST[NST], the ST rows and columns.

    Input, double AST[NST], the ST values.

    Input, char *TITLE, a title.
*/
{
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J       A\n" );
  printf ( "  ----  ----  ----  --------------\n" );
  printf ( "\n" );
  for ( k = 0; k < nst; k++ )
  {
    printf ( "  %4d  %4d  %4d  %16.8g\n", k, ist[k], jst[k], ast[k] );
  }

  return;
}
/******************************************************************************/

void timestamp ( )

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
