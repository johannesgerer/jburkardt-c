# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "st_io.h"

/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

void i4vec_dec ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_DEC decrements an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, int A[N], the vector to be decremented.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] - 1;
  }
  return;
}
/******************************************************************************/

void i4vec_inc ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_INC increments an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, int A[N], the vector to be incremented.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] + 1;
  }
  return;
}
/******************************************************************************/

int i4vec_max ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MAX returns the value of the maximum element in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int IVEC_MAX, the value of the maximum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
/******************************************************************************/

int i4vec_min ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MIN returns the minimum element in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int I4VEC_MIN, the value of the minimum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}
/******************************************************************************/

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

/******************************************************************************/
/*
  Purpose:

    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.

  Discussion:

    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.

    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2004

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the length of the input list.

    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.

    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.

    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
/*
  INDX = 0: This is the first call.
*/
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
/*
  INDX < 0: The user is returning the results of a comparison.
*/
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
/*
  0 < INDX: the user was asked to make an interchange.
*/
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
/******************************************************************************/

void st_data_read ( char *input_filename, int m, int n, int nst, int ist[], 
  int jst[], double ast[] )

/******************************************************************************/
/*
  Purpose:

    ST_DATA_READ reads the data of an ST file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the ST file.

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NST, the number of nonzeros.

    Output, int IST[NST], JST[NST], the row and column indices.

    Output, double AST[NST], the nonzero values.
*/
{
  double aij;
  int i;
  FILE *input;
  int j;
  int k;
  int num;

  input = fopen ( input_filename, "rt" );

  for ( k = 0; k < nst; k++ )
  {
    num = fscanf ( input, "%i%i%lf", &i, &j, &aij );

    if ( num != 3 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "ST_DATA_READ - Fatal error!\n" );
      fprintf ( stderr, "  I/O error reading data index %d\n", k );
      exit ( 1 );
    }

    ist[k] = i;
    jst[k] = j;
    ast[k] = aij;
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void st_header_print ( int i_min, int i_max, int j_min, int j_max, int m, 
  int n, int nst )

/******************************************************************************/
/*
  Purpose:

    ST_HEADER_PRINT prints the header of an ST file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int I_MIN, I_MAX, the minimum and maximum rows.

    Input, int J_MIN, J_MAX, the minimum and maximum columns.

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NST, the number of nonzeros.
*/
{
  printf ( "\n" );
  printf ( "  Sparse Triplet (ST) header information:\n" );
  printf ( "\n" );
  printf ( "  Minimum row index I_MIN = %d\n", i_min );
  printf ( "  Maximum row index I_MAX = %d\n", i_max );
  printf ( "  Minimum col index J_MIN = %d\n", j_min );
  printf ( "  Maximum col index J_MAX = %d\n", j_max );
  printf ( "  Number of rows        M = %d\n", m );
  printf ( "  Number of columns     N = %d\n", n );
  printf ( "  Number of nonzeros  NST = %d\n", nst );

  return;
}
/******************************************************************************/

void st_header_read ( char *input_filename, int *i_min, int *i_max, int *j_min, 
  int *j_max, int *m, int *n, int *nst )

/******************************************************************************/
/*
  Purpose:

    ST_HEADER_READ reads the header of an ST file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the ST file.

    Input, int *I_MIN, *I_MAX, the minimum and maximum rows.

    Input, int *J_MIN, *J_MAX, the minimum and maximum columns.

    Output, int *M, the number of rows.

    Output, int *N, the number of columns.

    Output, int *NST, the number of nonzeros.
*/
{
  double aij;
  int i;
  const int i4_huge = 2147483647;
  FILE *input;
  int j;
  int num;

  input = fopen ( input_filename, "rt" );

  *nst = 0;
  *i_min = + i4_huge;
  *i_max = - i4_huge;
  *j_min = + i4_huge;
  *j_max = - i4_huge;

  for ( ; ; )
  {
    num = fscanf ( input, "%i%i%lf", &i, &j, &aij );

    if ( num != 3 )
    {
      break;
    }

    *nst = *nst + 1;
    *i_min = i4_min ( *i_min, i );
    *i_max = i4_max ( *i_max, i );
    *j_min = i4_min ( *j_min, j );
    *j_max = i4_max ( *j_max, j );
  }

  fclose ( input );

  *m = *i_max - *i_min + 1;
  *n = *j_max - *j_min + 1;

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

void st_print_some ( int i_min, int i_max, int j_min, int j_max, int nst, 
  int  ist[], int jst[], double ast[], char *title )

/******************************************************************************/
/*
  Purpose:

    ST_PRINT_SOME prints some of a sparse matrix in ST format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int I_MIN, IMAX, the first and last rows to print.

    Input, int J_MIN, J_MAX, the first and last columns 
    to print.

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
    if ( i_min <= ist[k] && ist[k] <= i_max &&
         j_min <= jst[k] && jst[k] <= j_max )
    {
      printf ( "  %4d  %4d  %4d  %16.8g\n", k, ist[k], jst[k], ast[k] );
    }
  }

  return;
}
/******************************************************************************/

void st_sort_a ( int m, int n, int nst, int ist[], int jst[], double ast[] )

/******************************************************************************/
/*
  Purpose:

    ST_SORT_A sorts the entries of an ST matrix by column.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NST, the number of nonzeros.

    Input/output, int IST[NST], JST[NST], the row and column indices.

    Input/output, double AST[NST], the nonzero values.
*/
{
  double aij;
  int cij;
  int i;
  int indx;
  int isgn;
  int j;
  int rij;
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( nst, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      rij    = ist[i-1];
      ist[i-1] = ist[j-1];
      ist[j-1] = rij;

      cij    = jst[i-1];
      jst[i-1] = jst[j-1];
      jst[j-1] = cij;

      aij      = ast[i-1];
      ast[i-1] = ast[j-1];
      ast[j-1] = aij;
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      if ( jst[i-1] == jst[j-1] )
      {
        if ( ist[i-1] < ist[j-1] )
        {
          isgn = - 1;
        }
        else if ( ist[i-1] == ist[j-1] )
        {
          isgn = 0;
        }
        else if ( ist[j-1] < ist[i-1] )
        {
          isgn = + 1;
        }
      }
      else if ( jst[i-1] < jst[j-1] )
      {
        isgn = - 1;
      }
      else if ( jst[j-1] < jst[i-1] )
      {
        isgn = + 1;
      }
    }
    else if ( indx == 0 )
    {
      break;
    }
  }
  return;
}
/******************************************************************************/

void st_transpose ( int *m, int *n, int nst, int ist[], int jst[], 
  double ast[] )

/******************************************************************************/
/*
  Purpose:

    ST_TRANSPOSE transposes an ST matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2014

  Author:

    John Burkardt

  Parameters:

    Input/output, int *M, the number of rows.

    Input/output, int *N, the number of columns.

    Input, int NST, the number of nonzeros.

    Input/output, int IST[NST], JST[NST], the row and column indices.

    Input, double AST[NST], the nonzero values.
*/
{
  int k;
  int t;
  
  t  = *m;
  *m = *n;
  *n = t;

  for ( k = 0; k < nst; k++ )
  {
    t      = ist[k];
    ist[k] = jst[k];
    jst[k] = t;
  }

  return;
}
/******************************************************************************/

void st_write ( char *output_filename, int m, int n, int nst, int ist[], 
  int jst[], double ast[] )

/******************************************************************************/
/*
  Purpose:

    ST_WRITE writes an ST file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the name of the ST file.

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NST, the number of nonzeros.

    Input, int IST[NST], JST[NST], the row and column indices.

    Input, double AST[NST], the nonzero values.
*/
{
  FILE *output;
  int k;

  output = fopen ( output_filename, "wt" );

  for ( k = 0; k < nst; k++ )
  {
    fprintf ( output, "  %d  %d  %g\n", ist[k], jst[k], ast[k] );
  }

  fclose ( output );

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
