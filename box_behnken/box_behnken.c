# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

# include "box_behnken.h"

/******************************************************************************/

double *box_behnken ( int dim_num, int x_num, double range[] )

/******************************************************************************/
/*
  Purpose:

    BOX_BEHNKEN returns a Box-Behnken design for the given number of factors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2012

  Author:

    John Burkardt

  Reference:

    George Box, Donald Behnken,
    Some new three level designs for the study of quantitative variables,
    Technometrics,
    Volume 2, pages 455-475, 1960.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int X_NUM, the number of elements of the design.
    X_NUM should be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.

    Input, double RANGE[DIM_NUM*2], the minimum and maximum
    value for each component.

    Output, double BOX_BEHNKEN[DIM_NUM*X_NUM], the elements of the design.
*/
{
  int i;
  int i2;
  int j;
  int last_low;
  double *x;
/*
  Ensure that the range is legal.
*/
  for ( i = 0; i < dim_num; i++ )
  {
    if ( range[i+1*dim_num] <= range[i+0*dim_num] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "BOX_BEHNKEN - Fatal error!\n" );
      fprintf ( stderr, "  RANGE[%d,1] <= RANGE[%d,0].\n", i, i );
      exit ( 1 );
    }
  }

  x = ( double * ) malloc ( dim_num * x_num * sizeof ( double ) );
/*
  The first point is the center.
*/
  j = 0;

  for ( i = 0; i < dim_num; i++ )
  {
    x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
  }
/*
  For subsequent elements, one entry is fixed at the middle of the range.
  The others are set to either extreme.
*/
  for ( i = 0; i < dim_num; i++ )
  {
    j = j + 1;
    for ( i2 = 0; i2 < dim_num; i2++ )
    {
      x[i2+j*dim_num] = range[i2+0*dim_num];
    }
    x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
/*
  The next element is made by finding the last low value, making it
  high, and all subsequent high values low.
*/
    for ( ; ; )
    {
      last_low = -1;

      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        if ( x[i2+j*dim_num] == range[i2+0*dim_num] )
        {
          last_low = i2;
        }
      }

      if ( last_low == -1 )
      {
        break;
      }

      j = j + 1;
      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        x[i2+j*dim_num] = x[i2+(j-1)*dim_num];
      }
      x[last_low+j*dim_num] = range[last_low+1*dim_num];

      for ( i2 = last_low + 1; i2 < dim_num; i2++ )
      {
        if ( x[i2+j*dim_num] == range[i2+1*dim_num] )
        {
          x[i2+j*dim_num] = range[i2+0*dim_num];
        }
      }
    }
  }
  return x;
}
/******************************************************************************/

int box_behnken_size ( int dim_num )

/******************************************************************************/
/*
  Purpose:

    BOX_BEHNKEN_SIZE returns the size of a Box-Behnken design.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2012

  Author:

    John Burkardt

  Reference:

    George Box, Donald Behnken,
    Some new three level designs for the study of quantitative variables,
    Technometrics,
    Volume 2, pages 455-475, 1960.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Output, int X_NUM, the number of elements of the design.
    X_NUM will be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
*/
{
  int x_num;

  if ( 1 <= dim_num )
  {
    x_num = 1 + dim_num * i4_power ( 2, dim_num - 1 );
  }
  else
  {
    x_num = -1;
  }

  return x_num;
}
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

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14f", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

