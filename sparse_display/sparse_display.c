# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "sparse_display.h"

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

void spy_file ( char *header, char *data_filename )

/******************************************************************************/
/*
  Purpose:

    SPY_FILE plots a sparsity pattern stored in a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *HEADER, the name to be used for the
    title of the plot, and as part of the names of the command
    and plot files.

    Input, char *DATA_FILENAME, the name of the file
    containing the indices of nonzero matrix entries.
*/
{
  char command_filename[255];
  FILE *command_unit;
  FILE *data_unit;
  int i;
  const int i4_huge = 2147483647;
  int j;
  int m0;
  int m1;
  int n0;
  int n1;
  int nz_num;
  char png_filename[255];
  int status;

  n0 = + i4_huge;
  n1 = - i4_huge;
  m0 = + i4_huge;
  m1 = - i4_huge;
  nz_num = 0;

  data_unit = fopen ( data_filename, "rt" );

  for ( ; ; )
  {
    status = fscanf ( data_unit, "%d%d", &i, &j );

    if ( status != 2 )
    {
      break;
    }

    nz_num = nz_num + 1;
    m0 = i4_min ( m0, i );
    m1 = i4_max ( m1, i );
    n0 = i4_min ( n0, j );
    n1 = i4_max ( n1, j );
  }

  fclose ( data_unit );
/*
  Create command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set term png\n" );

  strcpy ( png_filename, header );
  strcat ( png_filename, ".png" );
  fprintf ( command_unit, "set output '%s'\n", png_filename );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set xlabel '<--- J --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- I --->'\n" );
  
  fprintf ( command_unit, "set title '%d nonzeros for \"%s\"'\n", nz_num, header );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, 
    "plot [y=%d:%d] [x=%d:%d] '%s' with points pt 5\n",
    m0, m1, n0, n1, data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

void spy_ge ( int m, int n, double a[], char *header )

/******************************************************************************/
/*
  Purpose:

    SPY_GE plots a sparsity pattern for a general storage (GE) matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    in the matrix.

    Input, double A[M*N], the matrix.

    Input, char *HEADER, the name to be used for the
    title of the plot, and as part of the names of the data, command
    and plot files.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  int j;
  int nz_num;
  char png_filename[255];
/*
  Create data file.
*/
  strcpy ( data_filename, header );
  strcat ( data_filename, "_data.txt" );
  data_unit = fopen ( data_filename, "wt" );
  nz_num = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        fprintf ( data_unit, "%d  %d\n", j, i );
        nz_num = nz_num + 1;
      }
    }
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created sparsity data file '%s'\n", data_filename );
/*
  Create command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set term png\n" );

  strcpy ( png_filename, header );
  strcat ( png_filename, ".png" );
  fprintf ( command_unit, "set output '%s'\n", png_filename );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set xlabel '<--- J --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- I --->'\n" );
  fprintf ( command_unit, "set title '%d nonzeros for \"%s\"'\n", 
    nz_num, header );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, 
    "plot [x=0:%d] [y=%d:0] '%s' with points pt 5\n",
    n-1, m-1, data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
