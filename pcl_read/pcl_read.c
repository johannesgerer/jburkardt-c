# include <ctype.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# define LINE_MAX_LEN 1024
char input[LINE_MAX_LEN];

int main ( int argc, char *argv[] );
int pcl_to_pure ( char *file_in_name, char *file_out_name );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    PCL_READ reads a PCL file for Arabidopsis and writes out a pure data file.

  Discussion:

    The PCL file has a (long) first line containing titles.  There follow
    N lines, each containing 3 labels, and 14 floating values, separated
    by TAB characters.

    The output file contains N lines, each containing 14 floating values,
    separated by spaces.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2003

  Author:

    John Burkardt
*/
{
  char file_in_name[81];
  char file_out_name[81];
  int iarg;
  int status;
  char *string;

  timestamp ( );
  printf ( "\n" );
  printf ( "PCL_READ:\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Read an Arabidopsis PCL file.\n" );
  printf ( "  Write a copy containing just the numeric data.\n" );

  iarg = 0;
/*
  Determination of the input file name.
*/
  if ( argc <= 1 )
  {
    printf ( "\n" );
    printf ( "Enter the input PCL file name.\n" );

    string = fgets ( input, LINE_MAX_LEN, stdin );

    if ( string == NULL )
    {
      return 1;
    }

    sscanf ( input, "%s", file_in_name );

  }
  else
  {
    iarg = iarg + 1;
    strcpy ( file_in_name, argv[iarg] );
  }
  printf ( "\n" );
  printf ( "Data will be read from %s.\n", file_in_name );
/*
  Determine the output file name.
*/
  if ( argc <= 2 )
  {
    printf ( "\n" );
    printf ( "Enter the output file name.\n" );

    string = fgets ( input, LINE_MAX_LEN, stdin );

    if ( string == NULL )
    {
      return 1;
    }

    sscanf ( input, "%s", file_out_name );

  }
  else
  {
    iarg = iarg + 1;
    strcpy ( file_out_name, argv[iarg] );
  }

  printf ( "\n" );
  printf ( "Data will be written to %s.\n", file_out_name );

  status = pcl_to_pure ( file_in_name, file_out_name );

  printf ( "\n" );
  printf ( "PCL_READ:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return status;
}
/******************************************************************************/

int pcl_to_pure ( char *file_in_name, char *file_out_name )

/******************************************************************************/
/*
  Purpose:

    PCL_TO_PURE reads data from a PCL file and writes the numeric data
    to another file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_IN_NAME, the name of the input file.

    Input, char *FILE_OUT_NAME, the name of the output file.

    Output, int PCL_TO_PURE, is 0 if there was no error, nonzero otherwise.
*/
{
  FILE *file_in;
  FILE *file_out;
  int i;
  int skip;
  int status;
  char *string;
  int text_num = 0;
  int width;
  float x[14];
/*
  Open the input file.
*/
  file_in = fopen ( file_in_name, "r" );

  if ( file_in == NULL )
  {
    printf ( "\n" );
    printf ( "PCL_READ - Fatal error!\n" );
    printf ( "  Could not open the input file '%s'!\n", file_in_name );
    return 1;
  }
/*
  Open the output file.
*/
  file_out = fopen ( file_out_name, "w" );

  if ( file_out == NULL )
  {
    printf ( "\n" );
    printf ( "PCL_READ - Fatal error!\n" );
    printf ( "  Could not open the output file '%s'!\n", file_out_name );
    return 1;
  }

  for ( ;; )
  {
/*
  Read a line from the file.
*/
    string = fgets ( input, LINE_MAX_LEN, file_in );

    if ( string == NULL )
    {
      break;
    }

    text_num = text_num + 1;
/*
  First line is the title.
*/
    if ( text_num == 1 )
    {
      printf ( "\n" );
      printf ( "PCL_TO_PURE: Title:\n" );
      printf ( "  %s\n", input );
      continue;
    }
/*
  All other lines have the form of 17 fields separated by tabs.
  The first 3 fields are identifiers.  The following 14 are numeric.
*/
    skip = 0;
    for ( i = 0; i < 3; i++ )
    {
      while ( *string != '\t' )
      {
        string = string + 1;
        skip = skip + 1;
      }

      string = string + 1;
      skip = skip + 1;
    }
/*
  Read the next word from the line.
*/
    for ( i = 0; i < 14; i++ )
    {
      sscanf ( string, "%f%n", &x[i], &width );
      string = string + width;
    }

    for ( i = 0; i < 14; i++ )
    {
      fprintf ( file_out, "%f ", x[i] );
    }
    fprintf ( file_out, "\n" );

  }

  status = fclose ( file_in );
  status = fclose ( file_out );

  printf ( "\n" );
  printf ( "PCL_TO_PURE:\n" );
  printf ( "  Number of input lines read was %d.\n", text_num );

  return 0;
}
/**********************************************************************/

void timestamp ( void )

/**********************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    May 31 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 August 2002

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 29

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  if ( len != 0 ) 
  {
    printf ( "%s\n", time_buffer );
  }

  return;
# undef TIME_SIZE
}
