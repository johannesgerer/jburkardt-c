# include "stdio.h"
# include "stdlib.h"
# include "string.h"
# include "time.h"

int main ( int argc, char *argv[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MRI_TO_ASCII.

  Usage:

    mri_to_ascii prefix

    where

    * prefix.dat is the binary MRI file to be read.
    * prefix.txt is the ASCII file to be created;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 December 2010

  Author:

    John Burkardt
*/
{
  int i;
  char input_filename[255];
  FILE *input_unit;
  char output_filename[255];
  FILE *output_unit;
  char prefix[255];
  unsigned long int value_long;
  unsigned short int value_max;
  unsigned short int value_min;
  unsigned short int value_short;
  int value_num;

  printf ( "\n" );
  timestamp ( );
  printf ( "\n" );
  printf ( "MRI_TO_ASCII:\n" );
  printf ( "  C version\n" );
  printf ( "  Convert MRI data to ASCII format.\n" );
/*
  Retrieve the common file prefix.
*/
  if ( 1 < argc )
  {
    strcpy ( prefix, argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the filename prefix:\n" );
    scanf ( "%s", prefix );
  }
/*
  Construct the filenames.
*/
  strcpy ( input_filename, prefix );
  strcat ( input_filename, ".dat" );
  strcpy ( output_filename, prefix );
  strcat ( output_filename, ".txt" );
/*
  Open the files.
*/
  input_unit = fopen ( input_filename, "rb" );

  if ( !input_unit )
  {
    printf ( "\n" );
    printf ( "MRI_TO_ASCII - Fatal error!\n" );
    printf ( "  Cannot open input file \"%s\".\n", input_filename );
    exit ( 1 );
  }

  output_unit = fopen ( output_filename, "wt" );

  if ( !output_unit )
  {
    printf ( "\n" );
    printf ( "MRI_TO_ASCII - Fatal error!\n" );
    printf ( "  Cannot open output file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  The conversion simply takes unsigned shorts and writes them to
  a text file as unsigned longs.  

  So the numerical limits here could be replaced by a WHILE statement.
*/
  value_num = 0;

  for ( i = 0; i < 64*64*26; i++ )
  {
    fread ( &value_short, sizeof ( unsigned short ), 1, input_unit );

    if ( i == 0 )
    {
      value_min = value_short;
      value_max = value_short;
    }
    else
    {
       if ( value_short < value_min )
       {
         value_min = value_short;
       }
       if ( value_max < value_short )
       {
         value_max = value_short;
       }
    }

    value_long = ( unsigned long ) value_short;

    fprintf ( output_unit, "%u\n", value_long );

    value_num = value_num + 1;
  }

  fclose ( input_unit );
  fclose ( output_unit );

  printf ( "\n" );
  printf ( "  Read the binary file \"%s\".\n", input_filename );
  printf ( "  Wrote the ASCII file \"%s\".\n", output_filename );
  printf ( "  Wrote %d values.\n", value_num );
  printf ( "  Minimum value was %u.\n", value_min );
  printf ( "  Maximum value was %u.\n", value_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MRI_TO_ASCII:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
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
