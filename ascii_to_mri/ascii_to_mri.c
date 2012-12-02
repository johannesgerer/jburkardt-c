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

    MAIN is the main program for ASCII_TO_MRI.

  Usage:

    ascii_to_mri prefix

    where

    * prefix.txt is the ASCII file to be read;
    * prefix.dat is the binary MRI file to be created.

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
  unsigned long int value_max;
  unsigned long int value_min;
  unsigned short int value_short;
  int value_num;

  printf ( "\n" );
  timestamp ( );
  printf ( "\n" );
  printf ( "ASCII_TO_MRI:\n" );
  printf ( "  C version\n" );
  printf ( "  Convert ASCII data to binary MRI format.\n" );
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
  strcat ( input_filename, ".txt" );
  strcpy ( output_filename, prefix );
  strcat ( output_filename, ".dat" );
/*
  Open the files.
*/
  input_unit = fopen ( input_filename, "r" );

  if ( !input_unit )
  {
    printf ( "\n" );
    printf ( "ASCII_TO_MRI - Fatal error!\n" );
    printf ( "  Cannot open input file \"%s\".\n", input_filename );
    exit ( 1 );
  }

  output_unit = fopen ( output_filename, "wb" );

  if ( !output_unit )
  {
    printf ( "\n" );
    printf ( "ASCII_TO_MRI - Fatal error!\n" );
    printf ( "  Cannot open output file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  The conversion simply takes unsigned longs and writes them to
  a binary file as unsigned shorts.  

  So the numerical limits here could be replaced by a WHILE statement.
*/
  value_num = 0;

  for ( i = 0; i < 64*64*26; i++ )
  {
    fscanf ( input_unit, "%u\n", &value_long );

    if ( i == 0 )
    {
      value_min = value_long;
      value_max = value_long;
    }
    else
    {
       if ( value_long < value_min )
       {
         value_min = value_long;
       }
       if ( value_max < value_long )
       {
         value_max = value_long;
       }
    }

    value_short = ( unsigned short ) value_long;

    fwrite ( &value_short, sizeof ( unsigned short ), 1, output_unit );

    value_num = value_num + 1;
  }

  fclose ( input_unit );
  fclose ( output_unit );

  printf ( "\n" );
  printf ( "  Read the ASCII file \"%s\".\n", input_filename );
  printf ( "  Wrote the binary file \"%s\".\n", output_filename );
  printf ( "  Wrote %d values.\n", value_num );
  printf ( "  Minimum value was %u.\n", value_min );
  printf ( "  Maximum value was %u.\n", value_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASCII_TO_MRI:\n" );
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
