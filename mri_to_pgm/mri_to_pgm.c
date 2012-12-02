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

    MAIN is the main program for MRI_TO_PGM.

  Usage:

    mri_to_pgm prefix

    where

    * prefix.dat is the binary MRI file to be read.
    * prefix_01.pgm through prefix_26.pgm are the PGM files created.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 December 2010

  Author:

    John Burkardt
*/
{
  int ascii = 1;
  int i;
  char input_filename[255];
  FILE *input_unit;
  int j;
  char output_filename[255];
  FILE *output_unit;
  char prefix[255];
  unsigned short int value;
  unsigned short int value_max;
  unsigned short int value_min;

  printf ( "\n" );
  timestamp ( );
  printf ( "\n" );
  printf ( "MRI_TO_PGM:\n" );
  printf ( "  C version\n" );
  printf ( "  Convert MRI data to PGM files.\n" );
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
  strcat ( output_filename, "_00.pgm" );
/*
  Determine the maximum value in the input file.
*/
  value_max = 0;

  input_unit = fopen ( input_filename, "rb" );

  if ( !input_unit )
  {
    printf ( "\n" );
    printf ( "MRI_TO_PGM - Fatal error!\n" );
    printf ( "  Cannot open input file \"%s\".\n", input_filename );
    exit ( 1 );
  }
  for ( j = 1; j <= 26; j++ )
  {
    for ( i = 0; i < 64*64; i++ )
    {
      fread ( &value, sizeof ( unsigned short ), 1, input_unit );

      if ( j == 1 && i == 0 )
      {
        value_min = value;
        value_max = value;
      }
      else 
      {
        if ( value < value_min )
        {
          value_min = value;
        }
        if ( value_max < value )
        {
          value_max = value;
        }
      }
    }
    fclose ( output_unit );
  }
  fclose ( input_unit );

  printf ( "\n" );
  printf ( "  Minimum value = %u\n", value_min );
  printf ( "  Maximum value = %d\n", value_max );
/*
  Open the input file.
*/
  input_unit = fopen ( input_filename, "rb" );

  if ( !input_unit )
  {
    printf ( "\n" );
    printf ( "MRI_TO_PGM - Fatal error!\n" );
    printf ( "  Cannot open input file \"%s\".\n", input_filename );
    exit ( 1 );
  }
/*
  Create the PGM files.
*/
  for ( j = 1; j <= 26; j++ )
  {
    sprintf ( output_filename, "%s_%02d.pgm", prefix, j );

    output_unit = fopen ( output_filename, "wt" );

    if ( !output_unit )
    {
      printf ( "\n" );
      printf ( "MRI_TO_PGM - Fatal error!\n" );
      printf ( "  Cannot open output file \"%s\".\n", output_filename );
      exit ( 1 );
    }

    if ( ascii )
    {
      fprintf ( output_unit, "P2\n" );
      fprintf ( output_unit, "64  64\n" );
      fprintf ( output_unit, "%u\n", value_max );
    }
    else
    {
      fprintf ( output_unit, "P5\n" );
      fprintf ( output_unit, "64  64\n" );
    }
    for ( i = 0; i < 64*64; i++ )
    {
      fread ( &value, sizeof ( unsigned short ), 1, input_unit );

      if ( ascii )
      {
        fprintf ( output_unit, "%u\n", value );
      }
      else
      {
        fwrite ( &value, sizeof ( unsigned short ), 1, output_unit );
      }
    }
    fclose ( output_unit );
    printf ( "%s\n", output_filename );
  }
  fclose ( input_unit );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MRI_TO_PGM:\n" );
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
