# include <stdlib.h>
# include <stdio.h>

# include "pbmlib.h"

int main ( long argc, char **argv );
int test01 ( char *file_name );
int test02 ( char *file_name );
int test03 ( char *file_name );
int test04 ( char *file_name );
int test05 ( char *file_name );
int test06 ( char *file_name );
int test07 ( char *file_name );
int test08 ( char *file_name );
int test09 ( char *file_name );
int test10 ( char *file_name );
int test11 ( char *file_name );
int test12 ( char *file_name );

/******************************************************************************/

int main ( long argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN calls the PBMLIB test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PBMLIB_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Tests for PBMLIB, portable bit map read/write routines.\n" );
  printf ( "\n" );

  test01 ( "pbmlib.pbma" );
  test02 ( "pbmlib.pbma" );

  test03 ( "pbmlib.pgma" );
  test04 ( "pbmlib.pgma" );

  test05 ( "pbmlib.ppma" );
  test06 ( "pbmlib.ppma" );

  test07 ( "pbmlib.pbmb" );
  test08 ( "pbmlib.pbmb" );

  test09 ( "pbmlib.pgmb" );
  test10 ( "pbmlib.pgmb" );

  test11 ( "pbmlib.ppmb" );
  test12 ( "pbmlib.ppmb" );

  printf ( "\n" );
  printf ( "PBMLIB_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int test01 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PBM_EXAMPLE, PBMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  int result;
  int xsize = 200;
  int ysize = 200;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PBM_EXAMPLE sets up PBM data.\n" );
  printf ( "  PBMA_WRITE writes an ASCII PBM file.\n" );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pbm_example ( xsize, ysize, b );

  printf ( "\n" );
  printf ( "  PBM_EXAMPLE has set up the data.\n" );

  result = pbma_write ( file_name, xsize, ysize, b );

  if ( result != 0 )
  {
    printf ( "  PBMA_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PBMA_WRITE was successful.\n" );
  }

  free ( b );

  return result;
}
/******************************************************************************/

int test02 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PBMA_READ_HEADER, PBMA_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  FILE *file_pointer;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PBMA_READ_HEADER reads the header of a PBMA file.\n" );
  printf ( "  PBMA_READ_DATA reads the data in a PBMA file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "r" );

  result = pbma_read_header ( file_pointer, &xsize, &ysize );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PBMA_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pbma_read_data ( file_pointer, xsize, ysize, b );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PBMA_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( b );

  return result;
}
/******************************************************************************/

int test03 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PGM_EXAMPLE, PGMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *g;
  int result;
  int xsize = 600;
  int ysize = 200;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  PGM_EXAMPLE sets up PGM data.\n" );
  printf ( "  PGMA_WRITE writes an ASCII PGM file.\n" );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pgm_example ( xsize, ysize, g );

  if ( result != 0 )
  {
    printf ( "  PGM_EXAMPLE failed!\n" );
    return result;
  }
  else
  {
    printf ( "  PGM_EXAMPLE has set up the data.\n" );
  }

  result = pgma_write ( file_name, xsize, ysize, g );

  if ( result != 0 )
  {
    printf ( "  PGMA_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PGMA_WRITE was successful.\n" );
  }

  free ( g );

  return result;
}
/******************************************************************************/

int test04 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests PGMA_READ_HEADER, PGMA_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  FILE *file_pointer;
  int *g;
  int maxgray;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  PGMA_READ_HEADER reads the header of a PGMA file.\n" );
  printf ( "  PGMA_READ_DATA reads the data in a PGMA file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "r" );

  result = pgma_read_header ( file_pointer, &xsize, &ysize, &maxgray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PGMA_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );
  printf ( "  Maximum data value =        %d\n", maxgray );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pgma_read_data ( file_pointer, xsize, ysize, g );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PGMA_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( g );

  return result;
}
/******************************************************************************/

int test05 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests PPM_EXAMPLE, PPMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  int *g;
  int *r;
  int result;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST05:\n" );
  printf ( "  PPM_EXAMPLE sets up PPM data.\n" );
  printf ( "  PPMA_WRITE writes an ASCII PPM file.\n" );

  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = ppm_example ( xsize, ysize, r, g, b );

  if ( result != 0 )
  {
    printf ( "  PPM_EXAMPLE failed!\n" );
    return result;
  }
  else
  {
    printf ( "  PPM_EXAMPLE has set up the data.\n" );
  }

  result = ppma_write ( file_name, xsize, ysize, r, g, b );

  if ( result != 0 )
  {
    printf ( "  PPMA_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PPMA_WRITE was successful.\n" );
  }

  free ( r );
  free ( g );
  free ( b );

  return result;
}
/******************************************************************************/

int test06 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests PPMA_READ_HEADER, PPMA_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  FILE *file_pointer;
  int *g;
  int maxrgb;
  int *r;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  PPMA_READ_HEADER reads the header of a PPMA file.\n" );
  printf ( "  PPMA_READ_DATA reads the data in a PPMA file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "r" );

  result = ppma_read_header ( file_pointer, &xsize, &ysize, &maxrgb );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PPMA_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );
  printf ( "  Maximum data value =        %d\n", maxrgb );

  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = ppma_read_data ( file_pointer, xsize, ysize, r, g, b );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PPMA_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( r );
  free ( g );
  free ( b );

  return result;
}
/******************************************************************************/

int test07 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests PBM_EXAMPLE, PBMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  int result;
  int xsize = 200;
  int ysize = 200;

  printf ( "\n" );
  printf ( "TEST07:\n" );
  printf ( "  PBM_EXAMPLE sets up PBM data.\n" );
  printf ( "  PBMB_WRITE writes a binary PBM file.\n" );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pbm_example ( xsize, ysize, b );

  printf ( "\n" );
  printf ( "  PBM_EXAMPLE has set up the data.\n" );

  result = pbmb_write ( file_name, xsize, ysize, b );

  if ( result != 0 )
  {
    printf ( "  PBMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PBMB_WRITE was successful.\n" );
  }

  free ( b );

  return result;
}
/******************************************************************************/

int test08 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests PBMB_READ_HEADER, PBMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  FILE *file_pointer;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  PBMB_READ_HEADER reads the header of a PBMB file.\n" );
  printf ( "  PBMB_READ_DATA reads the data in a PBMB file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "rb" );

  result = pbmb_read_header ( file_pointer, &xsize, &ysize );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PBMB_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pbmb_read_data ( file_pointer, xsize, ysize, b );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PBMB_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( b );

  return result;
}
/******************************************************************************/

int test09 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests PGM_EXAMPLE, PGMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *g;
  int result;
  int xsize = 600;
  int ysize = 200;

  printf ( "\n" );
  printf ( "TEST09:\n" );
  printf ( "  PGM_EXAMPLE sets up PGM data.\n" );
  printf ( "  PGMB_WRITE writes a binary PGM file.\n" );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pgm_example ( xsize, ysize, g );

  printf ( "\n" );
  printf ( "  PGM_EXAMPLE has set up the data.\n" );

  result = pgmb_write ( file_name, xsize, ysize, g );

  if ( result != 0 )
  {
    printf ( "  PGMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PGMB_WRITE was successful.\n" );
  }

  free ( g );

  return result;
}
/******************************************************************************/

int test10 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests PGMB_READ_HEADER, PBMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  FILE *file_pointer;
  int *g;
  int maxgray;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  PGMB_READ_HEADER reads the header of a PGMB file.\n" );
  printf ( "  PGMB_READ_DATA reads the data in a PGMB file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "rb" );

  result = pgmb_read_header ( file_pointer, &xsize, &ysize, &maxgray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PGMB_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );
  printf ( "  Maximum data value =        %d\n", maxgray );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = pgmb_read_data ( file_pointer, xsize, ysize, g );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PGMB_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( g );

  return result;
}
/******************************************************************************/

int test11 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests PPM_EXAMPLE, PPMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  int *g;
  int *r;
  int result;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST05:\n" );
  printf ( "  PPM_EXAMPLE sets up PPM data.\n" );
  printf ( "  PPMB_WRITE writes a binary PPM file.\n" );

  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = ppm_example ( xsize, ysize, r, g, b );

  printf ( "\n" );
  printf ( "  PPM_EXAMPLE has set up the data.\n" );

  result = ppmb_write ( file_name, xsize, ysize, r, g, b );

  if ( result != 0 )
  {
    printf ( "  PPMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PPMB_WRITE was successful.\n" );
  }

  free ( r );
  free ( g );
  free ( b );

  return result;
}
/******************************************************************************/

int test12 ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests PPMB_READ_HEADER, PPMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt
*/
{
  int *b;
  FILE *file_pointer;
  int *g;
  int maxrgb;
  int *r;
  int result;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  PPMB_READ_HEADER reads the header of a PPMB file.\n" );
  printf ( "  PPMB_READ_DATA reads the data in a PPMB file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_name );

  file_pointer = fopen ( file_name, "r" );

  result = ppmb_read_header ( file_pointer, &xsize, &ysize, &maxrgb );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PPMB_READ_HEADER failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );
  printf ( "  Maximum data value =        %d\n", maxrgb );

  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  result = ppmb_read_data ( file_pointer, xsize, ysize, r, g, b );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PPMB_READ_DATA failed!\n" );
    return result;
  }

  printf ( "\n" );
  printf ( "  The data was read successfully.\n" );

  free ( r );
  free ( g );
  free ( b );

  return result;
}
