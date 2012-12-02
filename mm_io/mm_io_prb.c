# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "mm_io.h"

int main ( void );
void test01 ( char *input_filename );
void test02 ( char *input_filename );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MM_IO_PRB calls the MM_IO tests.

  Modified:

    31 October 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MM_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MM_IO library.\n" );

  test01 ( "matrix_05_05_crg.txt" );
  test02 ( "matrix_05_05_crg.txt" );
  test02 ( "matrix_05_05_arg.txt" );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MM_IO_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests MM_READ_BANNER.

  Modified:

    02 November 2008

  Author:

    John Burkardt
*/
{
  FILE *input_file;
  MM_typecode matcode;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  MM_READ_BANNER reads the header line of\n" );
  printf ( "    a Matrix Market file.\n" );
  printf ( "\n" );
  printf ( "  Reading \"%s\".\n", input_filename );

  input_file = fopen ( input_filename, "r" );

  mm_read_banner ( input_file, &matcode );

  printf ( "\n" );
  printf ( "  MM_typecode[0] = %c\n", matcode[0] );
  printf ( "  MM_typecode[1] = %c\n", matcode[1] );
  printf ( "  MM_typecode[2] = %c\n", matcode[2] );
  printf ( "  MM_typecode[3] = %c\n", matcode[3] );

  return;
}
/******************************************************************************/

void test02 ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests MM_READ_MTX_ARRAY_SIZE and MM_READ_MTX_CRD_SIZE.

  Modified:

    03 November 2008

  Author:

    John Burkardt
*/
{
  int error;
  FILE *input_file;
  int m;
  MM_typecode matcode;
  int n;
  int nz;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  MM_READ_MTX_ARRAY_SIZE reads the sizes line of\n" );
  printf ( "    an MM matrix file;\n" );
  printf ( "  MM_READ_MTX_CRD_SIZE reads the sizes line of\n" );
  printf ( "    an MM coordinate file;\n" );
  printf ( "\n" );
  printf ( "  Reading \"%s\".\n", input_filename );

  input_file = fopen ( input_filename, "r" );

  mm_read_banner ( input_file, &matcode );

  if ( mm_is_array ( matcode ) )
  {
    error = mm_read_mtx_array_size ( input_file, &m, &n );

    printf ( "\n" );
    printf ( "  Array sizes:\n" );
    printf ( "    M = %d\n", m );
    printf ( "    N = %d\n", n );
  }
  else if ( mm_is_coordinate ( matcode ) )
  {
    error = mm_read_mtx_crd_size( input_file, &m, &n, &nz );

    printf ( "\n" );
    printf ( "  Coordinate sizes:\n" );
    printf ( "    M  = %d\n", m );
    printf ( "    N  = %d\n", n );
    printf ( "    NZ = %d\n", nz );
  }
  else
  {
    printf ( "\n" );
    printf ( "Warning:\n" );
    printf ( "  File does not seem to be of MATRIX or COORDINATE type.\n" );
  }

  return;
}
