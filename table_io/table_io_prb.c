# include <stdlib.h>
# include <stdio.h>

# include "table_io.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TABLE_IO_PRB.

  Discussion:

    TABLE_IO_PRB tests the TABLE_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TABLE_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TABLE_IO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TABLE_IO_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests R8MAT_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
# define M 5
# define N 20

  int i;
  int j;
  char *output_filename = "r8mat_05_00020.txt";
  double table[M*N];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  R8MAT_WRITE writes an R8MAT file.\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      table[i+j*M] = ( double ) ( 100 * ( j + 1 ) + ( i + 1 ) ) / 10.0;
    }
  }

  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", M );
  printf ( "  Number of points N  = %d\n", N );

  r8mat_print_some ( M, N, table, 1, 1, 5, 5, 
    "  5x5 portion of the data written to file:" );

  r8mat_transpose_print_some ( M, N, table, 1, 1, 5, 5, 
    "  5x5 portion of the TRANSPOSED data:" );

  r8mat_write ( output_filename, M, N, table );

  printf ( "\n" );
  printf ( "  Wrote the header and data for \"%s\".\n",
    output_filename );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8MAT_HEADER_READ and R8MAT_DATA_READ.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
  char *input_filename = "r8mat_05_00020.txt";
  int i;
  int j;
  int m;
  int n;
  double *table;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For an R8MAT file,\n" );
  printf ( "  R8MAT_HEADER_READ reads the header\n" );
  printf ( "  (Information about the dimension of the data)\n" );
  printf ( "  R8MAT_DATA_READ reads the data.\n" );

  r8mat_header_read ( input_filename, &m, &n );

  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", input_filename );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of points N  = %d\n", n );

  table = r8mat_data_read ( input_filename, m, n );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", input_filename );

  r8mat_print_some ( m, n, table, 1, 1, 5, 5, 
    "  5x5 portion of data read from file:" );

  free ( table );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests I4MAT_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
# define M 5
# define N 20

  int i;
  int j;
  char *output_filename = "i4mat_05_00020.txt";
  int table[M*N];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  I4MAT_WRITE writes an I4MAT file.\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      table[i+j*M] = ( 100 * ( j + 1 ) + ( i + 1 ) );
    }
  }

  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", M );
  printf ( "  Number of points N  = %d\n", N );

  i4mat_print_some ( M, N, table, 1, 1, 5, 5, 
    "  5 x 5 portion of data written to file:" );

  i4mat_write ( output_filename, M, N, table );

  printf ( "\n" );
  printf ( "  Wrote the header and data for \"%s\".\n",
    output_filename );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests I4MAT_HEADER_READ, I4MAT_DATA_READ.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
  char *input_filename = "i4mat_05_00020.txt";
  int m;
  int n;
  int *table;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For an I4MAT file,\n" );
  printf ( "  I4MAT_HEADER_READ reads the header\n" );
  printf ( "  (Information about the dimension of the data)\n" );
  printf ( "  I4MAT_DATA_READ reads the data.\n" );

  i4mat_header_read ( input_filename, &m, &n );

  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", input_filename );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of points N  = %d\n", n );

  table = i4mat_data_read (  input_filename, m, n );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", input_filename );

  i4mat_print_some ( m, n, table, 1, 1, 5, 5, 
    "  5x5 portion of data read from file:" );

  free ( table );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests R8MAT_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
# define M 2
# define N 10

  int seed = 123456789;
  double *table;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  R8MAT_UNIFORM_01 sets a random R8MAT.\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", M );
  printf ( "  Number of points N =  %d\n", N );

  table = r8mat_uniform_01 ( M, N, &seed );

  r8mat_print_some ( M, N, table, 1, 1, 5, 10, 
    "  5x10 portion of random real table dataset:" );

  free ( table );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests I4MAT_BORDER_CUT and I4MAT_BORDER_ADD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 August 2009

  Author:

    John Burkardt
*/
{
  int m = 6;
  int n = 4;
  int *table;
  int *table2;
  int *table3;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  I4MAT_BORDER_CUT cuts off the border;\n" );
  printf ( "  I4MAT_BORDER_ADD adds a zero border;\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of points N =  %d\n", n );

  table = i4mat_indicator_new ( m, n );

  i4mat_print ( m, n, table, "  Initial dataset:" );

  table2 = i4mat_border_cut ( m, n, table );

  i4mat_print ( m-2, n-2, table2, "  'Cut' dataset:" );

  table3 = i4mat_border_add ( m-2, n-2, table2 );

  i4mat_print ( m, n, table3, "  'Added' dataset:" );

  free ( table );
  free ( table2 );
  free ( table3 );

  return;
}

