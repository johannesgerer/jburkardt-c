# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( );
void filename_inc ( char *filename );
int *life_init ( double prob, int m, int n, int *seed );
void life_update ( int m, int n, int grid[] );
void life_write ( char *output_filename, int m, int n, int grid[] );
double r8_uniform_01 ( int *seed );
int s_len_trim ( char *s );
void timestamp ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LIFE_SERIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Reference:

    Martin Gardner,
    Mathematical Games:
    The Fantastic Combinations of John Conway's new solitaire game "Life",
    Scientific American,
    Volume 223, Number 4, October 1970, pages 120-123.
*/
{
  char filename[] = "life_000.txt";
  int it;
  int it_max;
  int m;
  int n;
  int *grid;
  double prob;
  int seed;

  timestamp ( );
  printf ( "\n" );
  printf ( "LIFE_SERIAL\n" );
  printf ( "  C version\n" );
  printf ( "  Carry out a few steps of John Conway's\n" );
  printf ( "  Game of Life.\n" );
  printf ( "\n" );

  it_max = 10;
  m = 10;
  n = 10;
  prob = 0.20;
  seed = 123456789;

  for ( it = 0; it <= it_max; it++ )
  {
    if ( it == 0 )
    {
      grid = life_init ( prob, m, n, &seed );
    }
    else
    {
      life_update ( m, n, grid );
    }
    life_write ( filename, m, n, grid );
    printf ( "  %s\n", filename );
    filename_inc ( filename );
  }
/*
  Free memory.
*/
  free ( grid );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LIFE_SERIAL\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void filename_inc ( char *filename )

/******************************************************************************/
/*
  Purpose:

    FILENAME_INC increments a partially numeric file name.

  Discussion:

    It is assumed that the digits in the name, whether scattered or
    connected, represent a number that is to be increased by 1 on
    each call.  If this number is all 9's on input, the output number
    is all 0's.  Non-numeric letters of the name are unaffected.

    If the name is empty, then the routine stops.

    If the name contains no digits, the empty string is returned.

  Example:

      Input            Output
      -----            ------
      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
      "a9to99.txt"     "a0to00.txt"  (wrap around)
      "cat.txt"        " "           (no digits to increment)
      " "              STOP!         (error)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 November 2011

  Author:

    John Burkardt

  Parameters:

    Input/output, char *FILENAME, the filename to be incremented.
*/
{
  char c;
  int change;
  int i;
  int n;
  char *t;

  n = s_len_trim ( filename );

  if ( n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILENAME_INC - Fatal error!\n" );
    fprintf ( stderr, "  The input string is empty.\n" );
    exit ( 1 );
  }

  change = 0;

  t = filename + n - 1;
  
  while ( 0 < n )
  {
    if ( '0' <= *t && *t <= '9' )
    {
      change = change + 1;

      if ( *t == '9' )
      {
        *t = '0';
      }
      else
      {
        *t = *t + 1;
        return;
      }
    }
    t--;
    n--;
  }
/*
  No digits were found.  Return blank.
*/
  if ( change == 0 )
  {
    n = s_len_trim ( filename );
    t = filename + n - 1;
    while ( 0 < n )
    {
      *t = ' ';
      t--;
      n--;
    }
  }

  return;
}
/******************************************************************************/

int *life_init ( double prob, int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    LIFE_INIT initializes the life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, double PROB, the probability that a grid cell
    should be alive.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, int LIFE_INIT[(1+M+1)*(1+N+1)], the initial grid.
*/
{
  int *grid;
  int i;
  int j;
  double r;

  grid = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );
  for ( j = 0; j <= n + 1; j++ )
  {
    for ( i = 0; i <= m + 1; i++ )
    {
      grid[i+j*(m+2)] = 0;
    }
  }

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      r = r8_uniform_01 ( seed );
      if ( r <= prob )
      {
        grid[i+j*(m+2)] = 1;
      }
    }
  }

  return grid;
}
/******************************************************************************/

void life_update ( int m, int n, int grid[] )

/******************************************************************************/
/*
  Purpose:

    LIFE_UPDATE updates a Life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
  int i;
  int j;
  int *s;

  s = ( int * ) malloc ( m * n * sizeof ( int ) );

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      s[i-1+(j-1)*m] = 
          grid[i-1+(j-1)*(m+2)] + grid[i-1+j*(m+2)] + grid[i-1+(j+1)*(m+2)]
        + grid[i  +(j-1)*(m+2)]                     + grid[i  +(j+1)*(m+2)]
        + grid[i+1+(j-1)*(m+2)] + grid[i+1+j*(m+2)] + grid[i+1+(j+1)*(m+2)];
    }
  }
/*
  Any dead cell with 3 live neighbors becomes alive.
  Any living cell with less than 2 or more than 3 neighbors dies.
*/
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      if ( grid[i+j*(m+2)] == 0 )
      {
        if ( s[i-1+(j-1)*m] == 3 )
        {
          grid[i+j*(m+2)] = 1;
        }
      }
      else if ( grid[i+j*(m+2)] == 1 )
      {
        if ( s[i-1+(j-1)*m] < 2 || 3 < s[i-1+(j-1)*m] )
        {
          grid[i+j*(m+2)] = 0;
        }
      }
    }
  }

  free ( s );

  return;
}
/******************************************************************************/

void life_write ( char *output_filename, int m, int n, int grid[] )

/******************************************************************************/
/*
  Purpose:

    LIFE_WRITE writes a grid to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output file name.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
  int i;
  int j;
  FILE *output_unit;
/*
  Open the file.
*/
  output_unit = fopen ( output_filename, "wt" );
/*
  Write the data.
*/
  for ( j = 0; j <= n + 1; j++ )
  {
    for ( i = 0; i <= m + 1; i++ )
    {
      fprintf ( output_unit, " %d", grid[i+j*(m+2)] );
    }
    fprintf ( output_unit, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output_unit );

  return;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a pseudorandom R8 scaled to [0,1].

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int i4_huge = 2147483647;
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
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

