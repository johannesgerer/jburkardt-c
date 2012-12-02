# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
void r8vla2_write ( char *output_filename, int m, int n, double a[m][n] );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    FD_PREDATOR_PREY solves a pair of predator-prey ODE's.

  Discussion:

    The physical system under consideration is a pair of animal populations.

    The PREY reproduce rapidly; for each animal alive at the beginning of the
    year, two more will be born by the end of the year.  The prey do not have
    a natural death rate; instead, they only die by being eaten by the predator.
    Every prey animal has 1 chance in 1000 of being eaten in a given year by
    a given predator.

    The PREDATORS only die of starvation, but this happens very quickly.
    If unfed, a predator will tend to starve in about 1/10 of a year.
    On the other hand, the predator reproduction rate is dependent on
    eating prey, and the chances of this depend on the number of available prey.

    The resulting differential equations can be written:

      PREY(0) = 5000
      PRED(0) =  100

      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)

    Here, the initial values (5000,100) are a somewhat arbitrary starting point.

    The pair of ordinary differential equations that result have an interesting
    behavior.  For certain choices of the interaction coefficients (such as
    those given here), the populations of predator and prey will tend to
    a periodic oscillation.  The two populations will be out of phase; the number
    of prey will rise, then after a delay, the predators will rise as the prey
    begins to fall, causing the predator population to crash again.

    In this program, the pair of ODE's is solved with a simple finite difference
    approximation using a fixed step size.  In general, this is NOT an efficient
    or reliable way of solving differential equations.  However, this program is
    intended to illustrate the ideas of finite difference approximation.

    In particular, if we choose a fixed time step size DT, then a derivative
    such as dPREY/dT is approximated by:

      d PREY / dT = approximately ( PREY(T+DT) - PREY(T) ) / DT

    which means that the first differential equation can be written as

      PREY(T+DT) = PREY(T) + DT * ( 2 * PREY(T) - 0.001 * PREY(T) * PRED(T) ).

    We can rewrite the equation for PRED as well.  Then, since we know the
    values of PRED and PREY at time 0, we can use these finite difference
    equations to estimate the values of PRED and PREY at time DT.  These values
    can be used to get estimates at time 2*DT, and so on.  To get from time
    T_START = 0 to time T_STOP = 5, we simply take STEP_NUM steps each of size
    DT = ( T_STOP - T_START ) / STEP_NUM.

    Because finite differences are only an approximation to derivatives, this
    process only produces estimates of the solution.  And these estimates tend
    to become more inaccurate for large values of time.  Usually, we can reduce
    this error by decreasing DT and taking more, smaller time steps.

    In this example, for instance, taking just 100 steps gives nonsensical
    answers.  Using STEP_NUM = 1000 gives an approximate solution that seems
    to have the right kind of oscillatory behavior, except that the amplitude
    of the waves increases with each repetition.  Using 10000 steps, the
    approximation begins to become accurate enough that we can see that the
    waves seem to have a fixed period and amplitude.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 February 2011

  Author:

    John Burkardt

  Reference:

    George Lindfield, John Penny,
    Numerical Methods Using MATLAB,
    Second Edition,
    Prentice Hall, 1999,
    ISBN: 0-13-012641-1,
    LC: QA297.P45.

  Parameters:

    Input, int STEP_NUM, the number of steps.
*/
{
# define LINE_MAX_LEN 80

  double dt;
  char filename[80];
  int i;
  int ierror;
  char input[LINE_MAX_LEN];
  int length;
  double pred_init = 100.0;
  double prey_init = 5000.0;
  char *s;
  int step_num;
  double t_start;
  double t_stop;

  timestamp ( );

  printf ( "\n" );
  printf ( "FD_PREDATOR_PREY\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  A finite difference approximate solution of a pair\n" );
  printf ( "  of ordinary differential equations for a population\n" );
  printf ( "  of predators and prey.\n" );
  printf ( "\n" );
  printf ( "  The exact solution shows wave behavior, with a fixed\n" );
  printf ( "  period and amplitude.  The finite difference approximation\n" );
  printf ( "  can provide a good estimate for this behavior if the stepsize\n" );
  printf ( "  DT is small enough.\n" );
/*
  STEP_NUM is an input argument or else read from the user interactively.
*/
  if ( 1 < argc )
  {
    step_num = atoi ( argv[1] );
/*
  s_to_i4 ( argv[1], &length, &ierror );
*/
  }
  else
  {
    printf ( "\n" );
    printf ( "FD_PREDATOR_PREY:\n" );
    printf ( "  Please enter the number of time steps:\n" );
    scanf ( "%d", &step_num );
  }

  t_start = 0.0;
  t_stop =  5.0;
  dt = ( t_stop - t_start ) / ( double ) ( step_num );

  printf ( "\n" );
  printf ( "  Initial time    = %f\n", t_start );
  printf ( "  Final time      = %f\n", t_stop );
  printf ( "  Initial prey    = %f\n", prey_init );
  printf ( "  Initial pred    = %f\n", pred_init );
  printf ( "  Number of steps = %d\n", step_num );
/*
  Declare TRF as a VLA (Variable Length Array).
*/
  double trf[3][step_num+1];

  trf[0][0] = t_start;
  trf[1][0] = prey_init;
  trf[2][0] = pred_init;

  for ( i = 0; i < step_num; i++ )
  {
    trf[0][i+1] = trf[0][i] + dt;

    trf[1][i+1] = trf[1][i] + dt
      * (    2.0 * trf[1][i] - 0.001 * trf[1][i] * trf[2][i] );

    trf[2][i+1] = trf[2][i] + dt
      * ( - 10.0 * trf[2][i] + 0.002 * trf[1][i] * trf[2][i] );
  }
/*
  Write data to files.
*/
  sprintf ( filename, "trf_%d.txt", step_num );

  r8vla2_write ( filename, 3, step_num + 1, trf );

  printf ( "  T, R, F values written to \"%s\".\n", filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FD_PREDATOR_PREY\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void r8vla2_write ( char *output_filename, int m, int n, double a[m][n] )

/******************************************************************************/
/*
  Purpose:

    R8VLA2_WRITE writes an R8VLA2 file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M][N], the table data.
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
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16e", a[i][j] );
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

int s_to_i4 ( char *s, int *last, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4 reads an I4 from a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a string to be examined.

    Output, int *LAST, the last character of S used to make IVAL.

    Output, int *ERROR is TRUE (1) if an error occurred and FALSE (0) otherwise.

    Output, int *S_TO_I4, the integer value read from the string.
    If the string is blank, then IVAL will be returned 0.
*/
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = 0;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s )
  {
    c = s[i];
    i = i + 1;
/*
  Haven't read anything.
*/
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read the sign, expecting digits.
*/
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read at least one digit, expecting more.
*/
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }
    }
  }
/*
  If we read all the characters in the string, see if we're OK.
*/
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = 1;
    *last = 0;
  }

  return ival;
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
