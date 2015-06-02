# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <time.h>

# include "backtrack_binary_rc.h"

/******************************************************************************/

void backbin_rc ( int n, bool reject, int *n2, int choice[] )

/******************************************************************************/
/*
  Purpose:

    BACKBIN_RC uses reverse communication for binary backtracking.

  Discussion:

    If this procedure returns a solution with N2 = N, which is acceptable
    to the user, then a full solution has been found.

    If this procedure returns N2 = -1, no more potential solutions are
    available to consider.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 January 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the length of the full solution.

    Input, bool REJECT, is TRUE if the proposed partial solution
    in the first N2 entries of CHOICE must be rejected.

    Input/output, int *N2, the length of the current
    partial solution.  On first call for a given problem, the user
    should set N2 to -1.  If the program has exhausted the search space,
    the value of N2 will be returned as -1.

    Input/output, int CHOICE[N], indicates the current
    partial solution in entries 1 through N2, which will contain 0 or 1.
*/
{
  int i;
/*
  N2 = -1 means an initialization call.
*/
  if ( *n2 == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      choice[i] = -1;
    }
    *n2 = 1;
    choice[*n2-1] = 1;
  }
/*
  1 <= FOCUS means we asked the user to evaluate CHOICE(1:N2).

  N2 = N means we returned a full prospective solution
  so in any case we must increment CHOICE.

  Returning REJECT = 1 means no solution begins this way
  so we must increment CHOICE.
*/
  else if ( *n2 == n || reject )
  {
    while ( 1 < *n2 )
    {
      if ( choice[*n2-1] == 1 )
      {
        choice[*n2-1] = 0;
        break;
      }
      choice[*n2-1] = -1;
      *n2 = *n2 - 1;
    }
/*
  Have we exhausted the solution space?
*/
    if ( *n2 == 1 )
    {
      if ( choice[*n2-1] == 1 )
      {
        choice[*n2-1] = 0;
      }
      else
      {
        choice[*n2-1] = -1;
        *n2 = -1;
      }
    }
  }
/*
  N2 < N and not REJECT means we can increment N2.
*/
  else
  {
    *n2 = *n2 + 1;
    choice[*n2-1] = 1;
  }

  return;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
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
