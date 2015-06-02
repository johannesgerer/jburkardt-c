# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "stochastic_rk.h"

int main ( void );
void test01 ( void );
double fi ( double x );
double gi ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for STOCHASTIC_RK_PRB.

  Discussion:

    STOCHASTIC_RK_PRB tests the STOCHASTIC_RK library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "STOCHASTIC_RK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the STOCHASTIC_RK library.\n" );
 
  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "STOCHASTIC_RK_PRB\n" );
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

    TEST01 tests RK1_TI_STEP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2010

  Author:

    John Burkardt
*/
{
  double h;
  int i;
  int n;
  double q;
  int seed;
  double t;
  double t0 = 0.0;
  double tn = 1.0;
  double *x;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST01\n" );
  fprintf ( stdout, "  RK1_TI_STEP uses a first order RK method\n" );
  fprintf ( stdout, "  for a problem whose right hand side does not\n" );
  fprintf ( stdout, "  depend explicitly on time.\n" );

  n = 10;
  x = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  h = ( tn - t0 ) / ( double ) ( n );
  q = 1.0;
  seed = 123456789;

  i = 0;
  t = t0;
  x[i] = 0.0;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "         I           T             X\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  %8d  %14f  %14f\n", i, t, x[i] );

  for ( i = 1; i <= n; i++ )
  {
    t = ( ( double ) ( n - i ) * t0   
        + ( double ) (     i ) * tn )
        / ( double ) ( n     );

    x[i] = rk1_ti_step ( x[i-1], t, h, q, fi, gi, &seed );

    fprintf ( stdout, "  %8d  %14f  %14f\n", i, t, x[i] );
  }

  free ( x );
  return;
}
/******************************************************************************/

double fi ( double x )

/******************************************************************************/
/*
  Purpose:

    FI is a time invariant deterministic right hand side.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double FI, the value.
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double gi ( double x )

/******************************************************************************/
/*
  Purpose:

    GI is a time invariant stochastic right hand side.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double GI, the value.
*/
{
  double value;

  value = 1.0;

  return value;
}
