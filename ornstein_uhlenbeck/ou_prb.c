# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "ou.h"

int main ( );
void ou_euler_test ( );
void ou_euler_maruyama_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for OU_PRB.

  Discussion:

    OU_PRB tests the OU library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "OU_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the OU library.\n" );

  ou_euler_test ( );
  ou_euler_maruyama_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "OU_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ou_euler_test ( )

/******************************************************************************/
/*
  Purpose:

    OU_EULER_TEST tests OU_EULER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2013

  Author:

    John Burkardt
*/
{
  double mu;
  int n;
  int seed;
  double sigma;
  double theta;
  double tmax;
  double x0;

  printf ( "\n" );
  printf ( "OU_EULER_TEST:\n" );
  printf ( "  Estimate a solution to the Ornstein-Uhlenbeck equation\n" );
  printf ( "  using the Euler method for stochastic differential equations.\n" );
  printf ( "\n" );

  theta = 2.0;
  printf ( "  Using decay rate THETA = %g\n", theta );
  mu = 1.0;
  printf ( "  Using mean MU = %g\n", mu );
  sigma = 0.15;
  printf ( "  Using variance SIGMA = %g\n", sigma );
  x0 = 2.0;
  printf ( "  Using initial value X0 = %g\n", x0 );
  tmax = 3.0;
  printf ( "  Using final time TMAX = %g\n", tmax );
  n = 10000;
  printf ( "  Using number of timesteps N = %d\n", n );
  seed = 123456789;
  printf ( "  Using value of random SEED = %d\n", seed );

  ou_euler ( theta, mu, sigma, x0, tmax, n, &seed );

  return;
}
/******************************************************************************/

void ou_euler_maruyama_test ( )

/******************************************************************************/
/*
  Purpose:

    OU_EULER_MARUYAMA_TEST tests OU_EULER_MARUYAMA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2013

  Author:

    John Burkardt
*/
{
  double mu;
  int n;
  int r;
  int seed;
  double sigma;
  double theta;
  double tmax;
  double x0;

  printf ( "\n" );
  printf ( "OU_EULER_MARUYAMA_TEST:\n" );
  printf ( "  Estimate a solution to the Ornstein-Uhlenbeck equation\n" );
  printf ( "  using the Euler-Maruyama method for stochastic \n" );
  printf ( "  differential equations.\n" );
  printf ( "\n" );

  theta = 2.0;
  printf ( "  Using decay rate THETA = %g\n", theta );
  mu = 1.0;
  printf ( "  Using mean MU = %g\n", mu );
  sigma = 0.15;
  printf ( "  Using variance SIGMA = %g\n", sigma );
  x0 = 2.0;
  printf ( "  Using initial value X0 = %g\n", x0 );
  tmax = 3.0;
  printf ( "  Using final time TMAX = %g\n", tmax );
  n = 10000;
  printf ( "  Using number of large timesteps N = %d\n", n );
  r = 16;
  printf ( "  Using number small time steps per one large time step R = %d\n", r );
  seed = 123456789;
  printf ( "  Using value of random SEED = %d\n", seed );

  ou_euler_maruyama ( theta, mu, sigma, x0, tmax, n, r, &seed );

  return;
}

