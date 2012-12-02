# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "burgers_solution.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BURGERS_SOLUTION_PRB.

  Discussion:

    BURGERS_SOLUTION_PRB tests the BURGERS_SOLUTION library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BURGERS_SOLUTION_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BURGERS_SOLUTION library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BURGERS_SOLUTION_PRB\n" );
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

    TEST01 tests sets up a small test case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2011

  Author:

    John Burkardt
*/
{
  char *filename = "burgers_test01.txt";
  double nu;
  double pi = 3.141592653589793;
  double thi;
  double tlo;
  double *vu;
  double *vt;
  int vtn = 11;
  double *vx;
  int vxn = 11;
  double xhi;
  double xlo;

  nu = 0.01 / pi;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compute an analytic solution to the Burgers equation.\n" );
  printf ( "\n" );
  printf ( "  Viscosity NU = %g\n", nu );
  printf ( "  NX = %d\n", vxn );
  printf ( "  NT = %d\n", vtn );

  xlo = -1.0;
  xhi = +1.0;
  vx = r8vec_even_new ( vxn, xlo, xhi );
  r8vec_print ( vxn, vx, "  X grid points:" );

  tlo = 0.0;
  thi = 3.0 / pi;
  vt = r8vec_even_new ( vtn, tlo, thi );
  r8vec_print ( vtn, vt, "  T grid points:" );

  vu = burgers_solution ( nu, vxn, vx, vtn, vt );

  r8mat_print ( vxn, vtn, vu, "  U(X,T) at grid points:" );

  r8mat_write ( filename, vxn, vtn, vu );

  printf ( "\n" );
  printf ( "  Data written to file \"%s\".\n", filename );

  free ( vt );
  free ( vu );
  free ( vx );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests sets up a finer test case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2011

  Author:

    John Burkardt
*/
{
  char *filename = "burgers_test02.txt";
  double nu;
  double pi = 3.141592653589793;
  double thi;
  double tlo;
  double *vu;
  double *vt;
  int vtn = 41;
  double *vx;
  int vxn = 41;
  double xhi;
  double xlo;

  nu = 0.01 / pi;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Compute an analytic solution to the Burgers equation.\n" );
  printf ( "\n" );
  printf ( "  Viscosity NU = %g\n", nu );
  printf ( "  NX = %d\n", vxn );
  printf ( "  NT = %d\n", vtn );

  xlo = -1.0;
  xhi = +1.0;
  vx = r8vec_even_new ( vxn, xlo, xhi );
  r8vec_print ( vxn, vx, "  X grid points:" );

  tlo = 0.0;
  thi = 3.0 / pi;
  vt = r8vec_even_new ( vtn, tlo, thi );
  r8vec_print ( vtn, vt, "  T grid points:" );

  vu = burgers_solution ( nu, vxn, vx, vtn, vt );

  r8mat_write ( filename, vxn, vtn, vu );

  printf ( "\n" );
  printf ( "  Data written to file \"%s\".\n", filename );

  free ( vt );
  free ( vu );
  free ( vx );

  return;
  return;
}
