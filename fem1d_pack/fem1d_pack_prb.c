# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem1d_pack.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM1D_PACK_PRB.

  Discussion:

    FEM1D_PACK_PRB tests the FEM1D_PACK library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_PACK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM1D_PACK library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_PACK_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 verifies LOCAL_BASIS_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 4

  double a;
  double b;
  int i;
  int j;
  int node_num = NODE_NUM;
  double node_x[NODE_NUM] = { 1.0, 2.0, 4.0, 4.5 };
  double *phi;
  double phi_matrix[NODE_NUM*NODE_NUM];
  double s;
  int seed;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  LOCAL_BASIS_1D evaluates the local basis functions\n" );
  printf ( "  for a 1D element.\n" );
  printf ( "\n" );
  printf ( "  Test that the basis functions, evaluated at the nodes,\n" );
  printf ( "  form the identity matrix.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", node_num );

  printf ( "\n" );
  printf ( "  Node coordinates:\n" );
  printf ( "\n" );
  for ( j = 0; j < node_num; j++ )
  {
    printf ( "  %8d  %7g\n", j, node_x[j] );
  }

  for ( j = 0; j < node_num; j++ )
  {
    x = node_x[j];
    phi = local_basis_1d ( node_num, node_x, x );
    for ( i = 0; i < node_num; i++ )
    {
      phi_matrix[i+j*node_num] = phi[i];
    }
    free ( phi );
  }

  r8mat_print ( node_num, node_num, phi_matrix, "  A(I,J) = PHI(I) at node (J):" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  The PHI functions should sum to 1 at random X values:\n" );
  printf ( "\n" );
  printf ( "       X        Sum ( PHI(:)(X) )\n" );
  printf ( "\n" );

  a = 1.0;
  b = 4.5;
  for ( j = 1; j <= 5; j++ )
  {
    x = r8_uniform_ab ( a, b, &seed );
    phi = local_basis_1d ( node_num, node_x, x );
    s = r8vec_sum ( node_num, phi );
    printf ( "  %14g  %14g\n", x, s );
    free ( phi );
  }

  return;
# undef NODE_NUM
}
