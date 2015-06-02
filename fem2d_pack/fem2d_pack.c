# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "fem2d_pack.h"

/******************************************************************************/

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.

  Discussion:

    The quantity computed here is the "geometric" bandwidth determined
    by the finite element mesh alone.

    If a single finite element variable is associated with each node
    of the mesh, and if the nodes and variables are numbered in the
    same way, then the geometric bandwidth is the same as the bandwidth
    of a typical finite element matrix.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where 

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    Because the finite element node adjacency relationship is symmetric,
    we are guaranteed that ML = MU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
/******************************************************************************/

void bandwidth_var ( int element_order, int element_num, int element_node[],
  int node_num, int var_node[], int var_num, int var[], int *ml, int *mu, 
  int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_VAR determines the bandwidth for finite element variables.

  Discussion:

    We assume that, attached to each node in the finite element mesh
    there are a (possibly zero) number of finite element variables.
    We wish to determine the bandwidth necessary to store the stiffness
    matrix associated with these variables.

    An entry K(I,J) of the stiffness matrix must be zero unless the
    variables I and J correspond to nodes N(I) and N(J) which are
    common to some element.

    In order to determine the bandwidth of the stiffness matrix, we
    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
    is maximized.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    We assume the finite element variable adjacency relationship is 
    symmetric, so we are guaranteed that ML = MU.

    Note that the user is free to number the variables in any way
    whatsoever, and to associate variables to nodes in any way,
    so that some nodes have no variables, some have one, and some
    have several.  

    The storage of the indices of the variables is fairly simple.
    In VAR, simply list all the variables associated with node 1, 
    then all those associated with node 2, and so on.  Then set up
    the pointer array VAR_NODE so that we can jump to the section of
    VAR where the list begins for any particular node.

    The routine does not check that each variable is only associated
    with a single node.  This would normally be the case in a finite
    element setting.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int node_global_i;
  int node_global_j;
  int node_local_i;
  int node_local_j;
  int var_global_i;
  int var_global_j;
  int var_local_i;
  int var_local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( node_local_i = 0; node_local_i < element_order; node_local_i++ )
    {
      node_global_i = element_node[node_local_i+element*element_order];

      for ( var_local_i = var_node[node_global_i-1]; 
            var_local_i <= var_node[node_global_i]-1; var_local_i++ )
      {
        var_global_i = var[var_local_i-1];

        for ( node_local_j = 0; node_local_j < element_order; node_local_j++ )
        {
          node_global_j = element_node[node_local_j+element*element_order];

          for ( var_local_j = var_node[node_global_j-1]; 
                var_local_j <= var_node[node_global_j]-1; var_local_j++ )
          {
            var_global_j = var[var_local_j-1];

            *mu = i4_max ( *mu, var_global_j - var_global_i );
            *ml = i4_max ( *ml, var_global_i - var_global_j );
          }
        }
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
/******************************************************************************/

void basis_11_t3 ( double t[2*3], int i, double p[2], double *qi, 
  double *dqidx, double *dqidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T3: one basis at one point for a T3 element.

  Discussion:

    The routine is given the coordinates of the nodes of a triangle.

           3
          / \
         /   \
        /     \
       1-------2

    It evaluates the linear basis function Q(I)(X,Y) associated with
    node I, which has the property that it is a linear function
    which is 1 at node I and zero at the other two nodes.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the coordinates of the nodes.

    Input, int I, the index of the desired basis function.
    I should be between 1 and 3.

    Input, double P[2], the coordinates of the point where 
    the basis function is to be evaluated.

    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
    and its X and Y derivatives.
*/
{
  double area;
  int ip1;
  int ip2;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Element has zero area.\n" );
    fprintf ( stderr, "  Area = %g\n", area );
    exit ( 1 );
  }

  if ( i < 1 || 3 < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Basis index I is not between 1 and 3.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }

  ip1 = i4_wrap ( i + 1, 1, 3 );
  ip2 = i4_wrap ( i + 2, 1, 3 );

  *qi = ( ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) 
        * ( p[1]           - t[1+(ip1-1)*2] ) 
        - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) 
        * ( p[0]           - t[0+(ip1-1)*2] ) ) / area;

  *dqidx = - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) / area;
  *dqidy =   ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) / area;

  return;
}
/******************************************************************************/

void basis_11_t3_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T3_TEST verifies BASIS_11_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 3

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  printf ( "\n" );
  printf ( "BASIS_11_T3_TEST:\n" );
  printf ( "  Verify basis functions for element T3.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      printf ( "  %10g", qj );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_11_t4 ( double t[2*4], int i, double p[], double *phi, 
  double *dphidx, double *dphidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4: one basis at one point for a T4 element.

  Discussion:

    The T4 element is the cubic bubble triangle.

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T4: 
       
            3
           / \
          /   \
         /  4  \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 August 2009

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*4], the coordinates of the vertices
    of the triangle, and the coordinates of the centroid.  
    It is common to list the first three points in counter clockwise
    order.

    Input, int I, the index of the basis function.

    Input, double P[2], the point where the basis function
    is to be evaluated.

    Output, double *PHI, the value of the basis function
    at the evaluation point.

    Output, double *DPHIDX, *DPHIDY, the value of the 
    derivatives at the evaluation point.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  double dpsidx[4];
  double dpsidy[4];
  int j;
  int n = 4;
  double psi[4];

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  psi[0] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1] - t[1+1*2] )     
                 - ( t[1+2*2] - t[1+1*2] ) * ( p[0] - t[0+1*2] ) );
  dpsidx[0] =    - ( t[1+2*2] - t[1+1*2] );
  dpsidy[0] =      ( t[0+2*2] - t[0+1*2] );

  psi[1] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1] - t[1+2*2] )     
                 - ( t[1+0*2] - t[1+2*2] ) * ( p[0] - t[0+2*2] ) );
  dpsidx[1] =    - ( t[1+0*2] - t[1+2*2] );
  dpsidy[1] =      ( t[0+0*2] - t[0+2*2] );

  psi[2] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1] - t[1+0*2] )     
                 - ( t[1+1*2] - t[1+0*2] ) * ( p[0] - t[0+0*2] ) );
  dpsidx[2] =    - ( t[1+1*2] - t[1+0*2] );
  dpsidy[2] =      ( t[0+1*2] - t[0+0*2] );
/*
  Normalize the first three functions.
*/
  for ( j = 0; j < 3; j++ )
  {
    psi[j]    = psi[j]    / area;
    dpsidx[j] = dpsidx[j] / area;
    dpsidy[j] = dpsidy[j] / area;
  }
/*
  Compute the cubic bubble function.
*/
  psi[3] = 27.0 * psi[0] * psi[1] * psi[2];

  dpsidx[3] = 27.0 * (
                  dpsidx[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidx[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidx[2] );

  dpsidy[3] = 27.0 * (
                  dpsidy[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidy[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidy[2] );
/*
  Subtract 1/3 of the cubic bubble function from each of the three linears.
*/
  for ( j = 0; j < 3; j++ )
  {
       psi[j] =    psi[j] -    psi[3] / 3.0;
    dpsidx[j] = dpsidx[j] - dpsidx[3] / 3.0;
    dpsidy[j] = dpsidy[j] - dpsidy[3] / 3.0;
  }

  *phi = psi[i-1];
  *dphidx = dpsidx[i-1];
  *dphidy = dpsidy[i-1];

  return;
}
/******************************************************************************/

void basis_11_t4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T4_TEST verifies BASIS_11_T4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2009

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 4

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    0.0, 0.0 };
/*
  The node associated with the fourth basis function is the centroid.
*/
  t[0+3*2] = ( t[0+0*2] + t[0+1*2] + t[0+2*2] ) / 3.0;
  t[1+3*2] = ( t[1+0*2] + t[1+1*2] + t[1+2*2] ) / 3.0;

  printf ( "\n" );
  printf ( "BASIS_11_T4_TEST:\n" );
  printf ( "  Verify basis functions for element T4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      printf ( "  %10g", qj );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_11_t6 ( double t[2*6], int i, double p[], double *bi, 
  double *dbidx, double *dbidy )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T6: one basis at one point for the T6 element.

  Discussion:

    The routine is given the coordinates of the nodes of a triangle. 
       
           3
          / \
         6   5
        /     \
       1---4---2

    It evaluates the quadratic basis function B(I)(X,Y) associated with
    node I, which has the property that it is a quadratic function
    which is 1 at node I and zero at the other five nodes.

    This routine assumes that the sides of the triangle are straight,
    so that the midside nodes fall on the line between two vertices.

    This routine relies on the fact that each basis function can be
    written as the product of two linear factors, which are easily
    computed and normalized.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*6], the coordinates of the nodes.

    Input, int I, the index of the desired basis function.
    I should be between 1 and 6.

    Input, double P[2], the coordinates of a point at which the basis
    function is to be evaluated.

    Output, double *BI, *DBIDX, *DBIDY, the values of the basis function
    and its X and Y derivatives.
*/
{
  double gf;
  double gn;
  double hf;
  double hn;
  int j1;
  int j2;
  int k1;
  int k2;

  if ( i < 1 || 6 < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_11_T6 - Fatal error!\n" );
    fprintf ( stderr, "  Basis index I is not between 1 and 6.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }
/*
  Determine the pairs of nodes.
*/
  if ( i <= 3 )
  {
    j1 = i4_wrap ( i + 1, 1, 3 );
    j2 = i4_wrap ( i + 2, 1, 3 );
    k1 = i + 3;
    k2 = i4_wrap ( i + 5, 4, 6 );
  }
  else
  {
    j1 = i - 3;
    j2 = i4_wrap ( i - 3 + 2, 1, 3 );
    k1 = i4_wrap ( i - 3 + 1, 1, 3 );
    k2 = i4_wrap ( i - 3 + 2, 1, 3 );
  }
/*
  For C indexing, it is helpful to knock the indices down by one.
*/
  i  = i  - 1;
  j1 = j1 - 1;
  j2 = j2 - 1;
  k1 = k1 - 1;
  k2 = k2 - 1;
/*
  Evaluate the two linear factors GF and HF, 
  and their normalizers GN and HN.
*/
  gf = ( p[0]      - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( p[1]      - t[1+j1*2] );

  gn = ( t[0+i*2]  - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( t[1+i*2]  - t[1+j1*2] );

  hf = ( p[0]      - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( p[1]      - t[1+k1*2] );

  hn = ( t[0+i*2]  - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( t[1+i*2]  - t[1+k1*2] );
/*
  Construct the basis function and its derivatives.
*/
  *bi =        ( gf                      / gn ) 
             * ( hf                      / hn );

  *dbidx =   ( ( t[1+j2*2] - t[1+j1*2] ) / gn ) 
           * (   hf                      / hn )
           + (   gf                      / gn ) 
           * ( ( t[1+k2*2] - t[1+k1*2] ) / hn );

  *dbidy = - ( ( t[0+j2*2] - t[0+j1*2] ) / gn ) 
           * (   hf                      / hn )
           - (   gf                      / gn ) 
           * ( ( t[0+k2*2] - t[0+k1*2] ) / hn );

  return;
}
/******************************************************************************/

void basis_11_t6_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_11_T6_TEST verifies BASIS_11_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double p[2];
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };
  double v1;
  double v2;
  double v3;

  printf ( "\n" );
  printf ( "BASIS_11_T6_TEST:\n" );
  printf ( "  Verify basis functions for element T6.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );
  for ( i = 1; i <= NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      p[0] = t[0+j*2];
      p[1] = t[1+j*2];

      basis_11_t6 ( t, i, p, &v1, &v2, &v3 );

      phi[i-1+j*NODE_NUM] = v1;
      dphidx[i-1+j*NODE_NUM] = v2;
      dphidy[i-1+j*NODE_NUM] = v3;
    }
  }
  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      printf ( "  %7g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_q4 ( double q[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_Q4: all bases at N points for a Q4 element.

  Discussion:

    The routine is given the coordinates of the vertices of a quadrilateral.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the element are presumed to lie along coordinate axes.

    The routine evaluates the basis functions associated with each corner,
    and their derivatives with respect to X and Y.

  Physical Element Q4:

    |
    |  4-----3
    |  |     |
    |  |     |
    Y  |     |
    |  |     |
    |  |     |
    |  1-----2
    |
    +-----X------>

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double Q[2*4], the coordinates of the vertices.
    It is common to list these points in counter clockwise order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the evaluation points.

    Output, double PHI[4*N], the bases at the evaluation points.

    Output, double DPHIDX[4*N], DPHIDY[4*N], the derivatives of the
    bases at the evaluation points.
*/
{
  double area;
  int i;
  int j;

  area =        ( q[0+2*2]         - q[0+0*2] ) 
              * ( q[1+2*2]         - q[1+0*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =      ( q[0+2*2] - p[0+j*2]             ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[1+j*4] =      (            p[0+j*2]  - q[0+0*2] ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[2+j*4] =      (            p[0+j*2]  - q[0+0*2] )
                    * (            p[1+j*2]  - q[1+0*2] );
    phi[3+j*4] =      ( q[0+2*2] - p[0+j*2]             )
                    * (            p[1+j*2]  - q[1+0*2] );

    dphidx[0+j*4] = - ( q[1+2*2] - p[1+j*2]             );
    dphidx[1+j*4] =   ( q[1+2*2] - p[1+j*2]             );
    dphidx[2+j*4] =   (            p[1+j*2]  - q[1+0*2] );
    dphidx[3+j*4] = - (            p[1+j*2]  - q[1+0*2] );

    dphidy[0+j*4] = - ( q[0+2*2] - p[0+j*2]             );
    dphidy[1+j*4] = - (            p[0+j*2]  - q[0+0*2] );
    dphidy[2+j*4] =   (            p[0+j*2]  - q[0+0*2] );
    dphidy[3+j*4] =   ( q[0+2*2] - p[0+j*2]             );
  }
/*
  Normalize.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
  return;
}
/******************************************************************************/

void basis_mn_q4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_Q4_TEST verifies BASIS_MN_Q4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double q[2*NODE_NUM] = {
    3.0, 1.0,
    5.0, 1.0,
    5.0, 4.0, 
    3.0, 4.0 };
  double sum_x;
  double sum_y;

  printf ( "\n" );
  printf ( "BASIS_MN_Q4_TEST:\n" );
  printf ( "  Verify basis functions for element Q4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( i = 0; i < NODE_NUM; i++ )
  {
    printf ( "  %3d  %10g  %10g\n", i, q[0+i*2], q[1+i*2] );
  }
 
  basis_mn_q4 ( q, NODE_NUM, q, phi, dphidx, dphidy );

  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
    printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T3: all bases at N points for a T3 element.

  Discussion:

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T3: 
       
            3
           / \
          /   \
         /     \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the coordinates of the vertices
    of the triangle.  It is common to list these points in counter clockwise
    order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the points where the basis functions 
    are to be evaluated.

    Output, double PHI[3*N], the value of the basis functions 
    at the evaluation points.

    Output, double DPHIDX[3*N], DPHIDY[3*N], the value of the 
    derivatives at the evaluation points.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  int i;
  int j;
  double temp;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BASIS_MN_T3 - Fatal error!\n" );
    fprintf ( stderr, "  Element has zero area.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*3] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*3] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*3] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*3] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*3] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*3] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*3] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*3] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*3] =      ( t[0+1*2] - t[0+0*2] );
  }
/*
  Normalize.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*3]    = phi[i+j*3]    / area;
      dphidx[i+j*3] = dphidx[i+j*3] / area;
      dphidy[i+j*3] = dphidy[i+j*3] / area;
    }
  }

  return;
}
/******************************************************************************/

void basis_mn_t3_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T3_TEST verifies BASIS_MN_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2006

  Author:

    John Burkardt

  Parameters:

    None.
*/
{
# define NODE_NUM 3

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T3_TEST:\n" );
  printf ( "  Verify basis functions for element T3.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t3 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t4 ( double t[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4: all bases at N points for a T4 element.

  Discussion:

    The T4 element is the cubic bubble triangle.

    The routine is given the coordinates of the vertices of a triangle.
    It works directly with these coordinates, and does not refer to a 
    reference element.

    The sides of the triangle DO NOT have to lie along a coordinate
    axis.

    The routine evaluates the basis functions associated with each vertex,
    and their derivatives with respect to X and Y.

  Physical Element T4: 
       
            3
           / \
          /   \
         /  4  \
        /       \
       1---------2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*4], the coordinates of the vertices
    of the triangle, and the coordinates of the centroid.  
    It is common to list the first three points in counter clockwise
    order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the points where the basis functions 
    are to be evaluated.

    Output, double PHI[4*N], the value of the basis functions 
    at the evaluation points.

    Output, double DPHIDX[4*N], DPHIDY[4*N], the value of the 
    derivatives at the evaluation points.

  Local parameters:

    Local, double AREA, is (twice) the area of the triangle.
*/
{
  double area;
  int i;
  int j;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*4] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*4] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*4] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*4] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*4] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*4] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*4] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*4] =      ( t[0+1*2] - t[0+0*2] );
  }
/*
  Normalize the first three functions.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
/*
  Compute the cubic bubble function.
*/
  for ( j = 0; j < n; j++ )
  {
       phi[3+j*4] = 27.0 * phi[0+j*4] * phi[1+j*4] * phi[2+j*4];

    dphidx[3+j*4] = 27.0 * (
                    dphidx[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidx[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidx[2+j*4] );

    dphidy[3+j*4] = 27.0 * (
                    dphidy[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidy[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidy[2+j*4] );
  }
/*
  Subtract 1/3 of the cubic bubble function from each of the three linears.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
         phi[i+j*4] =    phi[i+j*4] -    phi[3+j*4] / 3.0;
      dphidx[i+j*4] = dphidx[i+j*4] - dphidx[3+j*4] / 3.0;
      dphidy[i+j*4] = dphidy[i+j*4] - dphidy[3+j*4] / 3.0;
    }
  }

  return;
}
/******************************************************************************/

void basis_mn_t4_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T4_TEST verifies BASIS_MN_T4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2006

  Author:

    John Burkardt

  Parameters:
*/
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 2.0,
    0.0, 4.0,
    2.0, 2.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T4_TEST:\n" );
  printf ( "  Verify basis functions for element T4.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }

  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t4 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      printf ( "  %10g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }
  return;
# undef NODE_NUM
}
/******************************************************************************/

void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T6: all bases at N points for a T6 element.

  Discussion:

    The routine is given the coordinates of the vertices and midside
    nodes of a triangle.  It works directly with these coordinates, and does 
    not refer to a reference element.

    This routine requires that the midside nodes be "in line"
    with the vertices, that is, that the sides of the triangle be
    straight.  However, the midside nodes do not actually have to
    be halfway along the side of the triangle.  

  The physical element T6:

    This picture indicates the assumed ordering of the six nodes
    of the triangle.

    |
    |   
    |        3
    |       / \
    |      /   \
    Y     6     5
    |    /       \
    |   /         \
    |  1-----4-----2
    |
    +--------X-------->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 February 2006

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*6], the nodal oordinates of the element.
    It is common to list these points in counter clockwise order.

    Input, int N, the number of evaluation points.

    Input, double P[2*N], the coordinates of the point where
    the basis functions are to be evaluated.

    Output, double PHI[6*N], the value of the basis functions at P.

    Output, double DPHIDX[6*N], DPHIDY[6*N], the value of the X 
    and Y derivatives of the basis functions at P.
*/
{
  double gn;
  double gx;
  double hn;
  double hx;
  int j;

  for ( j = 0; j < n; j++ )
  {
/*
  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( p[1+j*2] - t[1+3*2] );

    hn = ( t[0+0*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( t[1+0*2] - t[1+3*2] );

    phi[0+j*6] =     ( gx * hx ) / ( gn * hn );
    dphidx[0+j*6] =  (      ( t[1+2*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+5*2] - t[1+3*2] ) ) / ( gn * hn );
    dphidy[0+j*6] = -(      ( t[0+2*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+5*2] - t[0+3*2] ) ) / ( gn * hn );
/*
  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
*/
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( p[1+j*2] - t[1+4*2] );

    hn = ( t[0+1*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( t[1+1*2] - t[1+4*2] );

    phi[1+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[1+j*6] =  (      ( t[1+2*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+3*2] - t[1+4*2] ) ) / ( gn * hn );
    dphidy[1+j*6] = -(      ( t[0+2*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+3*2] - t[0+4*2] ) ) / ( gn * hn );
/*
  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( p[1+j*2] - t[1+5*2] );

    hn = ( t[0+2*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( t[1+2*2] - t[1+5*2] );

    phi[2+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[2+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+4*2] - t[1+5*2] ) ) / ( gn * hn );
    dphidy[2+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+4*2] - t[0+5*2] ) ) / ( gn * hn );
/*
  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
*/
    gx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    gn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    hx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    hn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    phi[3+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[3+j*6] =  (      ( t[1+0*2] - t[1+2*2] ) * hx 
                     + gx * ( t[1+1*2] - t[1+2*2] ) ) / ( gn * hn );
    dphidy[3+j*6] = -(      ( t[0+0*2] - t[0+2*2] ) * hx 
                     + gx * ( t[0+1*2] - t[0+2*2] ) ) / ( gn * hn );
/*
  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
*/
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    hn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    phi[4+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[4+j*6] =  (      ( t[1+1*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+0*2] ) ) / ( gn * hn );
    dphidy[4+j*6] = -(      ( t[0+1*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+0*2] ) ) / ( gn * hn );
/*
  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
*/
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    hn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    phi[5+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[5+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+1*2] ) ) / ( gn * hn );
    dphidy[5+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+1*2] ) ) / ( gn * hn );
  }

  return;
}
/******************************************************************************/

void basis_mn_t6_test ( )

/******************************************************************************/
/*
  Purpose:

    BASIS_MN_T6_TEST verifies BASIS_MN_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 February 2006

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };

  printf ( "\n" );
  printf ( "BASIS_MN_T6_TEST:\n" );
  printf ( "  Verify basis functions for element T6.\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", NODE_NUM );

  printf ( "\n" );
  printf ( "  Physical Nodes:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    printf ( "  %3d  %10g  %10g\n", j, t[0+j*2], t[1+j*2] );
  }
 
  printf ( "\n" );
  printf ( "  The basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  basis_mn_t6 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      printf ( "  %7g", phi[i+j*NODE_NUM] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The X and Y derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dPhidX sum, dPhidY sum:\n" );
  printf ( "\n" );
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    printf ( "  %10g  %10g\n", sum_x, sum_y );
  }

  return;
# undef NODE_NUM
}
/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_CAP capitalizes a single character.

  Discussion:

    This routine should be equivalent to the library "toupper" function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 1998

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to capitalize.

    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
/******************************************************************************/

double degrees_to_radians ( double angle )

/******************************************************************************/
/*
  Purpose: 

    DEGREES_TO_RADIANS converts an angle from degrees to radians.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double ANGLE, an angle in degrees.

    Output, double DEGREES_TO_RADIANS, the equivalent angle
    in radians.
*/
{
# define PI 3.141592653589793

  return ( angle * PI / 180.0 );

# undef PI
}
/******************************************************************************/

void derivative_average_t3 ( int node_num, double node_xy[], int element_num,
  int element_node[], double c[], double dcdx[], double dcdy[] )

/******************************************************************************/
/*
  Purpose:

    DERIVATIVE_AVERAGE_T3 averages derivatives at the nodes of a T3 mesh.

  Discussion:

    This routine can be used to compute an averaged nodal value of any
    quantity associated with the finite element function.  At a node
    that is shared by several elements, the fundamental function
    U will be continuous, but its spatial derivatives, for instance,
    will generally be discontinuous.  This routine computes the
    value of the spatial derivatives in each element, and averages
    them, to make a reasonable assignment of a nodal value.

    Note that the ELEMENT_NODE array is assumed to be 1-based, rather
    than 0-based.  Thus, entries from this array must be decreased by
    1 before being used!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the element->node data.

    Input, double C[NODE_NUM], the finite element coefficient vector.

    Output, double DCDX[NODE_NUM], DCDY[NODE_NUM], the averaged
    values of dCdX and dCdY at the nodes.
*/
{
# define OFFSET 1

  int dim;
  double dphidx[3*3];
  double dphidy[3*3];
  int element;
  int j;
  int node;
  int node_count[node_num];
  int node_global1;
  int node_global2;
  int node_local1;
  int node_local2;
  double phi[3*3];
  double t[2*3];

  for ( node = 0; node < node_num; node++ )
  {
    node_count[node] = 0;
    dcdx[node] = 0.0;
    dcdy[node] = 0.0;
  }
/*
  Consider every element.
*/
  for ( element = 0; element < element_num; element++ )
  {
/*
  Get the coordinates of the nodes of the element.
*/
    for ( j = 0; j < 3; j++ )
    {
      for ( dim = 0; dim < 2; dim++ )
      {
        t[dim+2*j] = node_xy[dim+(element_node[j+element*3]-OFFSET)];
      }
    }
/*
  Evaluate the X and Y derivatives of the 3 basis functions at the
  3 nodes.
*/
    basis_mn_t3 ( t, 3, t, phi, dphidx, dphidy );
/*
  Evaluate dCdX and dCdY at each node in the element, and add
  them to the running totals.
*/
    for ( node_local1 = 0; node_local1 < 3; node_local1++ )
    {
      node_global1 = element_node[node_local1+element*3]-OFFSET;

      for ( node_local2 = 0; node_local2 < 3; node_local2++ )
      {
        node_global2 = element_node[node_local2+element*3]-OFFSET;

        dcdx[node_global1] = dcdx[node_global1]
          + c[node_global2] * dphidx[node_local2+node_local1*3];

        dcdy[node_global1] = dcdy[node_global1] 
          + c[node_global2] * dphidy[node_local2+node_local1*3];
      }
      node_count[node_global1] = node_count[node_global1] + 1;
    }
  }
/*
  Average the running totals.
*/
  for ( node = 0; node < node_num; node++ )
  {
    dcdx[node] = dcdx[node] / ( double ) node_count[node];
    dcdy[node] = dcdy[node] / ( double ) node_count[node];
  }
  return;
# undef OFFSET
}
/******************************************************************************/

void div_q4 ( int m, int n, double u[], double v[], double xlo, double xhi, 
  double ylo, double yhi, double div[], double vort[] )

/******************************************************************************/
/*
  Purpose: 

    DIV_Q4 estimates the divergence and vorticity of a discrete field.

  Discussion:

    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
    an array of points ( X(1:M), Y(1:N) ).

    The routine models the vector field over the interior of this region using
    a bilinear interpolant.  It then uses the interpolant to estimate the
    value of the divergence:

      DIV(X,Y) = dU/dX + dV/dY

    and the vorticity:

      VORT(X,Y) = dV/dX - dU/dY

    at the center point of each of the bilinear elements.

        |       |       |
      (3,1)---(3,2)---(3,3)---
        |       |       |
        | [2,1] | [2,2] |
        |       |       |
      (2,1)---(2,2)---(2,3)---
        |       |       |
        | [1,1] | [1,2] |
        |       |       |
      (1,1)---(1,2)---(1,3)---

    Here, the nodes labeled with parentheses represent the points at
    which the original (U,V) data is given, while the nodes labeled
    with square brackets represent the centers of the bilinear
    elements, where the approximations to the divergence and vorticity
    are made.

    The reason for evaluating the divergence and vorticity in this way
    is that the bilinear interpolant to the (U,V) data is not
    differentiable at the boundaries of the elements, nor especially at
    the nodes, but is an (infinitely differentiable) bilinear function
    in the interior of each element.  If a value at the original nodes
    is strongly desired, then the average at the four surrounding
    central nodes may be taken.

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of data rows.  M must be at least 2.

    Input, int N, the number of data columns.  N must be at least 2.

    Input, double U[M*N], V[M*N], the value of the components 
    of a vector quantity whose divergence and vorticity are desired. 
    A common example would be that U and V are the horizontal and 
    vertical velocity component of a flow field.

    Input, double XLO, XHI, the minimum and maximum X coordinates.

    Input, double YLO, YHI, the minimum and maximum Y coordinates.

    Output, double DIV[(M-1)*(N-1)], an estimate for
    the divergence in the bilinear element that lies between
    data rows I and I+1, and data columns J and J+1.

    Output, double VORT[(M-1)*(N-1)], an estimate for
    the vorticity in the bilinear element that lies between
    data rows I and I+1, and data columns J and J+1.
*/
{
  double dphidx[4];
  double dphidy[4];
  int i;
  int j;
  double p[2];
  double phi[4];
  double q[2*4];
  double xl;
  double xr;
  double yb;
  double yt;

  if ( m <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  M must be at least 2,\n" );
    fprintf ( stderr, "  but the input value of M is %d\n", m );
    exit ( 1 );
  }

  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  N must be at least 2,\n" );
    fprintf ( stderr, "  but the input value of N is %d\n", n );
    exit ( 1 );
  }

  if ( xhi == xlo )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  XHI and XLO must be distinct,\n" );
    fprintf ( stderr, "  but the input value of XLO is %g\n", xlo );
    fprintf ( stderr, "  and the input value of XHI is %g\n", xhi );
    exit ( 1 );
  }

  if ( yhi == ylo )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DIV_Q4 - Fatal error!\n" );
    fprintf ( stderr, "  YHI and YLO must be distinct,\n" );
    fprintf ( stderr, "  but the input value of YLO is %g\n", ylo );
    fprintf ( stderr, "  and the input value of YHI is %g\n", yhi );
    exit ( 1 );
  }

  for ( i = 1; i <= m-1; i++ )
  {
    yb = ( ( double ) ( 2 * m - 2 * i     ) * ylo   
         + ( double ) (         2 * i - 2 ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );
    p[1] = ( ( double ) ( 2 * m - 2 * i - 1 ) * ylo   
           + ( double ) (         2 * i - 1 ) * yhi ) 
           / ( double ) ( 2 * m         - 2 );
    yt = ( ( double ) ( 2 * m - 2 * i - 2 ) * ylo   
         + ( double ) (         2 * i     ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );

    q[1+0*2] = yb;
    q[1+1*2] = yb;
    q[1+2*2] = yt;
    q[1+3*2] = yt;

    for ( j = 1; j <= n-1; j++ )
    {
      xl = ( ( double ) ( 2 * n - 2 * j     ) * xlo   
           + ( double ) (         2 * j - 2 ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );
      p[0] = ( ( double ) ( 2 * n - 2 * j - 1 ) * xlo   
             + ( double ) (         2 * j - 1 ) * xhi ) 
             / ( double ) ( 2 * n         - 2 );
      xr = ( ( double ) ( 2 * n - 2 * j - 2 ) * xlo   
           + ( double ) (         2 * j     ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );

      q[0+0*2] = xl;
      q[0+1*2] = xr;
      q[0+2*2] = xr;
      q[0+3*2] = xl;
/*
  Evaluate the basis function and derivatives at the center of the element.
*/
      basis_mn_q4 ( q, 1, p, phi, dphidx, dphidy );
/*
  Note the following formula for the value of U and V at the same
  point that the divergence and vorticity are being evaluated.

         umid =  u(i  ,j  ) * phi[0] &
               + u(i  ,j+1) * phi[1] &
               + u(i+1,j+1) * phi[2] &
               + u(i+1,j  ) * phi[3] 

         vmid =  v(i  ,j  ) * phi[0] &
               + v(i  ,j+1) * phi[1] &
               + v(i+1,j+1) * phi[2] &
               + v(i+1,j  ) * phi[3] 
*/
      div[i-1+(j-1)*(m-1)] =  
                  u[i-1+(j-1)*m] * dphidx[0] + v[i-1+(j-1)*m] * dphidy[0] 
                + u[i-1+(j  )*m] * dphidx[1] + v[i-1+(j  )*m] * dphidy[1] 
                + u[i  +(j  )*m] * dphidx[2] + v[i  +(j  )*m] * dphidy[2] 
                + u[i  +(j-1)*m] * dphidx[3] + v[i  +(j-1)*m] * dphidy[3];

      vort[i-1+(j-1)*(m-1)] =  
                   v[i-1+(j-1)*m] * dphidx[0] - u[i-1+(j-1)*m] * dphidy[0] 
                 + v[i-1+(j  )*m] * dphidx[1] - u[i-1+(j  )*m] * dphidy[1] 
                 + v[i  +(j  )*m] * dphidx[2] - u[i  +(j  )*m] * dphidy[2] 
                 + v[i  +(j-1)*m] * dphidx[3] - u[i  +(j-1)*m] * dphidy[3]; 
    }
  }

  return;
}
/******************************************************************************/

char *element_code ( int i )

/******************************************************************************/
/*
  Purpose:

    ELEMENT_CODE returns the code for each element.

  List:

    I  ELEMENT_CODE   Definition
    -  ------------   ----------
    1  Q4             4 node linear Lagrange/serendipity quadrilateral;
    2  Q8             8 node quadratic serendipity quadrilateral;
    3  Q9             9 node quadratic Lagrange quadrilateral;
    4  Q12            12 node cubic serendipity quadrilateral;
    5  Q16            16 node cubic Lagrange quadrilateral;
    6  QL             6 node linear/quadratic quadrilateral;
    7  T3             3 node linear triangle;
    8  T4             4 node cubic bubble triangle
    9  T6             6 node quadratic triangle;
   10  T10            10 node cubic triangle.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number of the element.  

    Output, char *ELEMENT_CODE, the code for the element.
*/
{
  char *value;

  if ( i == 1 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q4" );
  }
  else if ( i == 2 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q8" );
  }
  else if ( i == 3 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "Q9" );
  }
  else if ( i == 4 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "Q12" );
  }
  else if ( i == 5 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "Q16" );
  }
  else if ( i == 6 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "QL" );
  }
  else if ( i == 7 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T3" );
  }
  else if ( i == 8 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T4" );
  }
  else if ( i == 9 )
  {
    value = ( char * ) malloc ( 3 * sizeof ( char ) );
    strcpy ( value, "T6" );
  }
  else if ( i == 10 )
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "T10" );
  }
  else
  {
    value = ( char * ) malloc ( 4 * sizeof ( char ) );
    strcpy ( value, "???" );
  }

  return value;
}
/******************************************************************************/

void elements_eps ( char *file_name, int node_num, double node_xy[], char *code, 
  int element_num, int element_mask[], int element_node[], int node_show, 
  int element_show )

/******************************************************************************/
/*
  Purpose: 

    ELEMENTS_EPS creates an EPS file image of the elements of a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to create.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.

    Input, char *CODE, the code for the element.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_MASK[ELEMENT_NUM], a mask for the elements.
    Only elements with a TRUE mask will be shown.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes making up
    each element.

    Input, int NODE_SHOW:
    0, do not show nodes;
    1, show nodes;
    2, show nodes and label them.

    Input, int TRIANGLE_SHOW:
    0, do not show triangles;
    1, show triangles;
    2, show triangles and label them.
*/
{
  double ave_x;
  double ave_y;
  int circle_size = 3;
  int delta;
  int e;
  int element;
  int element_order;
  FILE *file_unit;
  int i;
  int j;
  int local;
  int node;
  int *node_mask;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;

  element_order = order_code ( code );
/*
  Determine which nodes are visible, controlled by which elements are visible.
*/
  node_mask = ( int * ) malloc ( node_num * sizeof ( int ) );
  for ( node = 0; node < node_num; node++ )
  {
    node_mask[node] = 0;
  }
  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element*element_order]-1;
        node_mask[node] = 1;
      }
    }
  }
/*
  We need to do some figuring here, so that we can determine
  the range of the data, and hence the height and width
  of the piece of paper.
*/
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_mask[node] )
     {
       if ( x_max < node_xy[0+node*2] )
       {
         x_max = node_xy[0+node*2];
       }
    }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
       if ( node_xy[0+node*2] < x_min )
       {
         x_min = node_xy[0+node*2];
       }
    }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( y_max < node_xy[1+node*2] )
      {
        y_max = node_xy[1+node*2];
      }
    }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( node_xy[1+node*2] < y_min )
      {
        y_min = node_xy[1+node*2];
      }
    }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min ) 
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit = fopen ( file_name, "wt" );

  if ( !file_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ELEMENTS_EPS - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output EPS file.\n" );
    exit ( 1 );
  }

  fprintf ( file_unit, "%!PS-Adobe-3.0 EPSF-3.0\n" );
  fprintf ( file_unit, "%%Creator: elements_eps.C\n" );
  fprintf ( file_unit, "%%Title: %s\n", file_name );

  fprintf ( file_unit, "%%Pages: 1\n" );
  fprintf ( file_unit, "%%BoundingBox:  %d  %d  %d  %d\n", 
    x_ps_min, y_ps_min, x_ps_max, y_ps_max );
  fprintf ( file_unit, "%%Document-Fonts: Times-Roman\n" );
  fprintf ( file_unit, "%%LanguageLevel: 1\n" );
  fprintf ( file_unit, "%%EndComments\n" );
  fprintf ( file_unit, "%%BeginProlog\n" );
  fprintf ( file_unit, "/inch {72 mul} def\n" );
  fprintf ( file_unit, "%%EndProlog\n" );
  fprintf ( file_unit, "%%Page:      1     1\n" );
  fprintf ( file_unit, "save\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "% Set the RGB line color to very light gray.\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, " 0.9000 0.9000 0.9000 setrgbcolor\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "% Draw a gray border around the page.\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "stroke\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "% Set RGB line color to black.\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, " 0.0000 0.0000 0.0000 setrgbcolor\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Set the font and its size:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.50 inch scalefont\n" );
  fprintf ( file_unit, "setfont\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Print a title:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  210  702 moveto\n" );
  fprintf ( file_unit, "%(Pointset) show\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "% Define a clipping polygon\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "clip newpath\n" );
/*
  Draw the nodes.
*/
  if ( 1 <= node_show )
  {
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "%  Draw filled dots at each node:\n" );
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "%  Set the color to blue:\n" );
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "0.000  0.150  0.750  setrgbcolor\n" );
    fprintf ( file_unit, "%\n" );

    for ( node = 0; node < node_num; node++ )
    {
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "newpath  %d  %d  %d  0 360 arc closepath fill\n",
          x_ps, y_ps, circle_size );
      }
    }
  }
/*
  Label the nodes.
*/
  if ( 2 <= node_show )
  {
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "%  Label the nodes:\n" );
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "%  Set the color to darker blue:\n" );
    fprintf ( file_unit, "%\n" );
    fprintf ( file_unit, "0.000  0.250  0.850  setrgbcolor\n" );
    fprintf ( file_unit, "/Times-Roman findfont\n" );
    fprintf ( file_unit, "0.20 inch scalefont\n" );
    fprintf ( file_unit, "setfont\n" );

    fprintf ( file_unit, "%\n" );

    for ( node = 0; node < node_num; node++ )
    { 
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n",
          x_ps, y_ps + 5, node+1 );
      }
    }
  }
/*
  Draw the elements.
*/
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Draw the element sides:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, " 9.0000 0.0000 0.0000 setrgbcolor\n" );

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      local = 1;
      node = element_node[local-1+element*element_order] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max                     - y_min ) );

      fprintf ( file_unit, "newpath  %d  %d\n", x_ps, y_ps );

      for ( ; ; )
      {
        local = next_boundary_node ( local, code );
        node = element_node[local-1+element*element_order] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "  %d  %d  lineto", x_ps, y_ps ); 

        if ( local == 1 )
        {
          break;
        }
      }
      fprintf ( file_unit, "stroke\n" );
    }
  }
/*
  Label the elements.
*/
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Label the elements:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, " 1.0000 0.0000 0.0000 setrgbcolor\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.30 inch scalefont setfont\n" );

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element_order*element] - 1;

        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / ( double ) ( element_order );
      ave_y = ave_y / ( double ) ( element_order );

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )  
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )  
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max         - y_min ) );

      fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n", 
        x_ps, y_ps + 5, element+1 );
    }
  }
/*
  Finish up.
*/
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "restore  showpage\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  End of page.\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%%Trailer\n" );
  fprintf ( file_unit, "%%EOF\n" );

  close ( file_unit );

  free ( node_mask );

  return;
}
/******************************************************************************/

int *grid_element ( char *code, int element_order, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_ELEMENT returns the element grid associated with any available element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
    "T4", "T6" and "T10".

    Input, int ELEMENT_ORDER, the number of nodes per element.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY for quadrilaterals, or 2 * NELEMX * NELEMY for
    triangles.

    Output, int GRID_ELEMENT[ELEMENT_ORDER*ELEMENT_NUM], the nodes 
    that form each element.  
*/
{
  int *element_node;

  if ( s_eqi ( code, "Q4" ) )
  {
    element_node = grid_q4_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    element_node = grid_q8_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    element_node = grid_q9_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    element_node = grid_q12_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    element_node = grid_q16_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    element_node = grid_ql_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    element_node = grid_t3_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    element_node = grid_t4_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    element_node = grid_t6_element ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    element_node = grid_t10_element ( nelemx, nelemy );
  }
  else
  {
    element_node = NULL;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_ELEMENT - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    exit ( 1 );
  }

  return element_node;
}
/******************************************************************************/

int grid_element_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_ELEMENT_NUM returns the number of elements in a grid.

  Discussion:

    The number of elements generated will be NELEMX * NELEMY for
    quadrilaterals, or 2 * NELEMX * NELEMY for triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
    'T4', 'T6' and 'T10'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int GRID_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  if ( s_eqi ( code, "Q4" ) )
  {
    element_num = grid_q4_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    element_num = grid_q8_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    element_num = grid_q9_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    element_num = grid_q12_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    element_num = grid_q16_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    element_num = grid_ql_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    element_num = grid_t3_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    element_num = grid_t4_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    element_num = grid_t6_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    element_num = grid_t10_element_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_ELEMENT_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
/******************************************************************************/

int grid_node_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_NODE_NUM returns the number of nodes in a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
    'T4', 'T6' and 'T10'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int GRID_NODE_NUM, the number of elements in the grid.
*/
{
  int node_num;

  if ( s_eqi ( code, "Q4" ) )
  {
    node_num = grid_q4_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    node_num = grid_q8_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    node_num = grid_q9_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    node_num = grid_q12_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    node_num = grid_q16_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    node_num = grid_ql_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    node_num = grid_t3_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    node_num = grid_t4_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    node_num = grid_t6_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    node_num = grid_t10_node_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GRID_NODE_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
/******************************************************************************/

double *grid_nodes_01 ( int x_num, int y_num )

/******************************************************************************/
/*
  Purpose:

    GRID_NODES_01 returns an equally spaced rectangular grid of nodes in the unit square.

  Example:

    X_NUM = 5
    Y_NUM = 3

    NODE_XY =
    ( 0, 0.25, 0.5, 0.75, 1, 0,   0.25, 0.5, 0.75, 1,   0, 0.25, 0.5, 0.75, 1;
      0, 0,    0,   0,    0, 0.5, 0.5,  0.5, 0.5,  0.5, 1, 1.0,  1.0, 1.0,  1 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, Y_NUM, the number of nodes in the X and Y directions.

    Output, double GRID_NODES_01[2*X_NUM*Y_NUM], the coordinates of
    the nodes.
*/
{
  int i;
  int j;
  int k;
  int node_num;
  double *node_xy;
  double value;

  node_num = x_num * y_num;

  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );

  if ( x_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[0+k*2] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      value = ( double ) ( i ) / ( double ) ( x_num - 1 );
      for ( j = i; j < node_num; j = j + x_num )
      {
        node_xy[0+j*2] = value;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[1+k*2] = 0.5;
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      value = ( double ) ( j ) / ( double ) ( y_num - 1 );
      for ( i = j*x_num; i < ( j + 1 ) * x_num; i++ )
      {
        node_xy[1+i*2] = value;
      }
    }
  }

  return node_xy;
}
/******************************************************************************/

void grid_print ( int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose: 

    GRID_PRINT prints the elements that form a grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the number of nodes per element.

    Input, int ELEMENT_NUM, the number of elements.
 
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
    each element.
*/
{
  int element;
  int i;

  printf ( "\n" );
  printf ( "  GRID_PRINT: Element -> Node table.\n" );
  printf ( "\n" );
  printf ( "  Element order = %d\n", element_order );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "\n" );
  printf ( "    #   " );
  for ( i = 0; i < element_order; i++ )
  {
    printf ( "%3d", i );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( element = 0; element < element_num; element++ )
  {
    printf ( "  %3d   ", element );
    for ( i = 0; i < element_order; i++ )
    {
      printf ( "%3d", element_node[i+element*element_order] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

int *grid_q4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q4_ELEMENT produces a grid of 4 node quadrilaterals.

  Discussion:

    For each element, the nodes are listed in counter-clockwise order.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1, 2,  6,  5;
         2, 3,  7,  6;
         3, 4,  8,  7;
         5, 6, 10,  9;
         6, 7, 11, 10;
         7, 8, 12, 11.

  Grid:

    9---10---11---12
    |    |    |    |
    |    |    |    |
    |  4 |  5 |  6 |
    |    |    |    |
    5----6----7----8
    |    |    |    |
    |    |    |    |
    |  1 |  2 |  3 |
    |    |    |    |
    1----2----3----4

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q4[4*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW---NE
     |    |
    SW---SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );
  
      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q4_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q4_ELEMENT_NUM counts the elements in a grid of 4 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q4_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q4_NODE_NUM counts the nodes in a grid of 4 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q4_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_q8_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q8_ELEMENT produces a grid of 8 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE =
         1,  3, 14, 12,  2,  9, 13,  8;
         3,  5, 16, 14,  4, 10, 15,  9;
         5,  7, 18, 16,  6, 11, 17, 10;
        12, 14, 25, 23, 13, 20, 24, 19;
        14, 16, 27, 25, 15, 21, 26, 20;
        16, 18, 29, 27, 17, 22, 28, 21.

  Diagram:

   23---24---25---26---27---28---29
    |         |         |         |
    |         |         |         |
   19        20        21        22
    |         |         |         |
    | 4       | 5       | 6       |
   12---13---14---15---16---17---18
    |         |         |         |
    |         |         |         |
    8         9        10        11
    |         |         |         |
    | 1       | 2       | 3       |
    1----2----3----4----5----6----7

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q8[8*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int e;
  int element;
  int *element_node;
  int element_order = 8;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----N----NE
     |          |
     W   (C)    E
     |          |
    SW----S----SE
*/

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = ( j - 1 )  * ( 3 * nelemx + 2 ) + 2 * i - 1;
      w  = sw + 2 * nelemx + 2 - i;
      nw = sw + 3 * nelemx + 2;

      s =  sw + 1;
      n =  sw + ( 3 * nelemx + 2 ) + 1;

      se = sw + 2;
      e  = sw + 2 * nelemx + 2 - i + 1;
      ne = sw + ( 3 * nelemx + 2 ) + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q8_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q8_ELEMENT_NUM counts the elements in a grid of 8 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q8_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q8_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q8_NODE_NUM counts the nodes in a grid of 8 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q8_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 3 * nelemx * nelemy + 2 * nelemx + 2 * nelemy + 1;

  return node_num;
}
/******************************************************************************/

int *grid_q9_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q9_ELEMENT produces a grid of 9 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  3, 17, 15,  2, 10, 16,  8,  9;
         3,  5, 19, 17,  4, 12, 18, 10, 11;
         5,  7, 21, 19,  6, 14, 20, 12, 13;
        15, 17, 31, 29, 16, 24, 30, 22, 23;
        17, 19, 33, 31, 18, 26, 32, 24, 25;
        19, 21, 35, 33, 20, 28, 34, 26, 27.

  Grid:

   29---30---31---32---33---34---35
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
   22 . 23 . 24 . 25 . 26 . 27 . 28
    |    .    |    .    |    .    |
    | 4  .    | 5  .    | 6  .    |
   15---16---17---18---19---20---21
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
    8 .  9 . 10 . 11 . 12 . 13 . 14
    |    .    |    .    |    .    |
    | 1  .    | 2  .    | 3  .    |
    1----2----3----4----5----6----7

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q9[9*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 9;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----N----NE
     |          |
     W    C     E
     |          |
    SW----S----SE
*/
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;
      element_node[8+element*element_order] = c;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q9_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q9_ELEMENT_NUM counts the elements in a grid of 9 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q9_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q9_NODE_NUM counts the nodes in a grid of 9 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 March 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q9_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_q12_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q12_ELEMENT produces a grid of 12 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE =
         1,  2,  3,  4, 11, 12, 15, 16, 19, 20, 21, 22;
         4,  5,  6,  7, 12, 13, 16, 17, 22, 23, 24, 25;
         7,  8,  9, 10, 13, 14, 17, 18, 25, 26, 27, 28;
        19, 20, 21, 22, 29, 30, 33, 34, 37, 38, 39, 40;
        22, 23, 24, 25, 30, 31, 34, 35, 40, 41, 42, 43;
        25, 26, 27, 28, 31, 32, 35, 36, 43, 44, 45, 46.

  Grid:

   37-38-39-40-41-42-43-44-45-46
    |        |        |        |
   33       34       35       36
    |        |        |        |
   29       30       31       32
    | 4      | 5      | 6      |
   19-20-21-22-23-24-25-26-27-28
    |        |        |        |
   15       16       17       18
    |        |        |        |
   11       12       13       14
    | 1      | 2      | 3      |
    1--2--3--4--5--6--7--8--9-10

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q12[12*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 12;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 5 * nelemx + 3 ) + 1;

      element_node[ 0+element*element_order] =  base + ( i - 1 ) * 3;
      element_node[ 1+element*element_order] =  base + ( i - 1 ) * 3 + 1;
      element_node[ 2+element*element_order] =  base + ( i - 1 ) * 3 + 2;
      element_node[ 3+element*element_order] =  base + ( i - 1 ) * 3 + 3;

      element_node[ 4+element*element_order] =  base + 3 * nelemx + i;
      element_node[ 5+element*element_order] =  base + 3 * nelemx + i + 1;

      element_node[ 6+element*element_order] =  base + 4 * nelemx + i + 1;
      element_node[ 7+element*element_order] =  base + 4 * nelemx + i + 2;

      element_node[ 8+element*element_order] =  base + 5 * nelemx + 3 * i;
      element_node[ 9+element*element_order] = base + 5 * nelemx + 3 * i + 1;
      element_node[10+element*element_order] = base + 5 * nelemx + 3 * i + 2;
      element_node[11+element*element_order] = base + 5 * nelemx + 3 * i + 3;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q12_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q12_ELEMENT_NUM counts the elements in a grid of 12 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q12_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q12_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q12_NODE_NUM counts the nodes in a grid of 12 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q12_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 5 * nelemx * nelemy + 3 * nelemx + 3 * nelemy + 1;

  return node_num;
}
/******************************************************************************/

int *grid_q16_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q16_ELEMENT produces a grid of 16 node quadrilaterals.

  Example:

    Input:

      NELEMX = 2, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  4,  8,  9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25;
         4,  5,  6,  7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28;
        22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46;
        25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49. 
        
  Grid:

   43-44-45-46-47-48-49
    |        |        |
    |        |        |
   36 37 38 39 40 41 42
    |        |        |
    |        |        |
   29 30 31 32 33 34 35
    |        |        |
    | 3      | 4      |
   22-23-24-25-26-27-28
    |        |        |
    |        |        |
   15 16 17 18 19 20 21
    |        |        |
    |        |        |
    8  9 10 11 12 13 14
    |        |        |
    | 1      | 2      |
    1--2--3--4--5--6--7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q16[16*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 16;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[ 0+element*element_order] = base;
      element_node[ 1+element*element_order] = base                          + 1;
      element_node[ 2+element*element_order] = base                          + 2;
      element_node[ 3+element*element_order] = base                          + 3;
      element_node[ 4+element*element_order] = base +     ( 3 * nelemx + 1 );
      element_node[ 5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 1;
      element_node[ 6+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[ 7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[ 8+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[ 9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[10+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element_node[11+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[12+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[13+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[14+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[15+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_q16_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q16_ELEMENT_NUM counts the elements in a grid of 16 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_q16_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_Q16_NODE_NUM counts the nodes in a grid of 16 node quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_Q16_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_ql_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_QL_ELEMENT produces a grid of 6 node quadratics/linears.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  8,  9, 10;
         3,  4,  5, 10, 11, 12;
         5,  6,  7, 12, 13, 14;
         8,  9, 10, 15, 16, 17;
        10, 11, 12, 17, 18, 19;
        12, 13, 14, 19, 20, 21.

  Grid:

   15---16---17---18---19---20---21
    |         |         |         |
    |         |         |         |
    |    4    |    5    |    6    |
    |         |         |         |
    |         |         |         |
    8----9---10---11---12---13---14
    |         |         |         |
    |         |         |         |
    |    1    |    2    |    3    |
    |         |         |         |
    |         |         |         |
    1----2----3----4----5----6----7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.  X will the the "quadratic direction", and
    Y will be the "linear direction".

    Output, int GRID_QL[6*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * i - 1;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base + 1;
      element_node[2+element*element_order] = base + 2;
      element_node[3+element*element_order] = base + ( 2 * nelemx + 1 );
      element_node[4+element*element_order] = base + ( 2 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base + ( 2 * nelemx + 1 ) + 2;

      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_ql_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_QL_ELEMENT_NUM counts the elements in a grid of quadratic/linear quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int GRID_QL_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_ql_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_QL_NODE_NUM counts the nodes in a grid of quadratic/linear quadrilaterals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_QL_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 2 * nelemx * nelemy + 2 * nelemx + nelemy + 1;

  return node_num;
}
/******************************************************************************/

void grid_shape_2d ( int n, double a[], int *n1, int *n2 )

/******************************************************************************/
/*
  Purpose: 

    GRID_SHAPE_2D guesses the shape N1 by N2 of a vector of data.

  Discussion:

    The data vector A is assumed to contain N1 * N2 values, with
    where each of N2 values is repeated N1 times.

  Example:

    Input:

      A = ( 2, 2, 2, 7, 7, 7 )

    Output:

      N1 = 3, N2 = 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input, double A[N], the data, which should have the properties
    described above.

    Output, int *N1, *N2, the "shape" of the data in the array.
*/
{
  int i;
/*
  Make a guess for N1.
*/
  i = 1;
  *n1 = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[0] )
    {
      break;
    }
    *n1 = *n1 + 1;
  }
/*
  Guess that N2 = N / N1.
*/
  *n2 = n / (*n1);

  return;
}
/******************************************************************************/

int *grid_t3_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T3_ELEMENT produces a grid of pairs of 3 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  5;
         6,  5,  2;
         2,  3,  6;
         7,  6,  3;
         3,  4,  7;
         8,  7,  4;
         5,  6,  9;
        10,  9,  6;
         6,  7, 10;
        11, 10,  7;
         7,  8, 11;
        12, 11,  8.

  Grid:

    9---10---11---12
    |\ 8 |\10 |\12 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  7\|  9\| 11\|
    5----6----7----8
    |\ 2 |\ 4 |\ 6 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  1\|  3\|  5\|
    1----2----3----4

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T3[3*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW--NE
     |\ |
     | \|
    SW--SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t3_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T3_ELEMENT_NUM counts the elements in a grid of 3 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T3_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t3_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T3_NODE_NUM counts the nodes in a grid of 3 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T3_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_t4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T4_ELEMENT produces a grid of pairs of 4 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  11,  5;
        12, 11,   2,  8;
         2,  3,  12,  6;
        13, 12,   3,  9;
         3   4   13,  7;
        14, 13,   4,  10;
        11, 12,  21,  15;
        22, 21,  12,  18;
        12, 13,  22,  16;
        23, 22,  13,  19;
        13  14   23,  17;
        24, 23,  14,  20;

  Grid:

   21---22---23---24
    |\18 |\19 |\20 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    | 15\| 16\| 17\|
   11---12---13---14
    |\ 8 |\ 9 |\10 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    | 5 \|  6\|  7\|
    1----2----3----4

  Element T4:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  | 4 \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T4[4*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int nc;
  int ne;
  int nw;
  int sc;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----NE
     |\   |
     | \NC|
     |SC\ |
     |   \|
    SW---SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( 3 * nelemx + 1 );
      se = sw + 1;
      sc = sw +     nelemx + 1;
      nc = sw + 2 * nelemx + 1;
      nw = sw + 3 * nelemx + 1;
      ne = sw + 3 * nelemx + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = sc;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = nc;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t4_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T4_ELEMENT_NUM counts the elements in a grid of 4 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T4_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t4_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T4_NODE_NUM counts the nodes in a grid of 4 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T4_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 ) + 2 * nelemx * nelemy;

  return node_num;
}
/******************************************************************************/

int *grid_t6_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T6_ELEMENT produces a grid of pairs of 6 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  3, 15,  2,  9,  8;
        17, 15,  3, 16,  9, 10;
         3,  5, 17,  4, 11, 10;
        19, 17,  5, 18, 11, 12;
         5,  7, 19,  6, 13, 12;
        21, 19,  7, 20, 13, 14;
        15, 17, 29, 16, 23, 22;
        31, 29, 17, 30, 23, 24;
        17, 19, 31, 18, 25, 24;
        33, 31, 19, 32, 25, 26;
        19, 21, 33, 20, 27, 26;
        35, 33, 21, 34, 27, 28.

  Grid:

   29-30-31-32-33-34-35
    |\ 8  |\10  |\12  |
    | \   | \   | \   |
   22 23 24 25 26 27 28
    |   \ |   \ |   \ |
    |  7 \|  9 \| 11 \|
   15-16-17-18-19-20-21
    |\ 2  |\ 4  |\ 6  |
    | \   | \   | \   |
    8  9 10 11 12 13 14
    |   \ |   \ |   \ |
    |  1 \|  3 \|  5 \|
    1--2--3--4--5--6--7

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T6[6*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW---N--NE
     | \     |
     W   C   E
     |    \  |
    SW---S--SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = s;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = w;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = n;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = e;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t6_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T6_ELEMENT_NUM counts the elements in a grid of 6 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T6_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t6_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T6_NODE_NUM counts the nodes in a grid of 6 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T6_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

int *grid_t10_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_T10_ELEMENT produces a grid of pairs of 10 node triangles.

  Example:

    Input:

      NELEMX = 2, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  3,  4, 10, 16, 22, 15,  8,  9;
        25, 24, 23, 22, 16, 10,  4, 11, 18, 17;
         4,  5,  6,  7, 13, 19, 25, 18, 11, 12;
        28, 27, 26, 25, 19, 13,  7, 14, 21, 20;
        22, 23, 24, 25, 31, 37, 43, 36, 29, 30;
        46, 45, 44, 43, 37, 31, 25, 32, 39, 38;
        25, 26, 27, 28, 34, 40, 46, 39, 31, 33;
        49, 48, 47, 46, 40, 34, 28, 35, 42, 41.
        
  Grid:

   43-44-45-46-47-48-49
    |\     6 |\     8 |
    | \      | \      |
   36 37 38 39 40 41 42
    |   \    |   \    |
    |    \   |    \   |
   29 30 31 32 33 34 35
    |      \ |      \ |
    | 5     \| 7     \|
   22-23-24-25-26-27-28
    |\     2 |\     4 |
    | \      | \      |
   15 16 17 18 19 20 21
    |   \    |   \    |
    |    \   |    \   |
    8  9 10 11 12 13 14
    |      \ |      \ |
    | 1     \| 3     \|
    1--2--3--4--5--6--7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T10[10*2*NELEMX*NELEMY], the nodes that form
    each element.
*/
{
  int base;
  int element;
  int *element_node;
  int element_order = 10;
  int i;
  int j;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base                          + 1;
      element_node[2+element*element_order] = base                          + 2;
      element_node[3+element*element_order] = base                          + 3;
      element_node[4+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[5+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[6+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[7+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[8+element*element_order] = base +     ( 2 * nelemx + 1 ) + 2;
      element_node[9+element*element_order] = base +     ( 2 * nelemx + 1 ) + 3;
      element = element + 1;

      element_node[0+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;
      element_node[1+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[2+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[3+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[4+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[6+element*element_order] = base                          + 3;
      element_node[7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[8+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element = element + 1;
    }
  }

  return element_node;
}
/******************************************************************************/

int grid_t10_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T10_ELEMENT_NUM counts the elements in a grid of 10 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    2 * NELEMX * NELEMY.

    Output, int GRID_T10_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int grid_t10_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    GRID_T10_NODE_NUM counts the nodes in a grid of 10 node triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, int GRID_T10_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
/******************************************************************************/

void grid_test ( char *code )

/******************************************************************************/
/*
  Purpose: 

    GRID_TEST tests the grid routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, the code for the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.
*/
{
  int *element_node;
  int element_num;
  int element_order;
  int nelemx;
  int nelemy;
  int width;
/*
  NODE is defined as a vector rather than a two dimensional array,
  so that we can handle the various cases using a single array.
*/
  printf ( "\n" );
  printf ( "  GRID_TEST: Test the grid routine for element %s\n", code );

  nelemx = 3;
  nelemy = 2;

  if ( s_eqi ( code, "Q4" ) || 
       s_eqi ( code, "Q8" ) || 
       s_eqi ( code, "Q9" ) ||
       s_eqi ( code, "Q12" ) ||
       s_eqi ( code, "Q16" ) ||
       s_eqi ( code, "QL" ) )
  {
    element_num = nelemx * nelemy;
  }
  else if ( strcmp ( code, "T3" ) == 0 || 
            strcmp ( code, "T4" ) == 0 || 
            strcmp ( code, "T6" ) == 0 || 
            strcmp ( code, "T10" ) == 0 )
  {
    element_num = 2 * nelemx * nelemy;
  }

  element_order = order_code ( code );

  element_node = grid_element ( code, element_order, nelemx, nelemy );

  grid_print ( element_order, element_num, element_node );

  width = grid_width ( element_order, element_num, element_node );

  printf ( "\n" );
  printf ( "  Grid width is %d\n", width );

  free ( element_node );

  return;
}
/******************************************************************************/

int grid_width ( int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose: 

    GRID_WIDTH computes the width of a given grid.

  Definition:

    The grid width is defined to the maximum absolute
    difference of global indices of nodes in the same element.

  Example:

    For the following grid, the grid width is 13.

   23---24---25---26---27---28---29
    |         |         |         |
    |         |         |         |
   19        20        21        22
    |         |         |         |
    | 4       | 5       | 6       |
   12---13---14---15---16---17---18
    |         |         |         |
    |         |         |         |
    8         9        10        11
    |         |         |         |
    | 1       | 2       | 3       |
    1----2----3----4----5----6----7

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.
 
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
    each element.

    Output, int GRID_WIDTH, the grid width.
*/
{
  int element;
  int ip1;
  int ip2;
  int node1;
  int node2;
  int width;

  width = 0;
 
  for ( element = 0; element < element_num; element++ )
  {
    for ( node1 = 0; node1 < element_order; node1++ )
    {
      ip1 = element_node[node1+element*element_order];
      for ( node2 = 0; node2 < element_order; node2++ )
      {
        ip2 = element_node[node2+element*element_order];
        width = i4_max ( width, abs ( ip1 - ip2 ) );
      }
    }
  }
 
  return width;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_modp ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_MODP returns the nonnegative remainder of I4 division.

  Discussion:

    If
      NREM = I4_MODP ( I, J )
      NMULT = ( I - NREM ) / J
    then
      I = J * NMULT + NREM
    where NREM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, I4_MODP(A,360) is between 0 and 360, always.

  Example:

        I         J     MOD  I4_MODP   I4_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number to be divided.

    Input, int J, the number that divides I.

    Output, int I4_MODP, the nonnegative remainder when I is
    divided by J.
*/
{
  int value;

  if ( j == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_MODP - Fatal error!\n" );
    fprintf ( stderr, "  I4_MODP ( I, J ) called with J = %d\n", j );
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
/******************************************************************************/

int i4_wrap ( int ival, int ilo, int ihi )

/******************************************************************************/
/*
  Purpose:

    I4_WRAP forces an I4 to lie between given limits by wrapping.

  Example:

    ILO = 4, IHI = 8

    I   Value

    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 December 2012

  Author:

    John Burkardt

  Parameters:

    Input, int IVAL, an integer value.

    Input, int ILO, IHI, the desired bounds for the integer value.

    Output, int I4_WRAP, a "wrapped" version of IVAL.
*/
{
  int jhi;
  int jlo;
  int value;
  int wide;

  if ( ilo < ihi )
  {
   jlo = ilo;
   jhi = ihi;
  }
  else
  {
    jlo = ihi;
    jhi = ilo;
  }

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_WRITE writes an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the data.
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
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
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

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

void interp ( char *code, int element_order, double r, double s, 
  double ubase[], double *u, double *dudr, double *duds )

/******************************************************************************/
/*
  Purpose: 

    INTERP interpolates a quantity in an element from basis node values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T6' and 'T10'.

    Input, int ELEMENT_ORDER, order of the element.

    Input, double R, S, the reference coordinates of a point.

    Input, double UBASE[ELEMENT_ORDER], the value of the quantity 
    at the basis nodes.

    Output, double *U, *DUDR, *DUDS, the interpolated value of the
    quantity and its derivatives at the point (R,S).
*/
{
  double *dtdr;
  double *dtds;
  int i;
  double *t;

  dtdr = ( double * ) malloc ( element_order * sizeof ( double ) );
  dtds = ( double * ) malloc ( element_order * sizeof ( double ) );
  t = ( double * ) malloc ( element_order * sizeof ( double ) );

  shape ( code, r, s, t, dtdr, dtds );
 
  *u = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *u = *u + ubase[i] * t[i];
  }

  *dudr = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *dudr = *dudr + ubase[i] * dtdr[i];
  }

  *duds = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *duds = *duds + ubase[i] * dtds[i];
  }

  free ( dtdr );
  free ( dtds );
  free ( t );
 
  return;
}
/******************************************************************************/

void interp_test ( char *code )

/******************************************************************************/
/*
  Purpose:

    INTERP_TEST tests the interpolation property of an element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
    "T3", "T4", "T6" and "T10".
*/
{
  double area;
  double dudr;
  double dudr_exact;
  double duds;
  double duds_exact;
  int element_order;
  int i;
  int node;
  double *node_r;
  double *node_s;
  double *node_u;
  double r;
  double r_factor;
  int *rexp;
  double s;
  double s_factor;
  int *sexp;
  int seed;
  int test;
  int test_num = 5;
  double u;
  double u_exact;

  if ( s_eqi ( code, "T4" ) )
  {
    printf ( "\n" );
    printf ( "INTERP_TEST - Warning!\n" );
    printf ( "  Skipping test for element \"T4\".\n" );
    return;
  }

  printf ( "\n" );
  printf ( "INTERP_TEST for element \"%s\".\n", code );

  element_order = order_code ( code );

  printf ( "  Number of nodes = %d\n", element_order );

  node_r = ( double * ) malloc ( element_order * sizeof ( double ) );
  node_s = ( double * ) malloc ( element_order * sizeof ( double ) );
  node_u = ( double * ) malloc ( element_order * sizeof ( double ) );
  rexp = ( int * ) malloc ( element_order * sizeof ( int ) );
  sexp = ( int * ) malloc ( element_order * sizeof ( int ) );
/*
  Get the coordinates of the reference nodes.
*/
  node_reference ( code, node_r, node_s, &area );
/*
  Get the monomial exponents for which the element is exact.
*/
  poly ( code, rexp, sexp );

  seed = 123456789;

  for ( i = 0; i < element_order; i++ )
  {
    printf ( "\n" );
    printf ( "  Interpolate R^%d * S^%d\n", rexp[i], sexp[i] );
    printf ( "\n" );
/*
  Evaluate R^REXP(I) * S^SEXP(I) at the nodes.  This is our data.
*/
    for ( node = 0; node < element_order; node++ )
    {
      r = node_r[node];
      s = node_s[node];
      if ( rexp[i] == 0 )
      {
        r_factor = 1.0;
      }
      else
      {
        r_factor = pow ( r, rexp[i] );
      }
      if ( sexp[i] == 0 )
      {
        s_factor = 1.0;
      }
      else
      {
        s_factor = pow ( s, sexp[i] );
      }
      node_u[node] = r_factor * s_factor;
      printf ( "  (R,S,U):       %12g  %12g  %12g\n", r, s, node_u[node] );
    }
/*
  Now pick random points in the element, and compute the interpolated
  value of R^REXP(*) * S^SEXP(I) there.  Mathematically, these
  values should be exact.
*/
    for ( test = 1; test <= test_num; test++ )
    {
      reference_sample ( code, &seed, &r, &s );

      printf ( "\n" );
      printf ( "  (R,S):            %12g  %12g\n", r, s );

      u_exact = r8_power ( r, rexp[i] ) * r8_power ( s, sexp[i] );

      dudr_exact = ( double ) ( rexp[i] ) 
        * r8_power ( r, rexp[i] - 1 ) * r8_power ( s, sexp[i] );

      duds_exact = r8_power ( r, rexp[i] ) * ( double ) ( sexp[i] ) 
        * r8_power ( s, sexp[i] - 1 );

      interp ( code, element_order, r, s, node_u, &u, &dudr, &duds );

      printf ( "  (U ,U* ,Error):   %12g  %12g  %12g\n", 
        u_exact, u, fabs ( u_exact - u ) );

      printf ( "  (Ur,Ur*,Error):   %12g  %12g  %12g\n",
        dudr_exact, dudr, fabs ( dudr_exact - dudr ) );

      printf ( "  (Us,Us*,Error):   %12g  %12g  %12g\n", 
        duds_exact, duds, fabs ( duds_exact - duds ) );
    }
  }

  free ( node_r );
  free ( node_s );
  free ( node_u ); 
  free ( rexp );
  free ( sexp );

  return;
}
/******************************************************************************/

void legendre_compute_dr ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_COMPUTE_DR: Gauss-Legendre quadrature by Davis-Rabinowitz method.
  
  Discussion:
    
    The integral:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 August 2007
  
  Author:
  
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order.
    ORDER must be greater than 0.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_COMPUTE_DR - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) )
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
/*
  Initial approximation H:
*/
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
/*
  Refine H using one step of Newton's method:
*/
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
/*
  Shift the data up.
*/
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
/*
  Reflect values for the negative abscissas.
*/
  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
    weight[i-1] = weight[order-i];
  }

  return;
}
/******************************************************************************/

void legendre_set ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
  
  Discussion:
    
    The integral:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    Quadrature rule:
  
      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
  
    The quadrature rule is exact for polynomials through degree 2*N-1.
  
    The abscissas are the zeroes of the Legendre polynomial P(ORDER)(X).
  
    Mathematica can compute the abscissas and weights of a Gauss-Legendre
    rule of order N for the interval [A,B] with P digits of precision
    by the commands:

    Needs["NumericalDifferentialEquationAnalysis`"]
    GaussianQuadratureWeights [ n, a, b, p ]
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    20 April 2010
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Vladimir Krylov,
    Approximate Calculation of Integrals,
    Dover, 2006,
    ISBN: 0486445798,
    LC: QA311.K713.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int N, the order.
    N must be between 1 and 33 or 63/64/65, 127/128/129, 
    255/256/257.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  if ( n == 1 )
  {
    x[0] = 0.000000000000000000000000000000;

    w[0] = 2.000000000000000000000000000000;
  }
  else if ( n == 2 )
  {
    x[0] = -0.577350269189625764509148780502;
    x[1] = 0.577350269189625764509148780502;

    w[0] = 1.000000000000000000000000000000;
    w[1] = 1.000000000000000000000000000000;
  }
  else if ( n == 3 )
  {
    x[0] = -0.774596669241483377035853079956;
    x[1] = 0.000000000000000000000000000000;
    x[2] = 0.774596669241483377035853079956;

    w[0] = 0.555555555555555555555555555556;
    w[1] = 0.888888888888888888888888888889;
    w[2] = 0.555555555555555555555555555556;
  }
  else if ( n == 4 )
  {
    x[0] = -0.861136311594052575223946488893;
    x[1] = -0.339981043584856264802665759103;
    x[2] = 0.339981043584856264802665759103;
    x[3] = 0.861136311594052575223946488893;

    w[0] = 0.347854845137453857373063949222;
    w[1] = 0.652145154862546142626936050778;
    w[2] = 0.652145154862546142626936050778;
    w[3] = 0.347854845137453857373063949222;
  }
  else if ( n == 5 )
  {
    x[0] = -0.906179845938663992797626878299;
    x[1] = -0.538469310105683091036314420700;
    x[2] = 0.000000000000000000000000000000;
    x[3] = 0.538469310105683091036314420700;
    x[4] = 0.906179845938663992797626878299;

    w[0] = 0.236926885056189087514264040720;
    w[1] = 0.478628670499366468041291514836;
    w[2] = 0.568888888888888888888888888889;
    w[3] = 0.478628670499366468041291514836;
    w[4] = 0.236926885056189087514264040720;
  }
  else if ( n == 6 )
  {
    x[0] = -0.932469514203152027812301554494;
    x[1] = -0.661209386466264513661399595020;
    x[2] = -0.238619186083196908630501721681;
    x[3] = 0.238619186083196908630501721681;
    x[4] = 0.661209386466264513661399595020;
    x[5] = 0.932469514203152027812301554494;

    w[0] = 0.171324492379170345040296142173;
    w[1] = 0.360761573048138607569833513838;
    w[2] = 0.467913934572691047389870343990;
    w[3] = 0.467913934572691047389870343990;
    w[4] = 0.360761573048138607569833513838;
    w[5] = 0.171324492379170345040296142173;
  }
  else if ( n == 7 )
  {
    x[0] = -0.949107912342758524526189684048;
    x[1] = -0.741531185599394439863864773281;
    x[2] = -0.405845151377397166906606412077;
    x[3] = 0.000000000000000000000000000000;
    x[4] = 0.405845151377397166906606412077;
    x[5] = 0.741531185599394439863864773281;
    x[6] = 0.949107912342758524526189684048;

    w[0] = 0.129484966168869693270611432679;
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489;
    w[3] = 0.417959183673469387755102040816;
    w[4] = 0.381830050505118944950369775489;
    w[5] = 0.279705391489276667901467771424;
    w[6] = 0.129484966168869693270611432679;
  }
  else if ( n == 8 )
  {
    x[0] = -0.960289856497536231683560868569;
    x[1] = -0.796666477413626739591553936476;
    x[2] = -0.525532409916328985817739049189;
    x[3] = -0.183434642495649804939476142360;
    x[4] = 0.183434642495649804939476142360;
    x[5] = 0.525532409916328985817739049189;
    x[6] = 0.796666477413626739591553936476;
    x[7] = 0.960289856497536231683560868569;

    w[0] = 0.101228536290376259152531354310;
    w[1] = 0.222381034453374470544355994426;
    w[2] = 0.313706645877887287337962201987;
    w[3] = 0.362683783378361982965150449277;
    w[4] = 0.362683783378361982965150449277;
    w[5] = 0.313706645877887287337962201987;
    w[6] = 0.222381034453374470544355994426;
    w[7] = 0.101228536290376259152531354310;
  }
  else if ( n == 9 )
  {
    x[0] = -0.968160239507626089835576203;
    x[1] = -0.836031107326635794299429788;
    x[2] = -0.613371432700590397308702039;
    x[3] = -0.324253423403808929038538015;
    x[4] = 0.000000000000000000000000000;
    x[5] = 0.324253423403808929038538015;
    x[6] = 0.613371432700590397308702039;
    x[7] = 0.836031107326635794299429788;
    x[8] = 0.968160239507626089835576203;

    w[0] = 0.081274388361574411971892158111;
    w[1] = 0.18064816069485740405847203124;
    w[2] = 0.26061069640293546231874286942;
    w[3] = 0.31234707704000284006863040658;
    w[4] = 0.33023935500125976316452506929;
    w[5] = 0.31234707704000284006863040658;
    w[6] = 0.26061069640293546231874286942;
    w[7] = 0.18064816069485740405847203124;
    w[8] = 0.081274388361574411971892158111;
  }
  else if ( n == 10 )
  {
    x[0] = -0.973906528517171720077964012;
    x[1] = -0.865063366688984510732096688;
    x[2] = -0.679409568299024406234327365;
    x[3] = -0.433395394129247190799265943;
    x[4] = -0.148874338981631210884826001;
    x[5] = 0.148874338981631210884826001;
    x[6] = 0.433395394129247190799265943;
    x[7] = 0.679409568299024406234327365;
    x[8] = 0.865063366688984510732096688;
    x[9] = 0.973906528517171720077964012;

    w[0] = 0.066671344308688137593568809893;
    w[1] = 0.14945134915058059314577633966;
    w[2] = 0.21908636251598204399553493423;
    w[3] = 0.26926671930999635509122692157;
    w[4] = 0.29552422471475287017389299465;
    w[5] = 0.29552422471475287017389299465;
    w[6] = 0.26926671930999635509122692157;
    w[7] = 0.21908636251598204399553493423;
    w[8] = 0.14945134915058059314577633966;
    w[9] = 0.066671344308688137593568809893;
  }
  else if ( n == 11 )
  {
    x[0] = -0.978228658146056992803938001;
    x[1] = -0.887062599768095299075157769;
    x[2] = -0.730152005574049324093416252;
    x[3] = -0.519096129206811815925725669;
    x[4] = -0.269543155952344972331531985;
    x[5] = 0.000000000000000000000000000;
    x[6] = 0.269543155952344972331531985;
    x[7] = 0.519096129206811815925725669;
    x[8] = 0.730152005574049324093416252;
    x[9] = 0.887062599768095299075157769;
    x[10] = 0.978228658146056992803938001;

    w[0] = 0.055668567116173666482753720443;
    w[1] = 0.12558036946490462463469429922;
    w[2] = 0.18629021092773425142609764143;
    w[3] = 0.23319376459199047991852370484;
    w[4] = 0.26280454451024666218068886989;
    w[5] = 0.27292508677790063071448352834;
    w[6] = 0.26280454451024666218068886989;
    w[7] = 0.23319376459199047991852370484;
    w[8] = 0.18629021092773425142609764143;
    w[9] = 0.12558036946490462463469429922;
    w[10] = 0.055668567116173666482753720443;
  }
  else if ( n == 12 )
  {
    x[0] = -0.981560634246719250690549090;
    x[1] = -0.904117256370474856678465866;
    x[2] = -0.769902674194304687036893833;
    x[3] = -0.587317954286617447296702419;
    x[4] = -0.367831498998180193752691537;
    x[5] = -0.125233408511468915472441369;
    x[6] = 0.125233408511468915472441369;
    x[7] = 0.367831498998180193752691537;
    x[8] = 0.587317954286617447296702419;
    x[9] = 0.769902674194304687036893833;
    x[10] = 0.904117256370474856678465866;
    x[11] = 0.981560634246719250690549090;

    w[0] = 0.047175336386511827194615961485;
    w[1] = 0.10693932599531843096025471819;
    w[2] = 0.16007832854334622633465252954;
    w[3] = 0.20316742672306592174906445581;
    w[4] = 0.23349253653835480876084989892;
    w[5] = 0.24914704581340278500056243604;
    w[6] = 0.24914704581340278500056243604;
    w[7] = 0.23349253653835480876084989892;
    w[8] = 0.20316742672306592174906445581;
    w[9] = 0.16007832854334622633465252954;
    w[10] = 0.10693932599531843096025471819;
    w[11] = 0.047175336386511827194615961485;
  }
  else if ( n == 13 )
  {
    x[0] = -0.984183054718588149472829449;
    x[1] = -0.917598399222977965206547837;
    x[2] = -0.801578090733309912794206490;
    x[3] = -0.642349339440340220643984607;
    x[4] = -0.448492751036446852877912852;
    x[5] = -0.230458315955134794065528121;
    x[6] = 0.000000000000000000000000000;
    x[7] = 0.230458315955134794065528121;
    x[8] = 0.448492751036446852877912852;
    x[9] = 0.642349339440340220643984607;
    x[10] = 0.80157809073330991279420649;
    x[11] = 0.91759839922297796520654784;
    x[12] = 0.98418305471858814947282945;

    w[0] = 0.040484004765315879520021592201;
    w[1] = 0.092121499837728447914421775954;
    w[2] = 0.13887351021978723846360177687;
    w[3] = 0.17814598076194573828004669200;
    w[4] = 0.20781604753688850231252321931;
    w[5] = 0.22628318026289723841209018604;
    w[6] = 0.23255155323087391019458951527;
    w[7] = 0.22628318026289723841209018604;
    w[8] = 0.20781604753688850231252321931;
    w[9] = 0.17814598076194573828004669200;
    w[10] = 0.13887351021978723846360177687;
    w[11] = 0.092121499837728447914421775954;
    w[12] = 0.040484004765315879520021592201;
  }
  else if ( n == 14 )
  {
    x[0] = -0.986283808696812338841597267;
    x[1] = -0.928434883663573517336391139;
    x[2] = -0.827201315069764993189794743;
    x[3] = -0.687292904811685470148019803;
    x[4] = -0.515248636358154091965290719;
    x[5] = -0.319112368927889760435671824;
    x[6] = -0.108054948707343662066244650;
    x[7] = 0.108054948707343662066244650;
    x[8] = 0.31911236892788976043567182;
    x[9] = 0.51524863635815409196529072;
    x[10] = 0.68729290481168547014801980;
    x[11] = 0.82720131506976499318979474;
    x[12] = 0.92843488366357351733639114;
    x[13] = 0.98628380869681233884159727;

    w[0] = 0.035119460331751863031832876138;
    w[1] = 0.08015808715976020980563327706;
    w[2] = 0.12151857068790318468941480907;
    w[3] = 0.15720316715819353456960193862;
    w[4] = 0.18553839747793781374171659013;
    w[5] = 0.20519846372129560396592406566;
    w[6] = 0.21526385346315779019587644332;
    w[7] = 0.21526385346315779019587644332;
    w[8] = 0.20519846372129560396592406566;
    w[9] = 0.18553839747793781374171659013;
    w[10] = 0.15720316715819353456960193862;
    w[11] = 0.12151857068790318468941480907;
    w[12] = 0.08015808715976020980563327706;
    w[13] = 0.035119460331751863031832876138;
  }
  else if ( n == 15 )
  {
    x[0] = -0.987992518020485428489565719;
    x[1] = -0.937273392400705904307758948;
    x[2] = -0.848206583410427216200648321;
    x[3] = -0.724417731360170047416186055;
    x[4] = -0.570972172608538847537226737;
    x[5] = -0.394151347077563369897207371;
    x[6] = -0.201194093997434522300628303;
    x[7] = 0.00000000000000000000000000;
    x[8] = 0.20119409399743452230062830;
    x[9] = 0.39415134707756336989720737;
    x[10] = 0.57097217260853884753722674;
    x[11] = 0.72441773136017004741618605;
    x[12] = 0.84820658341042721620064832;
    x[13] = 0.93727339240070590430775895;
    x[14] = 0.98799251802048542848956572;

    w[0] = 0.030753241996117268354628393577;
    w[1] = 0.070366047488108124709267416451;
    w[2] = 0.107159220467171935011869546686;
    w[3] = 0.13957067792615431444780479451;
    w[4] = 0.16626920581699393355320086048;
    w[5] = 0.18616100001556221102680056187;
    w[6] = 0.19843148532711157645611832644;
    w[7] = 0.20257824192556127288062019997;
    w[8] = 0.19843148532711157645611832644;
    w[9] = 0.18616100001556221102680056187;
    w[10] = 0.16626920581699393355320086048;
    w[11] = 0.13957067792615431444780479451;
    w[12] = 0.107159220467171935011869546686;
    w[13] = 0.070366047488108124709267416451;
    w[14] = 0.030753241996117268354628393577;
  }
  else if ( n == 16 )
  {
    x[0] = -0.989400934991649932596154173;
    x[1] = -0.944575023073232576077988416;
    x[2] = -0.865631202387831743880467898;
    x[3] = -0.755404408355003033895101195;
    x[4] = -0.617876244402643748446671764;
    x[5] = -0.458016777657227386342419443;
    x[6] = -0.281603550779258913230460501;
    x[7] = -0.09501250983763744018531934;
    x[8] = 0.09501250983763744018531934;
    x[9] = 0.28160355077925891323046050;
    x[10] = 0.45801677765722738634241944;
    x[11] = 0.61787624440264374844667176;
    x[12] = 0.75540440835500303389510119;
    x[13] = 0.86563120238783174388046790;
    x[14] = 0.94457502307323257607798842;
    x[15] = 0.98940093499164993259615417;

    w[0] = 0.027152459411754094851780572456;
    w[1] = 0.062253523938647892862843836994;
    w[2] = 0.09515851168249278480992510760;
    w[3] = 0.12462897125553387205247628219;
    w[4] = 0.14959598881657673208150173055;
    w[5] = 0.16915651939500253818931207903;
    w[6] = 0.18260341504492358886676366797;
    w[7] = 0.18945061045506849628539672321;
    w[8] = 0.18945061045506849628539672321;
    w[9] = 0.18260341504492358886676366797;
    w[10] = 0.16915651939500253818931207903;
    w[11] = 0.14959598881657673208150173055;
    w[12] = 0.12462897125553387205247628219;
    w[13] = 0.09515851168249278480992510760;
    w[14] = 0.062253523938647892862843836994;
    w[15] = 0.027152459411754094851780572456;
  }
  else if ( n == 17 )
  {
    x[0] = -0.990575475314417335675434020;
    x[1] = -0.950675521768767761222716958;
    x[2] = -0.880239153726985902122955694;
    x[3] = -0.781514003896801406925230056;
    x[4] = -0.657671159216690765850302217;
    x[5] = -0.512690537086476967886246569;
    x[6] = -0.35123176345387631529718552;
    x[7] = -0.17848418149584785585067749;
    x[8] = 0.00000000000000000000000000;
    x[9] = 0.17848418149584785585067749;
    x[10] = 0.35123176345387631529718552;
    x[11] = 0.51269053708647696788624657;
    x[12] = 0.65767115921669076585030222;
    x[13] = 0.78151400389680140692523006;
    x[14] = 0.88023915372698590212295569;
    x[15] = 0.95067552176876776122271696;
    x[16] = 0.99057547531441733567543402;

    w[0] = 0.024148302868547931960110026288;
    w[1] = 0.055459529373987201129440165359;
    w[2] = 0.085036148317179180883535370191;
    w[3] = 0.111883847193403971094788385626;
    w[4] = 0.13513636846852547328631998170;
    w[5] = 0.15404576107681028808143159480;
    w[6] = 0.16800410215645004450997066379;
    w[7] = 0.17656270536699264632527099011;
    w[8] = 0.17944647035620652545826564426;
    w[9] = 0.17656270536699264632527099011;
    w[10] = 0.16800410215645004450997066379;
    w[11] = 0.15404576107681028808143159480;
    w[12] = 0.13513636846852547328631998170;
    w[13] = 0.111883847193403971094788385626;
    w[14] = 0.085036148317179180883535370191;
    w[15] = 0.055459529373987201129440165359;
    w[16] = 0.024148302868547931960110026288;
  }
  else if ( n == 18 )
  {
    x[0] = -0.991565168420930946730016005;
    x[1] = -0.955823949571397755181195893;
    x[2] = -0.892602466497555739206060591;
    x[3] = -0.803704958972523115682417455;
    x[4] = -0.691687043060353207874891081;
    x[5] = -0.55977083107394753460787155;
    x[6] = -0.41175116146284264603593179;
    x[7] = -0.25188622569150550958897285;
    x[8] = -0.08477501304173530124226185;
    x[9] = 0.08477501304173530124226185;
    x[10] = 0.25188622569150550958897285;
    x[11] = 0.41175116146284264603593179;
    x[12] = 0.55977083107394753460787155;
    x[13] = 0.69168704306035320787489108;
    x[14] = 0.80370495897252311568241746;
    x[15] = 0.89260246649755573920606059;
    x[16] = 0.95582394957139775518119589;
    x[17] = 0.99156516842093094673001600;

    w[0] = 0.021616013526483310313342710266;
    w[1] = 0.049714548894969796453334946203;
    w[2] = 0.07642573025488905652912967762;
    w[3] = 0.10094204410628716556281398492;
    w[4] = 0.12255520671147846018451912680;
    w[5] = 0.14064291467065065120473130375;
    w[6] = 0.15468467512626524492541800384;
    w[7] = 0.16427648374583272298605377647;
    w[8] = 0.16914238296314359184065647013;
    w[9] = 0.16914238296314359184065647013;
    w[10] = 0.16427648374583272298605377647;
    w[11] = 0.15468467512626524492541800384;
    w[12] = 0.14064291467065065120473130375;
    w[13] = 0.12255520671147846018451912680;
    w[14] = 0.10094204410628716556281398492;
    w[15] = 0.07642573025488905652912967762;
    w[16] = 0.049714548894969796453334946203;
    w[17] = 0.021616013526483310313342710266;
  }
  else if ( n == 19 )
  {
    x[0] = -0.992406843843584403189017670;
    x[1] = -0.960208152134830030852778841;
    x[2] = -0.903155903614817901642660929;
    x[3] = -0.822714656537142824978922487;
    x[4] = -0.72096617733522937861709586;
    x[5] = -0.60054530466168102346963816;
    x[6] = -0.46457074137596094571726715;
    x[7] = -0.31656409996362983199011733;
    x[8] = -0.16035864564022537586809612;
    x[9] = 0.00000000000000000000000000;
    x[10] = 0.16035864564022537586809612;
    x[11] = 0.31656409996362983199011733;
    x[12] = 0.46457074137596094571726715;
    x[13] = 0.60054530466168102346963816;
    x[14] = 0.72096617733522937861709586;
    x[15] = 0.82271465653714282497892249;
    x[16] = 0.90315590361481790164266093;
    x[17] = 0.96020815213483003085277884;
    x[18] = 0.99240684384358440318901767;

    w[0] = 0.019461788229726477036312041464;
    w[1] = 0.044814226765699600332838157402;
    w[2] = 0.069044542737641226580708258006;
    w[3] = 0.091490021622449999464462094124;
    w[4] = 0.111566645547333994716023901682;
    w[5] = 0.12875396253933622767551578486;
    w[6] = 0.14260670217360661177574610944;
    w[7] = 0.15276604206585966677885540090;
    w[8] = 0.15896884339395434764995643946;
    w[9] = 0.16105444984878369597916362532;
    w[10] = 0.15896884339395434764995643946;
    w[11] = 0.15276604206585966677885540090;
    w[12] = 0.14260670217360661177574610944;
    w[13] = 0.12875396253933622767551578486;
    w[14] = 0.111566645547333994716023901682;
    w[15] = 0.091490021622449999464462094124;
    w[16] = 0.069044542737641226580708258006;
    w[17] = 0.044814226765699600332838157402;
    w[18] = 0.019461788229726477036312041464;
  }
  else if ( n == 20 )
  {
    x[0] = -0.993128599185094924786122388;
    x[1] = -0.963971927277913791267666131;
    x[2] = -0.912234428251325905867752441;
    x[3] = -0.83911697182221882339452906;
    x[4] = -0.74633190646015079261430507;
    x[5] = -0.63605368072651502545283670;
    x[6] = -0.51086700195082709800436405;
    x[7] = -0.37370608871541956067254818;
    x[8] = -0.22778585114164507808049620;
    x[9] = -0.07652652113349733375464041;
    x[10] = 0.07652652113349733375464041;
    x[11] = 0.22778585114164507808049620;
    x[12] = 0.37370608871541956067254818;
    x[13] = 0.51086700195082709800436405;
    x[14] = 0.63605368072651502545283670;
    x[15] = 0.74633190646015079261430507;
    x[16] = 0.83911697182221882339452906;
    x[17] = 0.91223442825132590586775244;
    x[18] = 0.96397192727791379126766613;
    x[19] = 0.99312859918509492478612239;

    w[0] = 0.017614007139152118311861962352;
    w[1] = 0.040601429800386941331039952275;
    w[2] = 0.062672048334109063569506535187;
    w[3] = 0.08327674157670474872475814322;
    w[4] = 0.10193011981724043503675013548;
    w[5] = 0.11819453196151841731237737771;
    w[6] = 0.13168863844917662689849449975;
    w[7] = 0.14209610931838205132929832507;
    w[8] = 0.14917298647260374678782873700;
    w[9] = 0.15275338713072585069808433195;
    w[10] = 0.15275338713072585069808433195;
    w[11] = 0.14917298647260374678782873700;
    w[12] = 0.14209610931838205132929832507;
    w[13] = 0.13168863844917662689849449975;
    w[14] = 0.11819453196151841731237737771;
    w[15] = 0.10193011981724043503675013548;
    w[16] = 0.08327674157670474872475814322;
    w[17] = 0.062672048334109063569506535187;
    w[18] = 0.040601429800386941331039952275;
    w[19] = 0.017614007139152118311861962352;
  }
  else if ( n == 21 )
  {
    x[ 0] =  -0.99375217062038950026024204;
    x[ 1] =  -0.96722683856630629431662221;
    x[ 2] =  -0.92009933415040082879018713;
    x[ 3] =  -0.85336336458331728364725064;
    x[ 4] =  -0.76843996347567790861587785;
    x[ 5] =  -0.66713880419741231930596667;
    x[ 6] =  -0.55161883588721980705901880;
    x[ 7] =  -0.42434212020743878357366889;
    x[ 8] =  -0.28802131680240109660079252;
    x[9] =  -0.14556185416089509093703098;
    x[10] =   0.00000000000000000000000000;
    x[11] =  +0.14556185416089509093703098;
    x[12] =  +0.28802131680240109660079252;
    x[13] =  +0.42434212020743878357366889;
    x[14] =  +0.55161883588721980705901880;
    x[15] =  +0.66713880419741231930596667;
    x[16] =  +0.76843996347567790861587785;
    x[17] =  +0.85336336458331728364725064;
    x[18] =  +0.92009933415040082879018713;
    x[19] =  +0.96722683856630629431662221;
    x[20] =  +0.99375217062038950026024204;

    w[ 0] =   0.016017228257774333324224616858;
    w[ 1] =   0.036953789770852493799950668299; 
    w[ 2] =   0.057134425426857208283635826472;
    w[ 3] =   0.076100113628379302017051653300;
    w[ 4] =   0.093444423456033861553289741114;
    w[ 5] =   0.108797299167148377663474578070;
    w[ 6] =   0.12183141605372853419536717713;
    w[ 7] =   0.13226893863333746178105257450;
    w[ 8] =   0.13988739479107315472213342387;
    w[9] =   0.14452440398997005906382716655;
    w[10] =   0.14608113364969042719198514768;
    w[11] =   0.14452440398997005906382716655; 
    w[12] =   0.13988739479107315472213342387; 
    w[13] =   0.13226893863333746178105257450;
    w[14] =   0.12183141605372853419536717713;
    w[15] =   0.108797299167148377663474578070;
    w[16] =   0.093444423456033861553289741114;
    w[17] =   0.076100113628379302017051653300;
    w[18] =   0.057134425426857208283635826472;
    w[19] =   0.036953789770852493799950668299;
    w[20] =   0.016017228257774333324224616858;
  }
  else if ( n == 22 )
  {
    x[0] = -0.99429458548239929207303142;
    x[1] = -0.97006049783542872712395099;
    x[2] = -0.92695677218717400052069294;
    x[3] = -0.86581257772030013653642564;
    x[4] = -0.78781680597920816200427796;
    x[5] = -0.69448726318668278005068984;
    x[6] = -0.58764040350691159295887693;
    x[7] = -0.46935583798675702640633071;
    x[8] = -0.34193582089208422515814742;
    x[9] = -0.20786042668822128547884653;
    x[10] = -0.06973927331972222121384180;
    x[11] = 0.06973927331972222121384180;
    x[12] = 0.20786042668822128547884653;
    x[13] = 0.34193582089208422515814742;
    x[14] = 0.46935583798675702640633071;
    x[15] = 0.58764040350691159295887693;
    x[16] = 0.69448726318668278005068984;
    x[17] = 0.78781680597920816200427796;
    x[18] = 0.86581257772030013653642564;
    x[19] = 0.92695677218717400052069294;
    x[20] = 0.97006049783542872712395099;
    x[21] = 0.99429458548239929207303142;
 
    w[0] = 0.014627995298272200684991098047;
    w[1] = 0.033774901584814154793302246866;
    w[2] = 0.052293335152683285940312051273;
    w[3] = 0.06979646842452048809496141893;
    w[4] = 0.08594160621706772741444368137;
    w[5] = 0.10041414444288096493207883783;
    w[6] = 0.11293229608053921839340060742;
    w[7] = 0.12325237681051242428556098615;
    w[8] = 0.13117350478706237073296499253;
    w[9] = 0.13654149834601517135257383123;
    w[10] = 0.13925187285563199337541024834;
    w[11] = 0.13925187285563199337541024834;
    w[12] = 0.13654149834601517135257383123;
    w[13] = 0.13117350478706237073296499253;
    w[14] = 0.12325237681051242428556098615;
    w[15] = 0.11293229608053921839340060742;
    w[16] = 0.10041414444288096493207883783;
    w[17] = 0.08594160621706772741444368137;
    w[18] = 0.06979646842452048809496141893;
    w[19] = 0.052293335152683285940312051273;
    w[20] = 0.033774901584814154793302246866;
    w[21] = 0.014627995298272200684991098047;
  }
  else if ( n == 23 )
  {
    x[0] = -0.99476933499755212352392572;
    x[1] = -0.97254247121811523195602408;
    x[2] = -0.93297108682601610234919699;
    x[3] = -0.87675235827044166737815689;
    x[4] = -0.80488840161883989215111841;
    x[5] = -0.71866136313195019446162448;
    x[6] = -0.61960987576364615638509731;
    x[7] = -0.50950147784600754968979305;
    x[8] = -0.39030103803029083142148887;
    x[9] = -0.26413568097034493053386954;
    x[10] = -0.13325682429846611093174268;
    x[11] = 0.00000000000000000000000000;
    x[12] = 0.13325682429846611093174268;
    x[13] = 0.26413568097034493053386954;
    x[14] = 0.39030103803029083142148887;
    x[15] = 0.50950147784600754968979305;
    x[16] = 0.61960987576364615638509731;
    x[17] = 0.71866136313195019446162448;
    x[18] = 0.80488840161883989215111841;
    x[19] = 0.87675235827044166737815689;
    x[20] = 0.93297108682601610234919699;
    x[21] = 0.97254247121811523195602408;
    x[22] = 0.99476933499755212352392572;

    w[0] = 0.013411859487141772081309493459;
    w[1] = 0.030988005856979444310694219642;
    w[2] = 0.048037671731084668571641071632;
    w[3] = 0.064232421408525852127169615159;
    w[4] = 0.079281411776718954922892524742;
    w[5] = 0.092915766060035147477018617370;
    w[6] = 0.104892091464541410074086185015;
    w[7] = 0.11499664022241136494164351293;
    w[8] = 0.12304908430672953046757840067;
    w[9] = 0.12890572218808214997859533940;
    w[10] = 0.13246203940469661737164246470;
    w[11] = 0.13365457218610617535145711055;
    w[12] = 0.13246203940469661737164246470;
    w[13] = 0.12890572218808214997859533940;
    w[14] = 0.12304908430672953046757840067;
    w[15] = 0.11499664022241136494164351293;
    w[16] = 0.104892091464541410074086185015;
    w[17] = 0.092915766060035147477018617370;
    w[18] = 0.079281411776718954922892524742;
    w[19] = 0.064232421408525852127169615159;
    w[20] = 0.048037671731084668571641071632;
    w[21] = 0.030988005856979444310694219642;
    w[22] = 0.013411859487141772081309493459;
  }
  else if ( n == 24 )
  {
    x[0] = -0.99518721999702136017999741;
    x[1] = -0.97472855597130949819839199;
    x[2] = -0.93827455200273275852364900;
    x[3] = -0.88641552700440103421315434;
    x[4] = -0.82000198597390292195394987;
    x[5] = -0.74012419157855436424382810;
    x[6] = -0.64809365193697556925249579;
    x[7] = -0.54542147138883953565837562;
    x[8] = -0.43379350762604513848708423;
    x[9] = -0.31504267969616337438679329;
    x[10] = -0.19111886747361630915863982;
    x[11] = -0.06405689286260562608504308;
    x[12] = 0.06405689286260562608504308;
    x[13] = 0.19111886747361630915863982;
    x[14] = 0.31504267969616337438679329;
    x[15] = 0.43379350762604513848708423;
    x[16] = 0.54542147138883953565837562;
    x[17] = 0.64809365193697556925249579;
    x[18] = 0.74012419157855436424382810;
    x[19] = 0.82000198597390292195394987;
    x[20] = 0.88641552700440103421315434;
    x[21] = 0.93827455200273275852364900;
    x[22] = 0.97472855597130949819839199;
    x[23] = 0.99518721999702136017999741;

    w[0] = 0.012341229799987199546805667070;
    w[1] = 0.028531388628933663181307815952;
    w[2] = 0.044277438817419806168602748211;
    w[3] = 0.059298584915436780746367758500;
    w[4] = 0.07334648141108030573403361525;
    w[5] = 0.08619016153195327591718520298;
    w[6] = 0.09761865210411388826988066446;
    w[7] = 0.10744427011596563478257734245;
    w[8] = 0.11550566805372560135334448391;
    w[9] = 0.12167047292780339120446315348;
    w[10] = 0.12583745634682829612137538251;
    w[11] = 0.12793819534675215697405616522;
    w[12] = 0.12793819534675215697405616522;
    w[13] = 0.12583745634682829612137538251;
    w[14] = 0.12167047292780339120446315348;
    w[15] = 0.11550566805372560135334448391;
    w[16] = 0.10744427011596563478257734245;
    w[17] = 0.09761865210411388826988066446;
    w[18] = 0.08619016153195327591718520298;
    w[19] = 0.07334648141108030573403361525;
    w[20] = 0.059298584915436780746367758500;
    w[21] = 0.044277438817419806168602748211;
    w[22] = 0.028531388628933663181307815952;
    w[23] = 0.012341229799987199546805667070;
  }
  else if ( n == 25 )
  {
    x[0] = -0.99555696979049809790878495;
    x[1] = -0.97666392145951751149831539;
    x[2] = -0.94297457122897433941401117;
    x[3] = -0.89499199787827536885104201;
    x[4] = -0.83344262876083400142102111;
    x[5] = -0.75925926303735763057728287;
    x[6] = -0.67356636847346836448512063;
    x[7] = -0.57766293024122296772368984;
    x[8] = -0.47300273144571496052218212;
    x[9] = -0.36117230580938783773582173;
    x[10] = -0.24386688372098843204519036;
    x[11] = -0.12286469261071039638735982;
    x[12] = 0.00000000000000000000000000;
    x[13] = 0.12286469261071039638735982;
    x[14] = 0.24386688372098843204519036;
    x[15] = 0.36117230580938783773582173;
    x[16] = 0.47300273144571496052218212;
    x[17] = 0.57766293024122296772368984;
    x[18] = 0.67356636847346836448512063;
    x[19] = 0.75925926303735763057728287;
    x[20] = 0.83344262876083400142102111;
    x[21] = 0.89499199787827536885104201;
    x[22] = 0.94297457122897433941401117;
    x[23] = 0.97666392145951751149831539;
    x[24] = 0.99555696979049809790878495;

    w[0] = 0.0113937985010262879479029641132;
    w[1] = 0.026354986615032137261901815295;
    w[2] = 0.040939156701306312655623487712;
    w[3] = 0.054904695975835191925936891541;
    w[4] = 0.068038333812356917207187185657;
    w[5] = 0.080140700335001018013234959669;
    w[6] = 0.091028261982963649811497220703;
    w[7] = 0.100535949067050644202206890393;
    w[8] = 0.108519624474263653116093957050;
    w[9] = 0.11485825914571164833932554587;
    w[10] = 0.11945576353578477222817812651;
    w[11] = 0.12224244299031004168895951895;
    w[12] = 0.12317605372671545120390287308;
    w[13] = 0.12224244299031004168895951895;
    w[14] = 0.11945576353578477222817812651;
    w[15] = 0.11485825914571164833932554587;
    w[16] = 0.108519624474263653116093957050;
    w[17] = 0.100535949067050644202206890393;
    w[18] = 0.091028261982963649811497220703;
    w[19] = 0.080140700335001018013234959669;
    w[20] = 0.068038333812356917207187185657;
    w[21] = 0.054904695975835191925936891541;
    w[22] = 0.040939156701306312655623487712;
    w[23] = 0.026354986615032137261901815295;
    w[24] = 0.0113937985010262879479029641132;
  }
  else if ( n == 26 )
  {
    x[0] = -0.99588570114561692900321696;
    x[1] = -0.97838544595647099110058035;
    x[2] = -0.94715906666171425013591528;
    x[3] = -0.90263786198430707421766560;
    x[4] = -0.84544594278849801879750706;
    x[5] = -0.77638594882067885619296725;
    x[6] = -0.69642726041995726486381391;
    x[7] = -0.60669229301761806323197875;
    x[8] = -0.50844071482450571769570306;
    x[9] = -0.40305175512348630648107738;
    x[10] = -0.29200483948595689514283538;
    x[11] = -0.17685882035689018396905775;
    x[12] = -0.05923009342931320709371858;
    x[13] = 0.05923009342931320709371858;
    x[14] = 0.17685882035689018396905775;
    x[15] = 0.29200483948595689514283538;
    x[16] = 0.40305175512348630648107738;
    x[17] = 0.50844071482450571769570306;
    x[18] = 0.60669229301761806323197875;
    x[19] = 0.69642726041995726486381391;
    x[20] = 0.77638594882067885619296725;
    x[21] = 0.84544594278849801879750706;
    x[22] = 0.90263786198430707421766560;
    x[23] = 0.94715906666171425013591528;
    x[24] = 0.97838544595647099110058035;
    x[25] = 0.99588570114561692900321696;

    w[0] = 0.010551372617343007155651187685;
    w[1] = 0.024417851092631908789615827520;
    w[2] = 0.037962383294362763950303141249;
    w[3] = 0.050975825297147811998319900724;
    w[4] = 0.063274046329574835539453689907;
    w[5] = 0.07468414976565974588707579610;
    w[6] = 0.08504589431348523921044776508;
    w[7] = 0.09421380035591414846366488307;
    w[8] = 0.10205916109442542323841407025;
    w[9] = 0.10847184052857659065657942673;
    w[10] = 0.11336181654631966654944071844;
    w[11] = 0.11666044348529658204466250754;
    w[12] = 0.11832141527926227651637108570;
    w[13] = 0.11832141527926227651637108570;
    w[14] = 0.11666044348529658204466250754;
    w[15] = 0.11336181654631966654944071844;
    w[16] = 0.10847184052857659065657942673;
    w[17] = 0.10205916109442542323841407025;
    w[18] = 0.09421380035591414846366488307;
    w[19] = 0.08504589431348523921044776508;
    w[20] = 0.07468414976565974588707579610;
    w[21] = 0.063274046329574835539453689907;
    w[22] = 0.050975825297147811998319900724;
    w[23] = 0.037962383294362763950303141249;
    w[24] = 0.024417851092631908789615827520;
    w[25] = 0.010551372617343007155651187685;
  }
  else if ( n == 27 )
  {
    x[0] = -0.99617926288898856693888721;
    x[1] = -0.97992347596150122285587336;
    x[2] = -0.95090055781470500685190803;
    x[3] = -0.90948232067749110430064502;
    x[4] = -0.85620790801829449030273722;
    x[5] = -0.79177163907050822714439734;
    x[6] = -0.71701347373942369929481621;
    x[7] = -0.63290797194649514092773464;
    x[8] = -0.54055156457945689490030094;
    x[9] = -0.44114825175002688058597416;
    x[10] = -0.33599390363850889973031903;
    x[11] = -0.22645936543953685885723911;
    x[12] = -0.11397258560952996693289498;
    x[13] = 0.00000000000000000000000000;
    x[14] = 0.11397258560952996693289498;
    x[15] = 0.22645936543953685885723911;
    x[16] = 0.33599390363850889973031903;
    x[17] = 0.44114825175002688058597416;
    x[18] = 0.54055156457945689490030094;
    x[19] = 0.63290797194649514092773464;
    x[20] = 0.71701347373942369929481621;
    x[21] = 0.79177163907050822714439734;
    x[22] = 0.85620790801829449030273722;
    x[23] = 0.90948232067749110430064502;
    x[24] = 0.95090055781470500685190803;
    x[25] = 0.97992347596150122285587336;
    x[26] = 0.99617926288898856693888721;

    w[0] = 0.0097989960512943602611500550912;
    w[1] = 0.022686231596180623196034206447;
    w[2] = 0.035297053757419711022578289305;
    w[3] = 0.047449412520615062704096710114;
    w[4] = 0.058983536859833599110300833720;
    w[5] = 0.069748823766245592984322888357;
    w[6] = 0.079604867773057771263074959010;
    w[7] = 0.088423158543756950194322802854;
    w[8] = 0.096088727370028507565652646558;
    w[9] = 0.102501637817745798671247711533;
    w[10] = 0.107578285788533187212162984427;
    w[11] = 0.111252488356845192672163096043;
    w[12] = 0.113476346108965148620369948092;
    w[13] = 0.11422086737895698904504573690;
    w[14] = 0.113476346108965148620369948092;
    w[15] = 0.111252488356845192672163096043;
    w[16] = 0.107578285788533187212162984427;
    w[17] = 0.102501637817745798671247711533;
    w[18] = 0.096088727370028507565652646558;
    w[19] = 0.088423158543756950194322802854;
    w[20] = 0.079604867773057771263074959010;
    w[21] = 0.069748823766245592984322888357;
    w[22] = 0.058983536859833599110300833720;
    w[23] = 0.047449412520615062704096710114;
    w[24] = 0.035297053757419711022578289305;
    w[25] = 0.022686231596180623196034206447;
    w[26] = 0.0097989960512943602611500550912;
  }
  else if ( n == 28 )
  {
    x[0] = -0.99644249757395444995043639;
    x[1] = -0.98130316537087275369455995;
    x[2] = -0.95425928062893819725410184;
    x[3] = -0.91563302639213207386968942;
    x[4] = -0.86589252257439504894225457;
    x[5] = -0.80564137091717917144788596;
    x[6] = -0.73561087801363177202814451;
    x[7] = -0.65665109403886496121989818;
    x[8] = -0.56972047181140171930800328;
    x[9] = -0.47587422495511826103441185;
    x[10] = -0.37625151608907871022135721;
    x[11] = -0.27206162763517807767682636;
    x[12] = -0.16456928213338077128147178;
    x[13] = -0.05507928988403427042651653;
    x[14] = 0.05507928988403427042651653;
    x[15] = 0.16456928213338077128147178;
    x[16] = 0.27206162763517807767682636;
    x[17] = 0.37625151608907871022135721;
    x[18] = 0.47587422495511826103441185;
    x[19] = 0.56972047181140171930800328;
    x[20] = 0.65665109403886496121989818;
    x[21] = 0.73561087801363177202814451;
    x[22] = 0.80564137091717917144788596;
    x[23] = 0.86589252257439504894225457;
    x[24] = 0.91563302639213207386968942;
    x[25] = 0.95425928062893819725410184;
    x[26] = 0.98130316537087275369455995;
    x[27] = 0.99644249757395444995043639;

    w[0] = 0.009124282593094517738816153923;
    w[1] = 0.021132112592771259751500380993;
    w[2] = 0.032901427782304379977630819171;
    w[3] = 0.044272934759004227839587877653;
    w[4] = 0.055107345675716745431482918227;
    w[5] = 0.06527292396699959579339756678;
    w[6] = 0.07464621423456877902393188717;
    w[7] = 0.08311341722890121839039649824;
    w[8] = 0.09057174439303284094218603134;
    w[9] = 0.09693065799792991585048900610;
    w[10] = 0.10211296757806076981421663851;
    w[11] = 0.10605576592284641791041643700;
    w[12] = 0.10871119225829413525357151930;
    w[13] = 0.11004701301647519628237626560;
    w[14] = 0.11004701301647519628237626560;
    w[15] = 0.10871119225829413525357151930;
    w[16] = 0.10605576592284641791041643700;
    w[17] = 0.10211296757806076981421663851;
    w[18] = 0.09693065799792991585048900610;
    w[19] = 0.09057174439303284094218603134;
    w[20] = 0.08311341722890121839039649824;
    w[21] = 0.07464621423456877902393188717;
    w[22] = 0.06527292396699959579339756678;
    w[23] = 0.055107345675716745431482918227;
    w[24] = 0.044272934759004227839587877653;
    w[25] = 0.032901427782304379977630819171;
    w[26] = 0.021132112592771259751500380993;
    w[27] = 0.009124282593094517738816153923;
  }
  else if ( n == 29 )
  {
    x[0] = -0.99667944226059658616319153;
    x[1] = -0.98254550526141317487092602;
    x[2] = -0.95728559577808772579820804;
    x[3] = -0.92118023295305878509375344;
    x[4] = -0.87463780492010279041779342;
    x[5] = -0.81818548761525244498957221;
    x[6] = -0.75246285173447713391261008;
    x[7] = -0.67821453760268651515618501;
    x[8] = -0.59628179713822782037958621;
    x[9] = -0.50759295512422764210262792;
    x[10] = -0.41315288817400866389070659;
    x[11] = -0.31403163786763993494819592;
    x[12] = -0.21135228616600107450637573;
    x[13] = -0.10627823013267923017098239;
    x[14] = 0.00000000000000000000000000;
    x[15] = 0.10627823013267923017098239;
    x[16] = 0.21135228616600107450637573;
    x[17] = 0.31403163786763993494819592;
    x[18] = 0.41315288817400866389070659;
    x[19] = 0.50759295512422764210262792;
    x[20] = 0.59628179713822782037958621;
    x[21] = 0.67821453760268651515618501;
    x[22] = 0.75246285173447713391261008;
    x[23] = 0.81818548761525244498957221;
    x[24] = 0.87463780492010279041779342;
    x[25] = 0.92118023295305878509375344;
    x[26] = 0.95728559577808772579820804;
    x[27] = 0.98254550526141317487092602;
    x[28] = 0.99667944226059658616319153;

    w[0] = 0.0085169038787464096542638133022;
    w[1] = 0.019732085056122705983859801640;
    w[2] = 0.030740492202093622644408525375;
    w[3] = 0.041402062518682836104830010114;
    w[4] = 0.051594826902497923912594381180;
    w[5] = 0.061203090657079138542109848024;
    w[6] = 0.070117933255051278569581486949;
    w[7] = 0.078238327135763783828144888660;
    w[8] = 0.085472257366172527545344849297;
    w[9] = 0.091737757139258763347966411077;
    w[10] = 0.096963834094408606301900074883;
    w[11] = 0.101091273759914966121820546907;
    w[12] = 0.104073310077729373913328471285;
    w[13] = 0.105876155097320941406591327852;
    w[14] = 0.10647938171831424424651112691;
    w[15] = 0.105876155097320941406591327852;
    w[16] = 0.104073310077729373913328471285;
    w[17] = 0.101091273759914966121820546907;
    w[18] = 0.096963834094408606301900074883;
    w[19] = 0.091737757139258763347966411077;
    w[20] = 0.085472257366172527545344849297;
    w[21] = 0.078238327135763783828144888660;
    w[22] = 0.070117933255051278569581486949;
    w[23] = 0.061203090657079138542109848024;
    w[24] = 0.051594826902497923912594381180;
    w[25] = 0.041402062518682836104830010114;
    w[26] = 0.030740492202093622644408525375;
    w[27] = 0.019732085056122705983859801640;
    w[28] = 0.0085169038787464096542638133022;
  }
  else if ( n == 30 )
  {
    x[0] = -0.99689348407464954027163005;
    x[1] = -0.98366812327974720997003258;
    x[2] = -0.96002186496830751221687103;
    x[3] = -0.92620004742927432587932428;
    x[4] = -0.88256053579205268154311646;
    x[5] = -0.82956576238276839744289812;
    x[6] = -0.76777743210482619491797734;
    x[7] = -0.69785049479331579693229239;
    x[8] = -0.62052618298924286114047756;
    x[9] = -0.53662414814201989926416979;
    x[10] = -0.44703376953808917678060990;
    x[11] = -0.35270472553087811347103721;
    x[12] = -0.25463692616788984643980513;
    x[13] = -0.15386991360858354696379467;
    x[14] = -0.05147184255531769583302521;
    x[15] = 0.05147184255531769583302521;
    x[16] = 0.15386991360858354696379467;
    x[17] = 0.25463692616788984643980513;
    x[18] = 0.35270472553087811347103721;
    x[19] = 0.44703376953808917678060990;
    x[20] = 0.53662414814201989926416979;
    x[21] = 0.62052618298924286114047756;
    x[22] = 0.69785049479331579693229239;
    x[23] = 0.76777743210482619491797734;
    x[24] = 0.82956576238276839744289812;
    x[25] = 0.88256053579205268154311646;
    x[26] = 0.92620004742927432587932428;
    x[27] = 0.96002186496830751221687103;
    x[28] = 0.98366812327974720997003258;
    x[29] = 0.99689348407464954027163005;

    w[0] = 0.007968192496166605615465883475;
    w[1] = 0.018466468311090959142302131912;
    w[2] = 0.028784707883323369349719179611;
    w[3] = 0.038799192569627049596801936446;
    w[4] = 0.048402672830594052902938140423;
    w[5] = 0.057493156217619066481721689402;
    w[6] = 0.06597422988218049512812851512;
    w[7] = 0.07375597473770520626824385002;
    w[8] = 0.08075589522942021535469493846;
    w[9] = 0.08689978720108297980238753072;
    w[10] = 0.09212252223778612871763270709;
    w[11] = 0.09636873717464425963946862635;
    w[12] = 0.09959342058679526706278028210;
    w[13] = 0.10176238974840550459642895217;
    w[14] = 0.10285265289355884034128563671;
    w[15] = 0.10285265289355884034128563671;
    w[16] = 0.10176238974840550459642895217;
    w[17] = 0.09959342058679526706278028210;
    w[18] = 0.09636873717464425963946862635;
    w[19] = 0.09212252223778612871763270709;
    w[20] = 0.08689978720108297980238753072;
    w[21] = 0.08075589522942021535469493846;
    w[22] = 0.07375597473770520626824385002;
    w[23] = 0.06597422988218049512812851512;
    w[24] = 0.057493156217619066481721689402;
    w[25] = 0.048402672830594052902938140423;
    w[26] = 0.038799192569627049596801936446;
    w[27] = 0.028784707883323369349719179611;
    w[28] = 0.018466468311090959142302131912;
    w[29] = 0.007968192496166605615465883475;
  }
  else if ( n == 31 )
  {
    x[0] = -0.99708748181947707405562655;
    x[1] = -0.98468590966515248400246517;
    x[2] = -0.96250392509294966178905240;
    x[3] = -0.93075699789664816495694576;
    x[4] = -0.88976002994827104337419201;
    x[5] = -0.83992032014626734008690454;
    x[6] = -0.78173314841662494040636002;
    x[7] = -0.71577678458685328390597087;
    x[8] = -0.64270672292426034618441820;
    x[9] = -0.56324916140714926272094492;
    x[10] = -0.47819378204490248044059404;
    x[11] = -0.38838590160823294306135146;
    x[12] = -0.29471806998170161661790390;
    x[13] = -0.19812119933557062877241300;
    x[14] = -0.09955531215234152032517479;
    x[15] = 0.00000000000000000000000000;
    x[16] = 0.09955531215234152032517479;
    x[17] = 0.19812119933557062877241300;
    x[18] = 0.29471806998170161661790390;
    x[19] = 0.38838590160823294306135146;
    x[20] = 0.47819378204490248044059404;
    x[21] = 0.56324916140714926272094492;
    x[22] = 0.64270672292426034618441820;
    x[23] = 0.71577678458685328390597087;
    x[24] = 0.78173314841662494040636002;
    x[25] = 0.83992032014626734008690454;
    x[26] = 0.88976002994827104337419201;
    x[27] = 0.93075699789664816495694576;
    x[28] = 0.96250392509294966178905240;
    x[29] = 0.98468590966515248400246517;
    x[30] = 0.99708748181947707405562655;

    w[0] = 0.0074708315792487758586968750322;
    w[1] = 0.017318620790310582463157996087;
    w[2] = 0.027009019184979421800608708092;
    w[3] = 0.036432273912385464024392010468;
    w[4] = 0.045493707527201102902315857895;
    w[5] = 0.054103082424916853711666259087;
    w[6] = 0.062174786561028426910343543687;
    w[7] = 0.069628583235410366167756126255;
    w[8] = 0.076390386598776616426357674901;
    w[9] = 0.082392991761589263903823367432;
    w[10] = 0.087576740608477876126198069695;
    w[11] = 0.091890113893641478215362871607;
    w[12] = 0.095290242912319512807204197488;
    w[13] = 0.097743335386328725093474010979;
    w[14] = 0.099225011226672307874875514429;
    w[15] = 0.09972054479342645142753383373;
    w[16] = 0.099225011226672307874875514429;
    w[17] = 0.097743335386328725093474010979;
    w[18] = 0.095290242912319512807204197488;
    w[19] = 0.091890113893641478215362871607;
    w[20] = 0.087576740608477876126198069695;
    w[21] = 0.082392991761589263903823367432;
    w[22] = 0.076390386598776616426357674901;
    w[23] = 0.069628583235410366167756126255;
    w[24] = 0.062174786561028426910343543687;
    w[25] = 0.054103082424916853711666259087;
    w[26] = 0.045493707527201102902315857895;
    w[27] = 0.036432273912385464024392010468;
    w[28] = 0.027009019184979421800608708092;
    w[29] = 0.017318620790310582463157996087;
    w[30] = 0.0074708315792487758586968750322;
  }
  else if ( n == 32 )
  {
    x[0] = -0.99726386184948156354498113;
    x[1] = -0.98561151154526833540017504;
    x[2] = -0.96476225558750643077381193;
    x[3] = -0.93490607593773968917091913;
    x[4] = -0.89632115576605212396530724;
    x[5] = -0.84936761373256997013369300;
    x[6] = -0.79448379596794240696309730;
    x[7] = -0.73218211874028968038742667;
    x[8] = -0.66304426693021520097511517;
    x[9] = -0.58771575724076232904074548;
    x[10] = -0.50689990893222939002374747;
    x[11] = -0.42135127613063534536411944;
    x[12] = -0.33186860228212764977991681;
    x[13] = -0.23928736225213707454460321;
    x[14] = -0.14447196158279649348518637;
    x[15] = -0.04830766568773831623481257;
    x[16] = 0.04830766568773831623481257;
    x[17] = 0.14447196158279649348518637;
    x[18] = 0.23928736225213707454460321;
    x[19] = 0.33186860228212764977991681;
    x[20] = 0.42135127613063534536411944;
    x[21] = 0.50689990893222939002374747;
    x[22] = 0.58771575724076232904074548;
    x[23] = 0.66304426693021520097511517;
    x[24] = 0.73218211874028968038742667;
    x[25] = 0.79448379596794240696309730;
    x[26] = 0.84936761373256997013369300;
    x[27] = 0.89632115576605212396530724;
    x[28] = 0.93490607593773968917091913;
    x[29] = 0.96476225558750643077381193;
    x[30] = 0.98561151154526833540017504;
    x[31] = 0.99726386184948156354498113;

    w[0] = 0.007018610009470096600407063739;
    w[1] = 0.016274394730905670605170562206;
    w[2] = 0.025392065309262059455752589789;
    w[3] = 0.034273862913021433102687732252;
    w[4] = 0.042835898022226680656878646606;
    w[5] = 0.050998059262376176196163244690;
    w[6] = 0.058684093478535547145283637300;
    w[7] = 0.06582222277636184683765006371;
    w[8] = 0.07234579410884850622539935648;
    w[9] = 0.07819389578707030647174091883;
    w[10] = 0.08331192422694675522219907460;
    w[11] = 0.08765209300440381114277146275;
    w[12] = 0.09117387869576388471286857711;
    w[13] = 0.09384439908080456563918023767;
    w[14] = 0.09563872007927485941908200220;
    w[15] = 0.09654008851472780056676483006;
    w[16] = 0.09654008851472780056676483006;
    w[17] = 0.09563872007927485941908200220;
    w[18] = 0.09384439908080456563918023767;
    w[19] = 0.09117387869576388471286857711;
    w[20] = 0.08765209300440381114277146275;
    w[21] = 0.08331192422694675522219907460;
    w[22] = 0.07819389578707030647174091883;
    w[23] = 0.07234579410884850622539935648;
    w[24] = 0.06582222277636184683765006371;
    w[25] = 0.058684093478535547145283637300;
    w[26] = 0.050998059262376176196163244690;
    w[27] = 0.042835898022226680656878646606;
    w[28] = 0.034273862913021433102687732252;
    w[29] = 0.025392065309262059455752589789;
    w[30] = 0.016274394730905670605170562206;
    w[31] = 0.007018610009470096600407063739;
  }
  else if ( n == 33 )
  {
    x[0] = -0.99742469424645521726616802;
    x[1] = -0.98645572623064248811037570;
    x[2] = -0.96682290968999276892837771;
    x[3] = -0.93869437261116835035583512;
    x[4] = -0.90231676774343358304053133;
    x[5] = -0.85800965267650406464306148;
    x[6] = -0.80616235627416658979620087;
    x[7] = -0.74723049644956215785905512;
    x[8] = -0.68173195996974278626821595;
    x[9] = -0.61024234583637902730728751;
    x[10] = -0.53338990478634764354889426;
    x[11] = -0.45185001727245069572599328;
    x[12] = -0.36633925774807334107022062;
    x[13] = -0.27760909715249702940324807;
    x[14] = -0.18643929882799157233579876;
    x[15] = -0.09363106585473338567074292;
    x[16] = 0.00000000000000000000000000;
    x[17] = 0.09363106585473338567074292;
    x[18] = 0.18643929882799157233579876;
    x[19] = 0.27760909715249702940324807;
    x[20] = 0.36633925774807334107022062;
    x[21] = 0.45185001727245069572599328;
    x[22] = 0.53338990478634764354889426;
    x[23] = 0.61024234583637902730728751;
    x[24] = 0.68173195996974278626821595;
    x[25] = 0.74723049644956215785905512;
    x[26] = 0.80616235627416658979620087;
    x[27] = 0.85800965267650406464306148;
    x[28] = 0.90231676774343358304053133;
    x[29] = 0.93869437261116835035583512;
    x[30] = 0.96682290968999276892837771;
    x[31] = 0.98645572623064248811037570;
    x[32] = 0.99742469424645521726616802;

    w[0] = 0.0066062278475873780586492352085;
    w[1] = 0.015321701512934676127945768534;
    w[2] = 0.023915548101749480350533257529;
    w[3] = 0.032300358632328953281561447250;
    w[4] = 0.040401541331669591563409790527;
    w[5] = 0.048147742818711695670146880138;
    w[6] = 0.055470846631663561284944495439;
    w[7] = 0.062306482530317480031627725771;
    w[8] = 0.068594572818656712805955073015;
    w[9] = 0.074279854843954149342472175919;
    w[10] = 0.079312364794886738363908384942;
    w[11] = 0.083647876067038707613928014518;
    w[12] = 0.087248287618844337607281670945;
    w[13] = 0.090081958660638577239743705500;
    w[14] = 0.092123986643316846213240977717;
    w[15] = 0.093356426065596116160999126274;
    w[16] = 0.09376844616020999656730454155;
    w[17] = 0.093356426065596116160999126274;
    w[18] = 0.092123986643316846213240977717;
    w[19] = 0.090081958660638577239743705500;
    w[20] = 0.087248287618844337607281670945;
    w[21] = 0.083647876067038707613928014518;
    w[22] = 0.079312364794886738363908384942;
    w[23] = 0.074279854843954149342472175919;
    w[24] = 0.068594572818656712805955073015;
    w[25] = 0.062306482530317480031627725771;
    w[26] = 0.055470846631663561284944495439;
    w[27] = 0.048147742818711695670146880138;
    w[28] = 0.040401541331669591563409790527;
    w[29] = 0.032300358632328953281561447250;
    w[30] = 0.023915548101749480350533257529;
    w[31] = 0.015321701512934676127945768534;
    w[32] = 0.0066062278475873780586492352085;
  }
  else if ( n == 63 )
  {
    x[0] = -0.99928298402912378037893614;
    x[1] = -0.99622401277797010860219336;
    x[2] = -0.99072854689218946681089467;
    x[3] = -0.98280881059372723486251141;
    x[4] = -0.97248403469757002280196068;
    x[5] = -0.95977944975894192707035417;
    x[6] = -0.94472613404100980296637532;
    x[7] = -0.92736092062184320544703138;
    x[8] = -0.90772630277853155803695313;
    x[9] = -0.88587032850785342629029846;
    x[10] = -0.86184648236412371953961184;
    x[11] = -0.83571355431950284347180777;
    x[12] = -0.80753549577345676005146599;
    x[13] = -0.7773812629903723355633302;
    x[14] = -0.7453246483178474178293217;
    x[15] = -0.7114440995848458078514315;
    x[16] = -0.6758225281149860901311033;
    x[17] = -0.6385471058213653850003070;
    x[18] = -0.5997090518776252357390089;
    x[19] = -0.5594034094862850132676978;
    x[20] = -0.5177288132900332481244776;
    x[21] = -0.4747872479948043999222123;
    x[22] = -0.4306837987951116006620889;
    x[23] = -0.3855263942122478924776150;
    x[24] = -0.3394255419745844024688344;
    x[25] = -0.2924940585862514400361572;
    x[26] = -0.2448467932459533627484046;
    x[27] = -0.1966003467915066845576275;
    x[28] = -0.1478727863578719685698391;
    x[29] = -0.0987833564469452795297037;
    x[30] = -0.0494521871161596272342338;
    x[31] = 0.0000000000000000000000000;
    x[32] = 0.0494521871161596272342338;
    x[33] = 0.0987833564469452795297037;
    x[34] = 0.1478727863578719685698391;
    x[35] = 0.1966003467915066845576275;
    x[36] = 0.2448467932459533627484046;
    x[37] = 0.2924940585862514400361572;
    x[38] = 0.3394255419745844024688344;
    x[39] = 0.3855263942122478924776150;
    x[40] = 0.4306837987951116006620889;
    x[41] = 0.4747872479948043999222123;
    x[42] = 0.5177288132900332481244776;
    x[43] = 0.5594034094862850132676978;
    x[44] = 0.5997090518776252357390089;
    x[45] = 0.6385471058213653850003070;
    x[46] = 0.6758225281149860901311033;
    x[47] = 0.7114440995848458078514315;
    x[48] = 0.7453246483178474178293217;
    x[49] = 0.7773812629903723355633302;
    x[50] = 0.8075354957734567600514660;
    x[51] = 0.8357135543195028434718078;
    x[52] = 0.8618464823641237195396118;
    x[53] = 0.8858703285078534262902985;
    x[54] = 0.9077263027785315580369531;
    x[55] = 0.9273609206218432054470314;
    x[56] = 0.9447261340410098029663753;
    x[57] = 0.9597794497589419270703542;
    x[58] = 0.9724840346975700228019607;
    x[59] = 0.9828088105937272348625114;
    x[60] = 0.9907285468921894668108947;
    x[61] = 0.9962240127779701086021934;
    x[62] = 0.9992829840291237803789361;

    w[0] = 0.0018398745955770841170924455540;
    w[1] = 0.0042785083468637618660784110826;
    w[2] = 0.0067102917659601362519069307298;
    w[3] = 0.0091259686763266563540586454218;
    w[4] = 0.011519376076880041750750606149;
    w[5] = 0.013884612616115610824866086368;
    w[6] = 0.016215878410338338882283672975;
    w[7] = 0.018507464460161270409260545805;
    w[8] = 0.020753761258039090775341953421;
    w[9] = 0.022949271004889933148942319562;
    w[10] = 0.025088620553344986618630138068;
    w[11] = 0.027166574359097933225189839439;
    w[12] = 0.029178047208280526945551502154;
    w[13] = 0.031118116622219817508215988557;
    w[14] = 0.032982034883779341765683179672;
    w[15] = 0.034765240645355877697180504643;
    w[16] = 0.036463370085457289630452409788;
    w[17] = 0.038072267584349556763638324928;
    w[18] = 0.039587995891544093984807928149;
    w[19] = 0.041006845759666398635110037009;
    w[20] = 0.042325345020815822982505485403;
    w[21] = 0.043540267083027590798964315704;
    w[22] = 0.044648638825941395370332669517;
    w[23] = 0.045647747876292608685885992609;
    w[24] = 0.046535149245383696510395418747;
    w[25] = 0.047308671312268919080604988339;
    w[26] = 0.047966421137995131411052756195;
    w[27] = 0.048506789097883847864090099146;
    w[28] = 0.048928452820511989944709361549;
    w[29] = 0.049230380423747560785043116988;
    w[30] = 0.049411833039918178967039646117;
    w[31] = 0.04947236662393102088866936042;
    w[32] = 0.049411833039918178967039646117;
    w[33] = 0.049230380423747560785043116988;
    w[34] = 0.048928452820511989944709361549;
    w[35] = 0.048506789097883847864090099146;
    w[36] = 0.047966421137995131411052756195;
    w[37] = 0.047308671312268919080604988339;
    w[38] = 0.046535149245383696510395418747;
    w[39] = 0.045647747876292608685885992609;
    w[40] = 0.044648638825941395370332669517;
    w[41] = 0.043540267083027590798964315704;
    w[42] = 0.042325345020815822982505485403;
    w[43] = 0.041006845759666398635110037009;
    w[44] = 0.039587995891544093984807928149;
    w[45] = 0.038072267584349556763638324928;
    w[46] = 0.036463370085457289630452409788;
    w[47] = 0.034765240645355877697180504643;
    w[48] = 0.032982034883779341765683179672;
    w[49] = 0.031118116622219817508215988557;
    w[50] = 0.029178047208280526945551502154;
    w[51] = 0.027166574359097933225189839439;
    w[52] = 0.025088620553344986618630138068;
    w[53] = 0.022949271004889933148942319562;
    w[54] = 0.020753761258039090775341953421;
    w[55] = 0.018507464460161270409260545805;
    w[56] = 0.016215878410338338882283672975;
    w[57] = 0.013884612616115610824866086368;
    w[58] = 0.011519376076880041750750606149;
    w[59] = 0.0091259686763266563540586454218;
    w[60] = 0.0067102917659601362519069307298;
    w[61] = 0.0042785083468637618660784110826;
    w[62] = 0.0018398745955770841170924455540;
  }
  else if ( n == 64 )
  {
    x[0] = -0.99930504173577213945690562;
    x[1] = -0.99634011677195527934692450;
    x[2] = -0.99101337147674432073938238;
    x[3] = -0.98333625388462595693129930;
    x[4] = -0.97332682778991096374185351;
    x[5] = -0.96100879965205371891861412;
    x[6] = -0.94641137485840281606248149;
    x[7] = -0.92956917213193957582149015;
    x[8] = -0.91052213707850280575638067;
    x[9] = -0.88931544599511410585340404;
    x[10] = -0.86599939815409281976078339;
    x[11] = -0.8406292962525803627516915;
    x[12] = -0.8132653151227975597419233;
    x[13] = -0.7839723589433414076102205;
    x[14] = -0.7528199072605318966118638;
    x[15] = -0.7198818501716108268489402;
    x[16] = -0.6852363130542332425635584;
    x[17] = -0.6489654712546573398577612;
    x[18] = -0.6111553551723932502488530;
    x[19] = -0.5718956462026340342838781;
    x[20] = -0.5312794640198945456580139;
    x[21] = -0.4894031457070529574785263;
    x[22] = -0.4463660172534640879849477;
    x[23] = -0.4022701579639916036957668;
    x[24] = -0.3572201583376681159504426;
    x[25] = -0.3113228719902109561575127;
    x[26] = -0.2646871622087674163739642;
    x[27] = -0.2174236437400070841496487;
    x[28] = -0.1696444204239928180373136;
    x[29] = -0.1214628192961205544703765;
    x[30] = -0.0729931217877990394495429;
    x[31] = -0.0243502926634244325089558;
    x[32] = 0.0243502926634244325089558;
    x[33] = 0.0729931217877990394495429;
    x[34] = 0.1214628192961205544703765;
    x[35] = 0.1696444204239928180373136;
    x[36] = 0.2174236437400070841496487;
    x[37] = 0.2646871622087674163739642;
    x[38] = 0.3113228719902109561575127;
    x[39] = 0.3572201583376681159504426;
    x[40] = 0.4022701579639916036957668;
    x[41] = 0.4463660172534640879849477;
    x[42] = 0.4894031457070529574785263;
    x[43] = 0.5312794640198945456580139;
    x[44] = 0.5718956462026340342838781;
    x[45] = 0.6111553551723932502488530;
    x[46] = 0.6489654712546573398577612;
    x[47] = 0.6852363130542332425635584;
    x[48] = 0.7198818501716108268489402;
    x[49] = 0.7528199072605318966118638;
    x[50] = 0.7839723589433414076102205;
    x[51] = 0.8132653151227975597419233;
    x[52] = 0.8406292962525803627516915;
    x[53] = 0.8659993981540928197607834;
    x[54] = 0.8893154459951141058534040;
    x[55] = 0.9105221370785028057563807;
    x[56] = 0.9295691721319395758214902;
    x[57] = 0.9464113748584028160624815;
    x[58] = 0.9610087996520537189186141;
    x[59] = 0.9733268277899109637418535;
    x[60] = 0.9833362538846259569312993;
    x[61] = 0.9910133714767443207393824;
    x[62] = 0.9963401167719552793469245;
    x[63] = 0.9993050417357721394569056;

    w[0] = 0.0017832807216964329472960791450;
    w[1] = 0.0041470332605624676352875357286;
    w[2] = 0.006504457968978362856117360400;
    w[3] = 0.008846759826363947723030914660;
    w[4] = 0.011168139460131128818590493019;
    w[5] = 0.013463047896718642598060766686;
    w[6] = 0.015726030476024719321965995298;
    w[7] = 0.017951715775697343085045302001;
    w[8] = 0.020134823153530209372340316729;
    w[9] = 0.022270173808383254159298330384;
    w[10] = 0.024352702568710873338177550409;
    w[11] = 0.026377469715054658671691792625;
    w[12] = 0.028339672614259483227511305200;
    w[13] = 0.030234657072402478867974059820;
    w[14] = 0.032057928354851553585467504348;
    w[15] = 0.033805161837141609391565482111;
    w[16] = 0.035472213256882383810693146715;
    w[17] = 0.037055128540240046040415101810;
    w[18] = 0.038550153178615629128962496947;
    w[19] = 0.039953741132720341386656926128;
    w[20] = 0.041262563242623528610156297474;
    w[21] = 0.042473515123653589007339767909;
    w[22] = 0.043583724529323453376827860974;
    w[23] = 0.044590558163756563060134710031;
    w[24] = 0.045491627927418144479770996971;
    w[25] = 0.046284796581314417295953249232;
    w[26] = 0.046968182816210017325326285755;
    w[27] = 0.047540165714830308662282206944;
    w[28] = 0.04799938859645830772812617987;
    w[29] = 0.04834476223480295716976952716;
    w[30] = 0.04857546744150342693479906678;
    w[31] = 0.04869095700913972038336539073;
    w[32] = 0.04869095700913972038336539073;
    w[33] = 0.04857546744150342693479906678;
    w[34] = 0.04834476223480295716976952716;
    w[35] = 0.04799938859645830772812617987;
    w[36] = 0.047540165714830308662282206944;
    w[37] = 0.046968182816210017325326285755;
    w[38] = 0.046284796581314417295953249232;
    w[39] = 0.045491627927418144479770996971;
    w[40] = 0.044590558163756563060134710031;
    w[41] = 0.043583724529323453376827860974;
    w[42] = 0.042473515123653589007339767909;
    w[43] = 0.041262563242623528610156297474;
    w[44] = 0.039953741132720341386656926128;
    w[45] = 0.038550153178615629128962496947;
    w[46] = 0.037055128540240046040415101810;
    w[47] = 0.035472213256882383810693146715;
    w[48] = 0.033805161837141609391565482111;
    w[49] = 0.032057928354851553585467504348;
    w[50] = 0.030234657072402478867974059820;
    w[51] = 0.028339672614259483227511305200;
    w[52] = 0.026377469715054658671691792625;
    w[53] = 0.024352702568710873338177550409;
    w[54] = 0.022270173808383254159298330384;
    w[55] = 0.020134823153530209372340316729;
    w[56] = 0.017951715775697343085045302001;
    w[57] = 0.015726030476024719321965995298;
    w[58] = 0.013463047896718642598060766686;
    w[59] = 0.011168139460131128818590493019;
    w[60] = 0.008846759826363947723030914660;
    w[61] = 0.006504457968978362856117360400;
    w[62] = 0.0041470332605624676352875357286;
    w[63] = 0.0017832807216964329472960791450;
  }
  else if ( n == 65 )
  {
    x[0] = -0.99932609707541287726569361;
    x[1] = -0.99645094806184916305579494;
    x[2] = -0.99128527617680166872182118;
    x[3] = -0.98383981218703494137763778;
    x[4] = -0.97413153983355116907496789;
    x[5] = -0.96218275471805523771198375;
    x[6] = -0.94802092816840750637376974;
    x[7] = -0.93167862822874933796567699;
    x[8] = -0.91319344054284626173654692;
    x[9] = -0.89260788050473893142328554;
    x[10] = -0.8699692949264070361941320;
    x[11] = -0.8453297528999302839424500;
    x[12] = -0.8187459259226514534339191;
    x[13] = -0.7902789574921218430473804;
    x[14] = -0.7599943224419997868739828;
    x[15] = -0.7279616763294246790119737;
    x[16] = -0.6942546952139916335526225;
    x[17] = -0.6589509061936251330409408;
    x[18] = -0.6221315090854002415825996;
    x[19] = -0.5838811896604873133271545;
    x[20] = -0.5442879248622271385455725;
    x[21] = -0.5034427804550068823410431;
    x[22] = -0.4614397015691450576978341;
    x[23] = -0.4183752966234090092641990;
    x[24] = -0.3743486151220660120087939;
    x[25] = -0.3294609198374864076452867;
    x[26] = -0.2838154539022487306176554;
    x[27] = -0.2375172033464168065707124;
    x[28] = -0.1906726556261427697749124;
    x[29] = -0.1433895546989751711312496;
    x[30] = -0.0957766532091975056522186;
    x[31] = -0.0479434623531718575225298;
    x[32] = 0.0000000000000000000000000;
    x[33] = 0.0479434623531718575225298;
    x[34] = 0.0957766532091975056522186;
    x[35] = 0.1433895546989751711312496;
    x[36] = 0.1906726556261427697749124;
    x[37] = 0.2375172033464168065707124;
    x[38] = 0.2838154539022487306176554;
    x[39] = 0.3294609198374864076452867;
    x[40] = 0.3743486151220660120087939;
    x[41] = 0.4183752966234090092641990;
    x[42] = 0.4614397015691450576978341;
    x[43] = 0.5034427804550068823410431;
    x[44] = 0.5442879248622271385455725;
    x[45] = 0.5838811896604873133271545;
    x[46] = 0.6221315090854002415825996;
    x[47] = 0.6589509061936251330409408;
    x[48] = 0.6942546952139916335526225;
    x[49] = 0.7279616763294246790119737;
    x[50] = 0.7599943224419997868739828;
    x[51] = 0.7902789574921218430473804;
    x[52] = 0.8187459259226514534339191;
    x[53] = 0.8453297528999302839424500;
    x[54] = 0.8699692949264070361941320;
    x[55] = 0.8926078805047389314232855;
    x[56] = 0.9131934405428462617365469;
    x[57] = 0.9316786282287493379656770;
    x[58] = 0.9480209281684075063737697;
    x[59] = 0.9621827547180552377119837;
    x[60] = 0.9741315398335511690749679;
    x[61] = 0.9838398121870349413776378;
    x[62] = 0.9912852761768016687218212;
    x[63] = 0.9964509480618491630557949;
    x[64] = 0.9993260970754128772656936;

    w[0] = 0.0017292582513002508983395851463;
    w[1] = 0.0040215241720037363470786599528;
    w[2] = 0.0063079425789717545501888719039;
    w[3] = 0.0085801482668814598936358121592;
    w[4] = 0.0108326787895979686215140551272;
    w[5] = 0.013060311639994846336168342922;
    w[6] = 0.015257912146448310349265388145;
    w[7] = 0.017420421997670248495365759969;
    w[8] = 0.019542865836750062826837429313;
    w[9] = 0.021620361284934062841654274667;
    w[10] = 0.023648129691287236698780978994;
    w[11] = 0.025621506938037758214084978694;
    w[12] = 0.027535954088450343942499722327;
    w[13] = 0.029387067789310668062644859210;
    w[14] = 0.031170590380189142464431845777;
    w[15] = 0.032882419676368574984049638008;
    w[16] = 0.034518618398549058625221276859;
    w[17] = 0.036075423225565273932166270524;
    w[18] = 0.037549253448257709809772223198;
    w[19] = 0.038936719204051197616673806364;
    w[20] = 0.040234629273005533815446337743;
    w[21] = 0.041439998417240293022686299233;
    w[22] = 0.042550054246755802719217150803;
    w[23] = 0.043562243595800486532284821661;
    w[24] = 0.044474238395082974427323504000;
    w[25] = 0.045283941026300230657128240574;
    w[26] = 0.045989489146651696963893390818;
    w[27] = 0.046589259972233498302255136790;
    w[28] = 0.047081874010454522246006808290;
    w[29] = 0.047466198232885503152644458740;
    w[30] = 0.047741348681240621559038972227;
    w[31] = 0.047906692500495862031347289176;
    w[32] = 0.04796184939446661812070762137;
    w[33] = 0.047906692500495862031347289176;
    w[34] = 0.047741348681240621559038972227;
    w[35] = 0.047466198232885503152644458740;
    w[36] = 0.047081874010454522246006808290;
    w[37] = 0.046589259972233498302255136790;
    w[38] = 0.045989489146651696963893390818;
    w[39] = 0.045283941026300230657128240574;
    w[40] = 0.044474238395082974427323504000;
    w[41] = 0.043562243595800486532284821661;
    w[42] = 0.042550054246755802719217150803;
    w[43] = 0.041439998417240293022686299233;
    w[44] = 0.040234629273005533815446337743;
    w[45] = 0.038936719204051197616673806364;
    w[46] = 0.037549253448257709809772223198;
    w[47] = 0.036075423225565273932166270524;
    w[48] = 0.034518618398549058625221276859;
    w[49] = 0.032882419676368574984049638008;
    w[50] = 0.031170590380189142464431845777;
    w[51] = 0.029387067789310668062644859210;
    w[52] = 0.027535954088450343942499722327;
    w[53] = 0.025621506938037758214084978694;
    w[54] = 0.023648129691287236698780978994;
    w[55] = 0.021620361284934062841654274667;
    w[56] = 0.019542865836750062826837429313;
    w[57] = 0.017420421997670248495365759969;
    w[58] = 0.015257912146448310349265388145;
    w[59] = 0.013060311639994846336168342922;
    w[60] = 0.0108326787895979686215140551272;
    w[61] = 0.0085801482668814598936358121592;
    w[62] = 0.0063079425789717545501888719039;
    w[63] = 0.0040215241720037363470786599528;
    w[64] = 0.0017292582513002508983395851463;
  }
  else if ( n == 127 )
  {
    x[0] = -0.9998221304153061462673512;
    x[1] = -0.9990629343553118951383159;
    x[2] = -0.9976975661898046210744170;
    x[3] = -0.9957265513520272266354334;
    x[4] = -0.9931510492545171473611308;
    x[5] = -0.9899726145914841576077867;
    x[6] = -0.9861931740169316667104383;
    x[7] = -0.9818150208038141100334631;
    x[8] = -0.9768408123430703268174439;
    x[9] = -0.9712735681615291922889469;
    x[10] = -0.9651166679452921210908251;
    x[11] = -0.9583738494252387711491029;
    x[12] = -0.9510492060778803105479076;
    x[13] = -0.9431471846248148273454496;
    x[14] = -0.9346725823247379685736349;
    x[15] = -0.9256305440562338491274647;
    x[16] = -0.9160265591914658093130886;
    x[17] = -0.9058664582618213828024613;
    x[18] = -0.8951564094170837089690438;
    x[19] = -0.8839029146800265699452579;
    x[20] = -0.8721128059985607114196375;
    x[21] = -0.8597932410977408098120313;
    x[22] = -0.8469516991340975984533393;
    x[23] = -0.8335959761548995143795572;
    x[24] = -0.8197341803650786741551191;
    x[25] = -0.8053747272046802146665608;
    x[26] = -0.7905263342398137999454500;
    x[27] = -0.7751980158702023824449628;
    x[28] = -0.7593990778565366715566637;
    x[29] = -0.7431391116709545129205669;
    x[30] = -0.7264279886740726855356929;
    x[31] = -0.7092758541221045609994446;
    x[32] = -0.6916931210077006701564414;
    x[33] = -0.6736904637382504853466825;
    x[34] = -0.6552788116554826302767651;
    x[35] = -0.6364693424002972413476082;
    x[36] = -0.6172734751268582838576392;
    x[37] = -0.5977028635700652293844120;
    x[38] = -0.5777693889706125800032517;
    x[39] = -0.5574851528619322329218619;
    x[40] = -0.5368624697233975674581664;
    x[41] = -0.5159138595042493572772773;
    x[42] = -0.4946520400227821173949402;
    x[43] = -0.4730899192454052416450999;
    x[44] = -0.4512405874502662273318986;
    x[45] = -0.4291173092801933762625441;
    x[46] = -0.4067335156897825634086729;
    x[47] = -0.3841027957915169357790778;
    x[48] = -0.3612388886058697060709248;
    x[49] = -0.3381556747203985013760003;
    x[50] = -0.3148671678628949814860148;
    x[51] = -0.2913875063937056207945188;
    x[52] = -0.2677309447223886208883435;
    x[53] = -0.2439118446539178579707132;
    x[54] = -0.2199446666696875424545234;
    x[55] = -0.1958439611486108515042816;
    x[56] = -0.1716243595336421650083449;
    x[57] = -0.1473005654490856693893293;
    x[58] = -0.1228873457740829717260337;
    x[59] = -0.0983995216776989707510918;
    x[60] = -0.0738519596210485452734404;
    x[61] = -0.0492595623319266303153793;
    x[62] = -0.0246372597574209446148971;
    x[63] = 0.0000000000000000000000000;
    x[64] = 0.0246372597574209446148971;
    x[65] = 0.0492595623319266303153793;
    x[66] = 0.0738519596210485452734404;
    x[67] = 0.0983995216776989707510918;
    x[68] = 0.1228873457740829717260337;
    x[69] = 0.1473005654490856693893293;
    x[70] = 0.1716243595336421650083449;
    x[71] = 0.1958439611486108515042816;
    x[72] = 0.2199446666696875424545234;
    x[73] = 0.2439118446539178579707132;
    x[74] = 0.2677309447223886208883435;
    x[75] = 0.2913875063937056207945188;
    x[76] = 0.3148671678628949814860148;
    x[77] = 0.3381556747203985013760003;
    x[78] = 0.3612388886058697060709248;
    x[79] = 0.3841027957915169357790778;
    x[80] = 0.4067335156897825634086729;
    x[81] = 0.4291173092801933762625441;
    x[82] = 0.4512405874502662273318986;
    x[83] = 0.4730899192454052416450999;
    x[84] = 0.4946520400227821173949402;
    x[85] = 0.5159138595042493572772773;
    x[86] = 0.5368624697233975674581664;
    x[87] = 0.5574851528619322329218619;
    x[88] = 0.5777693889706125800032517;
    x[89] = 0.5977028635700652293844120;
    x[90] = 0.6172734751268582838576392;
    x[91] = 0.6364693424002972413476082;
    x[92] = 0.6552788116554826302767651;
    x[93] = 0.6736904637382504853466825;
    x[94] = 0.6916931210077006701564414;
    x[95] = 0.7092758541221045609994446;
    x[96] = 0.7264279886740726855356929;
    x[97] = 0.7431391116709545129205669;
    x[98] = 0.7593990778565366715566637;
    x[99] = 0.7751980158702023824449628;
    x[100] = 0.7905263342398137999454500;
    x[101] = 0.8053747272046802146665608;
    x[102] = 0.8197341803650786741551191;
    x[103] = 0.8335959761548995143795572;
    x[104] = 0.8469516991340975984533393;
    x[105] = 0.8597932410977408098120313;
    x[106] = 0.8721128059985607114196375;
    x[107] = 0.8839029146800265699452579;
    x[108] = 0.8951564094170837089690438;
    x[109] = 0.9058664582618213828024613;
    x[110] = 0.9160265591914658093130886;
    x[111] = 0.9256305440562338491274647;
    x[112] = 0.9346725823247379685736349;
    x[113] = 0.9431471846248148273454496;
    x[114] = 0.9510492060778803105479076;
    x[115] = 0.9583738494252387711491029;
    x[116] = 0.965116667945292121090825;
    x[117] = 0.971273568161529192288947;
    x[118] = 0.976840812343070326817444;
    x[119] = 0.981815020803814110033463;
    x[120] = 0.986193174016931666710438;
    x[121] = 0.989972614591484157607787;
    x[122] = 0.993151049254517147361131;
    x[123] = 0.995726551352027226635433;
    x[124] = 0.997697566189804621074417;
    x[125] = 0.999062934355311895138316;
    x[126] = 0.999822130415306146267351;

    w[0] = 0.00045645726109586662791936519265;
    w[1] = 0.00106227668695384869596523598532;
    w[2] = 0.0016683488125171936761028862915;
    w[3] = 0.0022734860707492547802810840776;
    w[4] = 0.0028772587656289004082883197514;
    w[5] = 0.0034792893810051465908910894100;
    w[6] = 0.0040792095178254605327114733457;
    w[7] = 0.0046766539777779034772638165663;
    w[8] = 0.0052712596565634400891303815906;
    w[9] = 0.0058626653903523901033648343751;
    w[10] = 0.0064505120486899171845442463869;
    w[11] = 0.0070344427036681608755685893033;
    w[12] = 0.0076141028256526859356393930849;
    w[13] = 0.0081891404887415730817235884719;
    w[14] = 0.0087592065795403145773316804234;
    w[15] = 0.0093239550065309714787536985834;
    w[16] = 0.0098830429087554914716648010900;
    w[17] = 0.0104361308631410052256731719977;
    w[18] = 0.0109828830900689757887996573761;
    w[19] = 0.011522967656921087154811609735;
    w[20] = 0.012056056679400848183529562145;
    w[21] = 0.012581826520465013101514365424;
    w[22] = 0.013099957986718627426172681913;
    w[23] = 0.013610136522139249906034237534;
    w[24] = 0.014112052399003395774044161634;
    w[25] = 0.014605400905893418351737288079;
    w[26] = 0.015089882532666922992635733981;
    w[27] = 0.015565203152273955098532590263;
    w[28] = 0.016031074199309941802254151843;
    w[29] = 0.016487212845194879399346060358;
    w[30] = 0.016933342169871654545878815295;
    w[31] = 0.017369191329918731922164721250;
    w[32] = 0.017794495722974774231027912900;
    w[33] = 0.018208997148375106468721469154;
    w[34] = 0.018612443963902310429440419899;
    w[35] = 0.019004591238555646611148901045;
    w[36] = 0.019385200901246454628112623489;
    w[37] = 0.019754041885329183081815217323;
    w[38] = 0.020110890268880247225644623956;
    w[39] = 0.020455529410639508279497065713;
    w[40] = 0.020787750081531811812652137291;
    w[41] = 0.021107350591688713643523847922;
    w[42] = 0.021414136912893259295449693234;
    w[43] = 0.021707922796373466052301324695;
    w[44] = 0.021988529885872983756478409759;
    w[45] = 0.022255787825930280235631416460;
    w[46] = 0.022509534365300608085694429903;
    w[47] = 0.022749615455457959852242553241;
    w[48] = 0.022975885344117206754377437839;
    w[49] = 0.023188206663719640249922582982;
    w[50] = 0.023386450514828194170722043497;
    w[51] = 0.023570496544381716050033676844;
    w[52] = 0.023740233018760777777714726703;
    w[53] = 0.023895556891620665983864481754;
    w[54] = 0.024036373866450369675132086026;
    w[55] = 0.024162598453819584716522917711;
    w[56] = 0.024274154023278979833195063937;
    w[57] = 0.024370972849882214952813561907;
    w[58] = 0.024452996155301467956140198472;
    w[59] = 0.024520174143511508275183033290;
    w[60] = 0.024572466031020653286354137335;
    w[61] = 0.024609840071630254092545634003;
    w[62] = 0.024632273575707679066033370218;
    w[63] = 0.02463975292396109441957941748;
    w[64] = 0.024632273575707679066033370218;
    w[65] = 0.024609840071630254092545634003;
    w[66] = 0.024572466031020653286354137335;
    w[67] = 0.024520174143511508275183033290;
    w[68] = 0.024452996155301467956140198472;
    w[69] = 0.024370972849882214952813561907;
    w[70] = 0.024274154023278979833195063937;
    w[71] = 0.024162598453819584716522917711;
    w[72] = 0.024036373866450369675132086026;
    w[73] = 0.023895556891620665983864481754;
    w[74] = 0.023740233018760777777714726703;
    w[75] = 0.023570496544381716050033676844;
    w[76] = 0.023386450514828194170722043497;
    w[77] = 0.023188206663719640249922582982;
    w[78] = 0.022975885344117206754377437839;
    w[79] = 0.022749615455457959852242553241;
    w[80] = 0.022509534365300608085694429903;
    w[81] = 0.022255787825930280235631416460;
    w[82] = 0.021988529885872983756478409759;
    w[83] = 0.021707922796373466052301324695;
    w[84] = 0.021414136912893259295449693234;
    w[85] = 0.021107350591688713643523847922;
    w[86] = 0.020787750081531811812652137291;
    w[87] = 0.020455529410639508279497065713;
    w[88] = 0.020110890268880247225644623956;
    w[89] = 0.019754041885329183081815217323;
    w[90] = 0.019385200901246454628112623489;
    w[91] = 0.019004591238555646611148901045;
    w[92] = 0.018612443963902310429440419899;
    w[93] = 0.018208997148375106468721469154;
    w[94] = 0.017794495722974774231027912900;
    w[95] = 0.017369191329918731922164721250;
    w[96] = 0.016933342169871654545878815295;
    w[97] = 0.016487212845194879399346060358;
    w[98] = 0.016031074199309941802254151843;
    w[99] = 0.015565203152273955098532590263;
    w[100] = 0.015089882532666922992635733981;
    w[101] = 0.014605400905893418351737288079;
    w[102] = 0.014112052399003395774044161634;
    w[103] = 0.013610136522139249906034237534;
    w[104] = 0.013099957986718627426172681913;
    w[105] = 0.012581826520465013101514365424;
    w[106] = 0.012056056679400848183529562145;
    w[107] = 0.011522967656921087154811609735;
    w[108] = 0.0109828830900689757887996573761;
    w[109] = 0.0104361308631410052256731719977;
    w[110] = 0.0098830429087554914716648010900;
    w[111] = 0.0093239550065309714787536985834;
    w[112] = 0.0087592065795403145773316804234;
    w[113] = 0.0081891404887415730817235884719;
    w[114] = 0.0076141028256526859356393930849;
    w[115] = 0.0070344427036681608755685893033;
    w[116] = 0.0064505120486899171845442463869;
    w[117] = 0.0058626653903523901033648343751;
    w[118] = 0.0052712596565634400891303815906;
    w[119] = 0.0046766539777779034772638165663;
    w[120] = 0.0040792095178254605327114733457;
    w[121] = 0.0034792893810051465908910894100;
    w[122] = 0.0028772587656289004082883197514;
    w[123] = 0.0022734860707492547802810840776;
    w[124] = 0.0016683488125171936761028862915;
    w[125] = 0.00106227668695384869596523598532;
    w[126] = 0.00045645726109586662791936519265;
  }
  else if ( n == 128 )
  {
    x[0] = -0.9998248879471319144736081;
    x[1] = -0.9990774599773758950119878;
    x[2] = -0.9977332486255140198821574;
    x[3] = -0.9957927585349811868641612;
    x[4] = -0.9932571129002129353034372;
    x[5] = -0.9901278184917343833379303;
    x[6] = -0.9864067427245862088712355;
    x[7] = -0.9820961084357185360247656;
    x[8] = -0.9771984914639073871653744;
    x[9] = -0.9717168187471365809043384;
    x[10] = -0.9656543664319652686458290;
    x[11] = -0.9590147578536999280989185;
    x[12] = -0.9518019613412643862177963;
    x[13] = -0.9440202878302201821211114;
    x[14] = -0.9356743882779163757831268;
    x[15] = -0.9267692508789478433346245;
    x[16] = -0.9173101980809605370364836;
    x[17] = -0.9073028834017568139214859;
    x[18] = -0.8967532880491581843864474;
    x[19] = -0.8856677173453972174082924;
    x[20] = -0.8740527969580317986954180;
    x[21] = -0.8619154689395484605906323;
    x[22] = -0.8492629875779689691636001;
    x[23] = -0.8361029150609068471168753;
    x[24] = -0.8224431169556438424645942;
    x[25] = -0.8082917575079136601196422;
    x[26] = -0.7936572947621932902433329;
    x[27] = -0.7785484755064119668504941;
    x[28] = -0.7629743300440947227797691;
    x[29] = -0.7469441667970619811698824;
    x[30] = -0.7304675667419088064717369;
    x[31] = -0.7135543776835874133438599;
    x[32] = -0.6962147083695143323850866;
    x[33] = -0.6784589224477192593677557;
    x[34] = -0.6602976322726460521059468;
    x[35] = -0.6417416925623075571535249;
    x[36] = -0.6228021939105849107615396;
    x[37] = -0.6034904561585486242035732;
    x[38] = -0.5838180216287630895500389;
    x[39] = -0.5637966482266180839144308;
    x[40] = -0.5434383024128103634441936;
    x[41] = -0.5227551520511754784539479;
    x[42] = -0.5017595591361444642896063;
    x[43] = -0.4804640724041720258582757;
    x[44] = -0.4588814198335521954490891;
    x[45] = -0.4370245010371041629370429;
    x[46] = -0.4149063795522750154922739;
    x[47] = -0.3925402750332674427356482;
    x[48] = -0.3699395553498590266165917;
    x[49] = -0.3471177285976355084261628;
    x[50] = -0.3240884350244133751832523;
    x[51] = -0.3008654388776772026671541;
    x[52] = -0.2774626201779044028062316;
    x[53] = -0.2538939664226943208556180;
    x[54] = -0.2301735642266599864109866;
    x[55] = -0.2063155909020792171540580;
    x[56] = -0.1823343059853371824103826;
    x[57] = -0.1582440427142249339974755;
    x[58] = -0.1340591994611877851175753;
    x[59] = -0.1097942311276437466729747;
    x[60] = -0.0854636405045154986364980;
    x[61] = -0.0610819696041395681037870;
    x[62] = -0.0366637909687334933302153;
    x[63] = -0.0122236989606157641980521;
    x[64] = 0.0122236989606157641980521;
    x[65] = 0.0366637909687334933302153;
    x[66] = 0.0610819696041395681037870;
    x[67] = 0.0854636405045154986364980;
    x[68] = 0.1097942311276437466729747;
    x[69] = 0.1340591994611877851175753;
    x[70] = 0.1582440427142249339974755;
    x[71] = 0.1823343059853371824103826;
    x[72] = 0.2063155909020792171540580;
    x[73] = 0.2301735642266599864109866;
    x[74] = 0.2538939664226943208556180;
    x[75] = 0.2774626201779044028062316;
    x[76] = 0.3008654388776772026671541;
    x[77] = 0.3240884350244133751832523;
    x[78] = 0.3471177285976355084261628;
    x[79] = 0.3699395553498590266165917;
    x[80] = 0.3925402750332674427356482;
    x[81] = 0.4149063795522750154922739;
    x[82] = 0.4370245010371041629370429;
    x[83] = 0.4588814198335521954490891;
    x[84] = 0.4804640724041720258582757;
    x[85] = 0.5017595591361444642896063;
    x[86] = 0.5227551520511754784539479;
    x[87] = 0.5434383024128103634441936;
    x[88] = 0.5637966482266180839144308;
    x[89] = 0.5838180216287630895500389;
    x[90] = 0.6034904561585486242035732;
    x[91] = 0.6228021939105849107615396;
    x[92] = 0.6417416925623075571535249;
    x[93] = 0.6602976322726460521059468;
    x[94] = 0.6784589224477192593677557;
    x[95] = 0.6962147083695143323850866;
    x[96] = 0.7135543776835874133438599;
    x[97] = 0.7304675667419088064717369;
    x[98] = 0.7469441667970619811698824;
    x[99] = 0.7629743300440947227797691;
    x[100] = 0.7785484755064119668504941;
    x[101] = 0.7936572947621932902433329;
    x[102] = 0.8082917575079136601196422;
    x[103] = 0.8224431169556438424645942;
    x[104] = 0.8361029150609068471168753;
    x[105] = 0.8492629875779689691636001;
    x[106] = 0.8619154689395484605906323;
    x[107] = 0.8740527969580317986954180;
    x[108] = 0.8856677173453972174082924;
    x[109] = 0.8967532880491581843864474;
    x[110] = 0.9073028834017568139214859;
    x[111] = 0.9173101980809605370364836;
    x[112] = 0.926769250878947843334625;
    x[113] = 0.935674388277916375783127;
    x[114] = 0.944020287830220182121111;
    x[115] = 0.951801961341264386217796;
    x[116] = 0.959014757853699928098919;
    x[117] = 0.965654366431965268645829;
    x[118] = 0.971716818747136580904338;
    x[119] = 0.977198491463907387165374;
    x[120] = 0.982096108435718536024766;
    x[121] = 0.986406742724586208871236;
    x[122] = 0.990127818491734383337930;
    x[123] = 0.993257112900212935303437;
    x[124] = 0.995792758534981186864161;
    x[125] = 0.997733248625514019882157;
    x[126] = 0.999077459977375895011988;
    x[127] = 0.999824887947131914473608;

    w[0] = 0.00044938096029209037639429223999;
    w[1] = 0.0010458126793403487793128516001;
    w[2] = 0.0016425030186690295387908755948;
    w[3] = 0.0022382884309626187436220542727;
    w[4] = 0.0028327514714579910952857346468;
    w[5] = 0.0034255260409102157743377846601;
    w[6] = 0.0040162549837386423131943434863;
    w[7] = 0.0046045842567029551182905419803;
    w[8] = 0.0051901618326763302050707671348;
    w[9] = 0.0057726375428656985893346176261;
    w[10] = 0.006351663161707188787214327826;
    w[11] = 0.006926892566898813563426670360;
    w[12] = 0.007497981925634728687671962688;
    w[13] = 0.008064589890486057972928598698;
    w[14] = 0.008626377798616749704978843782;
    w[15] = 0.009183009871660874334478743688;
    w[16] = 0.009734153415006805863548266094;
    w[17] = 0.010279479015832157133215340326;
    w[18] = 0.010818660739503076247659646277;
    w[19] = 0.011351376324080416693281668453;
    w[20] = 0.011877307372740279575891106926;
    w[21] = 0.012396139543950922968821728197;
    w[22] = 0.012907562739267347220442834004;
    w[23] = 0.013411271288616332314488951616;
    w[24] = 0.013906964132951985244288007396;
    w[25] = 0.014394345004166846176823892009;
    w[26] = 0.014873122602147314252385498520;
    w[27] = 0.015343010768865144085990853741;
    w[28] = 0.015803728659399346858965631687;
    w[29] = 0.016255000909785187051657456477;
    w[30] = 0.016696557801589204589091507954;
    w[31] = 0.017128135423111376830680987619;
    w[32] = 0.017549475827117704648706925634;
    w[33] = 0.017960327185008685940196927525;
    w[34] = 0.018360443937331343221289290991;
    w[35] = 0.018749586940544708650919548474;
    w[36] = 0.019127523609950945486518531668;
    w[37] = 0.019494028058706602823021918681;
    w[38] = 0.019848881232830862219944413265;
    w[39] = 0.020191871042130041180673158406;
    w[40] = 0.020522792486960069432284967788;
    w[41] = 0.020841447780751149113583948423;
    w[42] = 0.021147646468221348537019535180;
    w[43] = 0.021441205539208460137111853878;
    w[44] = 0.021721949538052075375260957768;
    w[45] = 0.021989710668460491434122106599;
    w[46] = 0.022244328893799765104629133607;
    w[47] = 0.022485652032744966871824603941;
    w[48] = 0.022713535850236461309712635923;
    w[49] = 0.022927844143686846920410987209;
    w[50] = 0.023128448824387027879297902403;
    w[51] = 0.023315229994062760122415671273;
    w[52] = 0.023488076016535913153025273282;
    w[53] = 0.023646883584447615143651392303;
    w[54] = 0.023791557781003400638780709885;
    w[55] = 0.023922012136703455672450408817;
    w[56] = 0.024038168681024052637587316820;
    w[57] = 0.024139957989019284997716653890;
    w[58] = 0.024227319222815248120093308442;
    w[59] = 0.024300200167971865323442606364;
    w[60] = 0.024358557264690625853268520246;
    w[61] = 0.024402355633849582093297989694;
    w[62] = 0.02443156909785004505484856143;
    w[63] = 0.02444618019626251821132585261;
    w[64] = 0.02444618019626251821132585261;
    w[65] = 0.02443156909785004505484856143;
    w[66] = 0.024402355633849582093297989694;
    w[67] = 0.024358557264690625853268520246;
    w[68] = 0.024300200167971865323442606364;
    w[69] = 0.024227319222815248120093308442;
    w[70] = 0.024139957989019284997716653890;
    w[71] = 0.024038168681024052637587316820;
    w[72] = 0.023922012136703455672450408817;
    w[73] = 0.023791557781003400638780709885;
    w[74] = 0.023646883584447615143651392303;
    w[75] = 0.023488076016535913153025273282;
    w[76] = 0.023315229994062760122415671273;
    w[77] = 0.023128448824387027879297902403;
    w[78] = 0.022927844143686846920410987209;
    w[79] = 0.022713535850236461309712635923;
    w[80] = 0.022485652032744966871824603941;
    w[81] = 0.022244328893799765104629133607;
    w[82] = 0.021989710668460491434122106599;
    w[83] = 0.021721949538052075375260957768;
    w[84] = 0.021441205539208460137111853878;
    w[85] = 0.021147646468221348537019535180;
    w[86] = 0.020841447780751149113583948423;
    w[87] = 0.020522792486960069432284967788;
    w[88] = 0.020191871042130041180673158406;
    w[89] = 0.019848881232830862219944413265;
    w[90] = 0.019494028058706602823021918681;
    w[91] = 0.019127523609950945486518531668;
    w[92] = 0.018749586940544708650919548474;
    w[93] = 0.018360443937331343221289290991;
    w[94] = 0.017960327185008685940196927525;
    w[95] = 0.017549475827117704648706925634;
    w[96] = 0.017128135423111376830680987619;
    w[97] = 0.016696557801589204589091507954;
    w[98] = 0.016255000909785187051657456477;
    w[99] = 0.015803728659399346858965631687;
    w[100] = 0.015343010768865144085990853741;
    w[101] = 0.014873122602147314252385498520;
    w[102] = 0.014394345004166846176823892009;
    w[103] = 0.013906964132951985244288007396;
    w[104] = 0.013411271288616332314488951616;
    w[105] = 0.012907562739267347220442834004;
    w[106] = 0.012396139543950922968821728197;
    w[107] = 0.011877307372740279575891106926;
    w[108] = 0.011351376324080416693281668453;
    w[109] = 0.010818660739503076247659646277;
    w[110] = 0.010279479015832157133215340326;
    w[111] = 0.009734153415006805863548266094;
    w[112] = 0.009183009871660874334478743688;
    w[113] = 0.008626377798616749704978843782;
    w[114] = 0.008064589890486057972928598698;
    w[115] = 0.007497981925634728687671962688;
    w[116] = 0.006926892566898813563426670360;
    w[117] = 0.006351663161707188787214327826;
    w[118] = 0.0057726375428656985893346176261;
    w[119] = 0.0051901618326763302050707671348;
    w[120] = 0.0046045842567029551182905419803;
    w[121] = 0.0040162549837386423131943434863;
    w[122] = 0.0034255260409102157743377846601;
    w[123] = 0.0028327514714579910952857346468;
    w[124] = 0.0022382884309626187436220542727;
    w[125] = 0.0016425030186690295387908755948;
    w[126] = 0.0010458126793403487793128516001;
    w[127] = 0.00044938096029209037639429223999;
  }
  else if ( n == 129 )
  {
    x[0] = -0.9998275818477487191077441;
    x[1] = -0.9990916504696409986514389;
    x[2] = -0.9977681080525852721429460;
    x[3] = -0.9958574393142831982149111;
    x[4] = -0.9933607326210712814854011;
    x[5] = -0.9902794486488178389207689;
    x[6] = -0.9866153978313475022005761;
    x[7] = -0.9823707352517413115507418;
    x[8] = -0.9775479582993672474447814;
    x[9] = -0.9721499048427034297274163;
    x[10] = -0.9661797514202097197778763;
    x[11] = -0.9596410113101918904168119;
    x[12] = -0.9525375324342090471027732;
    x[13] = -0.9448734950776734726784764;
    x[14] = -0.9366534094216514605284616;
    x[15] = -0.9278821128840036204317296;
    x[16] = -0.9185647672698286252225115;
    x[17] = -0.9087068557320696331245539;
    x[18] = -0.8983141795436338850435985;
    x[19] = -0.8873928546826803665034968;
    x[20] = -0.8759493082329433892035217;
    x[21] = -0.8639902746011257878940216;
    x[22] = -0.8515227915535356930243826;
    x[23] = -0.8385541960742664442975407;
    x[24] = -0.8250921200473358809210133;
    x[25] = -0.8111444857653120742087717;
    x[26] = -0.7967195012670592680339606;
    x[27] = -0.7818256555073413245387500;
    x[28] = -0.7664717133611208816717785;
    x[29] = -0.7506667104654910227632368;
    x[30] = -0.7344199479022727047791516;
    x[31] = -0.7177409867244055767721220;
    x[32] = -0.7006396423293521790044710;
    x[33] = -0.6831259786828258512462248;
    x[34] = -0.6652103023962409818802202;
    x[35] = -0.6469031566613704719753373;
    x[36] = -0.6282153150457794374886895;
    x[37] = -0.6091577751526861909563306;
    x[38] = -0.5897417521489813916767844;
    x[39] = -0.5699786721652138894754096;
    x[40] = -0.5498801655714271702189358;
    x[41] = -0.5294580601328034000099406;
    x[42] = -0.5087243740491428186199463;
    x[43] = -0.4876913088822746111853066;
    x[44] = -0.4663712423755613514331869;
    x[45] = -0.4447767211697226217818454;
    x[46] = -0.4229204534192644388475065;
    x[47] = -0.4008153013138596117693121;
    x[48] = -0.3784742735090801012801265;
    x[49] = -0.3559105174709357969672656;
    x[50] = -0.3331373117387248575049982;
    x[51] = -0.3101680581107488341147318;
    x[52] = -0.2870162737574911929568755;
    x[53] = -0.2636955832669005409666949;
    x[54] = -0.2402197106264598167721148;
    x[55] = -0.2166024711467599103221439;
    x[56] = -0.1928577633313305998663880;
    x[57] = -0.1689995606975133227390302;
    x[58] = -0.1450419035531891084328306;
    x[59] = -0.1209988907342009817690539;
    x[60] = -0.0968846713073332753086909;
    x[61] = -0.0727134362437305599118207;
    x[62] = -0.0484994100676562986191764;
    x[63] = -0.0242568424855058415749954;
    x[64] = 0.0000000000000000000000000;
    x[65] = 0.0242568424855058415749954;
    x[66] = 0.0484994100676562986191764;
    x[67] = 0.0727134362437305599118207;
    x[68] = 0.0968846713073332753086909;
    x[69] = 0.1209988907342009817690539;
    x[70] = 0.1450419035531891084328306;
    x[71] = 0.1689995606975133227390302;
    x[72] = 0.1928577633313305998663880;
    x[73] = 0.2166024711467599103221439;
    x[74] = 0.2402197106264598167721148;
    x[75] = 0.2636955832669005409666949;
    x[76] = 0.2870162737574911929568755;
    x[77] = 0.3101680581107488341147318;
    x[78] = 0.3331373117387248575049982;
    x[79] = 0.3559105174709357969672656;
    x[80] = 0.3784742735090801012801265;
    x[81] = 0.4008153013138596117693121;
    x[82] = 0.4229204534192644388475065;
    x[83] = 0.4447767211697226217818454;
    x[84] = 0.4663712423755613514331869;
    x[85] = 0.4876913088822746111853066;
    x[86] = 0.5087243740491428186199463;
    x[87] = 0.5294580601328034000099406;
    x[88] = 0.5498801655714271702189358;
    x[89] = 0.5699786721652138894754096;
    x[90] = 0.5897417521489813916767844;
    x[91] = 0.6091577751526861909563306;
    x[92] = 0.6282153150457794374886895;
    x[93] = 0.6469031566613704719753373;
    x[94] = 0.6652103023962409818802202;
    x[95] = 0.6831259786828258512462248;
    x[96] = 0.7006396423293521790044710;
    x[97] = 0.7177409867244055767721220;
    x[98] = 0.7344199479022727047791516;
    x[99] = 0.7506667104654910227632368;
    x[100] = 0.7664717133611208816717785;
    x[101] = 0.7818256555073413245387500;
    x[102] = 0.7967195012670592680339606;
    x[103] = 0.8111444857653120742087717;
    x[104] = 0.8250921200473358809210133;
    x[105] = 0.8385541960742664442975407;
    x[106] = 0.8515227915535356930243826;
    x[107] = 0.8639902746011257878940216;
    x[108] = 0.875949308232943389203522;
    x[109] = 0.887392854682680366503497;
    x[110] = 0.898314179543633885043599;
    x[111] = 0.908706855732069633124554;
    x[112] = 0.918564767269828625222511;
    x[113] = 0.927882112884003620431730;
    x[114] = 0.936653409421651460528462;
    x[115] = 0.944873495077673472678476;
    x[116] = 0.952537532434209047102773;
    x[117] = 0.959641011310191890416812;
    x[118] = 0.966179751420209719777876;
    x[119] = 0.972149904842703429727416;
    x[120] = 0.977547958299367247444781;
    x[121] = 0.982370735251741311550742;
    x[122] = 0.986615397831347502200576;
    x[123] = 0.990279448648817838920769;
    x[124] = 0.993360732621071281485401;
    x[125] = 0.995857439314283198214911;
    x[126] = 0.997768108052585272142946;
    x[127] = 0.999091650469640998651439;
    x[128] = 0.999827581847748719107744;

    w[0] = 0.00044246794182939296923668005717;
    w[1] = 0.00102972844619622394463273519315;
    w[2] = 0.0016172530556785534682413679271;
    w[3] = 0.0022039015180966937075786419741;
    w[4] = 0.0027892681877797554940944677057;
    w[5] = 0.0033729979506246246117755709288;
    w[6] = 0.0039547444682113562172392974765;
    w[7] = 0.0045341644298525434513226874954;
    w[8] = 0.0051109164669246267289761565766;
    w[9] = 0.0056846609912469045788016012203;
    w[10] = 0.0062550602724461408889348709586;
    w[11] = 0.0068217785893519121070498527769;
    w[12] = 0.0073844824072454014447165055698;
    w[13] = 0.0079428405646668029041114107832;
    w[14] = 0.0084965244635723279730542832506;
    w[15] = 0.0090452082602137316404219313819;
    w[16] = 0.0095885690555104190787301294510;
    w[17] = 0.0101262870842733548093160774580;
    w[18] = 0.0106580459029055185304204093001;
    w[19] = 0.0111835325753305049735380697538;
    w[20] = 0.011702437856964778185746436834;
    w[21] = 0.012214456376582979416221105914;
    w[22] = 0.012719286815944623465099036330;
    w[23] = 0.013216632087061724231482387345;
    w[24] = 0.013706199506993971244060563234;
    w[25] = 0.014187700970062900419317230938;
    w[26] = 0.014660853117380060971041027493;
    w[27] = 0.015125377503587024690403432771;
    w[28] = 0.015581000760707523415881287558;
    w[29] = 0.016027454759014214436403950465;
    w[30] = 0.016464476764814667467169189640;
    w[31] = 0.016891809595063204177526208819;
    w[32] = 0.017309201768707240731293596444;
    w[33] = 0.017716407654678809269702031810;
    w[34] = 0.018113187616443980503999783812;
    w[35] = 0.018499308153024985727791918518;
    w[36] = 0.018874542036411948181617592169;
    w[37] = 0.019238668445283284085199492202;
    w[38] = 0.019591473094956024580283987216;
    w[39] = 0.019932748363489542089706675388;
    w[40] = 0.020262293413868438317104423081;
    w[41] = 0.020579914312192665948185517085;
    w[42] = 0.020885424141805311409990024684;
    w[43] = 0.021178643113290860912881038703;
    w[44] = 0.021459398670279205389981598196;
    w[45] = 0.021727525590993110687305178710;
    w[46] = 0.021982866085479386179554968899;
    w[47] = 0.022225269888466526554736910919;
    w[48] = 0.022454594347794176432066564511;
    w[49] = 0.022670704508362374313093970958;
    w[50] = 0.022873473191551169638592083492;
    w[51] = 0.023062781070063872924670495006;
    w[52] = 0.023238516738149892544490435771;
    w[53] = 0.023400576777165831146714346635;
    w[54] = 0.023548865816436258377269094263;
    w[55] = 0.023683296589378342897341543485;
    w[56] = 0.023803789984857314051325299744;
    w[57] = 0.023910275093742530302367230296;
    w[58] = 0.024002689250636756075547029720;
    w[59] = 0.024080978070754089272959634041;
    w[60] = 0.024145095481924836783843156014;
    w[61] = 0.024195003751708503129818111597;
    w[62] = 0.024230673509598936275508460625;
    w[63] = 0.024252083764308562906498864071;
    w[64] = 0.02425922191612154143202867472;
    w[65] = 0.024252083764308562906498864071;
    w[66] = 0.024230673509598936275508460625;
    w[67] = 0.024195003751708503129818111597;
    w[68] = 0.024145095481924836783843156014;
    w[69] = 0.024080978070754089272959634041;
    w[70] = 0.024002689250636756075547029720;
    w[71] = 0.023910275093742530302367230296;
    w[72] = 0.023803789984857314051325299744;
    w[73] = 0.023683296589378342897341543485;
    w[74] = 0.023548865816436258377269094263;
    w[75] = 0.023400576777165831146714346635;
    w[76] = 0.023238516738149892544490435771;
    w[77] = 0.023062781070063872924670495006;
    w[78] = 0.022873473191551169638592083492;
    w[79] = 0.022670704508362374313093970958;
    w[80] = 0.022454594347794176432066564511;
    w[81] = 0.022225269888466526554736910919;
    w[82] = 0.021982866085479386179554968899;
    w[83] = 0.021727525590993110687305178710;
    w[84] = 0.021459398670279205389981598196;
    w[85] = 0.021178643113290860912881038703;
    w[86] = 0.020885424141805311409990024684;
    w[87] = 0.020579914312192665948185517085;
    w[88] = 0.020262293413868438317104423081;
    w[89] = 0.019932748363489542089706675388;
    w[90] = 0.019591473094956024580283987216;
    w[91] = 0.019238668445283284085199492202;
    w[92] = 0.018874542036411948181617592169;
    w[93] = 0.018499308153024985727791918518;
    w[94] = 0.018113187616443980503999783812;
    w[95] = 0.017716407654678809269702031810;
    w[96] = 0.017309201768707240731293596444;
    w[97] = 0.016891809595063204177526208819;
    w[98] = 0.016464476764814667467169189640;
    w[99] = 0.016027454759014214436403950465;
    w[100] = 0.015581000760707523415881287558;
    w[101] = 0.015125377503587024690403432771;
    w[102] = 0.014660853117380060971041027493;
    w[103] = 0.014187700970062900419317230938;
    w[104] = 0.013706199506993971244060563234;
    w[105] = 0.013216632087061724231482387345;
    w[106] = 0.012719286815944623465099036330;
    w[107] = 0.012214456376582979416221105914;
    w[108] = 0.011702437856964778185746436834;
    w[109] = 0.0111835325753305049735380697538;
    w[110] = 0.0106580459029055185304204093001;
    w[111] = 0.0101262870842733548093160774580;
    w[112] = 0.0095885690555104190787301294510;
    w[113] = 0.0090452082602137316404219313819;
    w[114] = 0.0084965244635723279730542832506;
    w[115] = 0.0079428405646668029041114107832;
    w[116] = 0.0073844824072454014447165055698;
    w[117] = 0.0068217785893519121070498527769;
    w[118] = 0.0062550602724461408889348709586;
    w[119] = 0.0056846609912469045788016012203;
    w[120] = 0.0051109164669246267289761565766;
    w[121] = 0.0045341644298525434513226874954;
    w[122] = 0.0039547444682113562172392974765;
    w[123] = 0.0033729979506246246117755709288;
    w[124] = 0.0027892681877797554940944677057;
    w[125] = 0.0022039015180966937075786419741;
    w[126] = 0.0016172530556785534682413679271;
    w[127] = 0.00102972844619622394463273519315;
    w[128] = 0.00044246794182939296923668005717;
  }
  else if ( n == 255 )
  {
    x[0] = -0.999955705317563751730191;
    x[1] = -0.999766621312000569367063;
    x[2] = -0.999426474680169959344386;
    x[3] = -0.998935241284654635142155;
    x[4] = -0.998292986136967889228248;
    x[5] = -0.997499804126615814044844;
    x[6] = -0.996555814435198617028738;
    x[7] = -0.995461159480026294089975;
    x[8] = -0.994216004616630164799381;
    x[9] = -0.992820538021989138984811;
    x[10] = -0.991274970630385567164523;
    x[11] = -0.989579536085920123498574;
    x[12] = -0.987734490699732356281248;
    x[13] = -0.985740113407419277752900;
    x[14] = -0.983596705724776358640192;
    x[15] = -0.981304591701017185126565;
    x[16] = -0.978864117869068155239121;
    x[17] = -0.976275653192735980815246;
    x[18] = -0.973539589010643617645393;
    x[19] = -0.970656338976880365477697;
    x[20] = -0.967626338998338798105523;
    x[21] = -0.964450047168726298761719;
    x[22] = -0.961127943699247839572910;
    x[23] = -0.957660530845962076295490;
    x[24] = -0.954048332833816317950921;
    x[25] = -0.950291895777368285733522;
    x[26] = -0.946391787598204251752103;
    x[27] = -0.942348597939064408301480;
    x[28] = -0.938162938074687317626793;
    x[29] = -0.933835440819386124349338;
    x[30] = -0.929366760431369935739045;
    x[31] = -0.924757572513824425220425;
    x[32] = -0.920008573912766315142721;
    x[33] = -0.915120482611686961035103;
    x[34] = -0.910094037623000801254172;
    x[35] = -0.904929998876314959753358;
    x[36] = -0.899629147103536800144342;
    x[37] = -0.894192283720836729335637;
    x[38] = -0.888620230707484040924981;
    x[39] = -0.882913830481574073645470;
    x[40] = -0.877073945772665439532627;
    x[41] = -0.871101459491346550796200;
    x[42] = -0.864997274595751144137121;
    x[43] = -0.858762313955042966785823;
    x[44] = -0.852397520209890250084237;
    x[45] = -0.845903855629951054143931;
    x[46] = -0.839282301968391021084600;
    x[47] = -0.832533860313455524647230;
    x[48] = -0.825659550937118650611534;
    x[49] = -0.818660413140831885432406;
    x[50] = -0.811537505098395829833580;
    x[51] = -0.804291903695978689734633;
    x[52] = -0.796924704369305728807154;
    x[53] = -0.789437020938044295117764;
    x[54] = -0.781829985437409458675147;
    x[55] = -0.774104747947015717207115;
    x[56] = -0.766262476417000644100858;
    x[57] = -0.758304356491446765092016;
    x[58] = -0.750231591329128358931528;
    x[59] = -0.742045401421610281838045;
    x[60] = -0.733747024408726316001889;
    x[61] = -0.725337714891464938687812;
    x[62] = -0.716818744242290800531501;
    x[63] = -0.708191400412930589382399;
    x[64] = -0.699456987739652339456557;
    x[65] = -0.690616826746067624571761;
    x[66] = -0.681672253943486448787259;
    x[67] = -0.672624621628855017806731;
    x[68] = -0.663475297680306939970658;
    x[69] = -0.654225665350358766508700;
    x[70] = -0.644877123056781136890077;
    x[71] = -0.635431084171177146547142;
    x[72] = -0.625888976805299900901619;
    x[73] = -0.616252243595141561442344;
    x[74] = -0.606522341482826526536576;
    x[75] = -0.596700741496341721653202;
    x[76] = -0.586788928527137300685706;
    x[77] = -0.576788401105631382036211;
    x[78] = -0.566700671174652760010815;
    x[79] = -0.556527263860855843833077;
    x[80] = -0.546269717244142383159817;
    x[81] = -0.535929582125124840335150;
    x[82] = -0.525508421790666565699453;
    x[83] = -0.515007811777534223035005;
    x[84] = -0.504429339634198197635551;
    x[85] = -0.493774604680816999489812;
    x[86] = -0.483045217767441948626854;
    x[87] = -0.472242801030478698742627;
    x[88] = -0.461368987647442418771401;
    x[89] = -0.450425421590043710043279;
    x[90] = -0.439413757375642589040685;
    x[91] = -0.428335659817108112494341;
    x[92] = -0.417192803771121462605751;
    x[93] = -0.405986873884960545511889;
    x[94] = -0.394719564341804385683361;
    x[95] = -0.383392578604595822734854;
    x[96] = -0.372007629158501235092510;
    x[97] = -0.360566437252006227074021;
    x[98] = -0.349070732636686422161576;
    x[99] = -0.337522253305692705554261;
    x[100] = -0.325922745230990453444769;
    x[101] = -0.314273962099392474845918;
    x[102] = -0.302577665047425574167140;
    x[103] = -0.290835622395070819082047;
    x[104] = -0.279049609378417768508970;
    x[105] = -0.267221407881273079721012;
    x[106] = -0.255352806165764071686080;
    x[107] = -0.243445598601977973686482;
    x[108] = -0.231501585396677734059116;
    x[109] = -0.219522572321135403508985;
    x[110] = -0.207510370438124240859625;
    x[111] = -0.195466795828110816293869;
    x[112] = -0.183393669314688508087976;
    x[113] = -0.171292816189293903533225;
    x[114] = -0.159166065935247723154292;
    x[115] = -0.147015251951161989456661;
    x[116] = -0.134842211273755257250625;
    x[117] = -0.122648784300117812092492;
    x[118] = -0.110436814509468826540991;
    x[119] = -0.098208148184447540736015;
    x[120] = -0.085964634131980604256000;
    x[121] = -0.073708123403767780288977;
    x[122] = -0.061440469016428270850728;
    x[123] = -0.049163525671349973093019;
    x[124] = -0.036879149474284021657652;
    x[125] = -0.024589197654727010541405;
    x[126] = -0.012295528285133320036860;
    x[127] = 0.000000000000000000000000;
    x[128] = 0.012295528285133320036860;
    x[129] = 0.024589197654727010541405;
    x[130] = 0.036879149474284021657652;
    x[131] = 0.049163525671349973093019;
    x[132] = 0.061440469016428270850728;
    x[133] = 0.073708123403767780288977;
    x[134] = 0.085964634131980604256000;
    x[135] = 0.098208148184447540736015;
    x[136] = 0.110436814509468826540991;
    x[137] = 0.122648784300117812092492;
    x[138] = 0.134842211273755257250625;
    x[139] = 0.147015251951161989456661;
    x[140] = 0.159166065935247723154292;
    x[141] = 0.171292816189293903533225;
    x[142] = 0.183393669314688508087976;
    x[143] = 0.195466795828110816293869;
    x[144] = 0.207510370438124240859625;
    x[145] = 0.219522572321135403508985;
    x[146] = 0.231501585396677734059116;
    x[147] = 0.243445598601977973686482;
    x[148] = 0.255352806165764071686080;
    x[149] = 0.267221407881273079721012;
    x[150] = 0.279049609378417768508970;
    x[151] = 0.290835622395070819082047;
    x[152] = 0.302577665047425574167140;
    x[153] = 0.314273962099392474845918;
    x[154] = 0.325922745230990453444769;
    x[155] = 0.337522253305692705554261;
    x[156] = 0.349070732636686422161576;
    x[157] = 0.360566437252006227074021;
    x[158] = 0.372007629158501235092510;
    x[159] = 0.383392578604595822734854;
    x[160] = 0.394719564341804385683361;
    x[161] = 0.405986873884960545511889;
    x[162] = 0.417192803771121462605751;
    x[163] = 0.428335659817108112494341;
    x[164] = 0.439413757375642589040685;
    x[165] = 0.450425421590043710043279;
    x[166] = 0.461368987647442418771401;
    x[167] = 0.472242801030478698742627;
    x[168] = 0.483045217767441948626854;
    x[169] = 0.493774604680816999489812;
    x[170] = 0.504429339634198197635551;
    x[171] = 0.515007811777534223035005;
    x[172] = 0.525508421790666565699453;
    x[173] = 0.535929582125124840335150;
    x[174] = 0.546269717244142383159817;
    x[175] = 0.556527263860855843833077;
    x[176] = 0.566700671174652760010815;
    x[177] = 0.576788401105631382036211;
    x[178] = 0.586788928527137300685706;
    x[179] = 0.596700741496341721653202;
    x[180] = 0.606522341482826526536576;
    x[181] = 0.616252243595141561442344;
    x[182] = 0.625888976805299900901619;
    x[183] = 0.635431084171177146547142;
    x[184] = 0.644877123056781136890077;
    x[185] = 0.654225665350358766508700;
    x[186] = 0.663475297680306939970658;
    x[187] = 0.672624621628855017806731;
    x[188] = 0.681672253943486448787259;
    x[189] = 0.690616826746067624571761;
    x[190] = 0.699456987739652339456557;
    x[191] = 0.708191400412930589382399;
    x[192] = 0.716818744242290800531501;
    x[193] = 0.725337714891464938687812;
    x[194] = 0.733747024408726316001889;
    x[195] = 0.742045401421610281838045;
    x[196] = 0.750231591329128358931528;
    x[197] = 0.758304356491446765092016;
    x[198] = 0.766262476417000644100858;
    x[199] = 0.774104747947015717207115;
    x[200] = 0.781829985437409458675147;
    x[201] = 0.789437020938044295117764;
    x[202] = 0.796924704369305728807154;
    x[203] = 0.804291903695978689734633;
    x[204] = 0.811537505098395829833580;
    x[205] = 0.818660413140831885432406;
    x[206] = 0.825659550937118650611534;
    x[207] = 0.832533860313455524647230;
    x[208] = 0.839282301968391021084600;
    x[209] = 0.845903855629951054143931;
    x[210] = 0.852397520209890250084237;
    x[211] = 0.858762313955042966785823;
    x[212] = 0.864997274595751144137121;
    x[213] = 0.871101459491346550796200;
    x[214] = 0.877073945772665439532627;
    x[215] = 0.882913830481574073645470;
    x[216] = 0.888620230707484040924981;
    x[217] = 0.894192283720836729335637;
    x[218] = 0.899629147103536800144342;
    x[219] = 0.904929998876314959753358;
    x[220] = 0.910094037623000801254172;
    x[221] = 0.915120482611686961035103;
    x[222] = 0.920008573912766315142721;
    x[223] = 0.924757572513824425220425;
    x[224] = 0.929366760431369935739045;
    x[225] = 0.933835440819386124349338;
    x[226] = 0.938162938074687317626793;
    x[227] = 0.942348597939064408301480;
    x[228] = 0.946391787598204251752103;
    x[229] = 0.950291895777368285733522;
    x[230] = 0.954048332833816317950921;
    x[231] = 0.957660530845962076295490;
    x[232] = 0.961127943699247839572910;
    x[233] = 0.964450047168726298761719;
    x[234] = 0.967626338998338798105523;
    x[235] = 0.970656338976880365477697;
    x[236] = 0.973539589010643617645393;
    x[237] = 0.976275653192735980815246;
    x[238] = 0.978864117869068155239121;
    x[239] = 0.981304591701017185126565;
    x[240] = 0.983596705724776358640192;
    x[241] = 0.985740113407419277752900;
    x[242] = 0.987734490699732356281248;
    x[243] = 0.989579536085920123498574;
    x[244] = 0.991274970630385567164523;
    x[245] = 0.992820538021989138984811;
    x[246] = 0.994216004616630164799381;
    x[247] = 0.995461159480026294089975;
    x[248] = 0.996555814435198617028738;
    x[249] = 0.997499804126615814044844;
    x[250] = 0.998292986136967889228248;
    x[251] = 0.998935241284654635142155;
    x[252] = 0.999426474680169959344386;
    x[253] = 0.999766621312000569367063;
    x[254] = 0.999955705317563751730191;

    w[0] = 0.00011367361999142272115645954414;
    w[1] = 0.00026459387119083065532790838855;
    w[2] = 0.00041569762526823913616284210066;
    w[3] = 0.00056675794564824918946626058353;
    w[4] = 0.00071773647800611087798371518325;
    w[5] = 0.00086860766611945667949717690640;
    w[6] = 0.00101934797642732530281229369360;
    w[7] = 0.0011699343729388079886897709773;
    w[8] = 0.0013203439900221692090523602144;
    w[9] = 0.0014705540427783843160097204304;
    w[10] = 0.0016205417990415653896921100325;
    w[11] = 0.0017702845706603213070421243905;
    w[12] = 0.0019197597117132050055085980675;
    w[13] = 0.0020689446195015801533643667413;
    w[14] = 0.0022178167367540171700373764020;
    w[15] = 0.0023663535543962867157201855305;
    w[16] = 0.0025145326145997073931298921370;
    w[17] = 0.0026623315139717112732749157331;
    w[18] = 0.0028097279068204407457332299361;
    w[19] = 0.0029566995084575002760043344138;
    w[20] = 0.0031032240985191112621977893133;
    w[21] = 0.0032492795242943133198690930777;
    w[22] = 0.0033948437040533928255056951665;
    w[23] = 0.0035398946303722552150296713510;
    w[24] = 0.0036844103734499176530742235517;
    w[25] = 0.0038283690844171626400743524999;
    w[26] = 0.0039717489986349171988699773906;
    w[27] = 0.0041145284389812475901826468094;
    w[28] = 0.0042566858191260658425395494472;
    w[29] = 0.0043981996467927779838546384780;
    w[30] = 0.0045390485270061921259394035112;
    w[31] = 0.0046792111653260640506279893190;
    w[32] = 0.0048186663710656988918572043815;
    w[33] = 0.0049573930604950563104281084148;
    w[34] = 0.0050953702600278273039420404117;
    w[35] = 0.0052325771093919661294970523234;
    w[36] = 0.0053689928647831724787741258653;
    w[37] = 0.0055045969020008281904902120813;
    w[38] = 0.0056393687195659001929970994675;
    w[39] = 0.0057732879418203275712033691864;
    w[40] = 0.0059063343220074160130475409466;
    w[41] = 0.0060384877453327676663371666884;
    w[42] = 0.0061697282320052788060812561217;
    w[43] = 0.0063000359402577418025981070425;
    w[44] = 0.0064293911693465917826140832500;
    w[45] = 0.0065577743625303421548456356354;
    w[46] = 0.0066851661100262568757892743568;
    w[47] = 0.0068115471519448109954345674817;
    w[48] = 0.0069368983812014946719507501243;
    w[49] = 0.0070612008464055194979848418291;
    w[50] = 0.0071844357547249896530757997058;
    w[51] = 0.0073065844747281040972736443146;
    w[52] = 0.0074276285391999597581348419714;
    w[53] = 0.0075475496479345294426435656724;
    w[54] = 0.0076663296705013920315933272426;
    w[55] = 0.0077839506489867963897419914623;
    w[56] = 0.0079003948007086443529587296692;
    w[57] = 0.0080156445209049821352946484008;
    w[58] = 0.0081296823853955935356080649925;
    w[59] = 0.0082424911532162924158504385939;
    w[60] = 0.0083540537692255160718568405530;
    w[61] = 0.0084643533666828253227353760036;
    w[62] = 0.0085733732697989214067758505840;
    w[63] = 0.0086810969962567940901133439612;
    w[64] = 0.0087875082597036197689825483144;
    w[65] = 0.0088925909722130327769834298578;
    w[66] = 0.0089963292467173975949700110383;
    w[67] = 0.0090987073994097142025303711406;
    w[68] = 0.0091997099521147934060534414075;
    w[69] = 0.0092993216346293436285393234867;
    w[70] = 0.0093975273870306153500305317074;
    w[71] = 0.0094943123619532541442165010292;
    w[72] = 0.0095896619268340180657610209655;
    w[73] = 0.0096835616661240200035669970076;
    w[74] = 0.0097759973834681605268499842249;
    w[75] = 0.0098669551038514217128483481814;
    w[76] = 0.0099564210757116974565448593910;
    w[77] = 0.0100443817730188408231888789497;
    w[78] = 0.0101308238973196141129538950955;
    w[79] = 0.0102157343797482324629939488415;
    w[80] = 0.0102991003830021970147153502911;
    w[81] = 0.0103809093032831189224876935085;
    w[82] = 0.0104611487722022407735015844669;
    w[83] = 0.0105398066586503673262517188088;
    w[84] = 0.0106168710706319228563864391054;
    w[85] = 0.0106923303570628578226139809571;
    w[86] = 0.0107661731095321330311788312990;
    w[87] = 0.0108383881640265149842990798832;
    w[88] = 0.0109089646026184216450603134401;
    w[89] = 0.0109778917551165634377595759712;
    w[90] = 0.0110451592006791299277436662993;
    w[91] = 0.0111107567693892782875426356195;
    w[92] = 0.0111746745437926853557086684962;
    w[93] = 0.0112369028603969308303734810332;
    w[94] = 0.0112974323111324849102690558722;
    w[95] = 0.0113562537447750795009464486204;
    w[96] = 0.011413358268329247942299599697;
    w[97] = 0.011468737248372824084374355981;
    w[98] = 0.011522382312362197440930930031;
    w[99] = 0.011574285349898127083439539046;
    w[100] = 0.011624438513951922901227922331;
    w[101] = 0.011672834222051808845465154244;
    w[102] = 0.011719465157429288794653489478;
    w[103] = 0.011764324270125341726399410909;
    w[104] = 0.011807404778056278953532930501;
    w[105] = 0.011848700168039102281222824051;
    w[106] = 0.011888204196776208064673282076;
    w[107] = 0.011925910891799288293359117699;
    w[108] = 0.011961814552372285996633285380;
    w[109] = 0.011995909750353268455989686823;
    w[110] = 0.012028191331015087920350431142;
    w[111] = 0.012058654413824705751531083631;
    w[112] = 0.012087294393181062176578184854;
    w[113] = 0.012114106939111380091025793650;
    w[114] = 0.012139087997925797641334635250;
    w[115] = 0.012162233792830230614908682534;
    w[116] = 0.012183540824497371981177306326;
    w[117] = 0.012203005871595742256331865516;
    w[118] = 0.012220625991276710706457005806;
    w[119] = 0.012236398519619413758040249691;
    w[120] = 0.012250321072033503350218104906;
    w[121] = 0.012262391543619664338660618398;
    w[122] = 0.012272608109487846445745237751;
    w[123] = 0.012280969225033162644659793962;
    w[124] = 0.012287473626169412265336919908;
    w[125] = 0.012292120329520193516690694701;
    w[126] = 0.012294908632567576531532225710;
    w[127] = 0.01229583811375831445681490730;
    w[128] = 0.012294908632567576531532225710;
    w[129] = 0.012292120329520193516690694701;
    w[130] = 0.012287473626169412265336919908;
    w[131] = 0.012280969225033162644659793962;
    w[132] = 0.012272608109487846445745237751;
    w[133] = 0.012262391543619664338660618398;
    w[134] = 0.012250321072033503350218104906;
    w[135] = 0.012236398519619413758040249691;
    w[136] = 0.012220625991276710706457005806;
    w[137] = 0.012203005871595742256331865516;
    w[138] = 0.012183540824497371981177306326;
    w[139] = 0.012162233792830230614908682534;
    w[140] = 0.012139087997925797641334635250;
    w[141] = 0.012114106939111380091025793650;
    w[142] = 0.012087294393181062176578184854;
    w[143] = 0.012058654413824705751531083631;
    w[144] = 0.012028191331015087920350431142;
    w[145] = 0.011995909750353268455989686823;
    w[146] = 0.011961814552372285996633285380;
    w[147] = 0.011925910891799288293359117699;
    w[148] = 0.011888204196776208064673282076;
    w[149] = 0.011848700168039102281222824051;
    w[150] = 0.011807404778056278953532930501;
    w[151] = 0.011764324270125341726399410909;
    w[152] = 0.011719465157429288794653489478;
    w[153] = 0.011672834222051808845465154244;
    w[154] = 0.011624438513951922901227922331;
    w[155] = 0.011574285349898127083439539046;
    w[156] = 0.011522382312362197440930930031;
    w[157] = 0.011468737248372824084374355981;
    w[158] = 0.011413358268329247942299599697;
    w[159] = 0.0113562537447750795009464486204;
    w[160] = 0.0112974323111324849102690558722;
    w[161] = 0.0112369028603969308303734810332;
    w[162] = 0.0111746745437926853557086684962;
    w[163] = 0.0111107567693892782875426356195;
    w[164] = 0.0110451592006791299277436662993;
    w[165] = 0.0109778917551165634377595759712;
    w[166] = 0.0109089646026184216450603134401;
    w[167] = 0.0108383881640265149842990798832;
    w[168] = 0.0107661731095321330311788312990;
    w[169] = 0.0106923303570628578226139809571;
    w[170] = 0.0106168710706319228563864391054;
    w[171] = 0.0105398066586503673262517188088;
    w[172] = 0.0104611487722022407735015844669;
    w[173] = 0.0103809093032831189224876935085;
    w[174] = 0.0102991003830021970147153502911;
    w[175] = 0.0102157343797482324629939488415;
    w[176] = 0.0101308238973196141129538950955;
    w[177] = 0.0100443817730188408231888789497;
    w[178] = 0.0099564210757116974565448593910;
    w[179] = 0.0098669551038514217128483481814;
    w[180] = 0.0097759973834681605268499842249;
    w[181] = 0.0096835616661240200035669970076;
    w[182] = 0.0095896619268340180657610209655;
    w[183] = 0.0094943123619532541442165010292;
    w[184] = 0.0093975273870306153500305317074;
    w[185] = 0.0092993216346293436285393234867;
    w[186] = 0.0091997099521147934060534414075;
    w[187] = 0.0090987073994097142025303711406;
    w[188] = 0.0089963292467173975949700110383;
    w[189] = 0.0088925909722130327769834298578;
    w[190] = 0.0087875082597036197689825483144;
    w[191] = 0.0086810969962567940901133439612;
    w[192] = 0.0085733732697989214067758505840;
    w[193] = 0.0084643533666828253227353760036;
    w[194] = 0.0083540537692255160718568405530;
    w[195] = 0.0082424911532162924158504385939;
    w[196] = 0.0081296823853955935356080649925;
    w[197] = 0.0080156445209049821352946484008;
    w[198] = 0.0079003948007086443529587296692;
    w[199] = 0.0077839506489867963897419914623;
    w[200] = 0.0076663296705013920315933272426;
    w[201] = 0.0075475496479345294426435656724;
    w[202] = 0.0074276285391999597581348419714;
    w[203] = 0.0073065844747281040972736443146;
    w[204] = 0.0071844357547249896530757997058;
    w[205] = 0.0070612008464055194979848418291;
    w[206] = 0.0069368983812014946719507501243;
    w[207] = 0.0068115471519448109954345674817;
    w[208] = 0.0066851661100262568757892743568;
    w[209] = 0.0065577743625303421548456356354;
    w[210] = 0.0064293911693465917826140832500;
    w[211] = 0.0063000359402577418025981070425;
    w[212] = 0.0061697282320052788060812561217;
    w[213] = 0.0060384877453327676663371666884;
    w[214] = 0.0059063343220074160130475409466;
    w[215] = 0.0057732879418203275712033691864;
    w[216] = 0.0056393687195659001929970994675;
    w[217] = 0.0055045969020008281904902120813;
    w[218] = 0.0053689928647831724787741258653;
    w[219] = 0.0052325771093919661294970523234;
    w[220] = 0.0050953702600278273039420404117;
    w[221] = 0.0049573930604950563104281084148;
    w[222] = 0.0048186663710656988918572043815;
    w[223] = 0.0046792111653260640506279893190;
    w[224] = 0.0045390485270061921259394035112;
    w[225] = 0.0043981996467927779838546384780;
    w[226] = 0.0042566858191260658425395494472;
    w[227] = 0.0041145284389812475901826468094;
    w[228] = 0.0039717489986349171988699773906;
    w[229] = 0.0038283690844171626400743524999;
    w[230] = 0.0036844103734499176530742235517;
    w[231] = 0.0035398946303722552150296713510;
    w[232] = 0.0033948437040533928255056951665;
    w[233] = 0.0032492795242943133198690930777;
    w[234] = 0.0031032240985191112621977893133;
    w[235] = 0.0029566995084575002760043344138;
    w[236] = 0.0028097279068204407457332299361;
    w[237] = 0.0026623315139717112732749157331;
    w[238] = 0.0025145326145997073931298921370;
    w[239] = 0.0023663535543962867157201855305;
    w[240] = 0.0022178167367540171700373764020;
    w[241] = 0.0020689446195015801533643667413;
    w[242] = 0.0019197597117132050055085980675;
    w[243] = 0.0017702845706603213070421243905;
    w[244] = 0.0016205417990415653896921100325;
    w[245] = 0.0014705540427783843160097204304;
    w[246] = 0.0013203439900221692090523602144;
    w[247] = 0.0011699343729388079886897709773;
    w[248] = 0.00101934797642732530281229369360;
    w[249] = 0.00086860766611945667949717690640;
    w[250] = 0.00071773647800611087798371518325;
    w[251] = 0.00056675794564824918946626058353;
    w[252] = 0.00041569762526823913616284210066;
    w[253] = 0.00026459387119083065532790838855;
    w[254] = 0.00011367361999142272115645954414;
  }
  else if ( n == 256 )
  {
    x[0] = -0.999956050018992230734801;
    x[1] = -0.999768437409263186104879;
    x[2] = -0.999430937466261408240854;
    x[3] = -0.998943525843408856555026;
    x[4] = -0.998306266473006444055500;
    x[5] = -0.997519252756720827563409;
    x[6] = -0.996582602023381540430504;
    x[7] = -0.995496454481096356592647;
    x[8] = -0.994260972922409664962878;
    x[9] = -0.992876342608822117143534;
    x[10] = -0.991342771207583086922189;
    x[11] = -0.989660488745065218319244;
    x[12] = -0.987829747564860608916488;
    x[13] = -0.985850822286125956479245;
    x[14] = -0.983724009760315496166686;
    x[15] = -0.981449629025464405769303;
    x[16] = -0.979028021257622038824238;
    x[17] = -0.976459549719234155621011;
    x[18] = -0.973744599704370405266079;
    x[19] = -0.970883578480743029320923;
    x[20] = -0.967876915228489454909004;
    x[21] = -0.964725060975706430932612;
    x[22] = -0.961428488530732144006407;
    x[23] = -0.957987692411178129365790;
    x[24] = -0.954403188769716241764448;
    x[25] = -0.950675515316628276363852;
    x[26] = -0.946805231239127481372052;
    x[27] = -0.942792917117462443183076;
    x[28] = -0.938639174837814804981926;
    x[29] = -0.934344627502003094292477;
    x[30] = -0.929909919334005641180246;
    x[31] = -0.925335715583316202872730;
    x[32] = -0.920622702425146495505047;
    x[33] = -0.915771586857490384526670;
    x[34] = -0.910783096595065011890907;
    x[35] = -0.905657979960144647082682;
    x[36] = -0.900397005770303544771620;
    x[37] = -0.895000963223084577441223;
    x[38] = -0.889470661777610888828677;
    x[39] = -0.883806931033158284859826;
    x[40] = -0.878010620604706543986435;
    x[41] = -0.872082599995488289130046;
    x[42] = -0.866023758466554519297515;
    x[43] = -0.859835004903376350696173;
    x[44] = -0.853517267679502965073036;
    x[45] = -0.847071494517296207187072;
    x[46] = -0.840498652345762713895068;
    x[47] = -0.833799727155504894348444;
    x[48] = -0.826975723850812514289093;
    x[49] = -0.820027666098917067403478;
    x[50] = -0.812956596176431543136410;
    x[51] = -0.805763574812998623257389;
    x[52] = -0.798449681032170758782543;
    x[53] = -0.791016011989545994546707;
    x[54] = -0.783463682808183820750670;
    x[55] = -0.775793826411325739132053;
    x[56] = -0.768007593352445635975891;
    x[57] = -0.760106151642655454941907;
    x[58] = -0.752090686575492059587530;
    x[59] = -0.743962400549111568455683;
    x[60] = -0.735722512885917834620373;
    x[61] = -0.727372259649652126586894;
    x[62] = -0.718912893459971448372640;
    x[63] = -0.710345683304543313394566;
    x[64] = -0.701671914348685159406084;
    x[65] = -0.692892887742576960105342;
    x[66] = -0.684009920426075953124877;
    x[67] = -0.675024344931162763855919;
    x[68] = -0.665937509182048559906408;
    x[69] = -0.656750776292973221887500;
    x[70] = -0.647465524363724862617016;
    x[71] = -0.638083146272911368668689;
    x[72] = -0.628605049469014975432210;
    x[73] = -0.619032655759261219430968;
    x[74] = -0.609367401096333939522311;
    x[75] = -0.599610735362968321730388;
    x[76] = -0.589764122154454300785786;
    x[77] = -0.579829038559082944921832;
    x[78] = -0.569806974936568759057668;
    x[79] = -0.559699434694481145136907;
    x[80] = -0.549507934062718557042427;
    x[81] = -0.539234001866059181127936;
    x[82] = -0.528879179294822261951476;
    x[83] = -0.518445019673674476221662;
    x[84] = -0.507933088228616036231925;
    x[85] = -0.497344961852181477119512;
    x[86] = -0.486682228866890350103621;
    x[87] = -0.475946488786983306390738;
    x[88] = -0.465139352078479313645570;
    x[89] = -0.454262439917589998774455;
    x[90] = -0.443317383947527357216926;
    x[91] = -0.432305826033741309953441;
    x[92] = -0.421229418017623824976812;
    x[93] = -0.410089821468716550006434;
    x[94] = -0.398888707435459127713463;
    x[95] = -0.387627756194515583637985;
    x[96] = -0.376308656998716390283056;
    x[97] = -0.364933107823654018533465;
    x[98] = -0.353502815112969989537790;
    x[99] = -0.342019493522371636480730;
    x[100] = -0.330484865662416976229187;
    x[101] = -0.318900661840106275631683;
    x[102] = -0.307268619799319076258610;
    x[103] = -0.295590484460135614563787;
    x[104] = -0.283868007657081741799766;
    x[105] = -0.272102947876336609505245;
    x[106] = -0.260297069991942541978561;
    x[107] = -0.248452145001056666833243;
    x[108] = -0.236569949758284018477508;
    x[109] = -0.224652266709131967147878;
    x[110] = -0.212700883622625957937040;
    x[111] = -0.200717593323126670068001;
    x[112] = -0.188704193421388826461504;
    x[113] = -0.176662486044901997403722;
    x[114] = -0.164594277567553849829285;
    x[115] = -0.152501378338656395374607;
    x[116] = -0.140385602411375885913025;
    x[117] = -0.128248767270607094742050;
    x[118] = -0.116092693560332804940735;
    x[119] = -0.103919204810509403639197;
    x[120] = -0.091730127163519552031146;
    x[121] = -0.079527289100232965903227;
    x[122] = -0.067312521165716400242290;
    x[123] = -0.055087655694633984104561;
    x[124] = -0.042854526536379098381242;
    x[125] = -0.030614968779979029366279;
    x[126] = -0.018370818478813665117926;
    x[127] = -0.006123912375189529501170;
    x[128] = 0.006123912375189529501170;
    x[129] = 0.018370818478813665117926;
    x[130] = 0.030614968779979029366279;
    x[131] = 0.042854526536379098381242;
    x[132] = 0.055087655694633984104561;
    x[133] = 0.067312521165716400242290;
    x[134] = 0.079527289100232965903227;
    x[135] = 0.091730127163519552031146;
    x[136] = 0.103919204810509403639197;
    x[137] = 0.116092693560332804940735;
    x[138] = 0.128248767270607094742050;
    x[139] = 0.140385602411375885913025;
    x[140] = 0.152501378338656395374607;
    x[141] = 0.164594277567553849829285;
    x[142] = 0.176662486044901997403722;
    x[143] = 0.188704193421388826461504;
    x[144] = 0.200717593323126670068001;
    x[145] = 0.212700883622625957937040;
    x[146] = 0.224652266709131967147878;
    x[147] = 0.236569949758284018477508;
    x[148] = 0.248452145001056666833243;
    x[149] = 0.260297069991942541978561;
    x[150] = 0.272102947876336609505245;
    x[151] = 0.283868007657081741799766;
    x[152] = 0.295590484460135614563787;
    x[153] = 0.307268619799319076258610;
    x[154] = 0.318900661840106275631683;
    x[155] = 0.330484865662416976229187;
    x[156] = 0.342019493522371636480730;
    x[157] = 0.353502815112969989537790;
    x[158] = 0.364933107823654018533465;
    x[159] = 0.376308656998716390283056;
    x[160] = 0.387627756194515583637985;
    x[161] = 0.398888707435459127713463;
    x[162] = 0.410089821468716550006434;
    x[163] = 0.421229418017623824976812;
    x[164] = 0.432305826033741309953441;
    x[165] = 0.443317383947527357216926;
    x[166] = 0.454262439917589998774455;
    x[167] = 0.465139352078479313645570;
    x[168] = 0.475946488786983306390738;
    x[169] = 0.486682228866890350103621;
    x[170] = 0.497344961852181477119512;
    x[171] = 0.507933088228616036231925;
    x[172] = 0.518445019673674476221662;
    x[173] = 0.528879179294822261951476;
    x[174] = 0.539234001866059181127936;
    x[175] = 0.549507934062718557042427;
    x[176] = 0.559699434694481145136907;
    x[177] = 0.569806974936568759057668;
    x[178] = 0.579829038559082944921832;
    x[179] = 0.589764122154454300785786;
    x[180] = 0.599610735362968321730388;
    x[181] = 0.609367401096333939522311;
    x[182] = 0.619032655759261219430968;
    x[183] = 0.628605049469014975432210;
    x[184] = 0.638083146272911368668689;
    x[185] = 0.647465524363724862617016;
    x[186] = 0.656750776292973221887500;
    x[187] = 0.665937509182048559906408;
    x[188] = 0.675024344931162763855919;
    x[189] = 0.684009920426075953124877;
    x[190] = 0.692892887742576960105342;
    x[191] = 0.701671914348685159406084;
    x[192] = 0.710345683304543313394566;
    x[193] = 0.718912893459971448372640;
    x[194] = 0.727372259649652126586894;
    x[195] = 0.735722512885917834620373;
    x[196] = 0.743962400549111568455683;
    x[197] = 0.752090686575492059587530;
    x[198] = 0.760106151642655454941907;
    x[199] = 0.768007593352445635975891;
    x[200] = 0.775793826411325739132053;
    x[201] = 0.783463682808183820750670;
    x[202] = 0.791016011989545994546707;
    x[203] = 0.798449681032170758782543;
    x[204] = 0.805763574812998623257389;
    x[205] = 0.812956596176431543136410;
    x[206] = 0.820027666098917067403478;
    x[207] = 0.826975723850812514289093;
    x[208] = 0.833799727155504894348444;
    x[209] = 0.840498652345762713895068;
    x[210] = 0.847071494517296207187072;
    x[211] = 0.853517267679502965073036;
    x[212] = 0.859835004903376350696173;
    x[213] = 0.866023758466554519297515;
    x[214] = 0.872082599995488289130046;
    x[215] = 0.878010620604706543986435;
    x[216] = 0.883806931033158284859826;
    x[217] = 0.889470661777610888828677;
    x[218] = 0.895000963223084577441223;
    x[219] = 0.900397005770303544771620;
    x[220] = 0.905657979960144647082682;
    x[221] = 0.910783096595065011890907;
    x[222] = 0.915771586857490384526670;
    x[223] = 0.920622702425146495505047;
    x[224] = 0.925335715583316202872730;
    x[225] = 0.929909919334005641180246;
    x[226] = 0.934344627502003094292477;
    x[227] = 0.938639174837814804981926;
    x[228] = 0.942792917117462443183076;
    x[229] = 0.946805231239127481372052;
    x[230] = 0.950675515316628276363852;
    x[231] = 0.954403188769716241764448;
    x[232] = 0.957987692411178129365790;
    x[233] = 0.961428488530732144006407;
    x[234] = 0.964725060975706430932612;
    x[235] = 0.967876915228489454909004;
    x[236] = 0.970883578480743029320923;
    x[237] = 0.973744599704370405266079;
    x[238] = 0.976459549719234155621011;
    x[239] = 0.979028021257622038824238;
    x[240] = 0.981449629025464405769303;
    x[241] = 0.983724009760315496166686;
    x[242] = 0.985850822286125956479245;
    x[243] = 0.987829747564860608916488;
    x[244] = 0.989660488745065218319244;
    x[245] = 0.991342771207583086922189;
    x[246] = 0.992876342608822117143534;
    x[247] = 0.994260972922409664962878;
    x[248] = 0.995496454481096356592647;
    x[249] = 0.996582602023381540430504;
    x[250] = 0.997519252756720827563409;
    x[251] = 0.998306266473006444055500;
    x[252] = 0.998943525843408856555026;
    x[253] = 0.999430937466261408240854;
    x[254] = 0.999768437409263186104879;
    x[255] = 0.999956050018992230734801;

    w[0] = 0.00011278901782227217551253887725;
    w[1] = 0.00026253494429644590628745756250;
    w[2] = 0.00041246325442617632843218583774;
    w[3] = 0.00056234895403140980281523674759;
    w[4] = 0.0007121541634733206669089891511;
    w[5] = 0.0008618537014200890378140934163;
    w[6] = 0.0010114243932084404526058128414;
    w[7] = 0.0011608435575677247239705981135;
    w[8] = 0.0013100886819025044578316804271;
    w[9] = 0.0014591373333107332010883864996;
    w[10] = 0.0016079671307493272424499395690;
    w[11] = 0.0017565557363307299936069145295;
    w[12] = 0.0019048808534997184044191411746;
    w[13] = 0.0020529202279661431745487818492;
    w[14] = 0.0022006516498399104996848834189;
    w[15] = 0.0023480529563273120170064609087;
    w[16] = 0.0024951020347037068508395354372;
    w[17] = 0.0026417768254274905641208292516;
    w[18] = 0.0027880553253277068805747610763;
    w[19] = 0.0029339155908297166460123254142;
    w[20] = 0.0030793357411993375832053528316;
    w[21] = 0.0032242939617941981570107134269;
    w[22] = 0.0033687685073155510120191062489;
    w[23] = 0.0035127377050563073309710549844;
    w[24] = 0.0036561799581425021693892413052;
    w[25] = 0.0037990737487662579981170192082;
    w[26] = 0.0039413976414088336277290349840;
    w[27] = 0.0040831302860526684085997759212;
    w[28] = 0.0042242504213815362723565049060;
    w[29] = 0.0043647368779680566815684200621;
    w[30] = 0.0045045685814478970686417923159;
    w[31] = 0.0046437245556800603139790923525;
    w[32] = 0.0047821839258926913729317340448;
    w[33] = 0.0049199259218138656695587765655;
    w[34] = 0.0050569298807868423875578160762;
    w[35] = 0.0051931752508692809303287536296;
    w[36] = 0.0053286415939159303170811114788;
    w[37] = 0.0054633085886443102775705318566;
    w[38] = 0.0055971560336829100775514452572;
    w[39] = 0.005730163850601437177384417555;
    w[40] = 0.005862312086922653060661598801;
    w[41] = 0.005993580919115338221127696870;
    w[42] = 0.006123950655567932542389081187;
    w[43] = 0.006253401739542401272063645975;
    w[44] = 0.006381914752107880570375164275;
    w[45] = 0.006509470415053660267809899951;
    w[46] = 0.006636049593781065044590038355;
    w[47] = 0.006761633300173798780927861108;
    w[48] = 0.006886202695446320346713323775;
    w[49] = 0.007009739092969822621234436194;
    w[50] = 0.007132223961075390071672422986;
    w[51] = 0.007253638925833913783829137214;
    w[52] = 0.007373965773812346437572440695;
    w[53] = 0.007493186454805883358599761133;
    w[54] = 0.007611283084545659461618719618;
    w[55] = 0.007728237947381555631110194958;
    w[56] = 0.007844033498939711866810316151;
    w[57] = 0.007958652368754348353613161227;
    w[58] = 0.008072077362873499500946974804;
    w[59] = 0.008184291466438269935619761004;
    w[60] = 0.008295277846235225425171412553;
    w[61] = 0.008405019853221535756180301698;
    w[62] = 0.008513501025022490693838354790;
    w[63] = 0.008620705088401014305368838410;
    w[64] = 0.008726615961698807140336632217;
    w[65] = 0.008831217757248750025318272685;
    w[66] = 0.008934494783758207548408417085;
    w[67] = 0.009036431548662873680227775572;
    w[68] = 0.009137012760450806402000472219;
    w[69] = 0.009236223330956302687378716714;
    w[70] = 0.009334048377623269712466014486;
    w[71] = 0.009430473225737752747352764482;
    w[72] = 0.009525483410629284811829685754;
    w[73] = 0.009619064679840727857162164401;
    w[74] = 0.009711202995266279964249670496;
    w[75] = 0.009801884535257327825498800250;
    w[76] = 0.009891095696695828602630683809;
    w[77] = 0.009978823097034910124733949495;
    w[78] = 0.010065053576306383309460978930;
    w[79] = 0.010149774199094865654634066042;
    w[80] = 0.010232972256478219656954857160;
    w[81] = 0.010314635267934015068260713997;
    w[82] = 0.010394750983211728997101725205;
    w[83] = 0.010473307384170403003569566927;
    w[84] = 0.010550292686581481517533575536;
    w[85] = 0.010625695341896561133961681801;
    w[86] = 0.010699504038979785603048200583;
    w[87] = 0.010771707705804626636653631927;
    w[88] = 0.010842295511114795995293477058;
    w[89] = 0.010911256866049039700796847788;
    w[90] = 0.010978581425729570637988203448;
    w[91] = 0.011044259090813901263517571044;
    w[92] = 0.011108280009009843630460815451;
    w[93] = 0.011170634576553449462710881938;
    w[94] = 0.011231313439649668572656802083;
    w[95] = 0.011290307495875509508367594121;
    w[96] = 0.011347607895545491941625714297;
    w[97] = 0.011403206043039185964847059552;
    w[98] = 0.011457093598090639152334392298;
    w[99] = 0.011509262477039497958586392439;
    w[100] = 0.011559704854043635772668656950;
    w[101] = 0.011608413162253105722084706677;
    w[102] = 0.011655380094945242121298939730;
    w[103] = 0.011700598606620740288189823359;
    w[104] = 0.011744061914060550305376732759;
    w[105] = 0.011785763497343426181690117627;
    w[106] = 0.011825697100823977771160737958;
    w[107] = 0.011863856734071078731904572908;
    w[108] = 0.011900236672766489754287204237;
    w[109] = 0.011934831459563562255873201696;
    w[110] = 0.011967635904905893729007282670;
    w[111] = 0.011998645087805811934536710071;
    w[112] = 0.012027854356582571161267533498;
    w[113] = 0.012055259329560149814347085327;
    w[114] = 0.012080855895724544655975183976;
    w[115] = 0.012104640215340463097757829736;
    w[116] = 0.012126608720527321034718492205;
    w[117] = 0.012146758115794459815559837664;
    w[118] = 0.012165085378535502061307291839;
    w[119] = 0.012181587759481772174047585032;
    w[120] = 0.012196262783114713518180974196;
    w[121] = 0.012209108248037240407514094371;
    w[122] = 0.012220122227303969191708737227;
    w[123] = 0.012229303068710278904146266083;
    w[124] = 0.012236649395040158109242574767;
    w[125] = 0.012242160104272800769728083260;
    w[126] = 0.012245834369747920142463857550;
    w[127] = 0.01224767164028975590407032649;
    w[128] = 0.01224767164028975590407032649;
    w[129] = 0.012245834369747920142463857550;
    w[130] = 0.012242160104272800769728083260;
    w[131] = 0.012236649395040158109242574767;
    w[132] = 0.012229303068710278904146266083;
    w[133] = 0.012220122227303969191708737227;
    w[134] = 0.012209108248037240407514094371;
    w[135] = 0.012196262783114713518180974196;
    w[136] = 0.012181587759481772174047585032;
    w[137] = 0.012165085378535502061307291839;
    w[138] = 0.012146758115794459815559837664;
    w[139] = 0.012126608720527321034718492205;
    w[140] = 0.012104640215340463097757829736;
    w[141] = 0.012080855895724544655975183976;
    w[142] = 0.012055259329560149814347085327;
    w[143] = 0.012027854356582571161267533498;
    w[144] = 0.011998645087805811934536710071;
    w[145] = 0.011967635904905893729007282670;
    w[146] = 0.011934831459563562255873201696;
    w[147] = 0.011900236672766489754287204237;
    w[148] = 0.011863856734071078731904572908;
    w[149] = 0.011825697100823977771160737958;
    w[150] = 0.011785763497343426181690117627;
    w[151] = 0.011744061914060550305376732759;
    w[152] = 0.011700598606620740288189823359;
    w[153] = 0.011655380094945242121298939730;
    w[154] = 0.011608413162253105722084706677;
    w[155] = 0.011559704854043635772668656950;
    w[156] = 0.011509262477039497958586392439;
    w[157] = 0.011457093598090639152334392298;
    w[158] = 0.011403206043039185964847059552;
    w[159] = 0.011347607895545491941625714297;
    w[160] = 0.011290307495875509508367594121;
    w[161] = 0.011231313439649668572656802083;
    w[162] = 0.011170634576553449462710881938;
    w[163] = 0.011108280009009843630460815451;
    w[164] = 0.011044259090813901263517571044;
    w[165] = 0.010978581425729570637988203448;
    w[166] = 0.010911256866049039700796847788;
    w[167] = 0.010842295511114795995293477058;
    w[168] = 0.010771707705804626636653631927;
    w[169] = 0.010699504038979785603048200583;
    w[170] = 0.010625695341896561133961681801;
    w[171] = 0.010550292686581481517533575536;
    w[172] = 0.010473307384170403003569566927;
    w[173] = 0.010394750983211728997101725205;
    w[174] = 0.010314635267934015068260713997;
    w[175] = 0.010232972256478219656954857160;
    w[176] = 0.010149774199094865654634066042;
    w[177] = 0.010065053576306383309460978930;
    w[178] = 0.009978823097034910124733949495;
    w[179] = 0.009891095696695828602630683809;
    w[180] = 0.009801884535257327825498800250;
    w[181] = 0.009711202995266279964249670496;
    w[182] = 0.009619064679840727857162164401;
    w[183] = 0.009525483410629284811829685754;
    w[184] = 0.009430473225737752747352764482;
    w[185] = 0.009334048377623269712466014486;
    w[186] = 0.009236223330956302687378716714;
    w[187] = 0.009137012760450806402000472219;
    w[188] = 0.009036431548662873680227775572;
    w[189] = 0.008934494783758207548408417085;
    w[190] = 0.008831217757248750025318272685;
    w[191] = 0.008726615961698807140336632217;
    w[192] = 0.008620705088401014305368838410;
    w[193] = 0.008513501025022490693838354790;
    w[194] = 0.008405019853221535756180301698;
    w[195] = 0.008295277846235225425171412553;
    w[196] = 0.008184291466438269935619761004;
    w[197] = 0.008072077362873499500946974804;
    w[198] = 0.007958652368754348353613161227;
    w[199] = 0.007844033498939711866810316151;
    w[200] = 0.007728237947381555631110194958;
    w[201] = 0.007611283084545659461618719618;
    w[202] = 0.007493186454805883358599761133;
    w[203] = 0.007373965773812346437572440695;
    w[204] = 0.007253638925833913783829137214;
    w[205] = 0.007132223961075390071672422986;
    w[206] = 0.007009739092969822621234436194;
    w[207] = 0.006886202695446320346713323775;
    w[208] = 0.006761633300173798780927861108;
    w[209] = 0.006636049593781065044590038355;
    w[210] = 0.006509470415053660267809899951;
    w[211] = 0.006381914752107880570375164275;
    w[212] = 0.006253401739542401272063645975;
    w[213] = 0.006123950655567932542389081187;
    w[214] = 0.005993580919115338221127696870;
    w[215] = 0.005862312086922653060661598801;
    w[216] = 0.005730163850601437177384417555;
    w[217] = 0.0055971560336829100775514452572;
    w[218] = 0.0054633085886443102775705318566;
    w[219] = 0.0053286415939159303170811114788;
    w[220] = 0.0051931752508692809303287536296;
    w[221] = 0.0050569298807868423875578160762;
    w[222] = 0.0049199259218138656695587765655;
    w[223] = 0.0047821839258926913729317340448;
    w[224] = 0.0046437245556800603139790923525;
    w[225] = 0.0045045685814478970686417923159;
    w[226] = 0.0043647368779680566815684200621;
    w[227] = 0.0042242504213815362723565049060;
    w[228] = 0.0040831302860526684085997759212;
    w[229] = 0.0039413976414088336277290349840;
    w[230] = 0.0037990737487662579981170192082;
    w[231] = 0.0036561799581425021693892413052;
    w[232] = 0.0035127377050563073309710549844;
    w[233] = 0.0033687685073155510120191062489;
    w[234] = 0.0032242939617941981570107134269;
    w[235] = 0.0030793357411993375832053528316;
    w[236] = 0.0029339155908297166460123254142;
    w[237] = 0.0027880553253277068805747610763;
    w[238] = 0.0026417768254274905641208292516;
    w[239] = 0.0024951020347037068508395354372;
    w[240] = 0.0023480529563273120170064609087;
    w[241] = 0.0022006516498399104996848834189;
    w[242] = 0.0020529202279661431745487818492;
    w[243] = 0.0019048808534997184044191411746;
    w[244] = 0.0017565557363307299936069145295;
    w[245] = 0.0016079671307493272424499395690;
    w[246] = 0.0014591373333107332010883864996;
    w[247] = 0.0013100886819025044578316804271;
    w[248] = 0.0011608435575677247239705981135;
    w[249] = 0.0010114243932084404526058128414;
    w[250] = 0.0008618537014200890378140934163;
    w[251] = 0.0007121541634733206669089891511;
    w[252] = 0.00056234895403140980281523674759;
    w[253] = 0.00041246325442617632843218583774;
    w[254] = 0.00026253494429644590628745756250;
    w[255] = 0.00011278901782227217551253887725;
  }
  else if ( n == 257 )
  {
    x[0] = -0.999956390712330402472857;
    x[1] = -0.999770232390338019056053;
    x[2] = -0.999435348366365078441838;
    x[3] = -0.998951714093223210129834;
    x[4] = -0.998319392445383847808766;
    x[5] = -0.997538475365520218731818;
    x[6] = -0.996609078365487004512326;
    x[7] = -0.995531339486830143483750;
    x[8] = -0.994305419008553630362377;
    x[9] = -0.992931499332908653172844;
    x[10] = -0.991409784923101705201254;
    x[11] = -0.989740502257507526030375;
    x[12] = -0.987923899788618253106809;
    x[13] = -0.985960247902290665366669;
    x[14] = -0.983849838875444644048531;
    x[15] = -0.981592986831381877693095;
    x[16] = -0.979190027692327124191591;
    x[17] = -0.976641319128992592610888;
    x[18] = -0.973947240507062326750976;
    x[19] = -0.971108192830542793021113;
    x[20] = -0.968124598681952354372943;
    x[21] = -0.964996902159337170373447;
    x[22] = -0.961725568810109767190665;
    x[23] = -0.958311085561711847074814;
    x[24] = -0.954753960649106318830855;
    x[25] = -0.951054723539105826691801;
    x[26] = -0.947213924851546682950881;
    x[27] = -0.943232136277318328151464;
    x[28] = -0.939109950493259404355123;
    x[29] = -0.934847981073932324370129;
    x[30] = -0.930446862400288909805510;
    x[31] = -0.925907249565240289235888;
    x[32] = -0.921229818276144817520964;
    x[33] = -0.916415264754228313295468;
    x[34] = -0.911464305630951423630955;
    x[35] = -0.906377677841339419411308;
    x[36] = -0.901156138514290206476301;
    x[37] = -0.895800464859876809085345;
    x[38] = -0.890311454053661045810287;
    x[39] = -0.884689923118035575018750;
    x[40] = -0.878936708800611938658765;
    x[41] = -0.873052667449672679799858;
    x[42] = -0.867038674886706051812473;
    x[43] = -0.860895626276042275514686;
    x[44] = -0.854624435991610735314055;
    x[45] = -0.848226037480837936478636;
    x[46] = -0.841701383125706473284556;
    x[47] = -0.835051444100995681967937;
    x[48] = -0.828277210229725073186687;
    x[49] = -0.821379689835822056081139;
    x[50] = -0.814359909594035880004229;
    x[51] = -0.807218914377120130552073;
    x[52] = -0.799957767100306523636066;
    x[53] = -0.792577548563093144962574;
    x[54] = -0.785079357288370682385816;
    x[55] = -0.777464309358910595129671;
    x[56] = -0.769733538251239556788216;
    x[57] = -0.761888194666924898264210;
    x[58] = -0.753929446361296162339238;
    x[59] = -0.745858477969628263337895;
    x[60] = -0.737676490830812123299244;
    x[61] = -0.729384702808539030149808;
    x[62] = -0.720984348110025333531072;
    x[63] = -0.712476677102304460118510;
    x[64] = -0.703862956126113592426171;
    x[65] = -0.695144467307402713168813;
    x[66] = -0.686322508366494071200553;
    x[67] = -0.677398392424920474813593;
    x[68] = -0.668373447809971163711735;
    x[69] = -0.659249017856974352220492;
    x[70] = -0.650026460709345873208532;
    x[71] = -0.640707149116433684724434;
    x[72] = -0.631292470229188329449219;
    x[73] = -0.621783825393689760680446;
    x[74] = -0.612182629942561267650033;
    x[75] = -0.602490312984301547488097;
    x[76] = -0.592708317190566281032495;
    x[77] = -0.582838098581430874902446;
    x[78] = -0.572881126308666332759406;
    x[79] = -0.562838882437060514424546;
    x[80] = -0.552712861723817332466074;
    x[81] = -0.542504571396066721967792;
    x[82] = -0.532215530926518500400434;
    x[83] = -0.521847271807293510797499;
    x[84] = -0.511401337321965712746629;
    x[85] = -0.500879282315849152005553;
    x[86] = -0.490282672964564000798817;
    x[87] = -0.479613086540916117008992;
    x[88] = -0.468872111180124821505728;
    x[89] = -0.458061345643433838720630;
    x[90] = -0.447182399080140586238810;
    x[91] = -0.436236890788079234603398;
    x[92] = -0.425226449972593188682213;
    x[93] = -0.414152715504032866791986;
    x[94] = -0.403017335673814873281489;
    x[95] = -0.391821967949078874408131;
    x[96] = -0.380568278725978696070941;
    x[97] = -0.369257943081644365255611;
    x[98] = -0.357892644524852014873858;
    x[99] = -0.346474074745438764010632;
    x[100] = -0.335003933362499872399782;
    x[101] = -0.323483927671405649204085;
    x[102] = -0.311915772389675771851948;
    x[103] = -0.300301189401748840754520;
    x[104] = -0.288641907502685160168097;
    x[105] = -0.276939662140840894253032;
    x[106] = -0.265196195159551900488370;
    x[107] = -0.253413254537865690008131;
    x[108] = -0.241592594130360106108882;
    x[109] = -0.229735973406087448117604;
    x[110] = -0.217845157186682897983880;
    x[111] = -0.205921915383676231351599;
    x[112] = -0.193968022735045913454182;
    x[113] = -0.181985258541054792946197;
    x[114] = -0.169975406399406713716337;
    x[115] = -0.157940253939763465806087;
    x[116] = -0.145881592557661591770148;
    x[117] = -0.133801217147868654144405;
    x[118] = -0.121700925837218653121859;
    x[119] = -0.109582519716966361063898;
    x[120] = -0.097447802574700412082119;
    x[121] = -0.085298580625855050603929;
    x[122] = -0.073136662244860502573600;
    x[123] = -0.060963857695971986730406;
    x[124] = -0.048781978863817431238958;
    x[125] = -0.036592838983704002816750;
    x[126] = -0.024398252371723591403953;
    x[127] = -0.012200034154697423345412;
    x[128] = 0.000000000000000000000000;
    x[129] = 0.012200034154697423345412;
    x[130] = 0.024398252371723591403953;
    x[131] = 0.036592838983704002816750;
    x[132] = 0.048781978863817431238958;
    x[133] = 0.060963857695971986730406;
    x[134] = 0.073136662244860502573600;
    x[135] = 0.085298580625855050603929;
    x[136] = 0.097447802574700412082119;
    x[137] = 0.109582519716966361063898;
    x[138] = 0.121700925837218653121859;
    x[139] = 0.133801217147868654144405;
    x[140] = 0.145881592557661591770148;
    x[141] = 0.157940253939763465806087;
    x[142] = 0.169975406399406713716337;
    x[143] = 0.181985258541054792946197;
    x[144] = 0.193968022735045913454182;
    x[145] = 0.205921915383676231351599;
    x[146] = 0.217845157186682897983880;
    x[147] = 0.229735973406087448117604;
    x[148] = 0.241592594130360106108882;
    x[149] = 0.253413254537865690008131;
    x[150] = 0.265196195159551900488370;
    x[151] = 0.276939662140840894253032;
    x[152] = 0.288641907502685160168097;
    x[153] = 0.300301189401748840754520;
    x[154] = 0.311915772389675771851948;
    x[155] = 0.323483927671405649204085;
    x[156] = 0.335003933362499872399782;
    x[157] = 0.346474074745438764010632;
    x[158] = 0.357892644524852014873858;
    x[159] = 0.369257943081644365255611;
    x[160] = 0.380568278725978696070941;
    x[161] = 0.391821967949078874408131;
    x[162] = 0.403017335673814873281489;
    x[163] = 0.414152715504032866791986;
    x[164] = 0.425226449972593188682213;
    x[165] = 0.436236890788079234603398;
    x[166] = 0.447182399080140586238810;
    x[167] = 0.458061345643433838720630;
    x[168] = 0.468872111180124821505728;
    x[169] = 0.479613086540916117008992;
    x[170] = 0.490282672964564000798817;
    x[171] = 0.500879282315849152005553;
    x[172] = 0.511401337321965712746629;
    x[173] = 0.521847271807293510797499;
    x[174] = 0.532215530926518500400434;
    x[175] = 0.542504571396066721967792;
    x[176] = 0.552712861723817332466074;
    x[177] = 0.562838882437060514424546;
    x[178] = 0.572881126308666332759406;
    x[179] = 0.582838098581430874902446;
    x[180] = 0.592708317190566281032495;
    x[181] = 0.602490312984301547488097;
    x[182] = 0.612182629942561267650033;
    x[183] = 0.621783825393689760680446;
    x[184] = 0.631292470229188329449219;
    x[185] = 0.640707149116433684724434;
    x[186] = 0.650026460709345873208532;
    x[187] = 0.659249017856974352220492;
    x[188] = 0.668373447809971163711735;
    x[189] = 0.677398392424920474813593;
    x[190] = 0.686322508366494071200553;
    x[191] = 0.695144467307402713168813;
    x[192] = 0.703862956126113592426171;
    x[193] = 0.712476677102304460118510;
    x[194] = 0.720984348110025333531072;
    x[195] = 0.729384702808539030149808;
    x[196] = 0.737676490830812123299244;
    x[197] = 0.745858477969628263337895;
    x[198] = 0.753929446361296162339238;
    x[199] = 0.761888194666924898264210;
    x[200] = 0.769733538251239556788216;
    x[201] = 0.777464309358910595129671;
    x[202] = 0.785079357288370682385816;
    x[203] = 0.792577548563093144962574;
    x[204] = 0.799957767100306523636066;
    x[205] = 0.807218914377120130552073;
    x[206] = 0.814359909594035880004229;
    x[207] = 0.821379689835822056081139;
    x[208] = 0.828277210229725073186687;
    x[209] = 0.835051444100995681967937;
    x[210] = 0.841701383125706473284556;
    x[211] = 0.848226037480837936478636;
    x[212] = 0.854624435991610735314055;
    x[213] = 0.860895626276042275514686;
    x[214] = 0.867038674886706051812473;
    x[215] = 0.873052667449672679799858;
    x[216] = 0.878936708800611938658765;
    x[217] = 0.884689923118035575018750;
    x[218] = 0.890311454053661045810287;
    x[219] = 0.895800464859876809085345;
    x[220] = 0.901156138514290206476301;
    x[221] = 0.906377677841339419411308;
    x[222] = 0.911464305630951423630955;
    x[223] = 0.916415264754228313295468;
    x[224] = 0.921229818276144817520964;
    x[225] = 0.925907249565240289235888;
    x[226] = 0.930446862400288909805510;
    x[227] = 0.934847981073932324370129;
    x[228] = 0.939109950493259404355123;
    x[229] = 0.943232136277318328151464;
    x[230] = 0.947213924851546682950881;
    x[231] = 0.951054723539105826691801;
    x[232] = 0.954753960649106318830855;
    x[233] = 0.958311085561711847074814;
    x[234] = 0.961725568810109767190665;
    x[235] = 0.964996902159337170373447;
    x[236] = 0.968124598681952354372943;
    x[237] = 0.971108192830542793021113;
    x[238] = 0.973947240507062326750976;
    x[239] = 0.976641319128992592610888;
    x[240] = 0.979190027692327124191591;
    x[241] = 0.981592986831381877693095;
    x[242] = 0.983849838875444644048531;
    x[243] = 0.985960247902290665366669;
    x[244] = 0.987923899788618253106809;
    x[245] = 0.989740502257507526030375;
    x[246] = 0.991409784923101705201254;
    x[247] = 0.992931499332908653172844;
    x[248] = 0.994305419008553630362377;
    x[249] = 0.995531339486830143483750;
    x[250] = 0.996609078365487004512326;
    x[251] = 0.997538475365520218731818;
    x[252] = 0.998319392445383847808766;
    x[253] = 0.998951714093223210129834;
    x[254] = 0.999435348366365078441838;
    x[255] = 0.999770232390338019056053;
    x[256] = 0.999956390712330402472857;

    w[0] = 0.00011191470145601756450862287886;
    w[1] = 0.00026049995580176964436806680831;
    w[2] = 0.00040926648283531339591138751432;
    w[3] = 0.00055799120546880640169677292533;
    w[4] = 0.00070663671051592291949335494247;
    w[5] = 0.00085517818446696565626595950963;
    w[6] = 0.00100359280467969441299468763292;
    w[7] = 0.0011518582377826677880963146741;
    w[8] = 0.0012999523174235227389668643832;
    w[9] = 0.0014478529559255120065233994722;
    w[10] = 0.0015955381166175133369701690235;
    w[11] = 0.0017429858051468299509941139300;
    w[12] = 0.0018901740676190104269878470891;
    w[13] = 0.0020370809914723626741694800322;
    w[14] = 0.0021836847075455253317921866057;
    w[15] = 0.0023299633927021828561308282641;
    w[16] = 0.0024758952727301488651840215879;
    w[17] = 0.0026214586253808109266552781372;
    w[18] = 0.0027666317834818283552560256501;
    w[19] = 0.0029113931380877846359302447381;
    w[20] = 0.0030557211416493711130936102459;
    w[21] = 0.0031995943111899437356540290142;
    w[22] = 0.0033429912314827618499065991316;
    w[23] = 0.0034858905582247143702551557840;
    w[24] = 0.0036282710212037760873102463983;
    w[25] = 0.0037701114274582873548537007645;
    w[26] = 0.0039113906644266662571543468015;
    w[27] = 0.0040520877030864825223229951262;
    w[28] = 0.0041921816010820254766367595011;
    w[29] = 0.0043316515058396297504806208252;
    w[30] = 0.0044704766576701092218388764046;
    w[31] = 0.0046086363928577081326523656522;
    w[32] = 0.0047461101467350184936945641585;
    w[33] = 0.0048828774567433411142588306018;
    w[34] = 0.0050189179654779878773297516544;
    w[35] = 0.0051542114237180378340642003713;
    w[36] = 0.0052887376934400710240953933529;
    w[37] = 0.0054224767508154127788846727083;
    w[38] = 0.0055554086891904284012033890901;
    w[39] = 0.0056875137220494140577838938236;
    w[40] = 0.0058187721859596348346566361185;
    w[41] = 0.0059491645434980654366600347567;
    w[42] = 0.0060786713861593931405204596709;
    w[43] = 0.0062072734372448464599330978665;
    w[44] = 0.0063349515547314166407936938524;
    w[45] = 0.0064616867341210426397202932350;
    w[46] = 0.0065874601112693336961737372300;
    w[47] = 0.0067122529651934070221351960200;
    w[48] = 0.0068360467208584215286561508406;
    w[49] = 0.0069588229519423919043121805236;
    w[50] = 0.0070805633835788707705149901066;
    w[51] = 0.0072012498950770900730828552207;
    w[52] = 0.0073208645226191563361371026044;
    w[53] = 0.0074393894619338979090297315972;
    w[54] = 0.0075568070709469658838993300454;
    w[55] = 0.0076730998724067939537782250476;
    w[56] = 0.0077882505564860261212726654404;
    w[57] = 0.0079022419833580248574070864277;
    w[58] = 0.0080150571857480760504667455353;
    w[59] = 0.0081266793714589108764118189068;
    w[60] = 0.0082370919258701685661946145361;
    w[61] = 0.0083462784144114279413811886655;
    w[62] = 0.0084542225850084395379670551258;
    w[63] = 0.0085609083705021941391459209280;
    w[64] = 0.0086663198910404675908861979240;
    w[65] = 0.0087704414564414858792445834744;
    w[66] = 0.0088732575685293586050755892934;
    w[67] = 0.0089747529234409331997949023068;
    w[68] = 0.0090749124139037264846862498962;
    w[69] = 0.0091737211314845944854270065178;
    w[70] = 0.0092711643688088057725325917169;
    w[71] = 0.0093672276217491880067391857021;
    w[72] = 0.0094618965915850218253881576301;
    w[73] = 0.0095551571871303607110514249099;
    w[74] = 0.0096469955268314600363329731559;
    w[75] = 0.0097373979408330030783691793250;
    w[76] = 0.0098263509730128164423854701706;
    w[77] = 0.0099138413829847720250916955489;
    w[78] = 0.0099998561480695773850435626986;
    w[79] = 0.0100843824652331611676814627839;
    w[80] = 0.0101674077529923650568895461852;
    w[81] = 0.0102489196532876585918958554047;
    w[82] = 0.0103289060333225980974485876288;
    w[83] = 0.0104073549873697559257355517893;
    w[84] = 0.0104842548385428511997370260353;
    w[85] = 0.0105595941405348182788823332058;
    w[86] = 0.0106333616793215542382761147904;
    w[87] = 0.0107055464748310917616231511294;
    w[88] = 0.0107761377825779489945556541150;
    w[89] = 0.0108451250952624130885928632830;
    w[90] = 0.0109124981443345193856719616965;
    w[91] = 0.0109782469015224934483083029166;
    w[92] = 0.0110423615803254284301924654946;
    w[93] = 0.0111048326374699756056269264803;
    w[94] = 0.0111656507743308312328559850485;
    w[95] = 0.0112248069383148083152535688671;
    w[96] = 0.0112822923242082872447042603128;
    w[97] = 0.0113380983754878447625379269120;
    w[98] = 0.011392216785593866154247619654;
    w[99] = 0.011444639499166951104119199270;
    w[100] = 0.011495358713246929174010288914;
    w[101] = 0.011544366878434306436012137033;
    w[102] = 0.011591656700013970380783131035;
    w[103] = 0.011637221139040985841125311445;
    w[104] = 0.011681053413388320313049670635;
    w[105] = 0.011723146998756342723302879656;
    w[106] = 0.011763495629643945382264331878;
    w[107] = 0.011802093300281144573421477037;
    w[108] = 0.011838934265523020964443424791;
    w[109] = 0.011874013041704866779344562066;
    w[110] = 0.011907324407458412445505183140;
    w[111] = 0.011938863404489011222535627643;
    w[112] = 0.011968625338313666131272065445;
    w[113] = 0.011996605778959789329711050159;
    w[114] = 0.012022800561624589927558893338;
    w[115] = 0.012047205787294992091420946532;
    w[116] = 0.012069817823327991167612855626;
    w[117] = 0.012090633303991361438266420912;
    w[118] = 0.012109649130964635027950450318;
    w[119] = 0.012126862473800277391553601370;
    w[120] = 0.012142270770344990738801546574;
    w[121] = 0.012155871727121082685623083829;
    w[122] = 0.012167663319667843366755737416;
    w[123] = 0.012177643792842880196606249581;
    w[124] = 0.012185811661083365425569178819;
    w[125] = 0.012192165708627157605870499188;
    w[126] = 0.012196704989693764053654538465;
    w[127] = 0.012199428828625117371582840212;
    w[128] = 0.01220033681998614507777289232;
    w[129] = 0.012199428828625117371582840212;
    w[130] = 0.012196704989693764053654538465;
    w[131] = 0.012192165708627157605870499188;
    w[132] = 0.012185811661083365425569178819;
    w[133] = 0.012177643792842880196606249581;
    w[134] = 0.012167663319667843366755737416;
    w[135] = 0.012155871727121082685623083829;
    w[136] = 0.012142270770344990738801546574;
    w[137] = 0.012126862473800277391553601370;
    w[138] = 0.012109649130964635027950450318;
    w[139] = 0.012090633303991361438266420912;
    w[140] = 0.012069817823327991167612855626;
    w[141] = 0.012047205787294992091420946532;
    w[142] = 0.012022800561624589927558893338;
    w[143] = 0.011996605778959789329711050159;
    w[144] = 0.011968625338313666131272065445;
    w[145] = 0.011938863404489011222535627643;
    w[146] = 0.011907324407458412445505183140;
    w[147] = 0.011874013041704866779344562066;
    w[148] = 0.011838934265523020964443424791;
    w[149] = 0.011802093300281144573421477037;
    w[150] = 0.011763495629643945382264331878;
    w[151] = 0.011723146998756342723302879656;
    w[152] = 0.011681053413388320313049670635;
    w[153] = 0.011637221139040985841125311445;
    w[154] = 0.011591656700013970380783131035;
    w[155] = 0.011544366878434306436012137033;
    w[156] = 0.011495358713246929174010288914;
    w[157] = 0.011444639499166951104119199270;
    w[158] = 0.011392216785593866154247619654;
    w[159] = 0.0113380983754878447625379269120;
    w[160] = 0.0112822923242082872447042603128;
    w[161] = 0.0112248069383148083152535688671;
    w[162] = 0.0111656507743308312328559850485;
    w[163] = 0.0111048326374699756056269264803;
    w[164] = 0.0110423615803254284301924654946;
    w[165] = 0.0109782469015224934483083029166;
    w[166] = 0.0109124981443345193856719616965;
    w[167] = 0.0108451250952624130885928632830;
    w[168] = 0.0107761377825779489945556541150;
    w[169] = 0.0107055464748310917616231511294;
    w[170] = 0.0106333616793215542382761147904;
    w[171] = 0.0105595941405348182788823332058;
    w[172] = 0.0104842548385428511997370260353;
    w[173] = 0.0104073549873697559257355517893;
    w[174] = 0.0103289060333225980974485876288;
    w[175] = 0.0102489196532876585918958554047;
    w[176] = 0.0101674077529923650568895461852;
    w[177] = 0.0100843824652331611676814627839;
    w[178] = 0.0099998561480695773850435626986;
    w[179] = 0.0099138413829847720250916955489;
    w[180] = 0.0098263509730128164423854701706;
    w[181] = 0.0097373979408330030783691793250;
    w[182] = 0.0096469955268314600363329731559;
    w[183] = 0.0095551571871303607110514249099;
    w[184] = 0.0094618965915850218253881576301;
    w[185] = 0.0093672276217491880067391857021;
    w[186] = 0.0092711643688088057725325917169;
    w[187] = 0.0091737211314845944854270065178;
    w[188] = 0.0090749124139037264846862498962;
    w[189] = 0.0089747529234409331997949023068;
    w[190] = 0.0088732575685293586050755892934;
    w[191] = 0.0087704414564414858792445834744;
    w[192] = 0.0086663198910404675908861979240;
    w[193] = 0.0085609083705021941391459209280;
    w[194] = 0.0084542225850084395379670551258;
    w[195] = 0.0083462784144114279413811886655;
    w[196] = 0.0082370919258701685661946145361;
    w[197] = 0.0081266793714589108764118189068;
    w[198] = 0.0080150571857480760504667455353;
    w[199] = 0.0079022419833580248574070864277;
    w[200] = 0.0077882505564860261212726654404;
    w[201] = 0.0076730998724067939537782250476;
    w[202] = 0.0075568070709469658838993300454;
    w[203] = 0.0074393894619338979090297315972;
    w[204] = 0.0073208645226191563361371026044;
    w[205] = 0.0072012498950770900730828552207;
    w[206] = 0.0070805633835788707705149901066;
    w[207] = 0.0069588229519423919043121805236;
    w[208] = 0.0068360467208584215286561508406;
    w[209] = 0.0067122529651934070221351960200;
    w[210] = 0.0065874601112693336961737372300;
    w[211] = 0.0064616867341210426397202932350;
    w[212] = 0.0063349515547314166407936938524;
    w[213] = 0.0062072734372448464599330978665;
    w[214] = 0.0060786713861593931405204596709;
    w[215] = 0.0059491645434980654366600347567;
    w[216] = 0.0058187721859596348346566361185;
    w[217] = 0.0056875137220494140577838938236;
    w[218] = 0.0055554086891904284012033890901;
    w[219] = 0.0054224767508154127788846727083;
    w[220] = 0.0052887376934400710240953933529;
    w[221] = 0.0051542114237180378340642003713;
    w[222] = 0.0050189179654779878773297516544;
    w[223] = 0.0048828774567433411142588306018;
    w[224] = 0.0047461101467350184936945641585;
    w[225] = 0.0046086363928577081326523656522;
    w[226] = 0.0044704766576701092218388764046;
    w[227] = 0.0043316515058396297504806208252;
    w[228] = 0.0041921816010820254766367595011;
    w[229] = 0.0040520877030864825223229951262;
    w[230] = 0.0039113906644266662571543468015;
    w[231] = 0.0037701114274582873548537007645;
    w[232] = 0.0036282710212037760873102463983;
    w[233] = 0.0034858905582247143702551557840;
    w[234] = 0.0033429912314827618499065991316;
    w[235] = 0.0031995943111899437356540290142;
    w[236] = 0.0030557211416493711130936102459;
    w[237] = 0.0029113931380877846359302447381;
    w[238] = 0.0027666317834818283552560256501;
    w[239] = 0.0026214586253808109266552781372;
    w[240] = 0.0024758952727301488651840215879;
    w[241] = 0.0023299633927021828561308282641;
    w[242] = 0.0021836847075455253317921866057;
    w[243] = 0.0020370809914723626741694800322;
    w[244] = 0.0018901740676190104269878470891;
    w[245] = 0.0017429858051468299509941139300;
    w[246] = 0.0015955381166175133369701690235;
    w[247] = 0.0014478529559255120065233994722;
    w[248] = 0.0012999523174235227389668643832;
    w[249] = 0.0011518582377826677880963146741;
    w[250] = 0.00100359280467969441299468763292;
    w[251] = 0.00085517818446696565626595950963;
    w[252] = 0.00070663671051592291949335494247;
    w[253] = 0.00055799120546880640169677292533;
    w[254] = 0.00040926648283531339591138751432;
    w[255] = 0.00026049995580176964436806680831;
    w[256] = 0.00011191470145601756450862287886;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of N = %d\n", n );
    fprintf ( stderr, "  Legal values are 1:33, 63/64/65, 127/128/129, 255/256/257\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

double *map ( char *code, int element_order )

/******************************************************************************/
/*
  Purpose: 

    MAP returns the interpolation matrix for any available element.

  Formula:

    For an element of order N, we suppose we are given N items of data 
    Q associated with the nodes.

   Let PHI(J)(R,S) be the Lagrange basis polynomial associated with 
   node J.  PHI(J)(R,S) is 1 at node J, and 0 at each of the other nodes.

   Let P(R,S) be the polynomial of N terms which interpolates the
   data Q, that is,

      P(R(J),S(J)) = Q(J)

   where the coordinates of node J are (R(J),S(J)).  Then we know
   that we can write

     P(R,S) = sum ( 1 <= J <= N ) Q(J) * PHI(J)(R,S)

   But P(R,S) also has a standard representation as

     P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I)

   where REXP(I) and SEXP(I) are the exponents of R and S and
   the A(I) are the appropriate coefficients.

   The interpolation matrix W allows us to immediately compute
   the standard basis coefficients A from the data Q to be interpolated
   using the formula:

      A(I) = sum ( 1 <= J <= N ) W(I,J) * Q(J)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.

    Input, int N, the order associated with the code.

    Output, double MAP[N*N], the interpolation matrix.
*/
{
  double area;
  int i;
  int info;
  int j;
  int *pivot;
  double *r;
  int *rexp;
  double rfact;
  double *s;
  int *sexp;
  double sfact;
  double *v;
  double *w;

  pivot = ( int * ) malloc ( element_order * sizeof ( int ) );
  r = ( double * ) malloc ( element_order * sizeof ( double ) );
  rexp = ( int * ) malloc ( element_order * sizeof ( int ) );
  s = ( double * ) malloc ( element_order * sizeof ( double ) );
  sexp = ( int * ) malloc ( element_order * sizeof ( int ) );
  v = ( double * ) malloc ( element_order * element_order * sizeof ( double ) );
/*
  Get the (R,S) location of the nodes.
*/
  node_reference ( code, r, s, &area );
/*
  Get the associated monomials.
*/
  poly ( code, rexp, sexp );
/*
  Set up the Vandermonde matrix.
  Factors of the form 0**0 are to be understood as 1.
*/
  for ( i = 0; i < element_order; i++ )
  {
    for ( j = 0; j < element_order; j++ )
    {
      if ( rexp[j] == 0 )
      {
        rfact = 1.0;
      }
      else
      {
        rfact = pow ( r[i], rexp[j] );
      }

      if ( sexp[j] == 0 )
      {
        sfact = 1.0;
      }
      else
      {
        sfact = pow ( s[i], sexp[j] );
      }
      v[i+j*element_order] = rfact * sfact;
    }
  }
/*
  Factor the Vandermonde matrix.
*/
  info = r8ge_fa ( element_order, v, pivot );

  if ( info != 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MAP - Fatal error!\n" );
    fprintf ( stderr, "  The Vandermonde matrix is singular.\n" );
    exit ( 1 );
  }
/*
  Invert the Vandermonde matrix.
*/
  w = r8ge_inverse ( element_order, v, pivot );

  free ( pivot );
  free ( r );
  free ( rexp );
  free ( s );
  free ( sexp );
  free ( v );

  return w;
}
/******************************************************************************/

void map_test ( char *code )

/******************************************************************************/
/*
  Purpose: 

    MAP_TEST tests the map routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, the code for the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.
*/
{
  int element_order;
  double *w;

  printf ( "\n" );
  printf ( "  MAP_TEST: The interpolation matrix for element %s\n", code );

  element_order = order_code ( code );

  w = map ( code, element_order );

  r8mat_print ( element_order, element_order, w, 
    "  The interpolation matrix:" );

  free ( w );

  return;
}
/******************************************************************************/

double *mass_matrix_t3 ( int node_num, int element_num, int element_node[], 
  double node_xy[] )

/******************************************************************************/
/*
  Purpose: 

    MASS_MATRIX_T3 computes the mass matrix, using 3-node triangles.

  Discussion:

    The mass matrix to be estimated has the form:

      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region

    where PHI(I) and PHI(J) are the shape functions associated with
    the I-th and J-th variables.

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the nodes that make up each element.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Output, double MASS_MATRIX_T3[NODE_NUM*NODE_NUM], the mass matrix.
*/
{
# define ELEMENT_ORDER 3

  double *a;
  double area;
  double dwdr[ELEMENT_ORDER];
  double dwds[ELEMENT_ORDER];
  int element;
  int global1;
  int global2;
  int local1;
  int local2;
  int norder;
  int p1;
  int p2;
  int p3;
  int quad;
  int quad_num;
  double r;
  double *rtab;
  int rule;
  double s;
  double *stab;
  double w[ELEMENT_ORDER];
  double *weight;

  a = ( double * ) malloc ( node_num * node_num * sizeof ( double ) );
/*
  Zero out the matrix.
*/
  for ( global2 = 0; global2 < node_num; global2++ )
  {
    for ( global1 = 0; global1 < node_num; global1++ )
    {
      a[global1+global2*node_num] = 0.0;
    }
  }
/*
  Get the weights and abscissas for a unit triangle.
*/
  rule = 4;
  quad_num = triangle_unit_size ( rule );

  rtab = ( double * ) malloc ( quad_num * sizeof ( double ) );
  stab = ( double * ) malloc ( quad_num * sizeof ( double ) );
  weight = ( double * ) malloc ( quad_num * sizeof ( double ) );

  triangle_unit_set ( rule, quad_num, rtab, stab, weight );
/*
  For each element.
*/
  for ( element = 0; element < element_num; element++ )
  {
    p1 = element_node[0+element*3] - 1;
    p2 = element_node[1+element*3] - 1;
    p3 = element_node[2+element*3] - 1;

    area = 0.5 * r8_abs ( 
        node_xy[0+p1*2] * ( node_xy[1+p2*2] - node_xy[1+p3*2] ) 
      + node_xy[0+p2*2] * ( node_xy[1+p3*2] - node_xy[1+p1*2] ) 
      + node_xy[0+p3*2] * ( node_xy[1+p1*2] - node_xy[1+p2*2] ) );

    if ( area == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "MASS_MATRIX_T3 - Fatal error!\n" );
      fprintf ( stderr, "  Zero area for element %d\n", element );
      fprintf ( stderr, "  Node 1 = %d\n", p1 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p1*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p1*2] );
      fprintf ( stderr, "  Node 2 = %d\n", p2 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p2*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p2*2] );
      fprintf ( stderr, "  Node 3 = %d\n", p3 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p3*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p3*2] );
      exit ( 1 );
    }
/*
  For each quadrature point in the element...
*/
    for ( quad = 0; quad < quad_num; quad++ )
    {
      r = rtab[quad];
      s = stab[quad];

      shape_t3 ( r, s, w, dwdr, dwds );
/*
  For each basis function PHI(I) associated with a node in the element,
*/
      for ( local1 = 0; local1 < 3; local1++ )
      {
        global1 = element_node[local1+element*3] - 1;
/*
  For each "neighbor" basis function PHI(J) associated with a node in
  the element.
*/
        for ( local2 = 0; local2 < 3; local2++ )
        {
          global2 = element_node[local2+element*3] - 1;

          a[global1+global2*node_num] = a[global1+global2*node_num] 
            + area * weight[quad] * w[local1] * w[local2];
        }
      }
    }
  }

  free ( rtab );
  free ( stab );
  free ( weight );

  return a;
# undef ELEMENT_ORDER
}
/******************************************************************************/

double *mass_matrix_t6 ( int node_num, int element_num, int element_node[], 
  double node_xy[] )

/******************************************************************************/
/*
  Purpose: 

    MASS_MATRIX_T6 computes the mass matrix, using 6-node triangles.

  Discussion:

    The mass matrix to be estimated has the form:

      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region

    where PHI(I) and PHI(J) are the shape functions associated with
    the I-th and J-th variables.

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2007

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that make up each element.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Output, double MASS_MATRIX_T6[NODE_NUM*NODE_NUM], the mass matrix.
*/
{
# define ELEMENT_ORDER 6

  double *a;
  double area;
  double dwdr[ELEMENT_ORDER];
  double dwds[ELEMENT_ORDER];
  int element;
  int global1;
  int global2;
  int local1;
  int local2;
  int norder;
  int p1;
  int p2;
  int p3;
  int quad;
  int quad_num;
  double r;
  double *rtab;
  int rule;
  double s;
  double *stab;
  double w[ELEMENT_ORDER];
  double *weight;

  a = ( double * ) malloc ( node_num * node_num * sizeof ( double ) );
/*
  Zero out the matrix.
*/
  for ( global2 = 0; global2 < node_num; global2++ )
  {
    for ( global1 = 0; global1 < node_num; global1++ )
    {
      a[global1+global2*node_num] = 0.0;
    }
  }
/*
  Get the weights and abscissas for a unit triangle.
*/
  rule = 12;
  quad_num = triangle_unit_size ( rule );

  rtab = ( double * ) malloc ( quad_num * sizeof ( double ) );
  stab = ( double * ) malloc ( quad_num * sizeof ( double ) );
  weight = ( double * ) malloc ( quad_num * sizeof ( double ) );

  triangle_unit_set ( rule, quad_num, rtab, stab, weight );
/*
  For each element.
*/
  for ( element = 0; element < element_num; element++ )
  {
    p1 = element_node[0+element*6] - 1;
    p2 = element_node[1+element*6] - 1;
    p3 = element_node[2+element*6] - 1;

    area = 0.5 * r8_abs ( 
        node_xy[0+p1*2] * ( node_xy[1+p2*2] - node_xy[1+p3*2] ) 
      + node_xy[0+p2*2] * ( node_xy[1+p3*2] - node_xy[1+p1*2] ) 
      + node_xy[0+p3*2] * ( node_xy[1+p1*2] - node_xy[1+p2*2] ) );

    if ( area == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "MASS_MATRIX_T6 - Fatal error!\n" );
      fprintf ( stderr, "  Zero area for element %d\n", element );
      fprintf ( stderr, "  Node 1 = %d\n", p1 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p1*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p1*2] );
      fprintf ( stderr, "  Node 2 = %d\n", p2 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p2*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p2*2] );
      fprintf ( stderr, "  Node 3 = %d\n", p3 );
      fprintf ( stderr, "  X = %g\n", node_xy[0+p3*2] );
      fprintf ( stderr, "  Y = %g\n", node_xy[1+p3*2] );
      exit ( 1 );
    }
/*
  For each quadrature point in the element...
*/
    for ( quad = 0; quad < quad_num; quad++ )
    {
      r = rtab[quad];
      s = stab[quad];

      shape_t6 ( r, s, w, dwdr, dwds );
/*
  For each basis function PHI(I) associated with a node in the element,
*/
      for ( local1 = 0; local1 < 6; local1++ )
      {
        global1 = element_node[local1+element*6] - 1;
/*
  For each "neighbor" basis function PHI(J) associated with a node in
  the element.
*/
        for ( local2 = 0; local2 < 6; local2++ )
        {
          global2 = element_node[local2+element*6] - 1;

          a[global1+global2*node_num] = a[global1+global2*node_num] 
            + area * weight[quad] * w[local1] * w[local2];
        }
      }
    }
  }

  free ( rtab );
  free ( stab );
  free ( weight );

  return a;
# undef ELEMENT_ORDER
}
/******************************************************************************/

int next_boundary_node ( int node, char *code )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE returns the next boundary node in any element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Input, char *CODE, identifies the element desired.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
    "T3", "T4", "T6" and "T10".

    Output, int NEXT_BOUNDARY_NODE, the index of the next edge node.
*/
{
  int value;

  if ( s_eqi ( code, "Q4" ) )
  {
    value = next_boundary_node_q4 ( node );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    value = next_boundary_node_q8 ( node );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    value = next_boundary_node_q9 ( node );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    value = next_boundary_node_q12 ( node );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    value = next_boundary_node_q16 ( node );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    value = next_boundary_node_ql ( node );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    value = next_boundary_node_t3 ( node );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    value = next_boundary_node_t4 ( node );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    value = next_boundary_node_t6 ( node );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    value = next_boundary_node_t10 ( node );
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_q4 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_Q4 returns the next boundary node in a Q4 element.

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_Q4, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_q8 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_Q8 returns the next boundary node in a Q8 element.

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_Q8, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_q9 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_Q9 returns the next boundary node in a Q9 element.

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_Q9, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else
  {
    value = 1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_q12 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_Q12 returns the next boundary node in a Q12 element.

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_Q12, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 12;
  }
  else if ( node == 12 )
  {
    value = 11;
  }
  else if ( node == 11 )
  {
    value = 10;
  }
  else if ( node == 10 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_q16 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_Q16 returns the next boundary node in a Q16 element.

  Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_Q16, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 12;
  }
  else if ( node == 12 )
  {
    value = 16;
  }
  else if ( node == 16 )
  {
    value = 15;
  }
  else if ( node == 15 )
  {
    value = 14;
  }
  else if ( node == 14 )
  {
    value = 13;
  }
  else if ( node == 13 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 5;
  }
  else if ( node == 5 ) 
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_ql ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_QL returns the next boundary node in a QL element.

  Element QL:

    |
    1  4---5---6
    |  |       |
    |  |       |
    S  |       |
    |  |       |
    |  |       |
    0  1---2---3
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_QL, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_t3 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_T3 returns the next boundary node in a T3 element.

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_T3, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_t4 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_T4 returns the next boundary node in a T4 element.

  Element T4:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  | 4 \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_T4, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_t6 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_T6 returns the next boundary node in a T6 element.

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_T6, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
/******************************************************************************/

int next_boundary_node_t10 ( int node )

/******************************************************************************/
/*
  Purpose: 

    NEXT_BOUNDARY_NODE_T10 returns the next boundary node in a T10 element.

  Element T10:

    |
    1  10
    |  |\
    |  | \
    |  8  9
    |  |   \
    S  |    \
    |  5  6  7
    |  |      \
    |  |       \
    0  1--2--3--4
    |
    +--0----R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NODE, the index of the current node.  An input
    value of 0 (or any "unusual" value") indicates that the
    first edge node is desired.

    Output, int NEXT_BOUNDARY_NODE_T10, the index of the next edge node.
*/
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 10;
  }
  else if ( node == 10 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 5;
  }
  else
  {
    value = 1;
  }

  return value;
}
/******************************************************************************/

void node_reference ( char *code, double r[], double s[], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE returns the basis nodes for any available element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.

    Output, double R[N], S[N], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  if ( s_eqi ( code, "Q4" ) )
  {
    node_reference_q4 ( r, s, area );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    node_reference_q8 ( r, s, area );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    node_reference_q9 ( r, s, area );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    node_reference_q12 ( r, s, area );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    node_reference_q16 ( r, s, area );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    node_reference_ql ( r, s, area );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    node_reference_t3 ( r, s, area );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    node_reference_t4 ( r, s, area );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    node_reference_t6 ( r, s, area );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    node_reference_t10 ( r, s, area );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "NODE_REFERENCE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = %s\n", code );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void node_reference_q4 ( double r[4], double s[4], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_Q4 returns the basis nodes for a 4 node quadrilateral.

  Element Q4:

    |
    1  4-------3
    |  |       |
    |  |       |
    S  |       |
    |  |       |
    |  |       |
    0  1-------2
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[4], S[4], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_q8 ( double r[8], double s[8], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_Q8 returns the basis nodes for an 8 node quadrilateral.

  Comment:

    This element is known as the quadratic "serendipity" element.

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[8], S[8], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0; 
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;
  r[6] = 0.5;
  r[7] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;
  s[4] = 0.0;
  s[5] = 0.5;
  s[6] = 1.0;
  s[7] = 0.5;

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_q9 ( double r[9], double s[9], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_Q9 returns the basis nodes for a 9 node quadrilateral.

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R[9], S[9], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;
  r[6] = 0.5;
  r[7] = 0.0;
  r[8] = 0.5;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;
  s[4] = 0.0;
  s[5] = 0.5;
  s[6] = 1.0;
  s[7] = 0.5;
  s[8] = 0.5;

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_q12 ( double r[12], double s[12], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_Q12 returns the basis nodes for a 12 node quadrilateral.

  Discussion:

    This element is known as the cubic "serendipity" element.

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[12], S[12], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  double a = 0.0;
  double b = 1.0 / 3.0;
  double c = 2.0 / 3.0;
  double d = 1.0;

  r[0] = a;
  r[1] = b;
  r[2] = c;
  r[3] = d;
  r[4] = a;
  r[5] = d;
  r[6] = a;
  r[7] = d;
  r[8] = a;
  r[9] = b;
  r[10] = c;
  r[11] = d;

  s[0] = a;
  s[1] = a;
  s[2] = a;
  s[3] = a;
  s[4] = b;
  s[5] = b;
  s[6] = c;
  s[7] = c;
  s[8] = d;
  s[9] = d;
  s[10] = d;
  s[11] = d;

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_q16 ( double r[16], double s[16], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_Q16 returns the basis nodes for a 16 node quadrilateral.

  Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[16], S[16], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  int i;
  int j;
  int k;

  k = 0;
  for ( i = 0; i <= 3; i++ )
  {
    for ( j = 0; j <= 3; j++ )
    {
      r[k] = ( double ) ( j ) / 3.0;
      s[k] = ( double ) ( i ) / 3.0;
      k = k + 1;
    }
  }

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_ql ( double r[6], double s[6], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_QL returns the basis nodes for a quadratic/linear.

  Element QL:

    |
    1  4---5---6
    |  |       |
    |  |       |
    S  |       |
    |  |       |
    |  |       |
    0  1---2---3
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[6], S[6], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 0.5;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 0.0;
  s[3] = 1.0;
  s[4] = 1.0;
  s[5] = 1.0;

  *area = 1.0;

  return;
}
/******************************************************************************/

void node_reference_t3 ( double r[3], double s[3], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_T3 returns the basis nodes for the 3 node triangle.

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[3], S[3], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;

  *area = 0.5;

  return;
}
/******************************************************************************/

void node_reference_t4 ( double r[4], double s[4], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_T4 returns the basis nodes for the 4 node triangle.

  Element T4:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  | 4 \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[4], S[4], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;
  r[3] = 1.0 / 3.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0 / 3.0;

  *area = 0.5;

  return;
}
/******************************************************************************/

void node_reference_t6 ( double r[6], double s[6], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_T6 returns the basis nodes for a 6 node triangle.

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[6], S[6], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;
  r[3] = 0.5;
  r[4] = 0.5;
  r[5] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 0.0;
  s[4] = 0.5;
  s[5] = 0.5;

  *area = 0.5;

  return;
}
/******************************************************************************/

void node_reference_t10 ( double r[10], double s[10], double *area )

/******************************************************************************/
/*
  Purpose: 

    NODE_REFERENCE_T10 returns the basis nodes for a 10 node triangle.

  Element T10:

    |
    1  10
    |  |\
    |  | \
    |  8  9
    |  |   \
    S  |    \
    |  5  6  7
    |  |      \
    |  |       \
    0  1--2--3--4
    |
    +--0----R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, double R[10], S[10], the coordinates of the basis nodes.

    Output, double *AREA, the area of the element.
*/
{
  r[0] = 0.0;
  s[0] = 0.0;

  r[1] = 1.0 / 3.0;
  s[1] = 0.0;

  r[2] = 2.0 / 3.0;
  s[2] = 0.0;

  r[3] = 1.0;
  s[3] = 0.0;

  r[4] = 0.0;
  s[4] = 1.0 / 3.0;

  r[5] = 1.0 / 3.0;
  s[5] = 1.0 / 3.0;

  r[6] = 2.0 / 3.0;
  s[6] = 1.0 / 3.0;

  r[7] = 0.0;
  s[7] = 2.0 / 3.0;

  r[8] = 1.0 / 3.0;
  s[8] = 2.0 / 3.0;

  r[9] = 0.0;
  s[9] = 1.0;

  *area = 0.5;

  return;
}
/******************************************************************************/

int ns_t6_var_count ( int element_num, int element_node[], int node_num, 
  int var_node[] )

/******************************************************************************/
/*
  Purpose:

    NS_T6_VAR_COUNT counts the Navier Stokes variables on a T6 grid.

  Discussion:

    We are given a mesh of T6 elements, and asked to count, in advance,
    the number of Navier-Stokes variables associated with the grid.
    In particular, every node has two velocity variables associated with
    it, but only a node that is a vertex of the element will also have
    an associated pressure variable.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Input, int NODE_NUM, the number of nodes.

    Output, int VAR_NODE[NODE_NUM+1], used to find the variables 
    associated with a given node, which are in VAR in locations 
    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
    this array points to the location just after the last location in VAR.

    Output, int NS_T6_VAR_COUNT, the number of variables.
*/
{
  int count;
  int element;
  int element_order = 6;
  int node;
  int num;
  int order;
  int var_num;
/*
  Our job is easy once we determine which nodes are vertices.
  So to begin with, let VAR_NODE count the number of variables
  associated with each node.
*/
  for ( node = 0; node < node_num; node++ )
  {
    var_node[node] = 2;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < 3; order++ )
    {
      node = element_node[order+element*element_order];
      var_node[node-1] = 3;
    }
  }
/*
  Count them.
*/
  var_num = 0;
  for ( node = 0; node < node_num; node++ )
  {
    var_num = var_num + var_node[node];
  }
/*
  Make pointers.
*/
  count = 1;

  for ( node = 0; node < node_num; node++ )
  {
    num = var_node[node];
    var_node[node] = count;
    count = count + num;
  }
  var_node[node_num] = count;

  return var_num;
}
/******************************************************************************/

int *ns_t6_var_set ( int element_num, int element_node[], int node_num, 
  int var_node[], int var_num )

/******************************************************************************/
/*
  Purpose:

    NS_T6_VAR_SET sets the Navier Stokes variables on a T6 grid.

  Discussion:

    We are given a mesh of T6 elements, and asked to create the natural
    list of indices for Navier-Stokes variables associated with each node.
    In particular, every node has two velocity variables associated with
    it, but only a node that is a vertex of the element will also have
    an associated pressure variable.

    The hard work has been done for us already, because the variables
    have been counted, and the pointers to the occurrence of the
    first variable associated with each node have been created.

    The indexing of the nodes can be arbitrary, although a bad
    indexing will result in a miserably large bandwidth (if band
    storage is being tried for the stiffness matrix).  Here, we
    simply try to natural ordering, that is, the variables are
    numbered in order of the node with which they are associated.

    For the Navier Stokes problem on a T6 grid, we take it as
    understood that each node has either 2 or 3 variables associated
    with it, that the first two are always the horizontal and
    then vertical velocity coefficients, and that the third, if
    present, is a pressure coefficient.

    In other settings, it might be necessary not merely to assign
    the variables an index, but also to identify them as to type.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Input, int NODE_NUM, the number of nodes.

    Input, int VAR_NODE[NODE_NUM+1], used to find the variables 
    associated with a given node, which are in VAR in locations 
    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
    this array points to the location just after the last location in VAR.

    Input, int VAR_NUM, the number of variables.

    Output, int NS_T6_VAR_SET[VAR_NUM], the indexes of the variables, which
    are simply 1, 2, 3, ..., VAR_NUM.
*/
{
  int element_order = 6;
  int i;
  int *var;

  var = ( int * ) malloc ( var_num * sizeof ( int ) );

  for ( i = 0; i < var_num; i++ )
  {
    var[i] = i + 1;
  }
  
  return var;
}
/******************************************************************************/

int order_code ( char *code )

/******************************************************************************/
/*
  Purpose: 

    ORDER_CODE returns the order for each element.

  List:

    CODE  Order  Definition
    ----  -----  ----------
    Q4     4     4 node linear Lagrange/serendipity quadrilateral;
    Q8     8     8 node quadratic serendipity quadrilateral;
    Q9     9     9 node quadratic Lagrange quadrilateral;
    Q12   12     12 node cubic serendipity quadrilateral;
    Q16   16     16 node cubic Lagrange quadrilateral;
    QL     6     6 node linear/quadratic quadrilateral;
    T3     3     3 node linear triangle;
    T4     4     4 node bubble triangle;
    T6     6     6 node quadratic triangle;
    T10   10     10 node cubic triangle.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, the code for the element.

    Output, int ORDER_CODE, the order of the element.
*/
{
  int value;

  if ( strcmp ( code, "Q4" ) == 0 )
  {
    value = 4;
  }
  else if ( strcmp ( code, "Q8" ) == 0 )
  {
    value = 8;
  }
  else if ( strcmp ( code, "Q9" ) == 0 )
  {
    value = 9;
  }
  else if ( strcmp ( code, "Q12" ) == 0 )
  {
    value = 12;
  }
  else if ( strcmp ( code, "Q16" ) == 0 )
  {
    value = 16;
  }
  else if ( strcmp ( code, "QL" ) == 0 )
  {
    value = 6;
  }
  else if ( strcmp ( code, "T3" ) == 0 )
  {
    value = 3;
  }
  else if ( strcmp ( code, "T4" ) == 0 )
  {
    value = 4;
  }
  else if ( strcmp ( code, "T6" ) == 0 )
  {
    value = 6;
  }
  else if ( strcmp ( code, "T10" ) == 0 )
  {
    value = 10;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ORDER_CODE - Fatal error!\n" );
    fprintf ( stderr, "  Did not recognize the code \"%s\".\n", code );
    exit ( 1 );
  }

  return value;
}
/******************************************************************************/

void physical_to_reference_t3 ( double t[], int n, double phy[], double ref[] )

/******************************************************************************/
/*
  Purpose:

    PHYSICAL_TO_REFERENCE_T3 maps physical points to reference points.

  Discussion:

    Given the vertices of an order 3 physical triangle and a point
    (X,Y) in the physical triangle, the routine computes the value
    of the corresponding image point (XSI,ETA) in reference space.

    Note that this routine may also be appropriate for an order 6
    triangle, if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image triangle are straight and that the "midside" nodes in the
    physical triangle are halfway along the sides of
    the physical triangle.

  Reference Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the X and Y coordinates
    of the vertices.  The vertices are assumed to be the images of
    (0,0), (1,0) and (0,1) respectively.

    Input, int N, the number of points to transform.

    Input, double PHY[2*N], the coordinates of physical points
    to be transformed.

    Output, double REF[2*N], the coordinates of the corresponding
    points in the reference space.
*/
{
  int j;

  for ( j = 0; j < n; j++ )
  {

    ref[0+j*2] = ( ( t[1+2*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] ) ) 
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );

    ref[1+j*2] = ( ( t[0+1*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] )   
                 - ( t[1+1*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] ) ) 
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );
  }
  return;
}
/******************************************************************************/

void points_plot ( char *file_name, int node_num, double node_xy[], 
  int node_label )

/******************************************************************************/
/*
  Purpose:

    POINTS_PLOT plots a pointset.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to create.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Input, int NODE_LABEL, is TRUE if the nodes are to be labeled.

  Local parameters:

    int CIRCLE_SIZE, controls the size of the circles depicting
    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
    barely visible.
*/
{
  int circle_size = 3;
  int delta;
  FILE *file_unit;
  int i;
  int node;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
/*
  We need to do some figuring here, so that we can determine
  the range of the data, and hence the height and width
  of the piece of paper.
*/
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min ) 
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit = fopen ( file_name, "wt" );

  if ( !file_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "POINTS_PLOT - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output EPS file.\n" );
    exit ( 1 );
  }

  fprintf ( file_unit, "%%!PS-Adobe-3.0 EPSF-3.0\n" );
  fprintf ( file_unit, "%%%%Creator: points_plot.C\n" );
  fprintf ( file_unit, "%%Title: %s\n", file_name );

  fprintf ( file_unit, "%%%%Pages: 1\n" );
  fprintf ( file_unit, "%%%%BoundingBox:  %d  %d  %d  %d\n", 
    x_ps_min, y_ps_min, x_ps_max, y_ps_max );
  fprintf ( file_unit, "%%%%Document-Fonts: Times-Roman\n" );
  fprintf ( file_unit, "%%%%LanguageLevel: 1\n" );
  fprintf ( file_unit, "%%%%EndComments\n" );
  fprintf ( file_unit, "%%%%BeginProlog\n" );
  fprintf ( file_unit, "/inch {72 mul} def\n" );
  fprintf ( file_unit, "%%%%EndProlog\n" );
  fprintf ( file_unit, "%%%%Page:      1     1\n" );
  fprintf ( file_unit, "save\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Set the RGB line color to very light gray.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 0.9000 0.9000 0.9000 setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Draw a gray border around the page.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_min  );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_max );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "stroke\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Set RGB line color to black.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 0.0000 0.0000 0.0000 setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Set the font and its size:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.50 inch scalefont\n" );
  fprintf ( file_unit, "setfont\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Print a title:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  210  702 moveto\n" );
  fprintf ( file_unit, "%%(Pointset) show\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Define a clipping polygon\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_min_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_max_clip );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "clip newpath\n" );
/*
  Draw the nodes.
*/
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Draw filled dots at each node:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Set the color to blue:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "0.000  0.150  0.750  setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );

  for ( node = 0; node < node_num; node++ )
  {
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
      / ( y_max                     - y_min ) );

    fprintf ( file_unit, "newpath  %d  %d  %d  0 360 arc closepath fill\n", 
      x_ps, y_ps, circle_size );
  }
/*
  Label the nodes.
*/
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Label the nodes:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%  Set the color to darker blue:\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "0.000  0.250  0.850  setrgbcolor\n" );
  fprintf ( file_unit, "/Times-Roman findfont\n" );
  fprintf ( file_unit, "0.20 inch scalefont\n" );
  fprintf ( file_unit, "setfont\n" );

  fprintf ( file_unit, "%%\n" );

  for ( node = 0; node < node_num; node++ )
  {
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
      / ( y_max                     - y_min ) );

    fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n", x_ps, y_ps + 5, node+1 );
  }

  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "restore showpage\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% End of page\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%Trailer\n" );
  fprintf ( file_unit, "%%EOF\n" );

  fclose ( file_unit );

  return;
}
/******************************************************************************/

void poly ( char *code, int rexp[], int sexp[] )

/******************************************************************************/
/*
  Purpose: 

    POLY returns the polynomial terms associated with any available element.

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
    'T4', 'T6' and 'T10'.

    Output, int REXP(N), SEXP(N), the powers of R and S associated
    with each polynomial.
*/
{
  if ( s_eqi ( code, "Q4" ) )
  {
    poly_q4 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    poly_q8 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    poly_q9 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    poly_q12 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    poly_q16 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    poly_ql ( rexp, sexp );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    poly_t3 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "POLY - Fatal error!\n" );
    fprintf ( stderr, "  The T4 element does not follow the pattern!\n" );
    exit ( 1 );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    poly_t6 ( rexp, sexp );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    poly_t10 ( rexp, sexp );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "POLY - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = %s\n", code );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void poly_q4 ( int rexp[4], int sexp[4] )

/******************************************************************************/
/*
  Purpose: 

    POLY_Q4 returns the monomials associated with a 4 node quadrilateral.

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[4], SEXP[4], the powers of R and S associated
    with each polynomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 1;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 1;

  return;
}
/******************************************************************************/

void poly_q8 ( int rexp[8], int sexp[8] )

/******************************************************************************/
/*
  Purpose: 

    POLY_Q8 returns the monomials associated with an 8 node quadrilateral.

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[8], SEXP[8], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 1;
  rexp[7] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 2;
  sexp[7] = 1;

  return;
}
/******************************************************************************/

void poly_q9 ( int rexp[9], int sexp[9] )

/******************************************************************************/
/*
  Purpose: 

    POLY_Q9 returns the monomials associated with a 9 node quadrilateral.

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[9], SEXP[9], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 1;
  rexp[7] = 2;
  rexp[8] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 2;
  sexp[7] = 1;
  sexp[8] = 2;

  return;
}
/******************************************************************************/

void poly_q12 ( int rexp[12], int sexp[12] )

/******************************************************************************/
/*
  Purpose: 

    POLY_Q12 returns the monomials associated with a 12 node quadrilateral.

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[12], SEXP[12], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;
  rexp[10] = 1;
  rexp[11] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;
  sexp[10] = 3;
  sexp[11] = 1;

  return;
}
/******************************************************************************/

void poly_q16 ( int rexp[16], int sexp[16] )

/******************************************************************************/
/*
  Purpose: 

    POLY_Q16 returns the monomials associated with a 16 node quadrilateral.

  Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[16], SEXP[16], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;
  rexp[10] = 1;
  rexp[11] = 2;
  rexp[12] = 3;
  rexp[13] = 2;
  rexp[14] = 3;
  rexp[15] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;
  sexp[10] = 3;
  sexp[11] = 2;
  sexp[12] = 1;
  sexp[13] = 3;
  sexp[14] = 2;
  sexp[15] = 3;

  return;
}
/******************************************************************************/

void poly_ql ( int rexp[6], int sexp[6] )

/******************************************************************************/
/*
  Purpose: 

    POLY_QL returns the monomials for a quadratic/linear quadrilateral.

  Element QL:

    |
    1  4---5---6
    |  |       |
    |  |       |
    S  |       |
    |  |       |
    |  |       |
    0  1---2---3
    |
    +--0---R---1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[6], SEXP[6], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 1;
  rexp[4] = 2;
  rexp[5] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 1;
  sexp[4] = 0;
  sexp[5] = 1;

  return;
}
/******************************************************************************/

void poly_t3 ( int rexp[3], int sexp[3] )

/******************************************************************************/
/*
  Purpose: 

    POLY_T3 returns the monomials associated with a 3 node triangle.

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[3], SEXP[3], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;

  return;
}
/******************************************************************************/

void poly_t6 ( int rexp[6], int sexp[6] )

/******************************************************************************/
/*
  Purpose: 

    POLY_T6 returns the monomials associated with a 6 node triangle.

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[6], SEXP[6], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;

  return;
}
/******************************************************************************/

void poly_t10 ( int rexp[10], int sexp[10] )

/******************************************************************************/
/*
  Purpose: 

    POLY_T10 returns the monomials associated with a 10 node triangle.

  Element T10:

    |
    1  10
    |  |\
    |  | \
    |  8  9
    |  |   \
    S  |    \
    |  5  6  7
    |  |      \
    |  |       \
    0  1--2--3--4
    |
    +--0----R---1-->

  Formula:

    Given coefficients A(I), the polynomial interpolant at (R,S) is

      P(R,S) = sum ( 1 <= I <= N ) A(I) * R^REXP(I) * S^SEXP(I) 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2013

  Author:

    John Burkardt

  Parameters:

    Output, int REXP[10], SEXP[10], the powers of R and S associated
    with each monomial.
*/
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;

  return;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_huge ( )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

int r8_nint ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_NINT returns the nearest integer to an R8.

  Example:

        X         R8_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value.

    Output, int R8_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

double r8_power ( double r, int p )

/******************************************************************************/
/*
  Purpose:

    R8_POWER computes an integer power of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, double R, the base.

    Input, int P, the power, which may be negative.

    Output, double R8_POWER, the value of R^P.
*/
{
  double value;
/*
  Special case.  R^0 = 1.
*/
  if ( p == 0 )
  {
    value = 1.0;
  }
/*
  Special case.  Positive powers of 0 are 0.
  We go ahead and compute negative powers, relying on the software to complain.
*/
  else if ( r == 0.0 )
  {
    if ( 0 < p )
    {
      value = 0.0;
    }
    else
    {
      value = pow ( r, p );
    }
  }
  else if ( 1 <= p )
  {
    value = pow ( r, p );
  }
  else
  {
    value = pow ( r, p );
  }

  return value;
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
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

int r8ge_fa ( int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_FA is a simplified version of the LINPACK routine SGEFA.

    The two dimensional array is stored by columns in a one dimensional
    array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N], the matrix to be factored.
    On output, A contains an upper triangular matrix and the multipliers
    which were used to obtain it.  The factorization can be written
    A = L * U, where L is a product of permutation and unit lower
    triangular matrices and U is upper triangular.

    Output, int PIVOT[N], a vector of pivot indices.

    Output, int R8GE_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int i;
  int j;
  int k;
  int l;
  double t;

  for ( k = 1; k <= n-1; k++ )
  {
/*
  Find L, the index of the pivot row.
*/
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*n] ) < r8_abs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
/*
  If the pivot index is zero, the algorithm has failed.
*/
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FA - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", k );
      exit ( 1 );
    }
/*
  Interchange rows L and K if necessary.
*/
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
/*
  Normalize the values that lie below the pivot entry A(K,K).
*/
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
/*
  Row elimination with column indexing.
*/
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8GE_FA - Fatal error!\n" );
    fprintf ( stderr, "  Zero pivot on step %d\n", n );
    exit ( 1 );
  }

  return 0;
}
/******************************************************************************/

double *r8ge_inverse ( int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
    SGEDI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N*N], the factor information computed by R8GE_FA.

    Input, int PIVOT(N), the pivot vector from R8GE_FA.

    Output, double R8GE_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );
/*
  Compute Inverse(U).
*/
  for ( k = 1; k <= n; k++ )
  {
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
    b[k-1+(k-1)*n] = 1.0 / a[k-1+(k-1)*n];

    for ( j = k+1; j <= n; j++ )
    {
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }
/*
  Multiply Inverse(U) by Inverse(L).
*/
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * a[j-1+(k-1)*n];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp = b[i-1+(k-1)*n];
        b[i-1+(k-1)*n] = b[i-1+(pivot[k-1]-1)*n];
        b[i-1+(pivot[k-1]-1)*n] = temp;
      }

    }

  }

  return b;
}
/******************************************************************************/

void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double C[N1*N3], the product matrix C = A * B.
*/
{
  double *d;
  int i;
  int j;
  int k;

  d = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      d[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        d[i+j*n1] = d[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

 for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = d[i+j*n1];
    }
  }

  free ( d );

  return;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
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
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
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

void reference_sample ( char *code, int *seed, double *r, double *s )

/******************************************************************************/
/*
  Purpose:

    REFERENCE_SAMPLE samples a reference element.

  Discussion:

    The routine either samples the unit triangle or the unit square.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
    "T4", "T6" and "T10".

    Input/output, int *SEED, a seed for the random number generator.

    Output, double *R, *S, a random point in the reference element.
*/
{
  if ( strcmp ( code, "Q4" ) == 0 || 
       strcmp ( code, "Q8" ) == 0 || 
       strcmp ( code, "Q9" ) == 0 ||
       strcmp ( code, "Q12" ) == 0 ||
       strcmp ( code, "Q16" ) == 0 ||
       strcmp ( code, "QL" ) == 0 )
  {
    *r = r8_uniform_01 ( seed );
    *s = r8_uniform_01 ( seed );
  }
  else if ( strcmp ( code, "T3" ) == 0 || 
            strcmp ( code, "T4" ) == 0 || 
            strcmp ( code, "T6" ) == 0 || 
            strcmp ( code, "T10" ) == 0 )
  {
    *r = r8_uniform_01 ( seed );
    *s = r8_uniform_01 ( seed );

    if ( 1.0 < *r + *s )
    {
      *r = 1.0 - *r;
      *s = 1.0 - *s;
    }
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "REFERENCE_SAMPLE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal code = \"%s\".\n", code );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void reference_to_physical_q4 ( double q4[2*4], int n, double rs[], 
  double xy[] )

/******************************************************************************/
/*
  Purpose:

    REFERENCE_TO_PHYSICAL_Q4 maps Q4 reference points to physical points.

  Discussion:

    XY(R,S) = XY(0,0) * (1-R) * (1-S)
            + XY(1,0) *    R  * (1-S)
            + XY(1,1) *    R  *    S
            + XY(0,1) * (1-R) *    S

  Reference Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double Q4[2*4], the coordinates of the vertices.
    The vertices are assumed to be the images of the reference vertices
    (0,0), (1,0), (1,1) and (0,1) respectively.

    Input, int N, the number of points to transform.

    Input, double RS[2*N], (R,S) points in the reference element.

    Output, double XY[2*N], (X,Y) points in the physical element.
*/
{
  int j;
  double *psi;

  psi = ( double * ) malloc ( 4 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    psi[0+j*2] = ( 1.0 - rs[0+j*2] ) * ( 1.0 - rs[1+j*2] );
    psi[1+j*2] =         rs[0+j*2]   * ( 1.0 - rs[1+j*2] );
    psi[2+j*2] =         rs[0+j*2]   *         rs[1+j*2];
    psi[3+j*2] = ( 1.0 - rs[0+j*2] ) *         rs[1+j*2];
  }

  r8mat_mm ( 2, 4, n, q4, psi, xy );

  free ( psi );

  return;
}
/******************************************************************************/

void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] )

/******************************************************************************/
/*
  Purpose:

    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.

  Discussion:

    Given the vertices of an order 3 physical triangle and a point
    (XSI,ETA) in the reference triangle, the routine computes the value
    of the corresponding image point (X,Y) in physical space.

    Note that this routine may also be appropriate for an order 6
    triangle, if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image triangle are straight and that the "midside" nodes in the
    physical triangle are halfway along the sides of
    the physical triangle.

  Reference Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*3], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0) and
    (0,1) respectively.

    Input, int N, the number of points to transform.

    Input, double REF[2*N], points in the reference triangle.

    Output, double PHY[2*N], corresponding points in the
    physical triangle.
*/
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] ) 
                 + t[i+1*2] *       + ref[0+j*2]                
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
/******************************************************************************/

void reference_to_physical_t6 ( double t[], int n, double ref[], double phy[] )

/******************************************************************************/
/*
  Purpose:

    REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.

  Discussion:

    Given the vertices of an order 6 physical triangle and a point
    (XSI,ETA) in the reference triangle, the routine computes the value
    of the corresponding image point (X,Y) in physical space.

    The mapping from (XSI,ETA) to (X,Y) has the form:

      X(ETA,XSI) = A1 * XSI^2 + B1 * XSI*ETA + C1 * ETA^2
                 + D1 * XSI   + E1 * ETA     + F1

      Y(ETA,XSI) = A2 * XSI^2 + B2 * XSI*ETA + C2 * ETA^2
                 + D2 * XSI   + E2 * ETA     + F2

  Reference Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T[2*6], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0),
    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.

    Input, int N, the number of points to transform.

    Input, double REF[2*N], points in the reference triangle.

    Output, double PHY[2*N], corresponding points in the
    physical triangle.
*/
{
  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double e[2];
  double f[2];
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    a[i] =   2.0 * t[i+0*2] + 2.0 * t[i+1*2]
           - 4.0 * t[i+3*2];

    b[i] =   4.0 * t[i+0*2] 
           - 4.0 * t[i+3*2] + 4.0 * t[i+4*2] - 4.0 * t[i+5*2];

    c[i] =   2.0 * t[i+0*2]                  + 2.0 * t[i+2*2] 
                                             - 4.0 * t[i+5*2];

    d[i] = - 3.0 * t[i+0*2] -       t[i+1*2]
           + 4.0 * t[i+3*2];

    e[i] = - 3.0 * t[i+0*2]                  -       t[i+2*2]
                                             + 4.0 * t[i+5*2];
    f[i] =         t[i+0*2];

  }

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = a[i] * ref[0+j*2] * ref[0+j*2] 
                 + b[i] * ref[0+j*2] * ref[1+j*2] 
                 + c[i] * ref[1+j*2] * ref[1+j*2]
                 + d[i] * ref[0+j*2]
                 + e[i] * ref[1+j*2] 
                 + f[i];
    }
  }

  return;
}
/******************************************************************************/

int s_eqi ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_EQI reports whether two strings are equal, ignoring case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, pointers to two strings.

    Output, int S_EQI, is true if the strings are equal.
*/
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  if ( nchar1 < nchar2 )
  {
    nchar = nchar1;
  }
  else
  {
    nchar = nchar2;
  }
/*
  The strings are not equal if they differ over their common length.
*/
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
/*
  The strings are not equal if the longer one includes nonblanks
  in the tail.
*/
  if ( nchar < nchar1 )
  {
    for ( i = nchar; i < nchar1; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return 0;
      }
    }
  }
  else if ( nchar < nchar2 )
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return 0;
      }
    }
  }

  return 1;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Discussion:

    It turns out that I also want to ignore the '\n' character!

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2014

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
    if ( *t != ' ' && *t != '\n' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

double serene ( char *type, double ve, double vn, double vne, double vnw, 
  double vs, double vse, double vsw, double vw )

/******************************************************************************/
/*
  Purpose: 

    SERENE interpolates data using a Q8 element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *TYPE, tells SERENE the geometry of the
    finite element that surrounds the point of interest.  The options
    are displayed in the following table, which suggests the meaning
    of each option by its position:

        |   |
     NW * N * NE
        |   |
     -*-*-*-*-*-
        |   |
      W * C * E
        |   |
     -*-*-*-*-*-
        |   |
     SW * S * SE
        |   |

    Input, double VE, VN, VNE, VNW, VS, VSE, VSW, VW,
    are the values of the function at the nodes to the east,
    north, northeast, northwest, south, southeast, southwest and
    west of the point of interest.  If the finite element is of
    type 'C', then all 8 values are needed.  However, if the
    finite element is of type 'SE', for instance, then only three
    values are needed, namely VE, VN, and VNW, since these are
    the only node positions defined in such a finite element.

    Output, double SERENE, the interpolated value of the
    function at the point of interest.
*/
{
  double eta;
  double pe;
  double pn;
  double pne;
  double pnw;
  double ps;
  double pse;
  double psw;
  double pw;
  double vterp;
  double xsi;
/*
  To make this routine more general, simply pass in the values of XSI
  and ETA at which the interpolated value is desired.

  By setting XSI = ETA = 0, we are asking for the interpolated value
  at the center of the finite element.
*/
  xsi = 0.0;
  eta = 0.0;
/*
  8 node center

  Polynomial space is spanned by:

         1
       x    y
    x^2  xy  y^2
      x^2y xy^2


    ^   1    4--7--3
    |        !     !
    E        !     !
    T   0    8  X  6
    A        !     !
    |        !     !
    V  -1    1--5--2

            -1  0  1

           <---XSI--->
*/
  if ( s_eqi ( type, "C" ) )
  {
    psw = - 0.25 * ( 1.0 - xsi ) * ( 1.0 - eta )
      * ( 1.0 + xsi + eta );
    pse = - 0.25 * ( 1.0 + xsi ) * ( 1.0 - eta )
      * ( 1.0 - xsi + eta );
    pne = - 0.25 * ( 1.0 + xsi ) * ( 1.0 + eta )
      * ( 1.0 - xsi - eta );
    pnw = - 0.25 * ( 1.0 - xsi ) * ( 1.0 + eta )
      * ( 1.0 + xsi - eta );
    ps =    0.50 * ( 1.0 - xsi ) * ( 1.0 + xsi )
      * ( 1.0 - eta );
    pe =    0.50 * ( 1.0 + xsi ) * ( 1.0 + eta ) 
      * ( 1.0 - eta );
    pn =    0.50 * ( 1.0 - xsi ) * ( 1.0 + xsi ) 
      * ( 1.0 + eta );
    pw =    0.50 * ( 1.0 - xsi ) * ( 1.0 + eta ) 
      * ( 1.0 - eta );

    vterp = vsw * psw + vse * pse + vne * pne + vnw * pnw 
      + vs * ps + ve * pe + vn * pn + vw * pw;
  }
/*
  5 node side

    ^   1
    |
    E
    T   0    8  X  6
    A        !     !
    |        !     !
    V  -1    1--5--2

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "N" ) )
  {
    psw =  0.5 * ( xsi - 1.0 ) * ( 1.0 + xsi + eta );
    pse = -0.5 * ( xsi + 1.0 ) * ( 1.0 - xsi + eta );
    ps =  -          ( xsi + 1.0 ) * ( xsi - 1.0 );
    pe =   0.5 * ( xsi + 1.0 ) * ( eta + 1.0 );
    pw =  -0.5 * ( xsi - 1.0 ) * ( eta + 1.0 );

    vterp = vsw * psw + vse * pse + vs * ps + ve * pe + vw * pw;
  }
/*
    ^   1    4--7
    |        !
    E        !
    T   0    8  X
    A        !
    |        !
    V  -1    1--5

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "E" ) )
  {
    pse =  0.5 * ( eta - 1.0 ) * ( 1.0 + xsi + eta );
    pne = -0.5 * ( eta + 1.0 ) * ( 1.0 + xsi - eta );
    ps =  -0.5 * ( xsi + 1.0 ) * ( eta - 1.0 );
    pn =   0.5 * ( xsi + 1.0 ) * ( eta + 1.0 );
    pw =  -          ( eta + 1.0 ) * ( eta - 1.0 );

    vterp = vse * pse + vne * pne + vs * ps + vn * pn + vw * pw;
  }
/*
  5 node side

    ^   1       7--3
    |              !
    E              !
    T   0       X  6
    A              !
    |              !
    V  -1       5--2

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "W" ) )
  {
    pse =  0.5 * ( eta - 1.0 ) * ( 1.0 - xsi + eta );
    pne = -0.5 * ( eta + 1.0 ) * ( 1.0 - xsi - eta );
    ps =   0.5 * ( xsi - 1.0 ) * ( eta - 1.0 );
    pe =  -          ( eta - 1.0 ) * ( eta + 1.0 );
    pn =  -0.5 * ( xsi - 1.0 ) * ( eta + 1.0 );

    vterp = vse * pse + vne * pne + vs * ps + ve * pe + vn * pn;
  }
/*
  5 node side

    ^   1    4--7--3
    |        !     !
    E        !     !
    T   0    8  X  6
    A
    |
    V  -1

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "S" ) )
  {
    pne = -0.5 * ( xsi + 1.0 ) * ( 1.0 - xsi - eta );
    pnw =  0.5 * ( xsi - 1.0 ) * ( 1.0 + xsi - eta );
    pe =  -0.5 * ( eta - 1.0 ) * ( xsi + 1.0 );
    pn =  -          ( xsi + 1.0 ) * ( xsi - 1.0 );
    pw =   0.5 * ( eta - 1.0 ) * ( xsi - 1.0 );

    vterp = vne * pne + vnw * pnw + ve * pe + vn * pn + vw * pw;
  }
/*
  3 node corner

  Polynomial space is spanned by:

         1
       x    y


    ^   1
    |
    E
    T   0    8  X
    A        !
    |        !
    V  -1    1--5

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "NE" ) )
  {
    psw = -1.0 - xsi - eta;
    ps =   1.0 + xsi;
    pw =   1.0       + eta;

    vterp = vsw * psw + vs * ps + vw * pw;
  }
/*
  3 node corner

  Polynomial space is spanned by:

         1
       x    y

    ^   1
    |
    E
    T   0       X  6
    A              !
    |              !
    V  -1       5--2

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "NW" ) )
  {
    pse = 1.0 + xsi - eta;
    ps =  1.0 - xsi;
    pe =  1.0       + eta;

    vterp = vse * pse + vs * ps + ve * pe;
  }
/*
  3 node corner

  Polynomial space is spanned by:
         1
       x    y


    ^   1    4--7
    |        !
    E        !
    T   0    8  X
    A
    |
    V  -1

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "SE" ) )
  {
    pnw = - 1.0 - xsi + eta;
    pn =    1.0 + xsi;
    pw =    1.0       - eta;

    vterp = vnw * pnw + vn * pn + vw * pw;
  }
/*
  3 node corner

  Polynomial space is spanned by:

         1
       x    y

    ^   1       7--3
    |              !
    E              !
    T   0       X  6
    A
    |
    V  -1

            -1  0  1

           <---XSI--->
*/
  else if ( s_eqi ( type, "SW" ) )
  {
    pne = - 1.0 + xsi + eta;
    pe =    1.0       - eta;
    pn =    1.0 - xsi;

    vterp = vne * pne + ve * pe + vn * pn;
  }

  return vterp;
}
/******************************************************************************/

void shape ( char *code, double r, double s, double t[], 
  double dtdr[], double dtds[] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE evaluates shape functions for any available element.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, string CODE, identifies the element.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T4', 'T6' and 'T10'.

    Input, double R, S, the reference coordinates of a point.

    Output, double T[N], the basis functions at the point.

    Output, double DTDR[N], the R basis derivatives at the point.

    Output, double DTDS[N], the S basis derivatives at the point.
*/
{
  if ( s_eqi ( code, "Q4" ) )
  {
    shape_q4 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    shape_q8 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    shape_q9 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    shape_q12 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    shape_q16 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    shape_ql ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    shape_t3 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    shape_t4 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    shape_t6 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    shape_t10 ( r, s, t, dtdr, dtds );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SHAPE - Fatal error!\n" );
    fprintf ( stderr, "  Unrecognized code = %s\n", code );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void shape_q4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_Q4 evaluates shape functions for a 4 node quadrilateral.

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[4], the basis functions at the point.

    Output, double DTDR[4], the R basis derivatives at the point.

    Output, double DTDS[4], the S basis derivatives at the point.
*/
{
  t[0] = ( 1.0 - r ) * ( 1.0 - s );
  t[1] =             r   * ( 1.0 - s );
  t[2] =             r   *             s;
  t[3] = ( 1.0 - r ) *             s;

  dtdr[0] = - 1.0 + s;
  dtdr[1] =   1.0 - s;   
  dtdr[2] =             s;
  dtdr[3] =           - s;

  dtds[0] = - 1.0 + r;
  dtds[1] =           - r;
  dtds[2] =             r;
  dtds[3] =   1.0 - r;

  return;
}
/******************************************************************************/

void shape_q8 ( double r, double s, double t[8], double dtdr[8], 
  double dtds[8] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_Q8 evaluates shape functions for an 8 node quadrilateral.

  Comment:

    This element is known as the "serendipity" element.

  Element Q8:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8     6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2005

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[8], the basis functions at the point.

    Output, double DTDR[8], the R basis derivatives at the point.

    Output, double DTDS[8], the S basis derivatives at the point.
*/
{
  t[0] =                 ( r - 1.0 )     * ( s - 1.0 ) 
    * ( 1.0 - 2.0 * r - 2.0 * s );
  t[1] =             r                       * ( s - 1.0 ) 
    * ( 1.0 - 2.0 * r + 2.0 * s );
  t[2] =             r                   * s                   
    * ( 2.0 * r + 2.0 * s - 3.0 );
  t[3] =                 ( r - 1.0 ) * s                   
    * ( 2.0 * r - 2.0 * s + 1.0 );
  t[4] =   4.0 * r * ( r - 1.0 )     * ( s - 1.0 );
  t[5] = - 4.0 * r                   * s * ( s - 1.0 );
  t[6] = - 4.0 * r * ( r - 1.0 ) * s;   
  t[7] =   4.0 *     ( r - 1.0 ) * s * ( s - 1.0 );

  dtdr[0] = ( s - 1.0 ) * ( - 4.0 * r - 2.0 * s + 3.0 );
  dtdr[1] = ( s - 1.0 ) * ( - 4.0 * r + 2.0 * s + 1.0 );
  dtdr[2] =   s         * (   4.0 * r + 2.0 * s - 3.0 );
  dtdr[3] =   s         * (   4.0 * r - 2.0 * s - 1.0 );
  dtdr[4] =   4.0 * ( 2.0 * r - 1.0 )     * ( s - 1.0 );
  dtdr[5] = - 4.0 *                     s * ( s - 1.0 );
  dtdr[6] = - 4.0 * ( 2.0 * r - 1.0 ) * s;
  dtdr[7] =   4.0 *                     s * ( s - 1.0 );

  dtds[0] = ( r - 1.0 ) * ( - 4.0 * s - 2.0 * r + 3.0 );
  dtds[1] =   r *       (   4.0 * s - 2.0 * r - 1.0 );
  dtds[2] =   r *       (   4.0 * s + 2.0 * r - 3.0 );
  dtds[3] = ( r - 1.0 ) * ( - 4.0 * s + 2.0 * r + 1.0 );
  dtds[4] =   4.0 * r * ( r - 1.0 );
  dtds[5] = - 4.0 * r               * ( 2.0 * s - 1.0 );
  dtds[6] = - 4.0 * r * ( r - 1.0 );
  dtds[7] =   4.0 *     ( r - 1.0 ) * ( 2.0 * s - 1.0 );

  return;
}
/******************************************************************************/

void shape_q9 ( double r, double s, double t[9], double dtdr[9], 
  double dtds[9] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_Q9 evaluates shape functions for a 9 node quadrilateral.

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[9], the basis functions at the point.

    Output, double DTDR[9], the R basis derivatives at the point.

    Output, double DTDS[9], the S basis derivatives at the point.
*/
{
  t[0] =    4.0 * ( r - 1.0 ) * ( r - 0.5 ) * ( s - 1.0 ) 
            * ( s - 0.5 );
  t[1] =    4.0 * r * ( r - 0.5 ) * ( s - 1.0 ) * ( s - 0.5 );
  t[2] =    4.0 * r * ( r - 0.5 ) * s * ( s - 0.5 );
  t[3] =    4.0 * ( r - 1.0 ) * ( r - 0.5 ) * s * ( s - 0.5 );
  t[4] = -  8.0 * r * ( r - 1.0 ) * ( s - 1.0 ) * ( s - 0.5 );
  t[5] = -  8.0 * r * ( r - 0.5 ) * s * ( s - 1.0 );
  t[6] = -  8.0 * r * ( r - 1.0 ) * s * ( s - 0.5 );
  t[7] = -  8.0 * ( r - 1.0 ) * ( r - 0.5 ) * s * ( s - 1.0 );
  t[8] =   16.0 * r * ( r - 1.0 ) * s * ( s - 1.0 );

  dtdr[0] =   4.0 * ( 2.0 * r - 1.5 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[1] =   4.0 * ( 2.0 * r - 0.5 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[2] =   4.0 * ( 2.0 * r - 0.5 ) * s * ( s - 0.5 );
  dtdr[3] =   4.0 * ( 2.0 * r - 1.5 ) * s * ( s - 0.5 );

  dtdr[4] = - 8.0 * ( 2.0 * r - 1.0 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[5] = - 8.0 * ( 2.0 * r - 0.5 ) * s * ( s - 1.0 );
  dtdr[6] = - 8.0 * ( 2.0 * r - 1.0 ) * s * ( s - 0.5 );
  dtdr[7] = - 8.0 * ( 2.0 * r - 1.5 ) * s * ( s - 1.0 );
  dtdr[8] =  16.0 * ( 2.0 * r - 1.0 ) * s * ( s - 1.0 );

  dtds[0] =   4.0 * ( r - 1.0 ) * ( r - 0.5 ) 
              * ( 2.0 * s - 1.5 );
  dtds[1] =   4.0 * r * ( r - 0.5 ) * ( 2.0 * s - 1.5 );
  dtds[2] =   4.0 * r * ( r - 0.5 ) * ( 2.0 * s - 0.5 );
  dtds[3] =   4.0 * ( r - 1.0 ) * ( r - 0.5 ) 
            * ( 2.0 * s - 0.5 );
  dtds[4] = - 8.0 * r * ( r - 1.0 ) * ( 2.0 * s - 1.5 );
  dtds[5] = - 8.0 * r * ( r - 0.5 ) * ( 2.0 * s - 1.0 );
  dtds[6] = - 8.0 * r * ( r - 1.0 ) * ( 2.0 * s - 0.5 );
  dtds[7] = - 8.0 * ( r - 1.0 ) * ( r - 0.5 ) 
            * ( 2.0 * s - 1.0 );
  dtds[8] =  16.0 * r * ( r - 1.0 ) * ( 2.0 * s - 1.0 );

  return;
}
/******************************************************************************/

void shape_q12 ( double r, double s, double t[12], double dtdr[12], 
  double dtds[12] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_Q12 evaluates shape functions for a 12 node quadrilateral.

  Element Q12:

    |
    1  9-10-11-12
    |  |        |
    |  7        8
    S  |        |
    |  5        6
    |  |        |
    0  1--2--3--4
    |
    +--0---R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[12], the basis functions at the point.

    Output, double DTDR[12], the R basis derivatives at the point.

    Output, double DTDS[12], the S basis derivatives at the point.
*/
{
  double a;
  double b;
  double c;
  double corner;
  double d;
  double dcdr;
  double dcds;

  a = 0.0;
  b = 1.0 / 3.0;
  c = 2.0 / 3.0;
  d = 1.0;

  corner = 9.0 * ( ( 2.0 * r - 1.0 ) * ( 2.0 * r - 1.0 )
                     + ( 2.0 * s - 1.0 ) * ( 2.0 * s - 1.0 ) ) 
    - 10.0;

  t[0] =     0.125  * ( r - d ) * ( s - d ) * corner;
  t[1] =  - 13.5    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - d );
  t[2] =    13.5    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - d );
  t[3] =   - 0.125  * ( r - a ) * ( s - d ) * corner;
  t[4] =  - 13.5    * ( r - d ) * ( s - a ) * ( s - c ) * ( s - d );
  t[5] =    13.5    * ( r - a ) * ( s - a ) * ( s - c ) * ( s - d );
  t[6] =    13.5    * ( r - d ) * ( s - a ) * ( s - b ) * ( s - d );
  t[7] =  - 13.5    * ( r - a ) * ( s - a ) * ( s - b ) * ( s - d );
  t[8] =   - 0.125  * ( r - d ) * ( s - a ) * corner;
  t[9] =   13.5    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - a );
  t[10] = - 13.5    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - a );
  t[11] =    0.125  * ( r - a ) * ( s - a ) * corner;
 
  dcdr = 36.0 * ( 2.0 * r - 1.0 );

  dtdr[0] =  0.125 * ( s - d ) * ( ( r - d ) * dcdr + corner );
  dtdr[1] =  - 13.5 * ( s - d ) * ( 3.0 * r * r 
    - 2.0 * ( a + c + d ) * r + a * c + c * d + d * a ); 
  dtdr[2] =    13.5 * ( s - d ) * ( 3.0 * r * r 
    - 2.0 * ( a + b + d ) * r + a * b + b * d + d * a );
  dtdr[3] = - 0.125 * ( s - d ) * ( ( r - a ) * dcdr + corner );
  dtdr[4] = - 13.5 * ( s - a ) * ( s - c ) * ( s - d );
  dtdr[5] =   13.5 * ( s - a ) * ( s - c ) * ( s - d );
  dtdr[6] =   13.5 * ( s - a ) * ( s - b ) * ( s - d );
  dtdr[7] = - 13.5 * ( s - a ) * ( s - b ) * ( s - d );
  dtdr[8] = - 0.125 * ( s - a ) * ( ( r - d ) * dcdr + corner );
  dtdr[9] =   13.5 * ( s - a ) * ( 3.0 * r * r 
    - 2.0 * ( a + c + d ) * r + a * c + c * d + d * a );
  dtdr[10] = - 13.5 * ( s - a ) * ( 3.0 * r * r
    - 2.0 * ( a + b + d ) * r + a * b + b * d + d * a );
  dtdr[11] = 0.125 * ( s - a ) * ( ( r - a ) * dcdr + corner );

  dcds = 36.0 * ( 2.0 * s - 1.0 );

  dtds[0] =  0.125 * ( r - d ) * ( corner + ( s - d ) * dcds );
  dtds[1] =  - 13.5 * ( r - a ) * ( r - c ) * ( r - d );
  dtds[2] =  13.5 * ( r - a ) * ( r - b ) * ( r - d );
  dtds[3] = - 0.125  * ( r - a ) * ( corner + ( s - d ) * dcds );
  dtds[4] =  - 13.5 * ( r - d ) * ( 3.0 * s * s 
    - 2.0 * ( a + c + d ) * s + a * c + c * d + d * a );
  dtds[5] =  13.5 * ( r - a ) * ( 3.0 * s * s 
    - 2.0 * ( a + c + d ) * s + a * c + c * d + d * a );
  dtds[6] =  13.5 * ( r - d ) * ( 3.0 * s * s 
    - 2.0 * ( a + b + d ) * s + a * b + b * d + d * a );
  dtds[7] =  - 13.5 * ( r - a ) * ( 3.0 * s * s 
    - 2.0 * ( a + b + d ) * s + a * b + b * d + d * a );
  dtds[8] =  - 0.125 * ( r - d ) * ( corner + ( s - a ) * dcds );
  dtds[9] = 13.5 * ( r - a ) * ( r - c ) * ( r - d ) ;
  dtds[10] = - 13.5 * ( r - a ) * ( r - b ) * ( r - d ) ;
  dtds[11] = 0.125 * ( r - a ) * ( corner + ( s - a ) * dcds );

  return;
}
/******************************************************************************/

void shape_q16 ( double r, double s, double t[16], double dtdr[16], 
  double dtds[16] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_Q16 evaluates shape functions for a 16 node quadrilateral.

  Diagram:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[16], the basis functions at the point.

    Output, double DTDR[16], the R basis derivatives at the point.

    Output, double DTDS[16], the S basis derivatives at the point.
*/
{
  double dabc;
  double dabd;
  double dacd;
  double dbcd;
  double ra;
  double rb;
  double rc;
  double rd;
  double sa;
  double sb;
  double sc;
  double sd;

  ra = r - 0.0;
  rb = r - 1.0 / 3.0;
  rc = r - 2.0 / 3.0;
  rd = r - 1.0;

  sa = s - 0.0;
  sb = s - 1.0 / 3.0;
  sc = s - 2.0 / 3.0;
  sd = s - 1.0;

  t[0]  =   (  81.0 / 4.0 ) * rb * rc * rd * sb * sc * sd;
  t[1]  = - ( 243.0 / 4.0 ) * ra * rc * rd * sb * sc * sd;
  t[2]  =   ( 243.0 / 4.0 ) * ra * rb * rd * sb * sc * sd;
  t[3]  = - (  81.0 / 4.0 ) * ra * rb * rc * sb * sc * sd;

  t[4]  = - ( 243.0 / 4.0 ) * rb * rc * rd * sa * sc * sd;
  t[5]  =   ( 729.0 / 4.0 ) * ra * rc * rd * sa * sc * sd;
  t[6]  = - ( 729.0 / 4.0 ) * ra * rb * rd * sa * sc * sd;
  t[7]  =   ( 243.0 / 4.0 ) * ra * rb * rc * sa * sc * sd;

  t[8]  =   ( 243.0 / 4.0 ) * rb * rc * rd * sa * sb * sd;
  t[9]  = - ( 729.0 / 4.0 ) * ra * rc * rd * sa * sb * sd;
  t[10] =   ( 729.0 / 4.0 ) * ra * rb * rd * sa * sb * sd;
  t[11] = - ( 243.0 / 4.0 ) * ra * rb * rc * sa * sb * sd;

  t[12] = - (  81.0 / 4.0 ) * rb * rc * rd * sa * sb * sc;
  t[13] =   ( 243.0 / 4.0 ) * ra * rc * rd * sa * sb * sc;
  t[14] = - ( 243.0 / 4.0 ) * ra * rb * rd * sa * sb * sc;
  t[15] =   (  81.0 / 4.0 ) * ra * rb * rc * sa * sb * sc;

  dbcd = 3.0 * r * r -  4.0 * r       + 11.0 / 9.0;
  dacd = 3.0 * r * r - 10.0 * r / 3.0 +  2.0 / 3.0;
  dabd = 3.0 * r * r -  8.0 * r / 3.0 +  1.0 / 3.0;
  dabc = 3.0 * r * r -  2.0 * r       +  2.0 / 9.0;

  dtdr[0]  =   (  81.0 / 4.0 ) * dbcd * sb * sc * sd;
  dtdr[1]  = - ( 243.0 / 4.0 ) * dacd * sb * sc * sd;
  dtdr[2]  =   ( 243.0 / 4.0 ) * dabd * sb * sc * sd;
  dtdr[3]  = - (  81.0 / 4.0 ) * dabc * sb * sc * sd;
  dtdr[4]  = - ( 243.0 / 4.0 ) * dbcd * sa * sc * sd;
  dtdr[5]  =   ( 729.0 / 4.0 ) * dacd * sa * sc * sd;
  dtdr[6]  = - ( 729.0 / 4.0 ) * dabd * sa * sc * sd;
  dtdr[7]  =   ( 243.0 / 4.0 ) * dabc * sa * sc * sd;
  dtdr[8]  =   ( 243.0 / 4.0 ) * dbcd * sa * sb * sd;
  dtdr[9]  = - ( 729.0 / 4.0 ) * dacd * sa * sb * sd;
  dtdr[10] =   ( 729.0 / 4.0 ) * dabd * sa * sb * sd;
  dtdr[11] = - ( 243.0 / 4.0 ) * dabc * sa * sb * sd;
  dtdr[12] = - (  81.0 / 4.0 ) * dbcd * sa * sb * sc;
  dtdr[13] =   ( 243.0 / 4.0 ) * dacd * sa * sb * sc;
  dtdr[14] = - ( 243.0 / 4.0 ) * dabd * sa * sb * sc;
  dtdr[15] =   (  81.0 / 4.0 ) * dabc * sa * sb * sc;

  dbcd = 3.0 * s * s -  4.0 * s       + 11.0 / 9.0;
  dacd = 3.0 * s * s - 10.0 * s / 3.0 +  2.0 / 3.0;
  dabd = 3.0 * s * s -  8.0 * s / 3.0 +  1.0 / 3.0;
  dabc = 3.0 * s * s -  2.0 * s       +  2.0 / 9.0;

  dtds[0]  =   (  81.0 / 4.0 ) * rb * rc * rd * dbcd;
  dtds[1]  = - ( 243.0 / 4.0 ) * ra * rc * rd * dbcd;
  dtds[2]  =   ( 243.0 / 4.0 ) * ra * rb * rd * dbcd;
  dtds[3]  = - (  81.0 / 4.0 ) * ra * rb * rc * dbcd;
  dtds[4]  = - ( 243.0 / 4.0 ) * rb * rc * rd * dacd;
  dtds[5]  =   ( 729.0 / 4.0 ) * ra * rc * rd * dacd;
  dtds[6]  = - ( 729.0 / 4.0 ) * ra * rb * rd * dacd;
  dtds[7]  =   ( 243.0 / 4.0 ) * ra * rb * rc * dacd;
  dtds[8]  =   ( 243.0 / 4.0 ) * rb * rc * rd * dabd;
  dtds[9]  = - ( 729.0 / 4.0 ) * ra * rc * rd * dabd;
  dtds[10] =   ( 729.0 / 4.0 ) * ra * rb * rd * dabd;
  dtds[11] = - ( 243.0 / 4.0 ) * ra * rb * rc * dabd;
  dtds[12] = - (  81.0 / 4.0 ) * rb * rc * rd * dabc;
  dtds[13] =   ( 243.0 / 4.0 ) * ra * rc * rd * dabc;
  dtds[14] = - ( 243.0 / 4.0 ) * ra * rb * rd * dabc;
  dtds[15] =   (  81.0 / 4.0 ) * ra * rb * rc * dabc;
  
  return;
}
/******************************************************************************/

void shape_ql ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_QL evaluates shape functions for a 6 node quadratic/linear.

  Diagram:

    |
    1  4--5--6
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1--2--3
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[6], the basis functions at the point.

    Output, double DTDR[6], the R basis derivatives at the point.

    Output, double DTDS[6], the S basis derivatives at the point.
*/
{
  t[0] = - 2.0 *     ( r - 0.5 ) * ( r - 1.0 )     * ( s - 1.0 );
  t[1] =   4.0 * r                   * ( r - 1.0 )     * ( s - 1.0 );
  t[2] = - 2.0 * r * ( r - 0.5 )                       * ( s - 1.0 );
  t[3] =   2.0 *     ( r - 0.5 ) * ( r - 1.0 ) * s;
  t[4] = - 4.0 * r                   * ( r - 1.0 ) * s;
  t[5] =   2.0 * r * ( r - 0.5 )                   * s;

  dtdr[0] = 2.0 * ( - 2.0 * r + 1.5 )     * ( s - 1.0 );
  dtdr[1] = 4.0 * (   2.0 * r - 1.0 )     * ( s - 1.0 );
  dtdr[2] = 2.0 * ( - 2.0 * r + 0.5 )     * ( s - 1.0 );
  dtdr[3] = 2.0 * (   2.0 * r - 1.5 ) * s;
  dtdr[4] = 4.0 * ( - 2.0 * r + 1.0 ) * s;
  dtdr[5] = 2.0 * (   2.0 * r - 0.5 ) * s;

  dtds[0] = - 2.0 *     ( r - 0.5 ) * ( r - 1.0 );
  dtds[1] =   4.0 * r                   * ( r - 1.0 );
  dtds[2] = - 2.0 * r * ( r - 0.5 );
  dtds[3] =   2.0 *     ( r - 0.5 ) * ( r - 1.0 );
  dtds[4] = - 4.0 * r                   * ( r - 1.0 );
  dtds[5] =   2.0 * r * ( r - 0.5 );

  return;
}
/******************************************************************************/

void shape_t3 ( double r, double s, double t[3], double dtdr[3], 
  double dtds[3] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_T3 evaluates shape functions for a 3 node triangle.

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[3], the basis functions at the point.

    Output, double DTDR[3], the R basis derivatives at the point.

    Output, double DTDS[3], the S basis derivatives at the point.
*/
{
  t[0] = 1.0 - r - s;
  t[1] =           r;
  t[2] =               s;

  dtdr[0] = -1.0;
  dtdr[1] =  1.0;
  dtdr[2] =  0.0;

  dtds[0] = -1.0;
  dtds[1] =  0.0;
  dtds[2] =  1.0;

  return;
}
/******************************************************************************/

void shape_t4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_T4 evaluates shape functions for a T4 triangle.

  Element T4:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  | 4 \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[4], the basis functions at the point.

    Output, double DTDR[3], the R basis derivatives at the point.

    Output, double DTDS[3], the S basis derivatives at the point.
*/
{
  t[0] = ( 1.0 - 9.0 * r * s ) * ( 1.0 - r - s );
  t[1] = r * ( 1.0 - 9.0 * ( 1.0 - r - s ) * s );
  t[2] = s * ( 1.0 - 9.0 * ( 1.0 - r - s ) * r );
  t[3] = 27.0            * ( 1.0 - r - s ) * r * s;

  dtdr[0] = -1.0 +  9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[1] =  1.0 +  9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[2] =         9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[3] =      - 27.0 * ( - s + 2.0 * r * s + s * s );

  dtds[0] = -1.0 +  9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[1] =         9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[2] =  1.0 +  9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[3] =      - 27.0 * ( - r + r * r + 2.0 * r * s );

  return;
}
/******************************************************************************/

void shape_t6 ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_T6 evaluates shape functions for a 6 node triangle.

  Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[6], the basis functions at the point.

    Output, double DTDR[6], the R basis derivatives at the point.

    Output, double DTDS[6], the S basis derivatives at the point.
*/
{
  t[0] = 2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
  t[1] = 2.0 * r * ( r - 0.5 );
  t[2] = 2.0 * s * ( s - 0.5 );
  t[3] = 4.0 * r * ( 1.0 - r - s );
  t[4] = 4.0 * r * s;
  t[5] = 4.0 * s * ( 1.0 - r - s );

  dtdr[0] = - 3.0 + 4.0 * r + 4.0 * s;
  dtdr[1] = - 1.0 + 4.0 * r;
  dtdr[2] =   0.0;
  dtdr[3] =   4.0 - 8.0 * r - 4.0 * s;
  dtdr[4] =                           4.0 * s;
  dtdr[5] =                         - 4.0 * s;

  dtds[0] = - 3.0 + 4.0 * r + 4.0 * s;
  dtds[1] =   0.0;
  dtds[2] = - 1.0               + 4.0 * s;
  dtds[3] =           - 4.0 * r;
  dtds[4] =             4.0 * r;
  dtds[5] =   4.0 - 4.0 * r - 8.0 * s;

  return;
}
/******************************************************************************/

void shape_t10 ( double r, double s, double t[10], double dtdr[10], 
  double dtds[10] )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_T10 evaluates shape functions for a 10 node triangle.

  Diagram:

    |
    1  10
    |  |\
    |  | \
    |  8  9
    |  |   \
    S  |    \
    |  5  6  7
    |  |      \
    |  |       \
    0  1--2--3--4
    |
    +--0----R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the reference coordinates of a point.

    Output, double T[10], the basis functions at the point.

    Output, double DTDR[10], the R basis derivatives at the point.

    Output, double DTDS[10], the S basis derivatives at the point.
*/
{
  double a;
  double b;
  double c;

  a = 1.0 / 3.0;
  b = 2.0 / 3.0;
  c = 1.0;

  t[0] = 4.5 * ( a - r - s ) * ( b - r - s ) * ( c - r - s );
  t[1] = 13.5 * r * ( b - r - s ) * ( c - r - s );
  t[2] = - 13.5 * r * ( a - r ) * ( c - r - s );
  t[3] = 4.5 * r * ( a - r ) * ( b - r );
  t[4] = 13.5 * s * ( b - r - s ) * ( c - r - s );
  t[5] = 27.0 * r * s * ( c - r - s );
  t[6] = 13.5 * r * s * ( r - a );
  t[7] = 13.5 * s * ( s - a ) * ( c - r - s );
  t[8] = 13.5 * r * s * ( s - a );
  t[9] = 4.5 * s * ( a - s ) * ( b - s );

  dtdr[0] = 4.5 * ( ( a - s ) * ( 2.0 * r - c - b + 2.0 * s ) 
    - ( s - b ) * ( s - c ) - 2.0 * ( 2.0 * s - b - c ) * r 
    - 3.0 * r * r );
  dtdr[1] = 13.5 * ( 
    ( s - b ) * ( s - c ) + 2.0 * ( 2.0 * s - b - c ) * r 
    + 3.0 * r * r );
  dtdr[2] = - 13.5 * ( a * ( c - s ) + 2.0 * ( s - a - c ) * r 
    + 3.0 * r * r );
  dtdr[3] = 4.5 * ( a * b - 2.0 * ( a + b ) * r + 3.0 * r * r );
  dtdr[4] = 13.5 * s * ( 2.0 * s - b - c + 2.0 * r );
  dtdr[5] = 27.0 * s * ( c - s - 2.0 * r );
  dtdr[6] = 13.5 * s * ( 2.0 * r - a );
  dtdr[7] = - 13.5 * s * ( s - a );
  dtdr[8] = 13.5 * s * ( s - a );
  dtdr[9] = 0.0;

  dtds[0] = 4.5 * ( ( a - r ) * ( 2.0 * s - c - b + 2.0 * r ) 
    - ( r - b ) * ( r - c ) - 2.0 * ( 2.0 * r - b - c ) * s 
    - 3.0 * s * s );
  dtds[1] = 13.5 * r * ( 2.0 * s + 2.0 * r - b - c );
  dtds[2] = 13.5 * r * ( a - r );
  dtds[3] = 0.0;
  dtds[4] = 13.5 * ( ( r - b ) * ( r - c ) + 
    2.0 * ( 2.0 * r - b - c ) * s + 3.0 * s * s );
  dtds[5] = 27.0 * r * ( c - r - 2.0 * s );
  dtds[6] = 13.5 * r * ( r - a );
  dtds[7] = - 13.5 * ( a * ( c - r ) + 2.0 * ( r - c - a ) * s 
    + 3.0 * s * s );
  dtds[8] = 13.5 * r * ( 2.0 * s - a );
  dtds[9] = 4.5 * ( a * b - 2.0 * ( a + b ) * s + 3.0 * s * s );

  return;
}
/******************************************************************************/

void shape_test ( char *code )

/******************************************************************************/
/*
  Purpose: 

    SHAPE_TEST verifies the shape function values at the basis nodes.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element to be used.
    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
    'T3', 'T6' and 'T10'.
*/
{
  double area;
  double *dtdr;
  double *dtds;
  int element_order;
  int i;
  int j;
  double *r;
  double rsum;
  double *s;
  double ssum;
  double *t;

  printf ( "\n" );
  printf ( "  SHAPE_TEST: Verify shape functions of type %s\n", code );

  element_order = order_code ( code );

  dtdr = ( double * ) malloc ( element_order * sizeof ( double ) );
  dtds = ( double * ) malloc ( element_order * sizeof ( double ) );
  r = ( double * ) malloc ( element_order * sizeof ( double ) );
  s = ( double * ) malloc ( element_order * sizeof ( double ) );
  t = ( double * ) malloc ( element_order * sizeof ( double ) );

  node_reference ( code, r, s, &area );

  printf ( "\n" );
  printf ( "  Element order = %d\n", element_order );
  printf ( "  Basis function values at basis nodes\n" );
  printf ( "  should form the identity matrix.\n" );
  printf ( "\n" );

  for ( i = 0; i < element_order; i++ )
  {
    shape ( code, r[i], s[i], t, dtdr, dtds );
    for ( j = 0; j < element_order; j++ )
    {
      printf ( "  %7f", t[j] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  The R and S derivatives should sum to 0.\n" );
  printf ( "\n" );
  printf ( "  dTdR sum, dTdS sum:\n" );
  printf ( "\n" );
  for ( i = 0; i < element_order; i++ )
  {
    shape ( code, r[i], s[i], t, dtdr, dtds );
    rsum = 0.0;
    for ( j = 0; j < element_order; j++ )
    {
      rsum = rsum + dtdr[j];
    }
    ssum = 0.0;
    for ( j = 0; j < element_order; j++ )
    {
      ssum = ssum + dtds[j];
    }
    printf ( "  %14g  %14g\n", rsum, ssum );
  }

  free ( dtdr );
  free ( dtds );
  free ( r );
  free ( s );
  free ( t );

  return;
}
/******************************************************************************/

int sphere_grid_element_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_ELEMENT_NUM returns the number of elements in a sphere grid.

  Discussion:

    The number of elements generated will be NELEMX * NELEMY for
    quadrilaterals, or 2 * NELEMX * ( NELEMY - 1 ) for triangles.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int SPHERE_GRID_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  if ( s_eqi ( code, "Q4" ) )
  {
    element_num = sphere_grid_q4_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    element_num = sphere_grid_q9_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    element_num = sphere_grid_q16_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    element_num = sphere_grid_t3_element_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    element_num = sphere_grid_t6_element_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPHERE_GRID_ELEMENT_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
/******************************************************************************/

int sphere_grid_node_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_NODE_NUM returns the number of nodes in a sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *CODE, identifies the element desired.
    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.

    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
    X and Y directions.  

    Output, int SPHERE_GRID_NODE_NUM, the number of elements in the grid.
*/
{
  int node_num;

  if ( s_eqi ( code, "Q" ) )
  {
    node_num = sphere_grid_q4_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    node_num = sphere_grid_q9_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    node_num = sphere_grid_q16_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    node_num = sphere_grid_t3_node_num ( nelemx, nelemy );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPHERE_GRID_NODE_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
/******************************************************************************/

int *sphere_grid_q4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q4_ELEMENT produces a Q4 sphere grid.

  Discussion:

    This would be the same as the grid in a plane, except that all the
    nodes along the bottom edge are identified (replaced by a single node
    that is the south pole) and similarly for the top edge, and the
    nodes on the extreme right edge are identified pairwise with those 
    on the extreme left edge.

  Example:

    Input:

      NELEMX = 3, NELEMY = 4

    Output:

      ELEMENT_NODE =
         1,  1,  3,  2;
         1,  1,  4,  3;
         1,  1,  2,  4;
         2,  3,  6,  5;
         3,  4,  7,  6;
         4,  2,  5,  7;
         5,  6,  9,  8;
         6,  7, 10,  9;
         7,  5,  8, 10;
         8,  9, 11, 11;
         9, 10, 11, 11;
        10,  8, 11, 11;

  Grid:

   11----11----11----11
    |     |     |     |
    | E10 | E11 | E12 |
    |     |     |     |
    8-----9----10-----8
    |     |     |     |
    | E7  | E8  | E9  |
    |     |     |     |
    5-----6-----7-----5
    |     |     |     |
    | E4  | E5  | E6  |
    |     |     |     |
    2-----3-----4-----2
    |     |     |     |
    | E1  | E2  | E3  |
    |     |     |     |
    1-----1-----1-----1

  Reference Element Q4:

    |
    1  4------3
    |  |      |
    S  |      |
    |  |      |
    |  |      |
    0  1------2
    |
    +--0--R---1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int SPHERE_GRID_Q4_ELEMENT[4*NELEMX*NELEMY], the nodes 
    that form each element.
*/
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 4;
  int i;
  int j;

  element_num = sphere_grid_q4_element_num ( nelemx, nelemy );

  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element = 0;

  for( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * nelemx + 2 - nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + i - 1;

      element_node[0+element*element_order] = base2;
      if ( i < nelemx )
      {
        element_node[1+element*element_order] = base2 + 1;
      }
      else
      {
        element_node[1+element*element_order] = base1;
      }
      element_node[2+element*element_order] = 
        element_node[1+element*element_order] + nelemx;
      element_node[3+element*element_order] = 
        element_node[0+element*element_order] + nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[2+element*element_order] = base1 + nelemx;
        element_node[3+element*element_order] = base1 + nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
/******************************************************************************/

int sphere_grid_q4_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q4_ELEMENT_NUM counts the elements in a Q4 sphere grid.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, int SPHERE_GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int sphere_grid_q4_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:
 
    SPHERE_GRID_Q4_NODE_NUM counts nodes in a Q4 sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_Q4_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = nelemx * ( nelemy - 1 ) + 2;

  return node_num;
}
/******************************************************************************/

double *sphere_grid_q4_node_xyz ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q4_NODE_XYZ produces node coordinates for a Q4 sphere grid.

  Discussion:

    The number of nodes to be generated is

      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, double SPHERE_GRID_Q4_NODE_XYZ[3*NODE_NUM], 
    the node coordinates.
*/
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;    

  for ( j = nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( nelemy );

    for ( i = 1; i <= nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;  
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
/******************************************************************************/

int *sphere_grid_q9_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q9_ELEMENT produces a Q9 sphere grid.

  Discussion:

    This would be the same as the grid in a plane, except that all the
    nodes along the bottom edge are identified (replaced by a single node
    that is the south pole) and similarly for the top edge, and the
    nodes on the extreme right edge are identified pairwise with those 
    on the extreme left edge.

  Example:

    Input:

      NELEMX = 3, NELEMY = 4

    Output:

      ELEMENT_NODE =
         1,  1, 10,  8,  1,  4,  9,  2,  3;
         1,  1, 12, 10,  1,  6, 11,  4,  5;
         1,  1,  8, 12,  1,  2, 13,  6,  7;
         8, 10, 22, 20,  9, 16, 21, 14, 15;
        10, 12, 24, 22, 11, 18, 23, 16, 17;
        12,  8, 20, 24, 13, 14, 25, 18, 19;
        20, 22, 34, 32, 21, 28, 33, 26, 27;
        22, 24, 36, 34, 23, 30, 35, 28, 29;
        24, 20, 32, 36, 25, 26, 37, 30, 31;
        32, 34, 44, 44, 33, 40, 44, 38, 39;
        34, 36, 44, 44, 35, 42, 44, 40, 41;
        36, 32, 44, 44, 37, 38, 44, 42, 43;

  Grid:

   44-44-44-44-44-44-44
    |     |     |     |
   38 39 40 41 42 43 38
    |     |     |     |
   32-33-34-35-36-37-32
    |     |     |     |
   26 27 28 29 30 31 26
    |     |     |     |
   20-21-22-23-24-25-20
    |     |     |     |
   14 15 16 17 18 19 14
    |     |     |     |
    8--9-10-11-12-13--8
    |     |     |     |
    2  3  4  5  6  7  2
    |     |     |     |
    1--1--1--1--1--1--1

  Reference Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_Q9_ELEMENT[9*NELEMX*NELEMY], 
    the nodes that form each element.
*/
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 9;
  int i;
  int j;

  element_num = sphere_grid_q9_element_num ( nelemx, nelemy );
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 2 * ( i - 1 );

      element_node[0+element*element_order] = base2;
      element_node[4+element*element_order] = base2 + 1;

      if ( i < nelemx )
      {
        element_node[1+element*element_order] = base2 + 2;
      }
      else
      {
        element_node[1+element*element_order] = base1;
      }

      element_node[7+element*element_order] = 
        element_node[0+element*element_order] + 2 * nelemx;
      element_node[8+element*element_order] = 
        element_node[4+element*element_order] + 2 * nelemx;
      element_node[5+element*element_order] = 
        element_node[1+element*element_order] + 2 * nelemx;

      element_node[3+element*element_order] = 
        element_node[7+element*element_order] + 2 * nelemx;
      element_node[6+element*element_order] = 
        element_node[8+element*element_order] + 2 * nelemx;
      element_node[2+element*element_order] = 
        element_node[5+element*element_order] + 2 * nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[4+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[3+element*element_order] = base1 + 4 * nelemx;
        element_node[6+element*element_order] = base1 + 4 * nelemx;
        element_node[2+element*element_order] = base1 + 4 * nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
/******************************************************************************/

int sphere_grid_q9_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q9_ELEMENT_NUM counts the elements in a Q9 sphere grid.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, int SPHERE_GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;
 
  return element_num;
}
/******************************************************************************/

int sphere_grid_q9_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q9_NODE_NUM counts nodes in a Q9 sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_Q9_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

  return node_num;
}
/******************************************************************************/

double *sphere_grid_q9_node_xyz ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q9_NODE_XYZ produces node coordinates for a Q9 sphere grid.

  Discussion:

    The number of nodes to be generated is

      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, double SPHERE_GRID_Q9_NODE_XYZ[3*NODE_NUM], 
    the node coordinates.
*/
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_q9_node_num ( nelemx, nelemy );
  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = 2 * nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 2 * nelemy );

    for ( i = 1; i <= 2 * nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 2 * nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
/******************************************************************************/

int *sphere_grid_q16_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q16_ELEMENT produces a Q16 sphere grid.

  Discussion:

    This would be the same as the grid in a plane, except that all the
    nodes along the bottom edge are identified (replaced by a single node
    that is the south pole) and similarly for the top edge, and the
    nodes on the extreme right edge are identified pairwise with those 
    on the extreme left edge.

  Example:

    Input:

      NELEMX = 2, NELEMY = 2

    Output:

      ELEMENT_NODE =
         1,  1,  1,  1,  2,  3,  4,  5,  8,  9, 10, 11, 14, 15, 16, 17;
         1,  1,  1,  1,  5,  6,  7,  2, 11, 12, 13,  8, 17, 18, 19, 14;
        14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29, 32, 32, 32, 32;
        17, 18, 19, 14, 23, 24, 25, 20, 29, 30, 31, 26, 32, 32, 32, 32.

  Grid:

   32-32-32-32-32-32-32
    |        |        |
    |        |        |
   26 27 28 29 30 31 26
    |        |        |
    |        |        |
   20 21 22 23 24 25 20
    |        |        |
    | E3     | E4     |
   14-15-16-17-18-19-14
    |        |        |
    |        |        |
    8  9 10 11 12 13  8
    |        |        |
    |        |        |
    2  3  4  5  6  7  2
    |        |        |
    | E1     | E2     |
    1--1--1--1--1--1--1

  Reference Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |
    0  1---2---3---4
    |
    +--0-----R-----1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  The number of elements generated will be
    NELEMX * NELEMY.

    Output, int SPHERE_GRID_Q16_ELEMENT[16*NELEMX*NELEMY], 
    the nodes that form each element.
*/
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 16;
  int i;
  int j;

  element_num = sphere_grid_q16_element_num ( nelemx, nelemy );
  element_node = ( int * ) malloc ( element_order * nelemx * nelemy * sizeof ( int ) );

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 3 * ( 3 * nelemx ) + 2 - 3 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 3 * ( i - 1 );

      element_node[0+element*element_order] = base2;
      element_node[1+element*element_order] = base2 + 1;
      element_node[2+element*element_order] = base2 + 2;

      if ( i < nelemx )
      {
        element_node[3+element*element_order] = base2 + 3;
      }
      else
      {
        element_node[3+element*element_order] = base1;
      }

      element_node[4+element*element_order] = 
        element_node[0+element*element_order] + 3 * nelemx;
      element_node[5+element*element_order] = 
        element_node[1+element*element_order] + 3 * nelemx;
      element_node[6+element*element_order] = 
        element_node[2+element*element_order] + 3 * nelemx;
      element_node[7+element*element_order] = 
        element_node[3+element*element_order] + 3 * nelemx;

      element_node[8+element*element_order] = 
        element_node[4+element*element_order] + 3 * nelemx;
      element_node[9+element*element_order] = 
        element_node[5+element*element_order] + 3 * nelemx;
      element_node[10+element*element_order] = 
        element_node[6+element*element_order] + 3 * nelemx;
      element_node[11+element*element_order] = 
        element_node[7+element*element_order] + 3 * nelemx;

      element_node[12+element*element_order] = 
        element_node[8+element*element_order] + 3 * nelemx;
      element_node[13+element*element_order] = 
        element_node[9+element*element_order] + 3 * nelemx;
      element_node[14+element*element_order] = 
        element_node[10+element*element_order] + 3 * nelemx;
      element_node[15+element*element_order] = 
        element_node[11+element*element_order] + 3 * nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
        element_node[2+element*element_order] = 1;
        element_node[3+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[12+element*element_order] = base1 + 9 * nelemx;
        element_node[13+element*element_order] = base1 + 9 * nelemx;
        element_node[14+element*element_order] = base1 + 9 * nelemx;
        element_node[15+element*element_order] = base1 + 9 * nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
/******************************************************************************/

int sphere_grid_q16_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q16_ELEMENT_NUM counts the elements in a Q16 sphere grid.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NUM = NELEMX * NELEMY = 6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, int SPHERE_GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
/******************************************************************************/

int sphere_grid_q16_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q16_NODE_NUM counts nodes in a Q16 sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_Q16_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 9 * nelemx * nelemy - 3 * nelemx + 2;

  return node_num;
}
/******************************************************************************/

double *sphere_grid_q16_node_xyz ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_Q16_NODE_XYZ produces node coordinates for a Q16 sphere grid.

  Discussion:

    The number of nodes to be generated is

      NODE_NUM = 9 * NELEMX * NELEMY - 3 * NELEMX + 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, double SPHERE_GRID_Q16_NODE_XYZ[3*NODE_NUM], the node coordinates.
*/
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_q16_node_num ( nelemx, nelemy );
  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  node = 0;

  for ( j = 3 * nelemy + 1; 1 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 3 * nelemy );

    if ( j == 1 )
    {
      node_xyz[0+node*3] =  0.0;
      node_xyz[1+node*3] =  0.0;
      node_xyz[2+node*3] =  1.0;
      node = node + 1;
    }
    else if ( j < 3 * nelemy + 1 )
    {
      for ( i = 1; i <= 3 * nelemx; i++ )
      {
        theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 3 * nelemx );

        node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
        node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
        node_xyz[2+node*3] =                 cos ( phi );
        node = node + 1;      
      }
    }
    else if ( j == 3 * nelemy + 1 )
    {
      node_xyz[0+node*3] =  0.0;
      node_xyz[1+node*3] =  0.0;
      node_xyz[2+node*3] = -1.0;
      node = node + 1;
    }
  }
  return node_xyz;
}
/******************************************************************************/

int *sphere_grid_t3_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T3_ELEMENT produces a T3 sphere grid.

  Discussion:

    This would be the same as the grid in a plane, except that all the
    nodes along the bottom edge are identified (replaced by a single node
    that is the south pole) and similarly for the top edge, and the
    nodes on the extreme right edge are identified pairwise with those 
    on the extreme left edge.

  Example:

    Input:

      NELEMX = 3, NELEMY = 4

    Output:

      ELEMENT_NODE =
         1,  3,  2;
         1,  4,  3;
         1,  2,  4;
         2,  3,  5
         6,  5,  3
         3,  4,  6
         7,  6,  4;
         4,  2,  7;
         5,  7,  2;
         5,  6,  8;
         9,  8,  6;
         6,  7,  9;
        10,  9,  7;
         7,  5, 10;
         8, 10,  5;
         8,  9, 11;
         9, 10, 11;
        10,  8, 11;

  Grid:

   11    11    11    11
    | \   | \   | \   |
    |  \  |  \  |  \  |
    |E16\ |E17 \|E18\ |
    8-----9----10-----8
    | \E11| \E13| \E15|
    |  \  |  \  |  \  |
    |E10\ |E12\ |E14\ |
    5-----6-----7-----5
    | \E5 | \E7 | \E9 |
    |  \  |  \  |  \  |
    |E4 \ |E6 \ |E8 \ |
    2-----3-----4-----2
      \E1 | \E2 | \E3 |
       \  |  \  |  \  |
        \ |   \ |   \ |
          1     1     1

  Reference Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_T3_ELEMENT[3*(2*NELEMX*(NELEMY-1)], 
    the nodes that form each element.
*/
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_num = sphere_grid_t3_element_num ( nelemx, nelemy );
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * nelemx + 2 - nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + i - 1;

      sw = base2;
      if ( i < nelemx )
      {
        se = base2 + 1;
      }
      else
      {
        se = base1;
      }
      nw = sw + nelemx;
      ne = se + nelemx;

      if ( j == 1 )
      {
        sw = 1;
        se = 1;
      }
      else if ( j == nelemx )
      {
        nw = base1 + nelemx;
        ne = base1 + nelemx;
      }

      if ( 1 < j )
      {
        element_node[0+element*element_order] = sw;
        element_node[1+element*element_order] = se;
        element_node[2+element*element_order] = nw;
        element = element + 1;
      }

      if ( j < nelemy )
      {
        element_node[0+element*element_order] = ne;
        element_node[1+element*element_order] = nw;
        element_node[2+element*element_order] = se;
        element = element + 1;
      }
    }
  }
  return element_node;
}
/******************************************************************************/

int sphere_grid_t3_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T3_ELEMENT_NUM counts the elements in a T3 sphere grid.

  Example:

    Input:

      NELEMX = 6, NELEMY = 4

    Output:

      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, int SPHERE_GRID_T3_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * ( nelemy - 1 );

  return element_num;
}
/******************************************************************************/

int sphere_grid_t3_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T3_NODE_NUM counts nodes in a T3 sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_T3_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = nelemx * ( nelemy - 1 ) + 2;

  return node_num;
}
/******************************************************************************/

double *sphere_grid_t3_node_xyz ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T3_NODE_XYZ produces node coordinates for a T3 sphere grid.

  Discussion:

    The number of nodes to be generated is

      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, double SPHERE_GRID_T3_NODE_XYZ[3*NODE_NUM], 
    the node coordinates.
*/
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t3_node_num ( nelemx, nelemy );
  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );

  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( nelemy );

    for ( i = 1; i <= nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
/******************************************************************************/

int *sphere_grid_t6_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T6_ELEMENT produces a T6 sphere grid.

  Discussion:

    This would be the same as the grid in a plane, except that all the
    nodes along the bottom edge are identified (replaced by a single node
    that is the south pole) and similarly for the top edge, and the
    nodes on the extreme right edge are identified pairwise with those 
    on the extreme left edge.

  Example:

    Input:

      NELEMX = 3, NELEMY = 4

    Output:

      ELEMENT_NODE =
        10,  8,  1,  9,  3,  4;
        12, 10,  1, 11,  5,  6;
         8, 12,  1, 13,  7,  2;
         8, 10, 20,  9, 15, 14;
        22, 20, 10, 21, 15, 16;
        10, 12, 22, 11, 17, 16;
        24, 22, 12, 23, 17, 18;
        12,  8, 24, 13, 19, 18;
        20, 24,  8, 25, 19, 14;
        20, 22, 32, 21, 27, 26;
        34, 32, 22, 33, 27, 28;
        22, 24, 34, 23, 29, 28;
        36, 34, 24, 35, 29, 30;
        24, 20, 36, 25, 31, 30;
        32, 36, 20, 37, 31, 26;
        32, 34, 44, 33, 39, 38;
        34, 36, 44, 35, 41, 40;
        36, 32, 44, 37, 43, 42;

  Grid:

   44    44    44
    |\    |\    |\
   38 39 40 41 42 43 
    |    \|    \|    \
   32-33-34-35-36-37-32
    |\    |\    |\    |
   26 27 28 29 30 31 26
    |    \|    \|    \|
   20-21-22-23-24-25-20
    |\    |\    |\    |
   14 15 16 17 18 19 14
    |    \|    \|    \|
    8--9-10-11-12-13--8
     \    |\    |\    |
       3  4  5  6  7  2
         \|    \|    \|
          1     1     1

  Reference Element T6:

    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_T6_ELEMENT[6*2*NELEMX*(NELEMY-1)], 
    the nodes that form each element.
*/
{
  int base1;
  int base2;
  int c;
  int e;
  int element;
  int *element_node;
  int element_num;
  int element_order = 6;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_num = sphere_grid_t6_element_num ( nelemx, nelemy );
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 2 * ( i - 1 );

      sw = base2;
      s = base2 + 1;
      if ( i < nelemx )
      {
        se = base2 + 2;
      }
      else
      {
        se = base1;
      }

      w = sw + 2 * nelemx;
      c = s  + 2 * nelemx;
      e = se + 2 * nelemx;

      nw = w + 2 * nelemx;
      n  = c + 2 * nelemx;
      ne = e + 2 * nelemx;

      if ( j == 1 )
      {
        sw = 1;
        s  = 1;
        se = 1;
      }
      else if ( j == nelemy )
      {
        nw = base1 + 4 * nelemx;
        n  = base1 + 4 * nelemx;
        ne = base1 + 4 * nelemx;
      }

      if ( 1 < j )
      {
        element_node[0+element*element_order] = sw;
        element_node[1+element*element_order] = se;
        element_node[2+element*element_order] = nw;
        element_node[3+element*element_order] = s;
        element_node[4+element*element_order] = c;
        element_node[5+element*element_order] = w;
        element = element + 1;
      }

      if ( j < nelemy )
      {
        element_node[0+element*element_order] = ne;
        element_node[1+element*element_order] = nw;
        element_node[2+element*element_order] = se;
        element_node[3+element*element_order] = n;
        element_node[4+element*element_order] = c;
        element_node[5+element*element_order] = e;
        element = element + 1;
      }
    }
  }
  return element_node;
}
/******************************************************************************/

int sphere_grid_t6_element_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T6_ELEMENT_NUM counts the elements in a T6 sphere grid.

  Example:

    Input:

      NELEMX = 6, NELEMY = 4

    Output:

      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions. 

    Output, int SPHERE_GRID_T6_ELEMENT_NUM, the number of elements in the grid.
*/
{
  int element_num;

  element_num = 2 * nelemx * ( nelemy - 1 );

  return element_num;
}
/******************************************************************************/

int sphere_grid_t6_node_num ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T6_NODE_NUM counts nodes in a T6 sphere grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.

    Output, int SPHERE_GRID_T6_NODE_NUM, the number of nodes in the grid.
*/
{
  int node_num;

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

  return node_num;
}
/******************************************************************************/

double *sphere_grid_t6_node_xyz ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    SPHERE_GRID_T6_NODE_XYZ produces node coordinates for a T6 sphere grid.

  Discussion:

    The number of nodes to be generated is

      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2006

  Author:

    John Burkardt

  Parameters:

    Input, int NELEMX, NELEMY, the number of elements along the
    X and Y directions.  

    Output, double SPHERE_GRID_T6_NODE_XYZ[3*NODE_NUM], 
    the node coordinates.
*/
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = 2 * nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 2 * nelemy );

    for ( i = 1; i <= 2 * nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 2 * nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
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
/******************************************************************************/

void triangle_unit_set ( int rule, int order, double xtab[], double ytab[], 
  double weight[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.

  Discussion:

    The user is responsible for determining the value of ORDER,
    and appropriately dimensioning the arrays XTAB, YTAB and
    WEIGHT so that they can accommodate the data.

    The value of ORDER for each rule can be found by invoking
    the function TRIANGLE_RULE_SIZE.

  Integration region:

      0 <= X,
    and
      0 <= Y, 
    and
      X + Y <= 1.

  Graph:

      ^
    1 | *
      | |\
    Y | | \
      | |  \
    0 | *---*
      +------->
        0 X 1

   The rules are accessed by an index number, RULE.  The indices,
   and the descriptions of the corresponding rules, are:

     1, ORDER =  1, precision 1, Zienkiewicz #1.
     2, ORDER =  2, precision 1, (the "vertex rule").
     3, ORDER =  3, precision 2, Strang and Fix formula #1.
     4, ORDER =  3, precision 2, Strang and Fix formula #2,
                                 Zienkiewicz #2.
     5, ORDER =  4, precision 3, Strang and Fix formula #3,
                                 Zienkiewicz #3.
     6, ORDER =  6, precision 3, Strang and Fix formula #4.
     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
     8, ORDER =  6, precision 4, Strang and Fix formula #5.
     9, ORDER =  7, precision 4, Strang and Fix formula #6.
    10, ORDER =  7, precision 5, Strang and Fix formula #7,
                                 Stroud formula T2:5-1, 
                                 Zienkiewicz #4, 
                                 Schwarz Table 2.2.
    11, ORDER =  9, precision 6, Strang and Fix formula #8.
    12, ORDER = 12, precision 6, Strang and Fix formula #9.
    13, ORDER = 13, precision 7, Strang and Fix formula #10.
        Note that there is a typographical error in Strang and Fix
        which lists the value of the XSI(3) component of the
        last generator point as 0.4869... when it should be 0.04869...
    14, ORDER =  7, precision 3.
    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
    16, ORDER = 64, precision 15, triangular product Gauss rule.
    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
    20, ORDER = 37, precision 13, from ACM TOMS #706.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 March 2008

  Author:

    John Burkardt

  Reference:

    Jarle Berntsen, Terje Espelid,
    Algorithm 706,
    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
    ACM Transactions on Mathematical Software,
    Volume 18, Number 3, September 1992, pages 329-342.

    Elise deDoncker, Ian Robinson,
    Algorithm 612:
    Integration over a Triangle Using Nonlinear Extrapolation,
    ACM Transactions on Mathematical Software,
    Volume 10, Number 1, March 1984, pages 17-22.

    Dirk Laurie,
    Algorithm 584,
    CUBTRI, Automatic Cubature Over a Triangle,
    ACM Transactions on Mathematical Software,
    Volume 8, Number 2, 1982, pages 210-218.

    James Lyness, Dennis Jespersen,
    Moderate Degree Symmetric Quadrature Rules for the Triangle,
    Journal of the Institute of Mathematics and its Applications,
    Volume 15, Number 1, February 1975, pages 19-32.

    Hans Rudolf Schwarz,
    Finite Element Methods,
    Academic Press, 1988,
    ISBN: 0126330107,
    LC: TA347.F5.S3313.

    Gilbert Strang, George Fix,
    An Analysis of the Finite Element Method,
    Cambridge, 1973,
    ISBN: 096140888X,
    LC: TA335.S77.

    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.

    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    LC: TA640.2.Z54

  Parameters:

    Input, int RULE, the index of the rule.

    Input, int ORDER, the order of the rule.

    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.

    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  int i;
  int j;
  int k;
  int order2;
  double p;
  double q;
  double r;
  double s;
  double t;
  double u;
  double v;
  double w;
  double w1;
  double w2;
  double w3;
  double w4;
  double w5;
  double w6;
  double w7;
  double w8;
  double w9;
  double weight1[8];
  double weight2[8];
  double wx;
  double x;
  double xtab1[8];
  double xtab2[8];
  double y;
  double z;
/*
  1 point, precision 1.
*/
  if ( rule == 1 )
  {
    xtab[0]   = 0.33333333333333333333;

    ytab[0]   = 0.33333333333333333333;

    weight[0] = 1.00000000000000000000;
  }
/*
  3 points, precision 1, the "vertex rule".
*/
  else if ( rule == 2 )
  {
    xtab[0] =   1.00000000000000000000;
    xtab[1] =   0.00000000000000000000;
    xtab[2] =   0.00000000000000000000;

    ytab[0] =   0.00000000000000000000;
    ytab[1] =   1.00000000000000000000;
    ytab[2] =   0.00000000000000000000;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
/*
  3 points, precision 2, Strang and Fix formula #1.
*/
  else if ( rule == 3 )
  {
    xtab[0]   = 0.66666666666666666667;
    xtab[1]   = 0.16666666666666666667;
    xtab[2]   = 0.16666666666666666667;

    ytab[0]   = 0.16666666666666666667;
    ytab[1]   = 0.66666666666666666667;
    ytab[2]   = 0.16666666666666666667;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
/*
  3 points, precision 2, Strang and Fix formula #2.
*/
  else if ( rule == 4 )
  {
    xtab[0]   = 0.50000000000000000000;
    xtab[1]   = 0.50000000000000000000;
    xtab[2]   = 0.00000000000000000000;

    ytab[0]   = 0.00000000000000000000;
    ytab[1]   = 0.50000000000000000000;
    ytab[2]   = 0.50000000000000000000;

    weight[0] = 0.33333333333333333333;
    weight[1] = 0.33333333333333333333;
    weight[2] = 0.33333333333333333333;
  }
/*
  4 points, precision 3, Strang and Fix formula #3.
*/
  else if ( rule == 5 )
  {
    a =   6.0 / 30.0;
    b =  10.0 / 30.0;
    c =  18.0 / 30.0;

    d =  25.0 / 48.0;
    e = -27.0 / 48.0;

    xtab[0] = b;
    xtab[1] = c;
    xtab[2] = a;
    xtab[3] = a;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = c;
    ytab[3] = a;

    weight[0] = e;
    weight[1] = d;
    weight[2] = d;
    weight[3] = d;
  }
/*
  6 points, precision 3, Strang and Fix formula #4.
*/
  else if ( rule == 6 )
  {
    a = 0.659027622374092;
    b = 0.231933368553031;
    c = 0.109039009072877;

    xtab[0] = a;
    xtab[1] = a;
    xtab[2] = b;
    xtab[3] = b;
    xtab[4] = c;
    xtab[5] = c;

    ytab[0] = b;
    ytab[1] = c;
    ytab[2] = a;
    ytab[3] = c;
    ytab[4] = a;
    ytab[5] = b;

    weight[0] = 0.16666666666666666667;
    weight[1] = 0.16666666666666666667;
    weight[2] = 0.16666666666666666667;
    weight[3] = 0.16666666666666666667;
    weight[4] = 0.16666666666666666667;
    weight[5] = 0.16666666666666666667;
  }
/*
  6 points, precision 3, Stroud T2:3-1.
*/
  else if ( rule == 7 )
  {
    a = 0.0;
    b = 0.5;
    c = 2.0 /  3.0;
    d = 1.0 /  6.0;
    v = 1.0 / 30.0;
    w = 3.0 / 10.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;

    weight[0] = v;
    weight[1] = v;
    weight[2] = v;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
  }
/*
  6 points, precision 4, Strang and Fix, formula #5.
*/
  else if ( rule == 8 )
  {
    a = 0.816847572980459;
    b = 0.091576213509771;
    c = 0.108103018168070;
    d = 0.445948490915965;
    v = 0.109951743655322;
    w = 0.223381589678011;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;

    weight[0] = v;
    weight[1] = v;
    weight[2] = v;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
  }
/*
  7 points, precision 4, Strang and Fix formula #6.
*/
  else if ( rule == 9 )
  {
    a = 1.0 / 3.0;
    c = 0.736712498968435;
    d = 0.237932366472434;
    e = 0.025355134551932;
    v = 0.375000000000000;
    w = 0.104166666666667;

    xtab[0] = a;
    xtab[1] = c;
    xtab[2] = c;
    xtab[3] = d;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;

    ytab[0] = a;
    ytab[1] = d;
    ytab[2] = e;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = c;
    ytab[6] = d;

    weight[0] = v;
    weight[1] = w;
    weight[2] = w;
    weight[3] = w;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
/*
  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
*/
  else if ( rule == 10 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    u = 0.225;
    v = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    w = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;

    weight[0] = u;
    weight[1] = v;
    weight[2] = v;
    weight[3] = v;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
/*
  9 points, precision 6, Strang and Fix formula #8.
*/
  else if ( rule == 11 )
  {
    a = 0.124949503233232;
    b = 0.437525248383384;
    c = 0.797112651860071;
    d = 0.165409927389841;
    e = 0.037477420750088;

    u = 0.205950504760887;
    v = 0.063691414286223;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = c;
    xtab[5] = d;
    xtab[6] = d;
    xtab[7] = e;
    xtab[8] = e;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = e;
    ytab[5] = c;
    ytab[6] = e;
    ytab[7] = c;
    ytab[8] = d;

    weight[0] = u;
    weight[1] = u;
    weight[2] = u;
    weight[3] = v;
    weight[4] = v;
    weight[5] = v;
    weight[6] = v;
    weight[7] = v;
    weight[8] = v;
  }
/*
  12 points, precision 6, Strang and Fix, formula #9.
*/
  else if ( rule == 12 )
  {
    a = 0.873821971016996;
    b = 0.063089014491502;
    c = 0.501426509658179;
    d = 0.249286745170910;
    e = 0.636502499121399;
    f = 0.310352451033785;
    g = 0.053145049844816;

    u = 0.050844906370207;
    v = 0.116786275726379;
    w = 0.082851075618374;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = b;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = d;
    xtab[6] = e;
    xtab[7] = e;
    xtab[8] = f;
    xtab[9] = f;
    xtab[10] = g;
    xtab[11] = g;

    ytab[0] = b;
    ytab[1] = a;
    ytab[2] = b;
    ytab[3] = d;
    ytab[4] = c;
    ytab[5] = d;
    ytab[6] = f;
    ytab[7] = g;
    ytab[8] = e;
    ytab[9] = g;
    ytab[10] = e;
    ytab[11] = f;

    weight[0] = u;
    weight[1] = u;
    weight[2] = u;
    weight[3] = v;
    weight[4] = v;
    weight[5] = v;
    weight[6] = w;
    weight[7] = w;
    weight[8] = w;
    weight[9] = w;
    weight[10] = w;
    weight[11] = w;
  }
/*
  13 points, precision 7, Strang and Fix, formula #10.

  Note that there is a typographical error in Strang and Fix
  which lists the value of the XSI[2] component of the
  last generator point as 0.4869... when it should be 0.04869...
*/
  else if ( rule == 13 )
  {
    h = 1.0 / 3.0;
    a = 0.479308067841923;
    b = 0.260345966079038;
    c = 0.869739794195568;
    d = 0.065130102902216;
    e = 0.638444188569809;
    f = 0.312865496004875;
    g = 0.048690315425316;

    w = -0.149570044467670;
    t =  0.175615257433204;
    u =  0.053347235608839;
    v =  0.077113760890257;

    xtab[0] = h;
    xtab[1] = a;
    xtab[2] = b;
    xtab[3] = b;
    xtab[4] = c;
    xtab[5] = d;
    xtab[6] = d;
    xtab[7] = e;
    xtab[8] = e;
    xtab[9] = f;
    xtab[10] = f;
    xtab[11] = g;
    xtab[12] = g;

    ytab[0] = h;
    ytab[1] = b;
    ytab[2] = a;
    ytab[3] = b;
    ytab[4] = d;
    ytab[5] = c;
    ytab[6] = d;
    ytab[7] = f;
    ytab[8] = g;
    ytab[9] = e;
    ytab[10] = g;
    ytab[11] = e;
    ytab[12] = f;

    weight[0] = w;
    weight[1] = t;
    weight[2] = t;
    weight[3] = t;
    weight[4] = u;
    weight[5] = u;
    weight[6] = u;
    weight[7] = v;
    weight[8] = v;
    weight[9] = v;
    weight[10] = v;
    weight[11] = v;
    weight[12] = v;
  }
/*
  7 points, precision 3.
*/
  else if ( rule == 14 )
  {
    a = 1.0 / 3.0;
    b = 1.0;
    c = 0.5;
    z = 0.0;

    u = 27.0 / 60.0;
    v =  3.0 / 60.0;
    w =  8.0 / 60.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = z;
    xtab[3] = z;
    xtab[4] = z;
    xtab[5] = c;
    xtab[6] = c;

    ytab[0] = a;
    ytab[1] = z;
    ytab[2] = b;
    ytab[3] = z;
    ytab[4] = c;
    ytab[5] = z;
    ytab[6] = c;

    weight[0] = u;
    weight[1] = v;
    weight[2] = v;
    weight[3] = v;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
  }
/*
  16 points, precision 5, Stroud T2:7-1.
*/
  else if ( rule == 15 )
  {
/*
  Legendre rule of order 4.
*/
    order2 = 4;

    xtab[0] = - 0.861136311594052575223946488893;
    xtab[1] = - 0.339981043584856264802665759103;
    xtab[2] =   0.339981043584856264802665759103;
    xtab[3] =   0.861136311594052575223946488893;

    weight1[0] = 0.347854845137453857373063949222;
    weight1[1] = 0.652145154862546142626936050778;
    weight1[2] = 0.652145154862546142626936050778;
    weight1[3] = 0.347854845137453857373063949222;

    for ( i = 0; i < order2; i++ )
    {
      xtab1[i] = 0.5 * ( xtab1[i] + 1.0 );
    }

    weight2[0] = 0.1355069134;
    weight2[1] = 0.2034645680;
    weight2[2] = 0.1298475476;
    weight2[3] = 0.0311809709;

    xtab2[0] = 0.0571041961;
    xtab2[1] = 0.2768430136;
    xtab2[2] = 0.5835904324;
    xtab2[3] = 0.8602401357;

    k = 0;
    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = xtab2[j];
        ytab[k] = xtab1[i] * ( 1.0 - xtab2[j] );
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
/*
  64 points, precision 15.
*/
  else if ( rule == 16 )
  {
/*
  Legendre rule of order 8.
*/
    order2 = 8;

    xtab1[0] = -0.960289856497536231683560868569;
    xtab1[1] = -0.796666477413626739591553936476;
    xtab1[2] = -0.525532409916328985817739049189;
    xtab1[3] = -0.183434642495649804939476142360;
    xtab1[4] =  0.183434642495649804939476142360;
    xtab1[5] =  0.525532409916328985817739049189;
    xtab1[6] =  0.796666477413626739591553936476;
    xtab1[7] =  0.960289856497536231683560868569;

    weight1[0] = 0.101228536290376259152531354310;
    weight1[1] = 0.222381034453374470544355994426;
    weight1[2] = 0.313706645877887287337962201987;
    weight1[3] = 0.362683783378361982965150449277;
    weight1[4] = 0.362683783378361982965150449277;
    weight1[5] = 0.313706645877887287337962201987;
    weight1[6] = 0.222381034453374470544355994426;
    weight1[7] = 0.101228536290376259152531354310;

    weight2[0] = 0.00329519144;
    weight2[1] = 0.01784290266;
    weight2[2] = 0.04543931950;
    weight2[3] = 0.07919959949;
    weight2[4] = 0.10604735944;
    weight2[5] = 0.11250579947;
    weight2[6] = 0.09111902364;
    weight2[7] = 0.04455080436;

    xtab2[0] = 0.04463395529;
    xtab2[1] = 0.14436625704;
    xtab2[2] = 0.28682475714;
    xtab2[3] = 0.45481331520;
    xtab2[4] = 0.62806783542;
    xtab2[5] = 0.78569152060;
    xtab2[6] = 0.90867639210;
    xtab2[7] = 0.98222008485;

    k = 0;
    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = 1.0 - xtab2[j];
        ytab[k] = 0.5 * ( 1.0 + xtab1[i] ) * xtab2[j];
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
/*
  19 points, precision 8, from CUBTRI.
*/
  else if ( rule == 17 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    f = ( 40.0 - 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
    g = ( 25.0 +  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
    p = ( 40.0 + 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
    q = ( 25.0 -  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
    r = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
    s = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      - sqrt ( 105.0 ) ) / 90.0;
    t = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      + sqrt ( 105.0 ) ) / 90.0;

    w1 = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
    w2 = -9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      + 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w2 = w2 / 3.0;
    w3 = -9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 
      + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      - 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w3 = w3 / 3.0;
    w4 = ( 102791225.0 - 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      + 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0;
    w4 = w4 / 3.0;
    w5 = ( 102791225.0 + 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      - 9914825 * sqrt ( 105.0 ) ) / 59157504.0;
    w5 = w5 / 3.0;
    w6 = ( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0;
    w6 = w6 / 6.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = r;
    xtab[15] = s;
    xtab[16] = s;
    xtab[17] = t;
    xtab[18] = t;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = t;
    ytab[15] = r;
    ytab[16] = t;
    ytab[17] = r;
    ytab[18] = s;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
/*
  19 points, precision 9.
  Lyness and Jesperson.
*/
  else if ( rule == 18 )
  {
    a = 1.0 / 3.0;
    b =  0.02063496160252593;
    c =  0.4896825191987370;
    d =  0.1258208170141290;
    e =  0.4370895914929355;
    f =  0.6235929287619356;
    g =  0.1882035356190322;
    r =  0.9105409732110941;
    s =  0.04472951339445297;
    t =  0.7411985987844980;
    u =  0.03683841205473626;
    v =  0.22196298916076574;

    w1 = 0.09713579628279610;
    w2 = 0.03133470022713983;
    w3 = 0.07782754100477543;
    w4 = 0.07964773892720910;
    w5 = 0.02557767565869810;
    w6 = 0.04328353937728940;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = r;
    xtab[11] = s;
    xtab[12] = s;
    xtab[13] = t;
    xtab[14] = t;
    xtab[15] = u;
    xtab[16] = u;
    xtab[17] = v;
    xtab[18] = v;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = s;
    ytab[10] = r;
    ytab[11] = s;
    ytab[12] = u;
    ytab[13] = v;
    ytab[14] = t;
    ytab[15] = t;
    ytab[16] = v;
    ytab[17] = t;
    ytab[18] = u;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
/*
  28 points, precision 11.
  Lyness and Jesperson.
*/
  else if ( rule == 19 )
  {
    a = 1.0 / 3.0;
    b = 0.9480217181434233;
    c = 0.02598914092828833;
    d = 0.8114249947041546;
    e = 0.09428750264792270;
    f = 0.01072644996557060;
    g = 0.4946367750172147;
    p = 0.5853132347709715;
    q = 0.2073433826145142;
    r = 0.1221843885990187;
    s = 0.4389078057004907;
    t = 0.6779376548825902;
    u = 0.04484167758913055;
    v = 0.27722066752827925;
    w = 0.8588702812826364;
    x = 0.0;
    y = 0.1411297187173636;

    w1 = 0.08797730116222190;
    w2 = 0.008744311553736190;
    w3 = 0.03808157199393533;
    w4 = 0.01885544805613125;
    w5 = 0.07215969754474100;
    w6 = 0.06932913870553720;
    w7 = 0.04105631542928860;
    w8 = 0.007362383783300573;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = s;
    xtab[15] = s;
    xtab[16] = t;
    xtab[17] = t;
    xtab[18] = u;
    xtab[19] = u;
    xtab[20] = v;
    xtab[21] = v;
    xtab[22] = w;
    xtab[23] = w;
    xtab[24] = x;
    xtab[25] = x;
    xtab[26] = y;
    xtab[27] = y;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = r;
    ytab[15] = s;
    ytab[16] = u;
    ytab[17] = v;
    ytab[18] = t;
    ytab[19] = v;
    ytab[20] = t;
    ytab[21] = u;
    ytab[22] = x;
    ytab[23] = y;
    ytab[24] = w;
    ytab[25] = y;
    ytab[26] = w;
    ytab[27] = x;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w7;
    weight[20] = w7;
    weight[21] = w7;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w8;
    weight[26] = w8;
    weight[27] = w8;
  }
/*
  37 points, precision 13.
*/
  else if ( rule == 20 )
  {
    a = 1.0 / 3.0;
    b = 0.950275662924105565450352089520;
    c = 0.024862168537947217274823955239;
    d = 0.171614914923835347556304795551;
    e = 0.414192542538082326221847602214;
    f = 0.539412243677190440263092985511;
    g = 0.230293878161404779868453507244;

    w1 = 0.051739766065744133555179145422;
    w2 = 0.008007799555564801597804123460;
    w3 = 0.046868898981821644823226732071;
    w4 = 0.046590940183976487960361770070;
    w5 = 0.031016943313796381407646220131;
    w6 = 0.010791612736631273623178240136;
    w7 = 0.032195534242431618819414482205;
    w8 = 0.015445834210701583817692900053;
    w9 = 0.017822989923178661888748319485;
    wx = 0.037038683681384627918546472190;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w8;
    weight[20] = w8;
    weight[21] = w8;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w9;
    weight[26] = w9;
    weight[27] = w9;
    weight[28] = w9;
    weight[29] = w9;
    weight[30] = w9;
    weight[31] = wx;
    weight[32] = wx;
    weight[33] = wx;
    weight[34] = wx;
    weight[35] = wx;
    weight[36] = wx;

    a = 0.772160036676532561750285570113;
    b = 0.113919981661733719124857214943;

    xtab[10] = a;
    ytab[10] = b;

    xtab[11] = b;
    ytab[11] = a;

    xtab[12] = b;
    ytab[12] = b;

    a = 0.009085399949835353883572964740;
    b = 0.495457300025082323058213517632;

    xtab[13] = a;
    ytab[13] = b;

    xtab[14] = b;
    ytab[14] = a;

    xtab[15] = b;
    ytab[15] = b;

    a = 0.062277290305886993497083640527;
    b = 0.468861354847056503251458179727;

    xtab[16] = a;
    ytab[16] = b;

    xtab[17] = b;
    ytab[17] = a;

    xtab[18] = b;
    ytab[18] = b;

    a = 0.022076289653624405142446876931;
    b = 0.851306504174348550389457672223;
    c = 0.126617206172027096933163647918263;

    xtab[19] = a;
    ytab[19] = b;

    xtab[20] = a;
    ytab[20] = c;

    xtab[21] = b;
    ytab[21] = a;

    xtab[22] = b;
    ytab[22] = c;

    xtab[23] = c;
    ytab[23] = a;

    xtab[24] = c;
    ytab[24] = b;

    a = 0.018620522802520968955913511549;
    b = 0.689441970728591295496647976487;
    c = 0.291937506468887771754472382212953;

    xtab[25] = a;
    ytab[25] = b;

    xtab[26] = a;
    ytab[26] = c;

    xtab[27] = b;
    ytab[27] = a;

    xtab[28] = b;
    ytab[28] = c;

    xtab[29] = c;
    ytab[29] = a;

    xtab[30] = c;
    ytab[30] = b;

    a = 0.096506481292159228736516560903;
    b = 0.635867859433872768286976979827;
    c = 0.267625659273967961282458816185681;

    xtab[31] = a;
    ytab[31] = b;

    xtab[32] = a;
    ytab[32] = c;

    xtab[33] = b;
    ytab[33] = a;

    xtab[34] = b;
    ytab[34] = c;

    xtab[35] = c;
    ytab[35] = a;

    xtab[36] = c;
    ytab[36] = b;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRIANGLE_UNIT_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of RULE = %d\n", rule );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

int triangle_unit_size ( int rule )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 March 2008

  Author:

    John Burkardt

  Reference:

    Jarle Berntsen, Terje Espelid,
    Algorithm 706,
    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
    ACM Transactions on Mathematical Software,
    Volume 18, Number 3, September 1992, pages 329-342.

    Elise deDoncker, Ian Robinson,
    Algorithm 612:
    Integration over a Triangle Using Nonlinear Extrapolation,
    ACM Transactions on Mathematical Software,
    Volume 10, Number 1, March 1984, pages 17-22.

    DP Laurie,
    Algorithm 584,
    CUBTRI, Automatic Cubature Over a Triangle,
    ACM Transactions on Mathematical Software,
    Volume 8, Number 2, 1982, pages 210-218.

    James Lyness, Dennis Jespersen,
    Moderate Degree Symmetric Quadrature Rules for the Triangle,
    Journal of the Institute of Mathematics and its Applications,
    Volume 15, Number 1, February 1975, pages 19-32.

    Hans Rudolf Schwarz,
    Methode der Finiten Elemente,
    Teubner Studienbuecher, 1980,
    ISBN: 3-519-02349-0.

    Gilbert Strang, George Fix,
    An Analysis of the Finite Element Method,
    Prentice Hall, 1973, page 184,
    ISBN: 096140888X,
    LC: TA335.S77.

    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.

    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    TA640.2.Z54

  Parameters:

    Input, int RULE, the index of the rule.
     1, ORDER =  1, precision 1, Zienkiewicz #1.
     2, ORDER =  2, precision 1, (the "vertex rule").
     3, ORDER =  3, precision 2, Strang and Fix formula #1.
     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
     6, ORDER =  6, precision 3, Strang and Fix formula #4.
     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
     8, ORDER =  6, precision 4, Strang and Fix formula #5.
     9, ORDER =  7, precision 4, Strang and Fix formula #6.
    10, ORDER =  7, precision 5, Strang and Fix formula #7,
        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
    11, ORDER =  9, precision 6, Strang and Fix formula #8.
    12, ORDER = 12, precision 6, Strang and Fix formula #9.
    13, ORDER = 13, precision 7, Strang and Fix formula #10.
    14, ORDER =  7, precision ?.
    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
    16, ORDER = 64, precision 15, triangular product Gauss rule.
    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
    20, ORDER = 37, precision 13, from ACM TOMS #706.

    Output, int TRIANGLE_UNIT_SIZE, the order of the rule.
*/
{
  int size;

  if ( rule == 1 )
  {
    size = 1;
  }
  else if ( rule == 2 )
  {
    size = 3;
  }
  else if ( rule == 3 )
  {
    size = 3;
  }
  else if ( rule == 4 )
  {
    size = 3;
  }
  else if ( rule == 5 )
  {
    size = 4;
  }
  else if ( rule == 6 )
  {
    size = 6;
  }
  else if ( rule == 7 )
  {
    size = 6;
  }
  else if ( rule == 8 )
  {
    size = 6;
  }
  else if ( rule == 9 )
  {
    size = 7;
  }
  else if ( rule == 10 )
  {
    size = 7;
  }
  else if ( rule == 11 )
  {
    size = 9;
  }
  else if ( rule == 12 )
  {
    size = 12;
  }
  else if ( rule == 13 )
  {
    size = 13;
  }
  else if ( rule == 14 )
  {
    size = 7;
  }
  else if ( rule == 15 )
  {
    size = 16;
  }
  else if ( rule == 16 )
  {
    size = 64;
  }
  else if ( rule == 17 )
  {
    size = 19;
  }
  else if ( rule == 18 )
  {
    size = 19;
  }
  else if ( rule == 19 )
  {
    size = 28;
  }
  else if ( rule == 20 )
  {
    size = 37;
  }
  else
  {
    size = -1;
  }

  return size;
}
