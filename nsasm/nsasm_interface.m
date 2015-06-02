function [ K, L ] = nsasm_interface ( p_file, t_file, e_file, np0, nu )

%*****************************************************************************80
%
%% NSASM_INTERFACE reads data and sends it to NSASM to compute the FEM system.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    25 January 2014
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Per-Olof Persson,
%    Implementation of Finite Element-Based Navier-Stokes Solver,
%    April 2002.
%
%  Parameters:
%
%    Input, string P_FILE, the name of the "point file", which contains 
%    2 rows and NP columns of (X,Y) node coordinates.
%
%    Input, string T_FILE, the name of the "triangle file", which contains 
%    6 rows and NT columns, with each column containing the (1-based) indices 
%    of nodes forming the quadratic triangle elements.  Note that the nodes
%    are ordered in a particular order.
%
%    Input, string E_FILE, the name of the "constraint file", which contains 
%    3 rows and NE columns, defining constraints on the data, including 
%    Dirichlet boundary values in particular.  Item #1 is a node index, 
%    #2 is a variable index (0 = horizontal velocity, 1 = vertical velocity, 
%    2 = pressure) and #3 is an associated value.
%
%    Input, integer NP0, the number of pressure nodes.
%
%    Input, integer NU, the kinematic viscosity.
%
%    Output, real sparse K(2*NP+NP0+NE,2*NP+NP0+NE), the stiffness matrix.
%
%    Output, real L(2*NP+NP0+NE), the residual vector.
%
%  Local Parameters:
%
%    Local, real E(3,NE), the constraint data:
%    Item #1 is a node index, 
%    Item #2 is 0 = horizontal velocity, 1 = vertical velocity, or 2 = pressure,
%    Item #3 is an associated value.
%
%    Local, sparse real K(NDOF,NDOF), the stiffness matrix, computed by NSASM.
%
%    Local, real L(NDOF), the residual for the finite element equations,
%    based on the current solution estimate U, as computed by NSASM.
%
%    Local, integer NDOF, the number of degrees of freedom.
%
%    Local, integer NE, the number of constraints.
%
%    Local, integer NP, the number of nodes.
%
%    Local, integer NT, the number of triangles.
%
%    Local, real P(2,NP), the coordinates of the nodes.
%
%    Local, integer T(6,NT), the indices of nodes that form triangular elements.
%
%    Local, real U0(1,NDOF), the solution estimate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'NSASM_INTERFACE:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Loading user node data from "%s".\n', p_file );

  p = load ( p_file );
  p = p';
  [ dim_num, np ] = size ( p );
  fprintf ( 1, '  Number of nodes = %d\n', np );
  fprintf ( 1, '  Number of pressure nodes NP0 = %d\n', np0 );

  fprintf ( 1, '  Loading user triangle data from "%s".\n', t_file );

  t = load ( t_file );
  t = t';
  [ triangle_order, nt ] = size ( t );
  fprintf ( 1, '  Number of triangles NT = %d\n', nt );

  fprintf ( 1, '  Loading user constraint data from "%s".\n', e_file );

  e = load ( e_file );
  [ e_order, ne ] = size ( e );
  fprintf ( 1, '  Number of constraints NE = %d\n', ne );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Viscosity NU = %f\n', nu );
%
%  U is the set of finite element coefficients.
%  We set this to zero.
%  In typical usage, NSASM is used inside a Newton iteration,
%  and we are trying to improve the estimate of U.
%
  ndof = 2 * np + np0 + ne;
  fprintf ( 1, '  Degrees of freedom NDOF =  %d\n', ndof );

  u0 = zeros ( 1, ndof );

  [ K, L ] = nsasm ( p, t, np0, e, u0, nu );

  fprintf ( 1, '\n' );
  fprintf ( 1, 'NSASM_INTERFACE:\n' );
  fprintf ( 1, '  Returning K, L data from NSASM.\n' );

  return
end
