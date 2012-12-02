function small_test ( )

%*****************************************************************************80
%
%% SMALL_TEST is a small sample program for NSASM.
%
%  Discussion:
%
%    The flow region uses a 5 x 5 grid of nodes indexed as follows:
%
%    1.0
%     ^       7  19   8  23   9
%     |      21  20  22  24  25
%     Y       4  10   5  15   6
%     |      12  11  13  16  17
%     v       1  14   2   18  3
%    0.0
%            0.0  <-- X -->  1.0
%
%    All 16 boundary nodes are constrained to have zero horizontal
%    and vertical velocities.
%
%    Node 1 is constrained to have a pressure of 1.
%
%    The 8 elements are quadratic six node triangles:
%
%      1: 1, 4, 5, 10, 11, 12.
%      2: 1, 2, 5, 13, 11, 14
%      3: 2, 5, 6, 15, 16, 13
%      4: 2, 3, 6, 17, 16, 18
%      5: 4, 7, 8, 19, 20, 21
%      6: 4, 5, 8, 22, 20, 10
%      7: 5, 8, 9, 23, 24, 22
%      8: 5, 6, 9, 25, 24, 15
%
%    The pressure nodes are the 9 nodes that are vertices of the triangles,
%    (occurring in any of the first three positions in the element definitions):
%
%      1, 2, 3, 4, 5, 6
%
%    The number of degrees of freedom "NDOF", is:
%      2 * NP (horizontal and vertical velocity at every node)
%      + NP0 (pressure at every interior node)
%      + NE (constraints)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 April 2012
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
%  Local Parameters:
%
%    Local, string E_FILE, the name of the constraint file, which contains 
%    3 rows and NE columns, defining constraints on the data, including Dirichlet 
%    boundary values in particular.
%    Item #1 is a node, #2 is a variable index (0 = horizontal velocity,
%    1 = vertical velocity, 2 = pressure) and #3 is an associated value.
%
%    Local, sparse real K(*), the stiffness matrix, in compressed column format.
%
%    Local, real L(2*NP+NP0+NE), the residual vector.
%
%    Local, integer NE, the number of constraints.
%
%    Local, integer NP, the number of nodes.
%
%    Local, integer NP0, the number of pressure nodes (midside nodes
%    of triangles.)
%
%    Local, integer NT, the number of triangles.
%
%    Local, real NU, the viscosity.
%
%    Local, string P_FILE, the name of the file that contains 2 rows and NP 
%    columns of (X,Y) node coordinates.
%
%    Local, string T_FILE, the name of the element file, which contains 6 rows 
%    and NT columns, with each column containing the (1-based) indices of nodes 
%    forming the triangles, in a particular order.
%
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SMALL_TEST:\n' );
  fprintf ( 1, '  Specify data files for a "small" problem;\n" );
  fprintf ( 1, '  Call NSASM_INTERFACE, which reads that data and\n' );
  fprintf ( 1, '  calls NSASM to get the stiffness matrix K and residual L.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  NSASM is a C library, and must be compiled with MATLAB''s MEX compiler.\n' );
%
%  Call NSASM_INTERFACE.
%
  p_file = 'small_nodes.txt';
  t_file = 'small_elements.txt';
  e_file = 'small_constraints.txt';
  np0 = 9;
  nu = 100.0;

  [ K, L ] = nsasm_interface ( 'small_nodes.txt', t_file, e_file, np0, nu );
%
%  Print some of L.
%
  l_norm = norm ( L, inf );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  L-Infinity norm of L = %f\n', l_norm );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Beginning and end of L vector:\n' );
  fprintf ( 1, '\n' );

  ndof = length ( L );

  for i = 1 : ndof

    if ( i < 10 || ndof - 10 < i )
      fprintf ( 1, '  %8d  %14f\n', i, L(i) );
    end

    if ( i == 10 )
      fprintf ( 1, '..(skipping some entries)...\n' );
    end

  end
%
%  Print some information from K.
%
  nz_num = nnz ( K );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  nz_num = nnz ( K )\n' );
  fprintf ( 1, '  Matrix nonzeros NZ_NUM = %d\n', nz_num );

  [ row, col, val ] = find ( K );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  [ row, col, val ] = find ( K )\n' );
  fprintf ( 1, '  Matrix sparse triplet representation:\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '       ROW       COL     VAL\n' );
  fprintf ( 1, '\n' );

  for i = 1 : nz_num

    if ( i < 10 || nz_num - 10 < i )
      fprintf ( 1, '  %8d  %8d  %14f\n', row(i), col(i), val(i) );
    end

    if ( i == 10 )
      fprintf ( 1, '..(skipping some entries)...\n' );
    end

  end

  r8sp_print_some ( ndof, ndof, nz_num, col, row, val, 1, 1, 10, ...
    10, '  Initial part of K as a matrix' );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SMALL_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
