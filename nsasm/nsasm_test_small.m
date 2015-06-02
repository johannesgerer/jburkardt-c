function nsasm_test_small ( )

%*****************************************************************************80
%
%% NSASM_TEST_SMALL is a small sample program for NSASM.
%
%  Discussion:
%
%    The flow region is the unit square [0,1]x[0,1] on which
%    we impose a 5 x 5 grid of nodes indexed as follows:
%
%    1.0
%     ^       7  19   8  23   9
%     |      21  20  22  24  25
%     Y       4  10   5  15   6
%     |      12  11  13  16  17
%     |       1  14   2   18  3
%     |
%     0-------------- X -->  1.0
%
%    The node coordinates are stored in the file 'small_nodes.txt'.
%    The pressure nodes, numbers 1 through 9, form a 3x3 subgrid.
%
%    The nodes are used to form 8 triangular elements:  
%
%     7--19---8--23---9
%     |     / |     / |
%    21  20  22  24  25
%     | /     | /     |
%     4--10---5--15---6
%     |     / |     / |
%    12  11  13  16  17
%     | /     | /     |
%     1--14---2--18---3
%
%    The 8 elements are quadratic six node triangles, with a peculiar
%    ordering for the nodes, and one which is not consistently
%    counterclockwise.
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
%    The node indices forming each element are stored in 'small_elements.txt'.
%
%    All boundary nodes are constrained to have zero horizontal velocity.
%
%    Boundary nodes 3, 17, 6, 25 and 9 are constrained to have a vertical
%    velocity of 1; the other boundary nodes are constrained to have a
%    vertical velocity of 0.
%
%    Node 1 is constrained to have a pressure of 1.
%
%    The node index, the variable being set (0=horizontal
%    velocity, 1=vertical, 2=pressure), 
%    and the value of the variable being constrained are stored as triples in
%    'small_constraints.txt'.
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
%    22 Septenber 2013
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
%    Local, integer NP0, the number of pressure nodes (vertices of triangles.)
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
  fprintf ( 1, 'NSASM_TEST_SMALL:\n' );
  fprintf ( 1, '  MATLAB version.\n' );
  fprintf ( 1, '  Specify data files for a "small" problem;\n' );
  fprintf ( 1, '  Call NSASM_INTERFACE, which reads that data and\n' );
  fprintf ( 1, '  calls NSASM to get the stiffness matrix K and residual L.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  NSASM is a C library,\n' );
  fprintf ( 1, '  It must be compiled with MATLAB''s MEX compiler.\n' );
%
%  Define the input to NSASM_INTERFACE:
%
%    P_FILE contains 2 rows and NP columns of (X,Y) node coordinates.
%    T_FILE contains 6 rows and NT columns of node indices forming triangles. 
%    E_FILE contains 3 rows and NE columns, defining constraints.
%    NP0 is the number of pressure nodes.
%    NU is the fluid viscosity.
%
  p_file = 'small_nodes.txt';
  t_file = 'small_elements.txt';
  e_file = 'small_constraints.txt';
  np0 = 9;
  nu = 100.0;
%
%  Call NSASM_INTERFACE to compute the system matrix and right hand side.
%
  [ K, L ] = nsasm_interface ( p_file, t_file, e_file, np0, nu );
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
%  Get the number of nonzero entries in K.
%
  nz_num = nnz ( K );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  nz_num = nnz ( K )\n' );
  fprintf ( 1, '  Matrix nonzeros NZ_NUM = %d\n', nz_num );
%
%  Get the sparse triplet representation of K as a set of row indices,
%  column indices, and values.  Notice how the data is sorted.
%
  [ row, col, val ] = find ( K );
%
%  Print (some of) the sparse triplet data.
%
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
%
%  Print (some of) the matrix.
%
  r8sp_print_some ( ndof, ndof, nz_num, col, row, val, 1, 1, 10, ...
    10, '  Initial part of K as a matrix' );
%
%  Write matrix to Harwell-Boeing file.
%
  filename = 'small.hb';
  title = 'nsasm_test_small';
  key = 'Key';
  type = 'RUA';
  format = 8;
  job = 3;

  addpath ( '../../m_src/msm_to_hb' );

  msm_to_hb ( filename, K, L, title, key, type, format, job );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Matrix data saved to Harwell Boeing file "%s"\n', filename );

  rmpath ( '../../m_src/msm_to_hb' );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'NSASM_TEST_SMALL:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
