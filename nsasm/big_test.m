function big_test ( )

%*****************************************************************************80
%
%% BIG_TEST is a sample program for NSASM using a "big" set of data.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 April 2012
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
  fprintf ( 1, 'BIG_TEST:\n' );
  fprintf ( 1, '  Specify data files for a "small" problem;\n' );
  fprintf ( 1, '  Call NSASM_INTERFACE, which reads that data and\n' );
  fprintf ( 1, '  calls NSASM to get the stiffness matrix K and residual L.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  NSASM is a C library, and must be compiled with MATLAB''s MEX compiler.\n' );
%
%  Call NSASM_INTERFACE.
%
  p_file = 'big_nodes.txt';
  t_file = 'big_elements.txt';
  e_file = 'big_constraints.txt';
  np0 = 545;
  nu = 500.0;

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
  fprintf ( 1, 'BIG_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
