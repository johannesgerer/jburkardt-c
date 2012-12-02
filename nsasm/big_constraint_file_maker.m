function big_constraint_file_maker ( )

%*****************************************************************************80
%
%% BIG_CONSTRAINT_FILE_MAKER creates a constraint file for the "big" data set.
%
%  Discussion:
%
%    Either I lost, or never had, the constraint file associated with the
%    big data set.  But I did have an output file describing a run of the
%    problem, and a picture of the nodes, from which it was clear that there
%    were 65 nodes on every side (counting both corners), and that the problem
%    was probably an example of the driven cavity.  For the driven cavity,
%    a common constraint is to set the velocity to zero on left, bottom, and
%    right boundaries.  Depending on how you treat the top corners (we constrain
%    them) you will have 65 (left) + 63 (strictly bottom) + 65 (right) constrained
%    boundary nodes, for which X = -1, Y = -1, or X = +1, respectively.
%
%    The constraint file simply lists the node indices, the associated variable,
%    and the associated value.  So we can use the FIND command to find the
%    constrained nodes, and then with a little work, we can set up the appropriate
%    constraint array, including one specified pressure at the end, and write
%    it out to a file.
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
  timestamp ( )
  fprintf ( 1, '\n' );
  fprintf ( 1, 'BIG_CONSTRAINT_FILE_MAKER:\n' );
  fprintf ( 1, '  Create the constraint file associated with the "big" problem.\n' );
%
%  Read the node coordinates.
%
  p = load ( 'big_nodes.txt' );
  p = p';
%
%  Identify nodes on the left, bottom, and right boundaries;
%  such nodes have |X| = 1 or Y = -1.
%
  i = find ( abs ( p(1,:) ) == 1.0 | p(2,:) == -1.0 );
  [ i_row, i_num ] = size ( i );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of constrained velocity nodes is %d\n', i_num );
%
%  Set space for the constraint array.
%
  bob = zeros ( 3, 2*i_num + 1 );
%
%  Set the horizontal and vertical velocities to zero on left, bottom, and right boundaries.
%
  bob(1,1:2:2*i_num-1) = i(1,:);
  bob(2,1:2:2*i_num-1) = 0;
  bob(3,1:2:2*i_num-1) = 0.0;

  bob(1,2:2:2*i_num) = i(1,:);
  bob(2,2:2:2*i_num) = 1;
  bob(3,2:2:2*i_num) = 0.0;
%
%  Set one pressure (say, at node 1) to an arbitrary value (say 1.0).
%
  bob(1,2*i_num+1) = 1;
  bob(2,2*i_num+1) = 2;
  bob(3,2*i_num+1) = 1.0;
%
%  Save the file as "big_constraints.txt".
%
  filename = 'big_constraints.txt';

  save ( filename, 'bob', '-ascii' )

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Constraint data saved in the file "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'BIG_CONSTRAINT_FILE_MAKER:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  fprintf ( 1, '\n' );
  timestamp ( );

  return
end

