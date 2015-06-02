# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "treepack.h"

int main ( );
void test005 ( );
void test006 ( );
void test01 ( );
void test02 ( );
void test025 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TREEPACK_PRB.

  Discussion:

    TREEPACK_PRB tests the TREEPACK library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TREEPACK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TREEPACK library.\n" );

  test005 ( );
  test006 ( );
  test01 ( );
  test02 ( );
  test025 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TREEPACK_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests CATALAN and CATALAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 July 2013

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  CATALAN computes Catalan numbers.\n" );
  printf ( "  CATALAN_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = catalan ( n );

    printf ( "  %4d  %6d  %6d\n", n, c, c2[n] );
    free ( c2 );

  }

  return;
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests CBT_TRAVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 July 2013

  Author:

    John Burkardt
*/
{
  int depth = 4;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  CBT_TRAVERSE traverses a complete binary tree.\n" );
  printf ( "\n" );
  printf ( "  For this demonstration, we simply print our path.\n" );
  printf ( "  The tree depth is %d\n", depth );
  printf ( "\n" );

  cbt_traverse ( depth );

  return;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PRUEFER_TO_TREE_ARC.

  Discussion:

    The tree is

          5
          |
    2-3-6-8-1-9
      |   |
      7   4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2013

  Author:

    John Burkardt
*/
{
  int code[7] = { 1, 3, 8, 8, 3, 6, 8 };
  int *inode;
  int *jnode;
  int nnode = 9;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  PRUEFER_TO_TREE_ARC computes a tree from its Pruefer code.\n" );
  printf ( "\n" );
  printf ( "          5\n" );
  printf ( "          |\n" );
  printf ( "    2-3-6-8-1-9\n" );
  printf ( "      |   |\n" );
  printf ( "      7   4\n" );

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );

  inode = ( int * ) malloc ( ( nnode - 1 ) * sizeof ( int ) );
  jnode = ( int * ) malloc ( ( nnode - 1 ) * sizeof ( int ) );

  pruefer_to_tree_arc ( nnode, code, inode, jnode );
 
  graph_arc_print ( nnode-1, inode, jnode, "  The graph:" );

  free ( inode );
  free ( jnode );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PRUEFER_TO_TREE_2_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2013

  Author:

    John Burkardt
*/
{
# define NNODE 9

  int code[NNODE] = { 1, 3, 8, 8, 3, 6, 8, 0, 0 };
  int *itree;
  int nnode = NNODE;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PRUEFER_TO_TREE_2_NEW produces a tree from its Pruefer code\n" );

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );

  itree = pruefer_to_tree_2_new ( nnode, code );
 
  i4vec_print ( nnode-1, itree, "  The edge list of the tree:" );
 
  free ( itree );

  return;
# undef NNODE
}
/******************************************************************************/

void test025 ( )

/******************************************************************************/
/*
  Purpose:

    TEST025 tests PRUEFER_TO_TREE_2.

  Discussion:

    This example is used to illustrate the Nijenhuis and Wilf algorithm
    LBLTRE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2013

  Author:

    John Burkardt
*/
{
# define NNODE 4

  int code[NNODE];
  int i;
  int itree[NNODE];
  int j;
  int nnode = NNODE;

  printf ( "\n" );
  printf ( "TEST025\n" );
  printf ( "  PRUEFER_TO_TREE_2 produces a tree from its Pruefer code\n" );
  printf ( "\n" );
  printf ( "   Code      Tree\n" );
  printf ( "\n" );
  for ( j = 1; j <= nnode; j++ )
  {
    code[1] = j;
    for ( i = 1; i <= nnode; i++ )
    {
      code[0] = i;
      pruefer_to_tree_2 ( nnode, code, itree );
      printf ( "  %2d  %2d    %2d  %2d  %2d\n", code[0], code[1], 
        itree[0], itree[1], itree[2] );
    }
  }

  return;
# undef NNODE
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests TREE_ARC_TO_PRUEFER.

  Discussion:

    The tree is

          5
          |
    2-3-6-8-1-9
      |   |
      7   4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2013

  Author:

    John Burkardt
*/
{
  int *code;
  int inode[8] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[8] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int nnode = 9;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  TREE_ARC_TO_PRUEFER: Tree => Pruefer code\n" );
  printf ( "\n" );
  printf ( "          5\n" );
  printf ( "          |\n" );
  printf ( "    2-3-6-8-1-9\n" );
  printf ( "      |   |\n" );
  printf ( "      7   4\n" );

  graph_arc_print ( nnode-1, inode, jnode, "  The graph:" );
 
  code = tree_arc_to_pruefer ( nnode, inode, jnode );

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );
 
  free ( code );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests TREE_ARC_CENTER.

  Discussion:

    The tree is

    2---3---6---8---1---9
       /       / \
      7       5   4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2013

  Author:

    John Burkardt
*/
{
# define NNODE 9

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[NNODE-1] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int nnode = NNODE;
  int parity;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  TREE_ARC_CENTER computes the center of a tree.\n" );

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, &eccent, &parity );

  printf ( "\n" );
  printf ( "  Parity = %d\n", parity );
  printf ( "  Eccentricity is %d\n", eccent );

  if ( parity == 0 )
  {
    printf ( "  No center node (degenerate case).\n" );
  }
  else if ( parity == 1 )
  {
    printf ( "  Center node: %d\n", center[0] );
  }
  else
  {
    printf ( "  Center nodes: %d  %d\n", center[0], center[1] );
  }

  return;
# undef NNODE
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests TREE_ARC_CENTER.

  Discussion:

    Compare:

    2--1--3

    1--2--3

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2013

  Author:

    John Burkardt
*/
{
# define NNODE 3

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1];
  int jnode[NNODE-1];
  int nnode = NNODE;
  int parity;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  TREE_ARC_CENTER computes the center of a tree.\n" );

  inode[0] = 1;
  inode[1] = 1;
  jnode[0] = 2;
  jnode[1] = 3;

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, &eccent, &parity );

  printf ( "\n" );
  printf ( "  Parity = %d\n", parity );
  printf ( "  Eccentricity is %d\n", eccent );

  if ( parity == 0 )
  {
    printf ( "  No center node (degenerate case).\n" );
  }
  else if ( parity == 1 )
  {
    printf ( "  Center node: %d\n", center[0] );
  }
  else
  {
    printf ( "  Center nodes: %d  %d\n", center[0], center[1] );
  }

  inode[0] = 2;
  inode[1] = 2;
  jnode[0] = 1;
  jnode[1] = 3;

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, &eccent, &parity );

  printf ( "\n" );
  printf ( "  Parity = %d\n", parity );
  printf ( "  Eccentricity is %d\n", eccent );

  if ( parity == 0 )
  {
    printf ( "  No center node (degenerate case).\n" );
  }
  else if ( parity == 1 )
  {
    printf ( "  Center node: %d\n", center[0] );
  }
  else
  {
    printf ( "  Center nodes: %d  %d\n", center[0], center[1] );
  }

  return;
# undef NNODE
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests TREE_ARC_CENTER.

  Discussion:

    The tree is

     1-----2-----3
    /|\   / \   /|\
   4 5 6 8  10 7 9 11

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2013

  Author:

    John Burkardt
*/
{
# define NNODE 11

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1] = { 1, 1, 1, 2,  2, 3, 3,  3, 1, 2 };
  int jnode[NNODE-1] = { 4, 5, 6, 8, 10, 7, 9, 11, 2, 3 };
  int nnode = NNODE;
  int parity;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  TREE_ARC_CENTER computes the center of a tree.\n" );

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, &eccent, &parity );

  printf ( "\n" );
  printf ( "  Parity = %d\n", parity );
  printf ( "  Eccentricity is %d\n", eccent );

  if ( parity == 0 )
  {
    printf ( "  No center node (degenerate case).\n" );
  }
  else if ( parity == 1 )
  {
    printf ( "  Center node: %d\n", center[0] );
  }
  else
  {
    printf ( "  Center nodes: %d  %d\n", center[0], center[1] );
  }

  return;
# undef NNODE
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests TREE_ARC_DIAM.

  Discussion:

    The tree is:

    2---3---6---8---1---9
       /       / \
      7       5   4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int diam;
  int inode[8] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[8] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int label[9];
  int nnode = 9;
  int nnode1;
  int nnode2;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  TREE_ARC_DIAM computes the diameter of a tree.\n" );

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_diam ( nnode, inode, jnode, &diam, label, &nnode1, &nnode2 );

  printf ( "\n" );
  printf ( "  This tree has a diameter of %d\n", diam );
  printf ( "  between nodes %d and %d\n", nnode1, nnode2 );

  i4vec_print ( nnode, label, "  Nodes and labels:" );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests TREE_ARC_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2013

  Author:

    John Burkardt
*/
{
  int i;
  int icode[2];
  int inode[3];
  int jnode[3];
  int nnode = 4;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  TREE_ARC_RANDOM produces a random labeled\n" );
  printf ( "  tree and its Pruefer code.\n" );
  printf ( "\n" );
 
  for ( i = 1; i <= 5; i++ )
  {
    tree_arc_random ( nnode, &seed, icode, inode, jnode );

    graph_arc_print ( nnode-1, inode, jnode, "  The random tree:" );

    i4vec_print ( nnode-2, icode, "  The Pruefer code:" );
  }
  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests TREE_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2013

  Author:

    John Burkardt
*/
{
  int nnode;
  int num;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  TREE_ENUM enumerates the labeled trees on a given\n" );
  printf ( "  number of nodes.\n" );
  printf ( "\n" );

  for ( nnode = 0; nnode <= 10; nnode++ )
  {
    num = tree_enum ( nnode );
    printf ( "  %8d  %10d\n", nnode, num );
  } 
  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests TREE_PARENT_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2013

  Author:

    John Burkardt
*/
{
# define NNODE 4

  int icode[NNODE];
  int itree[NNODE];
  int more;
  int nnode = NNODE;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  TREE_PARENT_NEXT finds all labeled trees of a given \n" );
  printf ( "  order, and their Pruefer codes.\n" );
  printf ( "\n" );
  printf ( "  Pruefer code     Tree\n" );
  printf ( "\n" );
 
  more = 0;
 
  for ( ; ; )
  {
    tree_parent_next ( nnode, icode, itree, &more );
 
    printf ( "  %2d  %2d              %2d  %2d  %2d\n",
      icode[0], icode[1], itree[0], itree[1], itree[2] );

    if ( ! more )
    {
      break;
    }
  }
  return;
# undef NNODE
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests TREE_RB_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2013

  Author:

    John Burkardt
*/
{
  int nnode;
  int num;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  TREE_RB_ENUM enumerates the rooted binary trees on a \n" );
  printf ( "  given number of nodes.\n" );
  printf ( "\n" );

  for ( nnode = 0; nnode <= 11; nnode++ )
  {
    num = tree_rb_enum ( nnode );

    printf ( "  %8d  %8d\n", nnode, num );
  }
  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests TREE_RB_LEX_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int a[11];
  int i;
  int j;
  int more;
  int n = 11;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  TREE_RB_LEX_NEXT produces all rooted binary trees with\n" );
  printf ( "  a given number of nodes, in lexicographic order, using\n" );
  printf ( "  the preorder traversal representation.\n" );
  printf ( "\n" );
  printf ( "  The number of nodes N = %d\n", n );
  printf ( "\n" );

  more = 0;
  i = 0;

  for ( ; ; )
  {
    tree_rb_lex_next ( n, a, &more );

    if ( ! more )
    {
      break;
    }

    i = i + 1;
    printf ( "  %2d  ", i );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%d", a[j] );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests TREE_RB_LEX_NEXT, TREE_RB_TO_PARENT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2009

  Author:

    John Burkardt
*/
{
  int a[11];
  int i;
  int j;
  int more;
  int n = 11;
  int *parent;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  TREE_RB_LEX_NEXT produces all rooted binary trees with\n" );
  printf ( "  a given number of nodes, in lexicographic order,\n" );
  printf ( "  using the preorder traversal representation.\n" );
  printf ( "  TREE_RB_TO_PARENT converts the preorder traversal form\n" );
  printf ( "  to the more comprehensible parent node representation.\n" );
  printf ( "\n" );
  printf ( "  The number of nodes N = %d\n", n );
  printf ( "\n" );

  more = 0;
  i = 0;

  for ( ; ; )
  {
    tree_rb_lex_next ( n, a, &more );

    if ( ! more )
    {
      break;
    }

    parent = tree_rb_to_parent ( n, a );

    i = i + 1;
    printf ( "  %2d  ", i );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%3d", parent[j] );
    }
    printf ( "\n" );

    free ( parent );
  }

  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests TREE_RB_YULE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int a[11];
  int i;
  int j;
  int n;
  int n_max = 11;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  TREE_RB_YULE carries out one step of the Yule model\n" );
  printf ( "  on a rooted binary tree stored in preorder traversal form.\n" );
  printf ( "\n" );
  printf ( "  Each call adds two children to an arbitary leaf.\n" );

  for ( i = 1; i <= 5; i++ )
  {
    printf ( "\n" );
    printf ( "  Simulation %d\n", i );
    printf ( "\n" );
    printf ( "  Nodes  Preorder code\n" );
    printf ( "\n" );

    n = 0;

    for ( ; ; )
    {
      tree_rb_yule ( &n, &seed, a );

      printf ( "  %2d  ", n );
      for ( j = 0; j < n; j++ )
      {
        printf ( "%d", a[j] );
      }
      printf ( "\n" );

      if ( n_max < n + 2 )
      {
        break;
      }
    }
  }
  return;
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests TREE_ROOTED_CODE.

  Discussion:

      1
      |\
      | \
      |  \
      2   3
     /|\  |\
    4 5 6 7 8
     /|  \
    9 10  11
      |
      12

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2013

  Author:

    John Burkardt
*/
{
  int *code;
  int nnode = 12;
  int parent[12] = { 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  TREE_ROOTED_CODE: code of a rooted tree.\n" );

  i4vec_print ( nnode, parent, "  Parent vector for tree:" );

  code = tree_rooted_code ( nnode, parent );

  i4vec_print ( 2*nnode, code, "  The tree code:" );

  free ( code );

  return;
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests TREE_ROOTED_DEPTH.

  Discussion:

      1
      |\
      | \
      |  \
      2   3
     /|\  |\
    4 5 6 7 8
     /|  \
    9 10  11
      |
      12

    Depths

    1  2  3  4  5  6  7  8  9 10 11 12
    0  1  1  2  2  2  2  2  3  3  3  4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2013

  Author:

    John Burkardt
*/
{
  int depth;
  int *depth_node;
  int nnode = 12;
  int parent[12] = { 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  TREE_ROOTED_DEPTH: depth of a rooted tree.\n" );

  i4vec_print ( nnode, parent, "  Parent vector for tree:" );

  depth_node = ( int * ) malloc ( nnode * sizeof ( int ) );

  tree_rooted_depth ( nnode, parent, &depth, depth_node );

  i4vec_print ( nnode, depth_node, "  Individual node depths:" );

  printf ( "\n" );
  printf ( "  Overall rooted tree depth: %d\n", depth );

  free ( depth_node );

  return;
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests TREE_ROOTED_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2013

  Author:

    John Burkardt
*/
{
  int nnode = 10;
  int *ntree;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  TREE_ROOTED_ENUM counts unlabeled rooted trees.\n" );

  ntree = tree_rooted_enum ( nnode );

  i4vec_print ( nnode, ntree, 
    "  Number of trees with given number of nodes:" );

  free ( ntree );

  return;
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests TREE_ROOTED_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2013

  Author:

    John Burkardt
*/
{
  int i;
  int *itree;
  int j;
  int nnode = 5;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  TREE_ROOTED_RANDOM: random unlabeled rooted trees.\n" );
  printf ( "\n" );
  printf ( "  Selecting random trees, rooted at 1\n" );
  printf ( "  Number of nodes is NNODE = %d\n", nnode );
  printf ( "\n" );
  printf ( "  Each tree is described by the nodes that\n" );
  printf ( "  connect nodes 2 through NNODE.\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    itree = tree_rooted_random ( nnode, &seed );

    printf ( "  " );
    for ( j = 1; j < nnode; j++ )
    {
      printf ( "%4d", itree[j] );
    }
    printf ( "\n" );

    free ( itree );
  }
  return;
}
