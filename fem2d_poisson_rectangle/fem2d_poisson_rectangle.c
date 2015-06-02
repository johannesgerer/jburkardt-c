# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( void );
void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] );
void assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int nq,
  double wq[], double xq[], double yq[], double element_area[], int indx[],
  int ib, int nunk, double a[], double f[] );
int bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num, int indx[] );
void boundary ( int nx, int ny, int node_num, double node_xy[], int indx[],
  int ib, int nunk, double a[], double f[] );
void compare ( int node_num, double node_xy[], int indx[], int nunk, double f[] );
void errors ( double element_area[], int element_node[], int indx[],
  double node_xy[], double f[], int element_num, int nnodes,
  int nunk, int node_num, double *el2, double *eh1 );
void exact ( double x, double y, double *u, double *dudx, double *dudy );
void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4mat_write ( char *output_filename, int m, int n, int table[] );
void indx_set ( int nx, int ny, int node_num, int indx[], int *nunk );
void nodes_plot ( char *file_name, int node_num, double node_xy[],
  int node_label );
void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *bb, double *bx, double *by );
void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] );
void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] );
double r8_abs ( double x );
double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
int r8gb_fa ( int n, int ml, int mu, double a[], int pivot[] );
void r8gb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, char *title );
double *r8gb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title );
double rhs ( double x, double y );
void solution_write ( double f[], int indx[], int node_num, int nunk,
  char *output_filename, double node_xy[] );
void timestamp ( void );
void triangulation_order6_plot ( char *file_name, int node_num, double node_xy[],
  int tri_num, int triangle_node[], int node_show, int triangle_show );
void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM2D_POISSON_RECTANGLE.

  Discussion:

    FEM2D_POISSON_RECTANGLE solves

      -Laplacian U(X,Y) = F(X,Y)

    in a rectangular region in the plane.  Along the boundary,
    Dirichlet boundary conditions are imposed.

      U(X,Y) = G(X,Y)

    The code uses continuous piecewise quadratic basis functions on
    triangles determined by a uniform grid of NX by NY points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Local parameters:

    Local, double A[(3*IB+1)*NUNK], the coefficient matrix.

    Local, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.

    Local, double C[NUNK], the finite element coefficients, solution of A * C = F.

    Local, double EH1, the H1 seminorm error.

    Local, double EL2, the L2 error.

    Local, int ELEMENT_NODE[ELEMENT_NUM*NNODES]; ELEMENT_NODE(I,J) is the
    global node index of the local node J in element I.

    Local, int ELEMENT_NUM, the number of elements.

    Local, double F[NUNK], the right hand side.

    Local, int IB, the half-bandwidth of the matrix.

    Local, int INDX[NODE_NUM], gives the index of the unknown quantity
    associated with the given node.

    Local, int NNODES, the number of nodes used to form one element.

    Local, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.

    Local, int NQ, the number of quadrature points used for assembly.

    Local, int NUNK, the number of unknowns.

    Local, int NX, the number of points in the X direction.

    Local, int NY, the number of points in the Y direction.

    Local, double WQ[NQ], quadrature weights.

    Local, double XL, XR, YB, YT, the X coordinates of
    the left and right sides of the rectangle, and the Y coordinates
    of the bottom and top of the rectangle.

    Local, double XQ[NQ*ELEMENT_NUM], YQ[NQ*ELEMENT_NUM], the X and Y
    coordinates of the quadrature points in each element.
*/
{
# define NNODES 6
# define NQ 3
# define NX 7
# define NY 7
# define ELEMENT_NUM ( NX - 1 ) * ( NY - 1 ) * 2
# define NODE_NUM ( 2 * NX - 1 ) * ( 2 * NY - 1 )

  double *a;
  double *c;
  double eh1;
  double el2;
  int element;
  double element_area[ELEMENT_NUM];
  int *element_mask;
  int element_node[NNODES*ELEMENT_NUM];
  double *f;
  int i;
  int ib;
  int ierr;
  int indx[NODE_NUM];
  int job;
  int local;
  int node;
  char *node_eps_file_name = "fem2d_poisson_rectangle_nodes.eps";
  char *node_txt_file_name = "fem2d_poisson_rectangle_nodes.txt";
  int node_label;
  int node_show;
  double node_xy[2*NODE_NUM];
  int nunk;
  int *pivot;
  char *solution_txt_file_name = "fem2d_poisson_rectangle_solution.txt";
  int triangle_show;
  char *triangulation_eps_file_name = "fem2d_poisson_rectangle_elements.eps";
  char *triangulation_txt_file_name = "fem2d_poisson_rectangle_elements.txt";
  double wq[NQ];
  double xl = 0.0;
  double xq[NQ*ELEMENT_NUM];
  double xr = 1.0;
  double yb = 0.0;
  double yq[NQ*ELEMENT_NUM];
  double yt = 1.0;

  timestamp ( );

  printf ( "\n" );
  printf ( "FEM2D_POISSON_RECTANGLE:\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Solution of the Poisson equation on a unit box\n" );
  printf ( "  in 2 dimensions.\n" );
  printf ( "\n" );
  printf ( "  - Uxx - Uyy = F(x,y) in the box\n" );
  printf ( "       U(x,y) = G(x,y) on the boundary.\n" );
  printf ( "\n" );
  printf ( "  The finite element method is used, with piecewise\n" );
  printf ( "  quadratic basis functions on 6 node triangular\n" );
  printf ( "  elements.\n" );
  printf ( "\n" );
  printf ( "  The corner nodes of the triangles are generated by an\n" );
  printf ( "  underlying grid whose dimensions are\n" );
  printf ( "\n" );
  printf ( "  NX =                 %d\n", NX );
  printf ( "  NY =                 %d\n", NY );
  printf ( "\n" );
  printf ( "  Number of nodes    = %d\n", NODE_NUM );
  printf ( "  Number of elements = %d\n", ELEMENT_NUM );
/*
  Set the coordinates of the nodes.
*/
  xy_set ( NX, NY, NODE_NUM, xl, xr, yb, yt, node_xy );
/*
  Organize the nodes into a grid of 6-node triangles.
*/
  grid_t6 ( NX, NY, NNODES, ELEMENT_NUM, element_node );
/*
  Set the quadrature rule for assembly.
*/
  quad_a ( node_xy, element_node, ELEMENT_NUM, NODE_NUM,
    NNODES, wq, xq, yq );
/*
  Determine the areas of the elements.
*/
  area_set ( NODE_NUM, node_xy, NNODES, ELEMENT_NUM,
    element_node, element_area );
/*
  Determine which nodes are boundary nodes and which have a
  finite element unknown.  Then set the boundary values.
*/
  indx_set ( NX, NY, NODE_NUM, indx, &nunk );

  printf ( "  Number of unknowns =       %d\n", nunk );
/*
  Determine the bandwidth of the coefficient matrix.
*/
  ib = bandwidth ( NNODES, ELEMENT_NUM, element_node, NODE_NUM, indx );

  printf ( "\n" );
  printf ( "  Total bandwidth is %d\n", 3 * ib + 1 );
/*
  Make an EPS picture of the nodes.
*/
  if ( NX <= 10 && NY <= 10 )
  {
    node_label = 1;
    nodes_plot ( node_eps_file_name, NODE_NUM, node_xy, node_label );

    printf ( "\n" );
    printf ( "FEM2D_POISSON_RECTANGLE:\n" );
    printf ( "  Wrote an EPS file\n" );
    printf ( "    \"%s\".\n", node_eps_file_name );
    printf ( "  containing a picture of the nodes.\n" );
  }
/*
  Write the nodes to an ASCII file that can be read into MATLAB.
*/
  r8mat_write ( node_txt_file_name, 2, NODE_NUM, node_xy );

  printf ( "\n" );
  printf ( "FEM2D_POISSON_RECTANGLE:\n" );
  printf ( "  Wrote an ASCII node file \"%s\",\n", node_txt_file_name );
  printf ( "  of the form\n" );
  printf ( "    X(I), Y(I)\n" );
  printf ( "  which can be used for plotting.\n" );
/*
  Make a picture of the elements.
*/
  if ( NX <= 10 && NY <= 10 )
  {
    node_show = 1;
    triangle_show = 2;

    triangulation_order6_plot ( triangulation_eps_file_name, NODE_NUM,
      node_xy, ELEMENT_NUM, element_node, node_show, triangle_show );

    printf ( "\n" );
    printf ( "FEM2D_POISSON_RECTANGLE:\n" );
    printf ( "  Wrote an EPS file \"%s\",\n", triangulation_eps_file_name );
    printf ( "  containing a picture of the elements.\n" );
  }
/*
  Write the elements to a file that can be read into MATLAB.
*/
  i4mat_write ( triangulation_txt_file_name, NNODES, ELEMENT_NUM, element_node );

  printf ( "\n" );
  printf ( "FEM2D_POISSON_RECTANGLE:\n" );
  printf ( "  Wrote an ASCII element file \"%s\"\n", triangulation_txt_file_name );
  printf ( "  of the form\n" );
  printf ( "    Node(1) Node(2) Node(3) Node(4) Node(5) Node(6)\n" );
  printf ( "  which can be used for plotting.\n" );
/*
  Allocate space for the coefficient matrix A and right hand side F.
*/
  a = ( double * ) malloc ( ( 3 * ib + 1 ) * nunk * sizeof ( double ) );
  f = ( double * ) malloc ( nunk * sizeof ( double ) );
  pivot = ( int * ) malloc ( nunk * sizeof ( int ) );
/*
  Assemble the coefficient matrix A and the right-hand side F of the
  finite element equations.
*/
  assemble ( NODE_NUM, node_xy, NNODES,
    ELEMENT_NUM, element_node, NQ,
    wq, xq, yq, element_area, indx, ib, nunk, a, f );
/*
  Print a tiny portion of the matrix.
*/
  r8gb_print_some ( nunk, nunk, ib, ib, a, 1, 1, 5, 5,
    "  Initial 5 x 5 block of coefficient matrix A:" );

  r8vec_print_some ( nunk, f, 1, 10, "  Part of the right hand side F:" );
/*
  Modify the coefficient matrix and right hand side to account for
  boundary conditions.
*/
  boundary ( NX, NY, NODE_NUM, node_xy, indx, ib, nunk, a, f );
/*
  Print a tiny portion of the matrix.
*/
  r8gb_print_some ( nunk, nunk, ib, ib, a, 1, 1, 5, 5,
    "  A after boundary adjustment:" );

  r8vec_print_some ( nunk, f, 1, 10, "  F after boundary adjustment:" );
/*
  Solve the linear system using a banded solver.
*/
  ierr = r8gb_fa ( nunk, ib, ib, a, pivot );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "FEM2D_POISSON_RECTANGLE - Error!\n" );
    printf ( "  R8GB_FA returned an error condition.\n" );
    printf ( "\n" );
    printf ( "  The linear system was not factored, and the\n" );
    printf ( "  algorithm cannot proceed.\n" );
    exit ( 1 );
  }

  job = 0;
  c = r8gb_sl ( nunk, ib, ib, a, pivot, f, job );

  r8vec_print_some ( nunk, c, 1, 10, "  Part of the solution vector:" );
/*
  Calculate error using 13 point quadrature rule.
*/
  errors ( element_area, element_node, indx, node_xy, c,
    ELEMENT_NUM, NNODES, nunk, NODE_NUM, &el2, &eh1 );
/*
  Compare the exact and computed solutions just at the nodes.
*/
  compare ( NODE_NUM, node_xy, indx, nunk, c );
/*
  Write an ASCII file that can be read into MATLAB.
*/
  solution_write ( c, indx, NODE_NUM, nunk, solution_txt_file_name,
    node_xy );

  printf ( "\n" );
  printf ( "FEM2D_POISSON_RECTANGLE:\n" );
  printf ( "  Wrote an ASCII solution file \"%s\",\n", solution_txt_file_name );
  printf ( "  of the form\n" );
  printf ( "    U( X(I), Y(I) )\n" );
  printf ( "  which can be used for plotting.\n" );
/*
  Free memory.
*/
  free ( a );
  free ( c );
  free ( f );
  free ( pivot );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM2D_POISSON_RECTANGLE:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
# undef NNODES
# undef NQ
# undef NX
# undef NY
# undef ELEMENT_NUM
# undef NODE_NUM
}
/******************************************************************************/

void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] )

/******************************************************************************/
/*
  Purpose:

    AREA_SET sets the area of each element.

  Discussion:

    The areas of the elements are needed in order to adjust
    the integral estimates produced by the quadrature formulas.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the
    coordinates of the nodes.

    Input, int NNODES, the number of local nodes per element.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
*/
{
  int element;
  int i1;
  int i2;
  int i3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  for ( element = 0; element < element_num; element++ )
  {
    i1 = element_node[0+element*nnodes];
    x1 = node_xy[0+(i1-1)*2];
    y1 = node_xy[1+(i1-1)*2];

    i2 = element_node[1+element*nnodes];
    x2 = node_xy[0+(i2-1)*2];
    y2 = node_xy[1+(i2-1)*2];

    i3 = element_node[2+element*nnodes];
    x3 = node_xy[0+(i3-1)*2];
    y3 = node_xy[1+(i3-1)*2];

    element_area[element] = 0.5 * r8_abs
      ( y1 * ( x2 - x3 )
      + y2 * ( x3 - x1 )
      + y3 * ( x1 - x2 ) );
  }

  return;
}
/******************************************************************************/

void assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int nq,
  double wq[], double xq[], double yq[], double element_area[], int indx[],
  int ib, int nunk, double a[], double f[] )

/******************************************************************************/
/*
  Purpose:

    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.

  Discussion:

    The matrix is known to be banded.  A special matrix storage format
    is used to reduce the space required.  Details of this format are
    discussed in the routine R8GB_FA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Output, double A[(3*IB+1)*NUNK], the NUNK by NUNK coefficient matrix,
    stored in a compressed format.

    Output, double F[NUNK], the right hand side.

    Input, int IB, the half-bandwidth of the matrix.

    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.

    Input, double XQ[NQ*ELEMENT_NUM], YQ[NQ*ELEMENT_NUM], the X and Y
    coordinates of the quadrature points in each element.

    Input, double WQ[NQ], quadrature weights.

    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    index of local node I in element J.

    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
    associated with the given node.

    Input, int NNODES, the number of nodes used to form one element.

    Input, int NUNK, the number of unknowns.

    Input, int NQ, the number of quadrature points used in assembly.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int NODE_NUM, the number of nodes.

  Local parameters:

    Local, double BB, BX, BY, the value of some basis function
    and its first derivatives at a quadrature point.

    Local, double BJ, DBJDX, DBJDY, the value of another basis
    function and its first derivatives at a quadrature point.

    Local, int NODE_NUM, the number of global nodes.
*/
{
  double aij;
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int ii;
  int ij;
  int ip;
  int ipp;
  int j;
  int quad;
  int test;
  double w;
  double x;
  double y;
/*
  Initialize the arrays to zero.
*/
  for ( i = 1; i <= nunk; i++ )
  {
    f[i-1] = 0.0;
  }

  for ( j = 1; j <= nunk; j++ )
  {
    for ( i = 1; i <= 3*ib + 1; i++ )
    {
      a[i-1+(j-1)*(3*ib+1)] = 0.0;
    }
  }
/*
  The actual values of A and F are determined by summing up
  contributions from all the elements.
*/
  for ( element = 1; element <= element_num; element++ )
  {
    for ( quad = 1; quad <= nq; quad++ )
    {
      x = xq[quad-1+(element-1)*nq];
      y = yq[quad-1+(element-1)*nq];
      w = element_area[element-1] * wq[quad-1];

      for ( test = 1; test <= nnodes; test++ )
      {
        ip = element_node[test-1+(element-1)*nnodes];
        i = indx[ip-1];

        qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        f[i-1] = f[i-1] + w * rhs ( x, y ) * bi;
/*
  We are about to compute a contribution associated with the
  I-th test function and the J-th basis function, and add this
  to the entry A(I,J).

  Because of the compressed storage of the matrix, the element
  will actually be stored in A(I-J+2*IB+1,J).

  Two extra complications: we are storing the array as a vector,
  and C uses 0-based indices rather than 1-based indices.

  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + (J-1) * (3*IB+1)];
*/
        for ( basis = 1; basis <= nnodes; basis++ )
        {
          ipp = element_node[basis-1+(element-1)*nnodes];
          j = indx[ipp-1];

          qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          aij = dbidx * dbjdx + dbidy * dbjdy;

          a[i-j+2*ib+(j-1)*(3*ib+1)] = a[i-j+2*ib+(j-1)*(3*ib+1)] + w * aij;
        }
      }
    }
  }

  return;
}
/******************************************************************************/

int bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num, int indx[] )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH determines the bandwidth of the coefficient matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NNODES, the number of local nodes per element.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    index of local node I in element J.

    Input, int NODE_NUM, the number of nodes.

    Input, int INDX[NODE_NUM], indicates how the value associated with the
    node is to be determined.  If INDX(I) is positive, then this is the
    index of the unknown in the finite element linear system.  The value
    at the node will be determined by solving the finite element system.
    If INDX(I) is negative, then the node is associated with a boundary
    condition; the value of the boundary condition is stored in the array
    UB, in location -INDX(I).

    Output, int BANDWIDTH, the half bandwidth of the matrix.
*/
{
  int element;
  int i;
  int iln;
  int in;
  int j;
  int jln;
  int jn;
  int nhba;

  nhba = 0;

  for ( element = 1; element <= element_num; element++ )
  {
    for ( iln = 1; iln <= nnodes; iln++ )
    {
      in = element_node[iln-1+(element-1)*nnodes];
      i = indx[in-1];
      if ( 0 < i )
      {
        for ( jln = 1; jln <= nnodes; jln++ )
        {
          jn = element_node[jln-1+(element-1)*nnodes];
          j = indx[jn-1];
          nhba = i4_max ( nhba, j - i );
        }
      }
    }
  }

  return nhba;
}
/******************************************************************************/

void boundary ( int nx, int ny, int node_num, double node_xy[], int indx[],
  int ib, int nunk, double a[], double f[] )

/******************************************************************************/
/*
  Purpose:

    BOUNDARY modifies the linear system for boundary conditions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, controls the number of elements along the
    X and Y directions.  The number of elements will be
    2 * ( NX - 1 ) * ( NY - 1 ).

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.

    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
    associated with the given node.

    Input, int IB, the half-bandwidth of the matrix.

    Input, int NUNK, the number of unknowns.

    Input/output, double A[(3*IB+1)*NUNK], the NUNK by NUNK
    coefficient matrix, stored in a compressed format.
    On output, A has been adjusted for boundary conditions.

    Input/output, double F[NUNK], the right hand side.
    On output, F has been adjusted for boundary conditions.
*/
{
  int col;
  double dudx;
  double dudy;
  int i;
  int j;
  int jhi;
  int jlo;
  int node;
  int row;
  double u;
  double x;
  double y;
/*
  Consider each node.
*/
  node = 0;

  for ( row = 1; row <= 2 * ny - 1; row++ )
  {
    for ( col = 1; col <= 2 * nx - 1; col++ )
    {
      node = node + 1;

      if ( row == 1 ||
           row == 2 * ny - 1 ||
           col == 1 ||
           col == 2 * nx - 1 )
      {
        i = indx[node-1];
        x = node_xy[0+(node-1)*2];
        y = node_xy[1+(node-1)*2];
        exact ( x, y, &u, &dudx, &dudy );

        jlo = i4_max ( i - ib, 1 );
        jhi = i4_min ( i + ib, nunk );

        for ( j = jlo; j <= jhi; j++ )
        {
          a[i-j+2*ib+(j-1)*(3*ib+1)] = 0.0;
        }

        a[i-i+2*ib+(i-1)*(3*ib+1)] = 1.0;

        f[i-1] = u;
      }
    }
  }

  return;
}
/******************************************************************************/

void compare ( int node_num, double node_xy[], int indx[], int nunk,
  double f[] )

/******************************************************************************/
/*
  Purpose:

    COMPARE compares the exact and computed solution at the nodes.

  Discussion:

    This is a rough comparison, done only at the nodes.  Such a pointwise
    comparison is easy, because the value of the finite element
    solution is exactly the value of the finite element coefficient
    associated with that node.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Input, int INDX[NODE_NUM], the index of the unknown in the finite
    element linear system.

    Input, int NUNK, the number of unknowns in the finite element system.

    Input, double F[NUNK], the solution vector of the finite
    element system.
*/
{
  double dudx;
  double dudy;
  int i;
  int node;
  double u;
  double uh;
  double x;
  double y;

  printf ( "\n" );
  printf ( "COMPARE:\n" );
  printf ( "  Compare computed and exact solutions at the nodes.\n" );
  printf ( "\n" );
  printf ( "         X           Y          U           U\n" );
  printf ( "                             computed     exact\n" );
  printf ( "\n" );

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    exact ( x, y, &u, &dudx, &dudy );

    i = indx[node];
    uh = f[i-1];

    printf ( "  %12g  %12g  %12g  %12g\n", x, y, uh, u );
  }

  return;
}
/******************************************************************************/

void errors ( double element_area[], int element_node[], int indx[],
  double node_xy[], double f[], int element_num, int nnodes,
  int nunk, int node_num, double *el2, double *eh1 )

/******************************************************************************/
/*
  Purpose:

    ERRORS calculates the error in the L2 and H1-seminorm.

  Discussion:

    This routine uses a 13 point quadrature rule in each element,
    in order to estimate the values of

      EL2 = Sqrt ( Integral ( U(x,y) - Uh(x,y) )**2 dx dy )

      EH1 = Sqrt ( Integral ( Ux(x,y) - Uhx(x,y) )**2 +
                            ( Uy(x,y) - Uhy(x,y) )**2 dx dy )

    Here U is the exact solution, and Ux and Uy its spatial derivatives,
    as evaluated by a user-supplied routine.

    Uh, Uhx and Uhy are the computed solution and its spatial derivatives,
    as specified by the computed finite element solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    index of local node I in element J.

    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
    associated with the given node.

    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.

    Input, double F[NUNK], the coefficients of the solution.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int NNODES, the number of nodes used to form one element.

    Input, int NUNK, the number of unknowns.

    Input, int NODE_NUM, the number of nodes.

    Output, double precision *EL2, the L2 error.

    Output, double precision *EH1, the H1 seminorm error.

  Local Parameters:

    Local, double AR, the weight for a given quadrature point
    in a given element.

    Local, double BI, DBIDX, DBIDY, a basis function and its first
    derivatives evaluated at a particular quadrature point.

    Local, double EH1, the H1 seminorm error.

    Local, double EL2, the L2 error.

    Local, int NQE, the number of points in the quadrature rule.
    This is actually fixed at 13.

    Local, double UEX, UEXX, UEXY, the exact solution and its first
    derivatives evaluated at a particular quadrature point.

    Local, double UH, UHX, UHY, the computed solution and its first
    derivatives evaluated at a particular quadrature point.

    Local, double WQE(NQE), stores the quadrature weights.

    Local, double X, Y, the coordinates of a particular
    quadrature point.

    Local, double XQE(NQE), YQE(NQE), stores the location
    of quadrature points in a given element.
*/
{
# define NQE 13

  double ar;
  double bi;
  double dbidx;
  double dbidy;
  double dudx;
  double dudxh;
  double dudy;
  double dudyh;
  int element;
  int i;
  int in1;
  int ip;
  int quad;
  double u;
  double uh;
  double wqe[NQE];
  double x;
  double x1;
  double xqe[NQE];
  double y;
  double y1;
  double yqe[NQE];

  *el2 = 0.0;
  *eh1 = 0.0;
/*
  For each element, retrieve the nodes, area, quadrature weights,
  and quadrature points.
*/
  for ( element = 1; element <= element_num; element++ )
  {
    quad_e ( node_xy, element_node, element, element_num,
      nnodes, node_num, NQE, wqe, xqe, yqe );
/*
  For each quadrature point, evaluate the computed solution and its X and
  Y derivatives.
*/
    for ( quad = 1; quad <= NQE; quad++ )
    {
      ar = element_area[element-1] * wqe[quad-1];
      x = xqe[quad-1];
      y = yqe[quad-1];

      uh = 0.0;
      dudxh = 0.0;
      dudyh = 0.0;

      for ( in1 = 1; in1 <= nnodes; in1++ )
      {
        ip = element_node[in1-1+(element-1)*nnodes];

        qbf (x, y, element, in1, node_xy,
          element_node, element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        x1 = node_xy[0+(ip-1)*2];
        y1 = node_xy[1+(ip-1)*2];
        i = indx[ip-1];

        uh    = uh    + bi    * f[i-1];
        dudxh = dudxh + dbidx * f[i-1];
        dudyh = dudyh + dbidy * f[i-1];
      }
/*
  Evaluate the exact solution and its X and Y derivatives.
*/
      exact ( x, y, &u, &dudx, &dudy );
/*
  Add the weighted value at this quadrature point to the quadrature sum.
*/
      *el2 = *el2 + ar *   pow ( ( uh  - u  ), 2 );

      *eh1 = *eh1 + ar * ( pow ( ( dudxh - dudx ), 2 )
                         + pow ( ( dudyh - dudy ), 2 ) );
    }
  }

  *el2 = sqrt ( *el2 );
  *eh1 = sqrt ( *eh1 );

  printf ( "\n" );
  printf ( "*********************************************\n" );
  printf ( "*                                           *\n" );
  printf ( "*  ERRORS:                                  *\n" );
  printf ( "*    L2 error =          %14g     *\n", *el2 );
  printf ( "*    H1-seminorm error = %14g     *\n", *eh1 );
  printf ( "*                                           *\n" );
  printf ( "*********************************************\n" );

  return;
# undef NQE
}
/******************************************************************************/

void exact ( double x, double y, double *u, double *dudx, double *dudy )

/******************************************************************************/
/*
  Purpose:

    EXACT calculates the exact solution and its first derivatives.

  Discussion:

    The function specified here depends on the problem being
    solved.  This is one of the routines that a user will
    normally want to change.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the coordinates of a point
    in the region, at which the right hand side of the
    differential equation is to be evaluated.

    Output, double *U, *DUDX, *DUDY, the value of
    the exact solution U and its derivatives dUdX
    and dUdY at the point (X,Y).
*/
{
# define PI 3.14159265358979323846264338327950288419716939937510

  *u    =      sin ( PI * x ) * sin ( PI * y ) + x;
  *dudx = PI * cos ( PI * x ) * sin ( PI * y ) + 1.0;
  *dudy = PI * sin ( PI * x ) * cos ( PI * y );

  return;
# undef PI
}
/******************************************************************************/

void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GRID_T6 produces a grid of pairs of 6 node triangles.

  Example:

    Input:

      NX = 4, NY = 3

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

  Diagram:

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

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, controls the number of elements along the
    X and Y directions.  The number of elements will be
    2 * ( NX - 1 ) * ( NY - 1 ).

    Input, int NNODES, the number of local nodes per element.

    Input, int ELEMENT_NUM, the number of elements.

    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
*/
{
  int c;
  int e;
  int element;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element = 0;

  for ( j = 1; j <= ny - 1; j++ )
  {
    for ( i = 1; i <= nx - 1; i++ )
    {
      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s  + 1;
      n  = s  + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element = element + 1;
      element_node[0+(element-1)*nnodes] = sw;
      element_node[1+(element-1)*nnodes] = se;
      element_node[2+(element-1)*nnodes] = nw;
      element_node[3+(element-1)*nnodes] = s;
      element_node[4+(element-1)*nnodes] = c;
      element_node[5+(element-1)*nnodes] = w;

      element = element + 1;
      element_node[0+(element-1)*nnodes] = ne;
      element_node[1+(element-1)*nnodes] = nw;
      element_node[2+(element-1)*nnodes] = se;
      element_node[3+(element-1)*nnodes] = n;
      element_node[4+(element-1)*nnodes] = c;
      element_node[5+(element-1)*nnodes] = e;
    }
  }

  return;
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

void indx_set ( int nx, int ny, int node_num, int indx[], int *nunk )

/******************************************************************************/
/*
  Purpose:

    INDX_SET assigns a boundary value index or unknown value index at each node.

  Discussion:

    Every finite element node will is assigned an index which
    indicates the finite element basis function and its coefficient
    which are associated with that node.

  Example:

    On a simple 5 by 5 grid, where the nodes are numbered starting
    at the lower left, and increasing in X first, we would have the
    following values of INDX:

       21  22  23  24  25
       16  17  18  19  20
       11  12  13  14  15
        6   7   8   9  10
        1   2   3   4   5

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of elements in the X and Y directions.

    Input, int NODE_NUM, the number of nodes.

    Output, int INDX[NODE_NUM], the index of the unknown in the finite
    element linear system.

    Output, int *NUNK, the number of unknowns.
*/
{
  int i;
  int in;
  int j;

  *nunk = 0;
  in = 0;

  for ( j = 1; j <= 2 * ny - 1; j++ )
  {
    for ( i = 1; i <= 2 * nx - 1; i++ )
    {
      in = in + 1;
      *nunk = *nunk + 1;
      indx[in-1] = *nunk;
    }
  }

  return;
}
/******************************************************************************/

void nodes_plot ( char *file_name, int node_num, double node_xy[],
  int node_label )

/******************************************************************************/
/*
  Purpose:

    NODES_PLOT plots a pointset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to create.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Input, int NODE_LABEL, is TRUE if the nodes are to be labeled.
*/
{
  int circle_size;
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
    printf ( "\n" );
    printf ( "POINTS_PLOT - Fatal error!\n" );
    printf ( "  Could not open the output EPS file.\n" );
    exit ( 1 );
  }

  fprintf ( file_unit, "%!PS-Adobe-3.0 EPSF-3.0\n" );
  fprintf ( file_unit, "%%Creator: nodes_plot.C\n" );
  fprintf ( file_unit, "%%Title: %s\n", file_name );

  fprintf ( file_unit, "%%Pages: 1\n" );
  fprintf ( file_unit, "%%BoundingBox:  %d  %d  %d  %d\n", x_ps_min, y_ps_min, x_ps_max, y_ps_max );
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
  fprintf ( file_unit, " %d  %d moveto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_max, y_ps_min );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_max, y_ps_max );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_min, y_ps_max );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_min, y_ps_min );
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
  fprintf ( file_unit, " %d  %d moveto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_max_clip, y_ps_min_clip );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_max_clip, y_ps_max_clip );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_min_clip, y_ps_max_clip );
  fprintf ( file_unit, " %d  %d lineto\n", x_ps_min_clip, y_ps_min_clip );
  fprintf ( file_unit, "clip newpath\n" );
/*
  Draw the nodes.
*/
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Draw filled dots at each node:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%  Set the color to blue:\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "0.000  0.150  0.750  setrgbcolor\n" );
  fprintf ( file_unit, "%\n" );

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

    fprintf ( file_unit, "newpath  %d  %d  %d  0 360 arc closepath fill\n", x_ps, y_ps, circle_size );
  }
/*
  Label the nodes.
*/
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
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
      / ( y_max                     - y_min ) );

    fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n", x_ps, y_ps + 5, node + 1 );
  }

  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "restore showpage\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "% End of page\n" );
  fprintf ( file_unit, "%\n" );
  fprintf ( file_unit, "%%Trailer\n" );
  fprintf ( file_unit, "%%EOF\n" );

  fclose ( file_unit );

  return;
}
/******************************************************************************/

void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *b, double *dbdx, double *dbdy )

/******************************************************************************/
/*
  Purpose:

    QBF evaluates the quadratic basis functions.

  Discussion:

    This routine assumes that the "midpoint" nodes are, in fact,
    exactly the average of the two extreme nodes.  This is NOT true
    for a general quadratic triangular element.

    Assuming this property of the midpoint nodes makes it easy to
    determine the values of (R,S) in the reference element that
    correspond to (X,Y) in the physical element.

    Once we know the (R,S) coordinates, it's easy to evaluate the
    basis functions and derivatives.

  The physical element T6:

    In this picture, we don't mean to suggest that the bottom of
    the physical triangle is horizontal.  However, we do assume that
    each of the sides is a straight line, and that the intermediate
    points are exactly halfway on each side.

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

  Reference element T6:

    In this picture of the reference element, we really do assume
    that one side is vertical, one horizontal, of length 1.

    |
    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-------->

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the (global) coordinates of the point
    at which the basis function is to be evaluated.

    Input, int ELEMENT, the index of the element which contains the point.

    Input, int INODE, the local index (between 1 and 6) that
    specifies which basis function is to be evaluated.

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int NNODES, the number of nodes used to form one element.

    Input, int NODE_NUM, the number of nodes.

    Output, double *B, *DBDX, *DBDY, the value of the basis function
    and its X and Y derivatives at (X,Y).
*/
{
  double dbdr;
  double dbds;
  double det;
  double drdx;
  double drdy;
  double dsdx;
  double dsdy;
  int i;
  double r;
  double s;
  double xn[6];
  double yn[6];

  for ( i = 0; i < 6; i++ )
  {
    xn[i] = node_xy[0+(element_node[i+(element-1)*nnodes]-1)*2];
    yn[i] = node_xy[1+(element_node[i+(element-1)*nnodes]-1)*2];
  }
/*
  Determine the (R,S) coordinates corresponding to (X,Y).

  What is happening here is that we are solving the linear system:

    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )

  by computing the inverse of the coefficient matrix and multiplying
  it by the right hand side to get R and S.

  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
  for R and S.
*/
  det =   ( xn[1] - xn[0] ) * ( yn[2] - yn[0] )
        - ( xn[2] - xn[0] ) * ( yn[1] - yn[0] );

  r = ( ( yn[2] - yn[0] ) * ( x     - xn[0] )
      + ( xn[0] - xn[2] ) * ( y     - yn[0] ) ) / det;

  drdx = ( yn[2] - yn[0] ) / det;
  drdy = ( xn[0] - xn[2] ) / det;

  s = ( ( yn[0] - yn[1] ) * ( x     - xn[0] )
      + ( xn[1] - xn[0] ) * ( y     - yn[0] ) ) / det;

  dsdx = ( yn[0] - yn[1] ) / det;
  dsdy = ( xn[1] - xn[0] ) / det;
/*
  The basis functions can now be evaluated in terms of the
  reference coordinates R and S.  It's also easy to determine
  the values of the derivatives with respect to R and S.
*/
  if ( inode == 1 )
  {
    *b   =   2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
    dbdr = - 3.0 + 4.0 * r + 4.0 * s;
    dbds = - 3.0 + 4.0 * r + 4.0 * s;
  }
  else if ( inode == 2 )
  {
    *b   =   2.0 * r * ( r - 0.5 );
    dbdr = - 1.0 + 4.0 * r;
    dbds =   0.0;
  }
  else if ( inode == 3 )
  {
    *b   =   2.0 * s * ( s - 0.5 );
    dbdr =   0.0;
    dbds = - 1.0               + 4.0 * s;
  }
  else if ( inode == 4 )
  {
    *b   =   4.0 * r * ( 1.0 - r - s );
    dbdr =   4.0 - 8.0 * r - 4.0 * s;
    dbds =           - 4.0 * r;
  }
  else if ( inode == 5 )
  {
    *b   =   4.0 * r * s;
    dbdr =                           4.0 * s;
    dbds =             4.0 * r;
  }
  else if ( inode == 6 )
  {
    *b   =   4.0 * s * ( 1.0 - r - s );
    dbdr =                         - 4.0 * s;
    dbds =   4.0 - 4.0 * r - 8.0 * s;
  }
  else
  {
    printf ( "\n" );
    printf ( "QBF - Fatal error!\n" );
    printf ( "  Request for local basis function INODE = %d\n", inode );
    exit ( 1 );
  }
/*
  We need to convert the derivative information from (R(X,Y),S(X,Y))
  to (X,Y) using the chain rule.
*/
  *dbdx = dbdr * drdx + dbds * dsdx;
  *dbdy = dbdr * drdy + dbds * dsdy;

  return;
}
/******************************************************************************/

void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] )

/******************************************************************************/
/*
  Purpose:

    QUAD_A sets the quadrature rule for assembly.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double NODE_XY[2*NODE_NUM], the nodes.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int NODE_NUM, the number of nodes.

    Input, int NNODES, the number of nodes used to form one element.

    Output, double WQ[3], quadrature weights.

    Output, double XQ[3*ELEMENT_NUM], YQ[3*ELEMENT_NUM], the
    coordinates of the quadrature points in each element.
*/
{
  int element;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  wq[0] = 1.0 / 3.0;
  wq[1] = wq[0];
  wq[2] = wq[0];

  for ( element = 1; element <= element_num; element++ )
  {
    ip1 = element_node[0+(element-1)*nnodes];
    ip2 = element_node[1+(element-1)*nnodes];
    ip3 = element_node[2+(element-1)*nnodes];

    x1 = node_xy[0+(ip1-1)*2];
    x2 = node_xy[0+(ip2-1)*2];
    x3 = node_xy[0+(ip3-1)*2];

    y1 = node_xy[1+(ip1-1)*2];
    y2 = node_xy[1+(ip2-1)*2];
    y3 = node_xy[1+(ip3-1)*2];

    xq[0+(element-1)*3] = 0.5 * ( x1 + x2 );
    xq[1+(element-1)*3] = 0.5 * ( x2 + x3 );
    xq[2+(element-1)*3] = 0.5 * ( x1 + x3 );

    yq[0+(element-1)*3] = 0.5 * ( y1 + y2 );
    yq[1+(element-1)*3] = 0.5 * ( y2 + y3 );
    yq[2+(element-1)*3] = 0.5 * ( y1 + y3 );
  }

  return;
}
/******************************************************************************/

void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] )

/******************************************************************************/
/*
  Purpose:

    QUAD_E sets a quadrature rule for the error calculation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.

    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    index of local node I in element J.

    Input, int ELEMENT, the index of the element for which the quadrature
    points are to be computed.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int NNODES, the number of nodes used to form one element.

    Input, int NODE_NUM, the number of nodes.

    Input, int NQE, the number of points in the quadrature rule.
    This is actually fixed at 13.

    Output, double WQE[NQE], the quadrature weights.

    Output, double XQE[NQE], YQE[NQE], the X and Y coordinates
    of the quadrature points.
*/
{
  int i;
  int ii;
  int iii;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;
  double z4;
  double z5;
  double z6;
  double z7;

  for ( i = 1; i <= 3; i++ )
  {
    wqe[i-1] = 0.175615257433204;
    ii = i + 3;
    wqe[ii-1] = 0.053347235608839;
    ii = i + 6;
    iii = ii + 3;
    wqe[ii-1] = 0.077113760890257;
    wqe[iii-1] = wqe[ii-1];
  }

  wqe[13-1] = -0.14957004446767;

  z1 = 0.479308067841923;
  z2 = 0.260345966079038;
  z3 = 0.869739794195568;
  z4 = 0.065130102902216;
  z5 = 0.638444188569809;
  z6 = 0.312865496004875;
  z7 = 0.048690315425316;

  ip1 = element_node[0+(element-1)*nnodes];
  ip2 = element_node[1+(element-1)*nnodes];
  ip3 = element_node[2+(element-1)*nnodes];

  x1 = node_xy[0+(ip1-1)*2];
  x2 = node_xy[0+(ip2-1)*2];
  x3 = node_xy[0+(ip3-1)*2];

  y1 = node_xy[1+(ip1-1)*2];
  y2 = node_xy[1+(ip2-1)*2];
  y3 = node_xy[1+(ip3-1)*2];

  xqe[ 1-1] = z1 * x1 + z2 * x2 + z2 * x3;
  yqe[ 1-1] = z1 * y1 + z2 * y2 + z2 * y3;
  xqe[ 2-1] = z2 * x1 + z1 * x2 + z2 * x3;
  yqe[ 2-1] = z2 * y1 + z1 * y2 + z2 * y3;
  xqe[ 3-1] = z2 * x1 + z2 * x2 + z1 * x3;
  yqe[ 3-1] = z2 * y1 + z2 * y2 + z1 * y3;
  xqe[ 4-1] = z3 * x1 + z4 * x2 + z4 * x3;
  yqe[ 4-1] = z3 * y1 + z4 * y2 + z4 * y3;
  xqe[ 5-1] = z4 * x1 + z3 * x2 + z4 * x3;
  yqe[ 5-1] = z4 * y1 + z3 * y2 + z4 * y3;
  xqe[ 6-1] = z4 * x1 + z4 * x2 + z3 * x3;
  yqe[ 6-1] = z4 * y1 + z4 * y2 + z3 * y3;
  xqe[ 7-1] = z5 * x1 + z6 * x2 + z7 * x3;
  yqe[ 7-1] = z5 * y1 + z6 * y2 + z7 * y3;
  xqe[ 8-1] = z5 * x1 + z7 * x2 + z6 * x3;
  yqe[ 8-1] = z5 * y1 + z7 * y2 + z6 * y3;
  xqe[ 9-1] = z6 * x1 + z5 * x2 + z7 * x3;
  yqe[ 9-1] = z6 * y1 + z5 * y2 + z7 * y3;
  xqe[10-1] = z6 * x1 + z7 * x2 + z5 * x3;
  yqe[10-1] = z6 * y1 + z7 * y2 + z5 * y3;
  xqe[11-1] = z7 * x1 + z5 * x2 + z6 * x3;
  yqe[11-1] = z7 * y1 + z5 * y2 + z6 * y3;
  xqe[12-1] = z7 * x1 + z6 * x2 + z5 * x3;
  yqe[12-1] = z7 * y1 + z6 * y2 + z5 * y3;
  xqe[13-1] = ( x1 + x2 + x3 ) / 3.0;
  yqe[13-1] = ( y1 + y2 + y3 ) / 3.0;

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

double r8_huge ( void )

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

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
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
  value = s * ( int ) ( r8_abs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

int r8gb_fa ( int n, int ml, int mu, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GB_FA performs a LINPACK-style PLU factorization of a R8GB matrix.

  Discussion:

    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
    which may be required to store nonzero entries generated during Gaussian 
    elimination.

    The original M by N matrix is "collapsed" downward, so that diagonals
    become rows of the storage array, while columns are preserved.  The
    collapsed array is logically 2*ML+MU+1 by N.  

    The two dimensional array can be further reduced to a one dimensional
    array, stored by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2012

  Author:

    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, int ML, MU, the lower and upper bandwidths.
    ML and MU must be nonnegative, and no greater than N-1.

    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.  
    On output, A has been overwritten by the LU factors.

    Output, int PIVOT[N], the pivot vector.

    Output, int R8GB_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int col = 2 * ml + mu + 1;
  int i;
  int i0;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  double t;

  m = ml + mu + 1;
/*
  Zero out the initial fill-in columns.
*/
  j0 = mu + 2;
  j1 = i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      a[i-1+(jz-1)*col] = 0.0;
    }
  }

  jz = j1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
/*
  Zero out the next fill-in column.
*/
    jz = jz + 1;
    if ( jz <= n ) 
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(jz-1)*col] = 0.0;
      }
    }
/*
  Find L = pivot index.
*/
    lm = i4_min ( ml, n-k );
    l = m;

    for ( j = m+1; j <= m + lm; j++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*col] ) < r8_abs ( a[j-1+(k-1)*col] ) )
      {
        l = j;
      }
    }

    pivot[k-1] = l + k - m;
/*
  Zero pivot implies this column already triangularized.
*/
    if ( a[l-1+(k-1)*col] == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GB_FA - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", k );
      exit ( 1 );
    }
/*
  Interchange if necessary.
*/
    t                = a[l-1+(k-1)*col];
    a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
    a[m-1+(k-1)*col] = t;
/*
  Compute multipliers.
*/
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
    }
/*
  Row elimination with column indexing.
*/
    ju = i4_max ( ju, mu + pivot[k-1] );
    ju = i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;

      if ( l != mm )
      {
        t                 = a[l-1+(j-1)*col];
        a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
        a[mm-1+(j-1)*col] = t;
      }
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col] 
          + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
      }
    }
  }

  pivot[n-1] = n;

  if ( a[m-1+(n-1)*col] == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8GB_FA - Fatal error!\n" );
    fprintf ( stderr, "  Zero pivot on step %d\n", n );
    exit ( 1 );
  }

  return 0;
}
/******************************************************************************/

void r8gb_print_some ( int m, int n, int ml, int mu, double a[], int ilo, 
  int jlo, int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8GB_PRINT_SOME prints some of a R8GB matrix.

  Discussion:

    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
    which may be required to store nonzero entries generated during Gaussian 
    elimination.

    The original M by N matrix is "collapsed" downward, so that diagonals
    become rows of the storage array, while columns are preserved.  The
    collapsed array is logically 2*ML+MU+1 by N.  

    The two dimensional array can be further reduced to a one dimensional
    array, stored by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 February 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int ML, MU, the lower and upper bandwidths.
    ML and MU must be nonnegative, and no greater than min(M,N)-1..

    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int col = 2 * ml + mu + 1;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    printf ( "\n" );
    printf ( "  Col: " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%7d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu - ml );

    i2hi = i4_min ( ihi, m );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%6d  ", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i < j - mu - ml || j + ml < i )
        {
          printf ( "            " );
        }
        else
        {
          printf ( "%10g  ", a[i-j+ml+mu+(j-1)*col] );
        }
      }
      printf ( "\n" );
    }
  }

  printf ( "\n" );

  return;
# undef INCX
}
/******************************************************************************/

double *r8gb_sl ( int n, int ml, int mu, double a_lu[], int pivot[], 
  double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GB_SL solves a system factored by R8GB_FA.

  Discussion:

    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
    which may be required to store nonzero entries generated during Gaussian 
    elimination.

    The original M by N matrix is "collapsed" downward, so that diagonals
    become rows of the storage array, while columns are preserved.  The
    collapsed array is logically 2*ML+MU+1 by N.  

    The two dimensional array can be further reduced to a one dimensional
    array, stored by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 February 2012

  Author:

    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, int ML, MU, the lower and upper bandwidths.
    ML and MU must be nonnegative, and no greater than N-1.

    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA.

    Input, int PIVOT[N], the pivot vector from R8GB_FA.

    Input, double B[N], the right hand side vector.

    Input, int JOB.
    0, solve A * x = b.
    nonzero, solve A' * x = b.

    Output, double R8GB_SL[N], the solution.
*/
{
  int col = 2 * ml + mu + 1;
  int i;
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  m = mu + ml + 1;
/*
  Solve A * x = b.
*/
  if ( job == 0 )
  {
/*
  Solve L * Y = B.
*/
    if ( 1 <= ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
        for ( i = 1; i <= lm; i++ )
        {
          x[k+i-1] = x[k+i-1] + x[k-1] * a_lu[m+i-1+(k-1)*col];
        }
      }
    }
/*
  Solve U * X = Y.
*/
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a_lu[la+i-1+(k-1)*col];
      }
    }
  }
/*
  Solve A' * X = B.
*/
  else
  {
/*
  Solve U' * Y = B.
*/
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - x[lb+i-1] * a_lu[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
    }
/*
  Solve L' * X = Y.
*/
    if ( 1 <= ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 1; i <= lm; i++ )
        {
          x[k-1] = x[k-1] + x[k+i-1] * a_lu[m+i-1+(k-1)*col];
        }
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
      }
    }
  }

  return x;
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

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_SOME prints "some" of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N], the vector to be printed.

    Input, integer I_LO, I_HI, the first and last indices to print.
    The routine expects 1 <= I_LO <= I_HI <= N.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i-1] );
  }

  return;
}
/******************************************************************************/

double rhs ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    RHS gives the right-hand side of the differential equation.

  Discussion:

    The function specified here depends on the problem being
    solved.  This is one of the routines that a user will
    normally want to change.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the coordinates of a point
    in the region, at which the right hand side of the
    differential equation is to be evaluated.

    Output, double RHS, the value of the right
    hand side of the differential equation at (X,Y).
*/
{
# define PI 3.14159265358979323846264338327950288419716939937510

  double value;

  value = 2.0 * PI * PI * sin ( PI * x ) * sin ( PI * y );

  return value;
# undef PI
}
/******************************************************************************/

void solution_write ( double f[], int indx[], int node_num, int nunk,
  char *output_filename, double node_xy[] )

/******************************************************************************/
/*
  Purpose:

    SOLUTION_WRITE writes the solution to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, double F[NUNK], the coefficients of the solution.

    Input, int INDX[NODE_NUM], gives the index of the unknown quantity
    associated with the given node.

    Input, int NODE_NUM, the number of nodes.

    Input, int NUNK, the number of unknowns.

    Input, char *OUTPUT_FILENAME, the name of the file
    in which the data should be stored.

    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
*/
{
  double dudx;
  double dudy;
  int node;
  FILE *output;
  double u;
  double x;
  double y;

  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "SOLUTION_WRITE - Warning!\n" );
    printf ( "  Could not write the solution file.\n" );
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    if ( 0 < indx[node] )
    {
      u = f[indx[node]-1];
    }
    else
    {
      exact ( x, y, &u, &dudx, &dudy );
    }

    fprintf ( output, "%14g\n", u );
  }

  fclose ( output );

  return;
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

void triangulation_order6_plot ( char *file_name, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show )

/******************************************************************************/
/*
  Purpose:

    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.

  Discussion:

    The triangulation is most usually a Delaunay triangulation,
    but this is not necessary.

    This routine has been specialized to deal correctly ONLY with
    a mesh of 6 node elements, with the property that starting
    from local node 1 and traversing the edges of the element will
    result in encountering local nodes 1, 4, 2, 5, 3, 6 in that order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to create.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.

    Input, int TRIANGLE_NUM, the number of triangles.

    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists, for each triangle,
    the indices of the nodes that form the vertices and midsides
    of the triangle.

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
  int circle_size;
  int delta;
  int e;
  FILE *file_unit;
  int i;
  int ip1;
  int node;
  int order[6] = { 1, 4, 2, 5, 3, 6 };
  int triangle;
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
  x_max = - r8_huge ( );
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
    printf ( "\n" );
    printf ( "TRIANGULATION_ORDER6_PLOT - Fatal error!\n" );
    printf ( "  Could not open the output EPS file.\n" );
    exit ( 1 );
  }

  fprintf ( file_unit, "%%!PS-Adobe-3.0 EPSF-3.0\n" );
  fprintf ( file_unit, "%%%%Creator: triangulation_order6_plot.C\n" );
  fprintf ( file_unit, "%%%%Title: %s\n", file_name );

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
  fprintf ( file_unit, "%%  Increase line width from default 0.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "2 setlinewidth\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Set the RGB line color to very light gray.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, " 0.9000 0.9000 0.9000 setrgbcolor\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% Draw a gray border around the page.\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "newpath\n" );
  fprintf ( file_unit, "%d  %d  moveto\n", x_ps_min, y_ps_min );
  fprintf ( file_unit, "%d  %d  lineto\n", x_ps_max, y_ps_min );
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
  fprintf ( file_unit, "  210  702 moveto\n" );
  fprintf ( file_unit, "(Pointset) show\n" );
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
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
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
  }
/*
  Label the nodes.
*/
  if ( 2 <= node_show )
  {
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Label the nodes:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Set the color to darker blue:\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "0.000  0.250  0.850  setrgbcolor\n" );
    fprintf ( file_unit, "/Times-Roman findfont\n" );
    fprintf ( file_unit, "0.20 inch scalefont\n" );
    fprintf ( file_unit, "setfont\n" );

    fprintf ( file_unit, "%\n" );

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

      fprintf ( file_unit, "newpath  %d  %d  moveto (%d) show\n", x_ps, y_ps + 5, node + 1 );
    }
  }
/*
  Draw the triangles.
*/
  if ( 1 <= triangle_show )
  {
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Set the RGB color to red.\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "0.900  0.200  0.100 setrgbcolor\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Draw the triangles.\n" );
    fprintf ( file_unit, "%%\n" );

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      node = triangle_node[order[0]-1+triangle*6] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      fprintf ( file_unit, "newpath  %d  %d  moveto\n", x_ps, y_ps );

      for ( i = 1; i <= 6; i++ )
      {
        ip1 = ( i % 6 ) + 1;
        node = triangle_node[order[ip1-1]-1+triangle*6] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        fprintf ( file_unit, "%d  %d  lineto\n", x_ps, y_ps );
      }
      fprintf ( file_unit, "stroke\n" );
    }
  }
/*
  Label the triangles.
*/
  if ( 2 <= triangle_show )
  {
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Label the triangles.\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "%%  Set the RGB color to darker red.\n" );
    fprintf ( file_unit, "%%\n" );
    fprintf ( file_unit, "0.950  0.250  0.150 setrgbcolor\n" );
    fprintf ( file_unit, "/Times-Roman findfont\n" );
    fprintf ( file_unit, "0.20 inch scalefont\n" );
    fprintf ( file_unit, "setfont\n" );
    fprintf ( file_unit, "%%\n" );

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 0; i < 6; i++ )
      {
        node = triangle_node[i+triangle*6] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }

      ave_x = ave_x / 6.0;
      ave_y = ave_y / 6.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      fprintf ( file_unit, "%d  %d  moveto (%d) show\n", x_ps, y_ps, triangle + 1 );
    }
  }

  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "restore showpage\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%% End of page\n" );
  fprintf ( file_unit, "%%\n" );
  fprintf ( file_unit, "%%%%Trailer\n" );
  fprintf ( file_unit, "%%%%EOF\n" );

  fclose ( file_unit );

  return;
}
/******************************************************************************/

void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] )

/******************************************************************************/
/*
  Purpose:

    XY_SET sets the XY coordinates of the nodes.

  Discussion:

    The nodes are laid out in an evenly spaced grid, in the unit square.

    The first node is at the origin.  More nodes are created to the
    right until the value of X = 1 is reached, at which point
    the next layer is generated starting back at X = 0, and an
    increased value of Y.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of elements in the X and
    Y direction.

    Input, int NODE_NUM, the number of nodes.

    Input, double XL, XR, YB, YT, the X coordinates of
    the left and right sides of the rectangle, and the Y coordinates
    of the bottom and top of the rectangle.

    Output, double NODE_XY[2*NODE_NUM], the nodes.
*/
{
  int i;
  int j;

  for ( j = 1; j <= 2*ny-1; j++ )
  {
    for ( i = 1; i <= 2*nx - 1; i++ )
    {
      node_xy[0+(i-1+(j-1)*(2*nx-1))*2] =
        ( ( double ) ( 2 * nx - i - 1 ) * xl
        + ( double ) (          i - 1 ) * xr )
        / ( double ) ( 2 * nx     - 2 );

      node_xy[1+(i-1+(j-1)*(2*nx-1))*2] =
        ( ( double ) ( 2 * ny - j - 1 ) * yb
        + ( double ) (          j - 1 ) * yt )
        / ( double ) ( 2 * ny     - 2 );

    }
  }

  return;
}
