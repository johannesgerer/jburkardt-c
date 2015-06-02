/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 1.  It is not written to be comprehensible without the
explanation in that book.

Input: 2n integer coordinates for vertices of a simple polygon,
       in counterclockwise order.  NB: the code will not work
       for points in clockwise order!
Output: the diagonals of a triangulation, in PostScript.

Written by Joseph O'Rourke, with contributions by Min Xu.
Last modified: October 1997
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/
# include <stdio.h>
# include <stdlib.h>
# include <stdbool.h>
# include <math.h>

# define  X 0
# define  Y 1
/*
  typedef  enum { FALSE, TRUE } bool;
*/
#define  DIM 2               /* Dimension of points */
typedef  int tPointi[DIM];   /* Type integer point */

typedef  struct tVertexStructure tsVertex;
typedef  tsVertex *tVertex;
struct  tVertexStructure {
   int    vnum;    /* Index */
   tPointi  v;    /* Coordinates */
   bool   ear;    /* true iff an ear */
   tVertex   next,prev;
};

/*
  Global variable definitions
*/
tVertex  vertices  = NULL;  /* "Head" of circular list. */
/*
  Functions:
*/
int main ( );
int area_poly2 ( );
int area_sign ( tPointi a, tPointi b, tPointi c );
int area2 ( tPointi a, tPointi b, tPointi c );
bool between ( tPointi a, tPointi b, tPointi c );
bool collinear ( tPointi a, tPointi b, tPointi c );
bool diagonal ( tVertex a, tVertex b );
bool diagonalie ( tVertex a, tVertex b );
void ear_init ( );
bool in_cone ( tVertex a, tVertex b );
bool intersect ( tPointi a, tPointi b, tPointi c, tPointi d );
bool intersect_prop ( tPointi a, tPointi b, tPointi c, tPointi d );
bool left ( tPointi a, tPointi b, tPointi c );
bool left_on ( tPointi a, tPointi b, tPointi c );
tVertex make_null_vertex ( );
void print_poly ( );
void print_vertices ( int nvertices, int xmin, int xmax, int ymin, int ymax, 
  int scale );
void read_vertices ( int *nvertices );
void scale_data ( int *xmin, int *xmax, int *ymin, int *ymax, int *scale );
void triangulate ( int nvertices, int xmin, int xmax, int ymin, int ymax, 
  int scale );
bool xor ( bool x, bool y );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGULATE.

  Modified:

    02 January 2013

  Author:

    Joseph O'Rourke

  Reference:

    Joseph ORourke,
    Computational Geometry,
    Second Edition,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
*/
{
  tVertex a;
  double area;
  int nvertices;
  tVertex p;
  int scale;
  int xmax;
  int xmin;
  int ymax;
  int ymin;

  read_vertices ( &nvertices );

  scale_data ( &xmin, &xmax, &ymin, &ymax, &scale );

  print_vertices ( nvertices, xmin, xmax, ymin, ymax, scale );

  area = 0.5 * ( double ) area_poly2 ( );

  printf ( "%%Area of polygon = %g\n", area );

  if ( area == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRIANGULATE - Fatal error!\n" );
    fprintf ( stderr, "  Computed area of polygon is zero.\n" );
    exit ( 1 );
  }

  if ( area < 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRIANGULATE - Fatal error!\n" );
    fprintf ( stderr, "  Computed area of polygon is negative.\n" );
    fprintf ( stderr, "  The vertices are probably listed in clockwise order.\n" );
    fprintf ( stderr, "  Revise the input file by reversing the vertex order.\n" );
    exit ( 1 );
  }
/*
  Refuse to accept input if two consecutive vertices are equal.
*/
  p = vertices;
  a = p->next;
  do
  {
    if ( p->v[X] == a->v[X] && p->v[Y] == a->v[Y] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "TRIANGULATE - Fatal error!\n" );
      fprintf ( stderr, "  Vertices %d and %d are equal.\n", p->vnum, a->vnum );
      exit ( 1 );
    }
    p = a;
    a = p->next;
  } while ( a->next != vertices );
/*
  Compute the triangulation.
*/
  triangulate ( nvertices, xmin, xmax, ymin, ymax, scale );

  printf ( "showpage\n%%%%EOF\n" );
/*
  Terminate.
*/
  return 0;
}
/******************************************************************************/

int area_poly2 ( )

/******************************************************************************/
/*
  Purpose:

    AREA_POLY2 returns the area of a polygon.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Output, int AREA_POLY2, the area of the polygon.
*/
{
  tVertex a;
  tVertex p;
  int sum = 0;
/*
  Keep P fixed, but allow A to move.
*/
  p = vertices;
  a = p->next;

  do
  {
    sum = sum + area2 ( p->v, a->v, a->next->v );
    a = a->next;
  } while ( a->next != vertices );

  return sum;
}
/******************************************************************************/

int area_sign ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    AREA_SIGN returns the sign of the area defined by three points.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points that define a triangle.

    Output, int AREA_SIGN, the sign of the area of the triangle.
*/
{
  double area;

  area = ( b[0] - a[0] ) * ( double ) ( c[1] - a[1] ) -
         ( c[0] - a[0] ) * ( double ) ( b[1] - a[1] );
/*
  The area should be an integer.
*/
  if ( 0.5 < area )
  {
    return  1;
  }
  else if ( area < -0.5 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}
/******************************************************************************/

int area2 ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    AREA2 returns twice the signed area of a triangle.

  Discussion:

    The area is positive if points A, B, and C are oriented counter clockwise,
    negative if clockwise, and zero if the points are collinear.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points that define a triangle.

    Output, int AREA2, twice the signed area of the triangle.
*/
{
  int value;

  value = ( b[X] - a[X] ) * ( c[Y] - a[Y] ) -
          ( c[X] - a[X] ) * ( b[Y] - a[Y] );

  return value;
}
/******************************************************************************/

bool between ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    BETWEEN returns TRUE if and only if point C lies on the closed segement AB.

  Discussion:

    The function first checks whether C is collinear with A and B.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points to be tested.

    Output, bool BETWEEN is TRUE if point C lies on the closed
    segment between A and B.
*/
{
  bool value;

  if ( ! collinear ( a, b, c ) )
  {
    value = false;
    return value;
  }
/*
  If AB not vertical, check betweenness on x; else on y.
*/
  if ( a[X] != b[X] )
  {
    value = ( (a[X] <= c[X] ) && ( c[X] <= b[X] ) ) ||
            ( (a[X] >= c[X] ) && ( c[X] >= b[X] ) );
  }
  else
  {
    value = ( ( a[Y] <= c[Y] ) && ( c[Y] <= b[Y] ) ) ||
            ( ( a[Y] >= c[Y] ) && ( c[Y] >= b[Y] ) );
  }

  return value;
}
/******************************************************************************/

bool collinear ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    COLLINEAR is TRUE if the points A, B and C are collinear.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points to be tested.

    Output, bool COLLINEAR, is TRUE if the points are collinear.
*/
{
  bool value;

  value = ( area_sign ( a, b, c ) == 0 );

  return value;
}
/******************************************************************************/

bool diagonal ( tVertex a, tVertex b )

/******************************************************************************/
/*
  Purpose:

    DIAGONAL returns TRUE iff (A,B) is a proper internal diagonal of the polygon.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tVertex A, B, two vertices of the polygon.

    Output, bool DIAGONAL, is TRUE if the line connecting A and B is a
    proper internal diagonal of the polygon.
*/
{
  bool value;

  value = in_cone ( a, b ) && 
          in_cone ( b, a ) && 
          diagonalie ( a, b );

  return value;
}
/******************************************************************************/

bool diagonalie ( tVertex a, tVertex b )

/******************************************************************************/
/*
  Purpose:

    DIAGONALIE returns TRUE iff (A,B) is a proper diagonal of a polygon.

  Discussion:

    (A,B) may be an internal or external diagonal of the polygon, ignoring edges
    incident to A and B.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tVertex A, B, two vertices of the polygon.

    Output, bool DIAGONALIE, is TRUE if the line connecting A and B is a
    proper diagonal of the polygon.
*/
{
  tVertex c;
  tVertex c1;
/*
  For each edge (C,C1) of P.
*/
  c = vertices;
  do
  {
    c1 = c->next;
/*
  Skip edges incident to A or B.
*/
    if ( ( c != a ) && ( c1 != a ) && ( c != b ) && ( c1 != b ) &&
           intersect( a->v, b->v, c->v, c1->v ) )
    {
      return false;
    }
    c = c->next;
  } while ( c != vertices );

  return true;
}
/******************************************************************************/

void ear_init ( )

/******************************************************************************/
/*
  Purpose:

    EAR_INIT initializes the data structures, and calls Triangulate2 to clip ears.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Local Parameters:

    Local, tVertex V0, V1, V2, three consecutive vertices of the polygon.
*/
{
  tVertex v0;
  tVertex v1;
  tVertex v2;
/*
  Initialize v1->ear for all vertices.
*/
  v1 = vertices;
  do
  {
    v2 = v1->next;
    v0 = v1->prev;
    v1->ear = diagonal ( v0, v2 );
    v1 = v1->next;
  } while ( v1 != vertices );

  return;
}
/******************************************************************************/

bool in_cone ( tVertex a, tVertex b )

/******************************************************************************/
/*
  Purpose:

    IN_CONE returns TRUE iff the diagonal (A,B) is strictly internal.

  Discussion:

    More correctly, the diagonal (A,B) must be strictly internal to the
    polygon in the neighborhood of the A endpoint.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tVertex A, B, two vertices of the polygon.

    Output, bool IN_CONE, is TRUE if the line connecting A and B is
    strictly internal to the polygon in the neighborhood of A.

  Local Parameters:

    Local, tVertex A0, A1, are the vertices before and after A.
*/
{
  tVertex a0;
  tVertex a1;
  bool value;

  a1 = a->next;
  a0 = a->prev;
/*
  If A is a convex vertex ...
*/
  if ( left_on ( a->v, a1->v, a0->v ) )
  {
    value = left ( a->v, b->v, a0->v ) &&
            left ( b->v, a->v, a1->v );
  }
/*
  Else A is reflex vertex:
*/
  else
  {
    value = !( left_on ( a->v, b->v, a1->v ) &&
               left_on ( b->v, a->v, a0->v ) );
  }
  return value;
}
/******************************************************************************/

bool intersect ( tPointi a, tPointi b, tPointi c, tPointi d )

/******************************************************************************/
/*
  Purpose:

    INTERSECT returns TRUE iff segments AB and CD intersect.

  Discussion:

    The intersection may be proper or improper.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, D, points that define the segments.

    Output, bool INTERSECT, is TRUE if segments AB and CD intersect.
*/
{
  if ( intersect_prop ( a, b, c, d ) )
  {
    return true;
  }
  else if ( between ( a, b, c )
         || between ( a, b, d )
         || between ( c, d, a )
         || between ( c, d, b ) )
  {
    return true;
  }
  else
  {
    return false;
  }
}
/******************************************************************************/

bool intersect_prop ( tPointi a, tPointi b, tPointi c, tPointi d )

/******************************************************************************/
/*
  Purpose:

    INTERSECT_PROP returns true if and only if AB properly intersects CD.

  Discussion:

    AB and CD must share a point interior to both segments.  The properness of
    the intersection is ensured by using strict leftness.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, D, points that define segments AB and CD.

    Output, bool INTERSECT_PROP, is TRUE if AB properly intersects CD.
*/
{
  bool value;
/*
  Eliminate improper cases.
*/
  if ( collinear ( a, b, c ) ||
       collinear ( a, b, d ) ||
       collinear ( c, d, a ) ||
       collinear ( c, d, b ) )
  {
    value = false;
  }
  else
  {
    value = xor ( left ( a, b, c ), left ( a, b, d ) ) &&
            xor ( left ( c, d, a ), left ( c, d, b ) );
  }
  return value;
}
/******************************************************************************/

bool left ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    LEFT is TRUE if C is on the left side of the line from A to B.

  Discussion:

    More correctly, the function returns true if and only if C is strictly
    to the left of the directed line through A to B.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points to be tested.

    Output, bool LEFT, is TRUE if C is strictly to the left of the directed
    line from A to B.
*/
{
  bool value;

  value = ( area_sign ( a, b, c ) > 0 );

  return value;
}
/******************************************************************************/

bool left_on ( tPointi a, tPointi b, tPointi c )

/******************************************************************************/
/*
  Purpose:

    LEFT_ON is TRUE if C is to the left side, or on, the line from A to B.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, tPointi A, B, C, three points to be tested.

    Output, bool LEFT_ON, is TRUE if C is strictly to the left of, or on,
    the directed line from A to B.
*/
{
  bool value;

  value = ( area_sign ( a, b, c ) >= 0 );

  return value;
}
/******************************************************************************/

tVertex make_null_vertex ( )

/******************************************************************************/
/*
  Purpose:

    MAKE_NULL_VERTEX makes a vertex.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Output, tVertex MAKE_NULL_VERTEX, the new vertex.
*/
{
  tVertex v;

  v = ( tsVertex * ) malloc ( sizeof ( tsVertex ) );

  if ( v == NULL )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MAKE_NULL_VERTEX - Fatal error!\n" );
    fprintf ( stderr, "  Out of memory!\n" );
    exit ( 1 );
  }

  if ( vertices )
  {
    v->next = vertices;
    v->prev = vertices->prev;
    vertices->prev = v;
    v->prev->next = v;
  }
  else
  {
    vertices = v;
    vertices->next = v;
    vertices->prev = v;
  }

  return v;
}
/******************************************************************************/

void print_poly ( )

/******************************************************************************/
/*
  Purpose:

    PRINT_POLY prints the polygon data.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke
*/
{
  tVertex  v;

  printf ( "%%\n" );
  printf ( "%%Polygon circular list:\n" );
  printf ( "%%\n" );

  v = vertices;
  do
  {
    printf( "%% vnum=%5d:  ear=%d\n", v->vnum, v->ear );
    v = v->next;
  } while ( v != vertices );

  printf ( "%%\n" );

  return;
}
/******************************************************************************/

void print_vertices ( int nvertices, int xmin, int xmax, int ymin, int ymax, 
  int scale )

/******************************************************************************/
/*
  Purpose:

    PRINT_VERTICES prints the vertices.

  Discussion:

    This function uses the VNUM indices corresponding to the order in which
    the vertices were input.  The output is in PostScript format.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, int NVERTICES, the number of vertices.

    Input, int XMIN, XMAX, YMIN, YMAX, the minimum and maximum
    X and Y values of the coordinates of the vertices of the polygon.

    Input, int SCALE, an appropriate scaling for the data.
*/
{
/*
  Pointers to vertices, edges, faces.
*/
  tVertex v;
  int x;
  int y;
/*
  PostScript header
*/
  printf ( "%%!PS\n" );
  printf ( "%%%%Creator: triangulate.c (Joseph O'Rourke)\n" );
  printf ( "%%%%BoundingBox: %d %d %d %d\n",
    0, 0, 72 + scale * ( xmax - xmin ), 72 + scale * ( ymax - ymin ) );
  printf ( "%%%%EndComments\n" );
  printf ( "1 1 setlinewidth\n" );
/*
  Output vertex information as a PostScript comment.
*/
  printf ( "\n" );
  printf ( "%% number of vertices = %d\n", nvertices );

  v = vertices;
  do
  {
     printf ( "%% vnum=%5d:\tx=%5d\ty=%5d\n", v->vnum, v->v[X], v->v[Y] );
     v = v->next;
  } while ( v != vertices );
/*
  Draw the polygon.
*/
  printf ( "\n%%Polygon:\n" );
  printf ( "newpath\n" );
  v = vertices;
  x = 36 + scale * ( v->v[X] - xmin );
  y = 36 + scale * ( v->v[Y] - ymin );
  printf ( "%d\t%d\tmoveto\n", x, y );
  v = v->next;
  do
  {
    x = 36 + scale * ( v->v[X] - xmin );
    y = 36 + scale * ( v->v[Y] - ymin );
    printf ( "%d\t%d\tlineto\n", x, y );
    v = v->next;
  } while ( v != vertices );

  printf ( "closepath stroke\n" );

  return;
}
/******************************************************************************/

void read_vertices ( int *nvertices )

/******************************************************************************/
/*
  Purpose:

    READ_VERTICES reads the polygon vertice from standard input.

  Discussion:

    After reading the vertices, the function links them into a circular
    list with MAKE_NULL_VERTEX.  There is no need for the # of vertices
    to be the first line: the function looks for EOF instead.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Output, int *NVERTICES, the number of vertices.
*/
{
  tVertex v;
  int vnum = 0;
  int x;
  int y;

  while ( scanf ( "%d %d", &x, &y ) != EOF )
  {
    v = make_null_vertex ( );
    v->v[X] = x;
    v->v[Y] = y;
    v->vnum = vnum++;
  }

  if ( vnum < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "READ_VERTICES - Fatal error!\n" ),
    fprintf ( stderr, "  Number of vertices = %d < 3\n", vnum ),
    exit ( 1 );
  }

  *nvertices = vnum;

  return;
}
/******************************************************************************/

void scale_data ( int *xmin, int *xmax, int *ymin, int *ymax, int *scale )

/******************************************************************************/
/*
  Purpose:

    SCALE_DATA determines the scale for the polygonal data.

  Discussion:

    The integer scale factor SCALE is determined so that
      36 + SCALE * ( X - XMIN )
    and
      36 + SCALE * ( Y - YMIN )
    should fall between 36 and 540 for all reasonable values of X and Y .

    This makes the PostScript image appear in a 7.5 inch by 7.5 inch
    square on the paper.

  Modified:

    30 April 2007

  Author:

    John Burkardt

  Parameters:

    Output, int *XMIN, *XMAX, *YMIN, *YMAX, the minimum and maximum
    X and Y values of the coordinates of the vertices of the polygon.

    Output, int *SCALE, an appropriate scaling for the data.
*/
{
  tVertex v;
  int range;
/*
  Compute bounding box for Encapsulated PostScript.
*/
  v = vertices;

  *xmin = v->v[X];
  *xmax = v->v[X];
  *ymin = v->v[Y];
  *ymax = v->v[Y];

  do
  {
    if      ( v->v[X] > *xmax )
    {
      *xmax = v->v[X];
    }
    else if ( v->v[X] < *xmin )
    {
      *xmin = v->v[X];
    }

    if      ( v->v[Y] > *ymax )
    {
      *ymax = v->v[Y];
    }
    else if ( v->v[Y] < *ymin )
    {
      *ymin = v->v[Y];
    }

    v = v->next;
  } while ( v != vertices );

  range = *xmax - *xmin;
  if ( range < *ymax - *ymin )
  {
    range = *ymax - *ymin;
  }

  *scale = ( int ) ( ( double ) 540 / ( double ) range );

  return;
}
/******************************************************************************/

void triangulate ( int nvertices, int xmin, int xmax, int ymin, int ymax, 
  int scale )

/******************************************************************************/
/*
  Purpose:

    TRIANGULATE prints N-3 diagonals which triangulate the polygon.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, int NVERTICES, the number of vertices.

    Input, int XMIN, XMAX, YMIN, YMAX, the minimum and maximum
    X and Y values of the coordinates of the vertices of the polygon.

    Input, int SCALE, an appropriate scaling for the data.
*/
{
  tVertex v0, v1, v2, v3, v4;  /* five consecutive vertices */
  int x;
  int y;
  int n = nvertices;
/*
  N is the "current" number of vertices;
  It starts at NVERTICES, but shrinks to 3.
*/
  ear_init ( );

  printf ( "\n" );
  printf ( "newpath\n" );
/*
  Each step of outer loop removes one ear.
*/
  while ( 3 < n )
  {
/*
  Inner loop searches for an ear.
*/
    v2 = vertices;
    do
    {
      if ( v2->ear )
      {
/*
  Ear found. Fill variables.
*/
        v3 = v2->next;
        v4 = v3->next;
        v1 = v2->prev;
        v0 = v1->prev;
/*
  (v1,v3) is a diagonal
*/
        printf ( "%%Diagonal: (%d,%d)\n", v1->vnum, v3->vnum );
        x = 36 + scale * ( v1->v[X] - xmin );
        y = 36 + scale * ( v1->v[Y] - ymin );
        printf ( "%d\t%d\tmoveto\n", x, y );
        x = 36 + scale * ( v3->v[X] - xmin );
        y = 36 + scale * ( v3->v[Y] - ymin );
        printf ( "%d\t%d\tlineto\n", x, y );
/*
  Update earity of diagonal endpoints.
*/
        v1->ear = diagonal ( v0, v3 );
        v3->ear = diagonal ( v1, v4 );
/*
  Cut off the ear v2.
*/
        v1->next = v3;
        v3->prev = v1;
        vertices = v3;
/*
  In case the head was v2.
*/
        n--;
        break;
      }
/*
  End if ear found.
*/
      v2 = v2->next;
    } while ( v2 != vertices );
  }

  printf ( "closepath stroke\n" );
  printf ( "\n" );

  return;
}
/******************************************************************************/

bool xor ( bool x, bool y )

/******************************************************************************/
/*
  Purpose:

    XOR returns the exclusive OR of two values.

  Discussion:

    The function is TRUE if exactly one input argument is TRUE.
    The arguments are negated to ensure that they are 0/1 values.
    Idea due to Michael Baldwin.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke

  Parameters:

    Input, bool X, Y, the arguments to be tested.

    Output, bool XOR, is TRUE if exactly one of X and Y is true.
*/
{
  return ( !x ^ !y );
}
