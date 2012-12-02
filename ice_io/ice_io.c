# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "ice_io.h"
# include "netcdf.h"

/******************************************************************************/

void cyl248_data ( int dim, int vertices, int edges, int triangles, 
  int quadrilaterals, int tetrahedrons, int hexahedrons, 
  double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
  int edge_label[], int triangle_vertex[], int triangle_label[], 
  int quadrilateral_vertex[], int quadrilateral_label[], 
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[], 
  int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    CYL248_DATA defines the data for a 3D tetrahedral mesh.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
    of each vertex.

    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Output, int EDGE_LABEL[EDGES], a label for each edge.

    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for 
    each quadrilateral.

    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for 
    each tetrahedron.

    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
# define TETRAHEDRONS_SAVE 248
# define TRIANGLES_SAVE 154
# define VERTICES_SAVE 92

  int i;

  int tetrahedron_vertex_save[4*TETRAHEDRONS_SAVE] = { 
    23, 1, 9, 8, 
    27, 9, 23, 1, 
    26, 8, 23, 9, 
    26, 9, 7, 8, 
    2, 9, 27, 1, 
    26, 9, 10, 7, 
    26, 28, 7, 10, 
    11, 29, 3, 2, 
    7, 6, 10, 28, 
    10, 6, 31, 28, 
    11, 29, 30, 3, 
    11, 30, 4, 3, 
    11, 30, 32, 4, 
    10, 6, 5, 31, 
    11, 5, 4, 32, 
    19, 33, 34, 20, 
    39, 22, 40, 16, 
    39, 17, 36, 22, 
    39, 22, 16, 17, 
    40, 22, 15, 16, 
    12, 19, 20, 33, 
    19, 20, 34, 18, 
    12, 33, 20, 35, 
    38, 37, 14, 21, 
    36, 22, 17, 18, 
    38, 14, 15, 21, 
    13, 14, 37, 21, 
    12, 20, 13, 35, 
    80, 32, 11, 30, 
    80, 28, 10, 31, 
    80, 31, 59, 28, 
    80, 58, 57, 26, 
    80, 28, 58, 26, 
    80, 59, 58, 28, 
    80, 28, 26, 10, 
    80, 10, 26, 9, 
    80, 9, 11, 10, 
    80, 9, 26, 23, 
    80, 23, 26, 57, 
    80, 23, 27, 9, 
    80, 23, 56, 27, 
    80, 30, 11, 29, 
    80, 5, 10, 11, 
    80, 5, 11, 32, 
    80, 5, 32, 31, 
    80, 31, 10, 5, 
    80, 2, 11, 9, 
    80, 29, 11, 2, 
    80, 2, 9, 27, 
    80, 27, 29, 2, 
    81, 40, 39, 22, 
    81, 22, 39, 36, 
    81, 18, 36, 34, 
    81, 34, 20, 18, 
    81, 22, 36, 18, 
    81, 20, 22, 18, 
    81, 37, 38, 21, 
    81, 20, 33, 35, 
    81, 13, 21, 20, 
    81, 13, 20, 35, 
    81, 13, 37, 21, 
    81, 35, 37, 13, 
    81, 20, 21, 22, 
    81, 34, 33, 20, 
    81, 21, 38, 15, 
    81, 38, 40, 15, 
    81, 22, 21, 15, 
    81, 15, 40, 22, 
    82, 60, 74, 59, 
    82, 74, 25, 59, 
    82, 73, 72, 58, 
    82, 25, 73, 58, 
    82, 59, 25, 58, 
    82, 58, 72, 57, 
    82, 57, 80, 58, 
    82, 58, 80, 59, 
    83, 71, 79, 70, 
    83, 70, 76, 78, 
    83, 79, 76, 70, 
    83, 79, 60, 76, 
    83, 82, 60, 74, 
    84, 54, 64, 55, 
    84, 64, 65, 55, 
    84, 65, 63, 55, 
    84, 65, 71, 63, 
    85, 29, 62, 30, 
    85, 80, 29, 30, 
    85, 29, 61, 62, 
    85, 78, 83, 76, 
    85, 78, 76, 30, 
    85, 62, 78, 30, 
    85, 76, 83, 60, 
    85, 76, 32, 30, 
    85, 32, 80, 30, 
    85, 32, 76, 60, 
    85, 27, 61, 29, 
    85, 80, 27, 29, 
    85, 83, 82, 60, 
    85, 77, 78, 62, 
    85, 60, 82, 59, 
    85, 59, 82, 80, 
    85, 32, 60, 31, 
    85, 80, 32, 31, 
    85, 60, 59, 31, 
    85, 59, 80, 31, 
    86, 51, 68, 52, 
    86, 69, 68, 51, 
    86, 68, 67, 52, 
    86, 52, 67, 53, 
    86, 67, 66, 53, 
    86, 53, 66, 54, 
    87, 50, 70, 49, 
    87, 71, 70, 50, 
    87, 63, 71, 50, 
    87, 63, 84, 71, 
    87, 70, 69, 49, 
    87, 71, 83, 70, 
    87, 49, 69, 51, 
    87, 69, 86, 51, 
    88, 64, 66, 73, 
    88, 72, 73, 66, 
    88, 72, 82, 73, 
    88, 24, 72, 66, 
    88, 64, 73, 25, 
    88, 73, 82, 25, 
    88, 66, 64, 54, 
    88, 84, 54, 64, 
    88, 87, 86, 84, 
    88, 67, 24, 66, 
    88, 66, 86, 67, 
    88, 64, 25, 65, 
    88, 65, 84, 64, 
    88, 25, 74, 65, 
    88, 25, 82, 74, 
    88, 83, 87, 71, 
    88, 71, 87, 84, 
    88, 82, 83, 74, 
    88, 74, 83, 71, 
    88, 65, 74, 71, 
    88, 71, 84, 65, 
    89, 86, 87, 84, 
    89, 39, 48, 44, 
    89, 44, 49, 43, 
    89, 44, 43, 36, 
    89, 44, 48, 50, 
    89, 48, 63, 50, 
    89, 86, 84, 54, 
    89, 51, 87, 86, 
    89, 44, 50, 49, 
    89, 50, 87, 49, 
    89, 43, 49, 51, 
    89, 49, 87, 51, 
    89, 39, 44, 36, 
    89, 36, 81, 39, 
    89, 63, 48, 47, 
    89, 47, 48, 40, 
    89, 46, 55, 47, 
    89, 38, 46, 47, 
    89, 55, 63, 47, 
    89, 55, 84, 63, 
    89, 43, 42, 34, 
    89, 43, 51, 42, 
    89, 45, 53, 54, 
    89, 53, 86, 54, 
    89, 45, 54, 46, 
    89, 42, 52, 41, 
    89, 41, 52, 53, 
    89, 52, 86, 53, 
    89, 42, 51, 52, 
    89, 51, 86, 52, 
    89, 46, 54, 55, 
    89, 54, 84, 55, 
    90, 56, 75, 61, 
    90, 24, 75, 56, 
    90, 27, 56, 61, 
    90, 61, 85, 27, 
    90, 75, 77, 61, 
    90, 80, 82, 57, 
    90, 85, 82, 80, 
    90, 57, 24, 56, 
    90, 72, 24, 57, 
    90, 57, 82, 72, 
    90, 80, 56, 27, 
    90, 85, 80, 27, 
    91, 85, 90, 77, 
    91, 86, 87, 69, 
    91, 78, 77, 69, 
    91, 83, 88, 82, 
    91, 90, 82, 88, 
    91, 67, 88, 86, 
    91, 88, 87, 86, 
    91, 87, 88, 83, 
    91, 83, 85, 78, 
    91, 78, 85, 77, 
    91, 77, 75, 68, 
    91, 77, 90, 75, 
    91, 69, 77, 68, 
    91, 68, 86, 69, 
    91, 68, 75, 67, 
    91, 67, 86, 68, 
    91, 24, 88, 67, 
    91, 90, 88, 24, 
    91, 69, 87, 70, 
    91, 87, 83, 70, 
    91, 75, 24, 67, 
    91, 75, 90, 24, 
    92, 89, 46, 45, 
    92, 41, 53, 45, 
    92, 89, 45, 53, 
    92, 89, 53, 41, 
    92, 89, 41, 42, 
    92, 35, 41, 45, 
    92, 33, 41, 35, 
    92, 35, 81, 33, 
    92, 35, 45, 37, 
    92, 81, 35, 37, 
    92, 34, 89, 42, 
    92, 81, 89, 34, 
    92, 33, 42, 41, 
    92, 37, 45, 46, 
    92, 37, 46, 38, 
    92, 81, 37, 38, 
    92, 33, 34, 42, 
    92, 33, 81, 34, 
    83, 74, 60, 71, 
    83, 60, 79, 71, 
    89, 39, 40, 48, 
    89, 39, 81, 40, 
    89, 36, 43, 34, 
    89, 34, 81, 36, 
    89, 63, 87, 50, 
    89, 84, 87, 63, 
    54, 88, 66, 86, 
    54, 88, 86, 84, 
    90, 72, 88, 24, 
    90, 82, 88, 72, 
    38, 47, 89, 40, 
    38, 89, 81, 40, 
    92, 46, 89, 38, 
    92, 89, 81, 38, 
    80, 23, 57, 56, 
    80, 57, 90, 56, 
    61, 85, 62, 77, 
    61, 90, 85, 77, 
    82, 85, 91, 83, 
    82, 90, 91, 85, 
    70, 91, 78, 83, 
    70, 78, 91, 69 };
  int triangle_label_save[TRIANGLES_SAVE] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
    2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3 };
  int triangle_vertex_save[3*TRIANGLES_SAVE] = {
    12, 20, 19, 
    12, 13, 20, 
    19, 20, 18, 
    20, 22, 18, 
    22, 17, 18, 
    22, 16, 17, 
    13, 21, 20, 
    13, 14, 21, 
    14, 15, 21, 
    22, 15, 16, 
    22, 21, 15, 
    20, 21, 22, 
    1, 9, 8, 
    2, 9, 1, 
    9, 7, 8, 
    2, 11, 9, 
    11, 2, 3, 
    11, 3, 4, 
    9, 10, 7, 
    7, 10, 6, 
    10, 5, 6, 
    11, 4, 5, 
    5, 10, 11, 
    9, 11, 10, 
    23, 1, 8, 
    26, 23, 8, 
    26, 8, 7, 
    27, 1, 23, 
    2, 1, 27, 
    26, 7, 28, 
    7, 6, 28, 
    27, 29, 2, 
    29, 3, 2, 
    29, 30, 3, 
    30, 4, 3, 
    6, 31, 28, 
    6, 5, 31, 
    5, 32, 31, 
    5, 4, 32, 
    12, 19, 33, 
    19, 34, 33, 
    19, 18, 34, 
    12, 33, 35, 
    12, 35, 13, 
    18, 36, 34, 
    36, 18, 17, 
    35, 37, 13, 
    13, 37, 14, 
    38, 14, 37, 
    38, 15, 14, 
    39, 36, 17, 
    39, 17, 16, 
    38, 40, 15, 
    40, 16, 15, 
    39, 16, 40, 
    33, 41, 35, 
    33, 42, 41, 
    33, 34, 42, 
    36, 43, 34, 
    43, 42, 34, 
    39, 44, 36, 
    44, 43, 36, 
    35, 45, 37, 
    35, 41, 45, 
    37, 46, 38, 
    37, 45, 46, 
    38, 47, 40, 
    38, 46, 47, 
    39, 48, 44, 
    39, 40, 48, 
    47, 48, 40, 
    44, 49, 43, 
    44, 50, 49, 
    44, 48, 50, 
    43, 51, 42, 
    43, 49, 51, 
    42, 52, 41, 
    42, 51, 52, 
    41, 53, 45, 
    41, 52, 53, 
    45, 54, 46, 
    45, 53, 54, 
    46, 55, 47, 
    46, 54, 55, 
    30, 32, 4, 
    23, 56, 27, 
    23, 57, 56, 
    23, 26, 57, 
    28, 58, 26, 
    58, 57, 26, 
    31, 59, 28, 
    59, 58, 28, 
    32, 60, 31, 
    60, 59, 31, 
    27, 61, 29, 
    27, 56, 61, 
    29, 62, 30, 
    29, 61, 62, 
    55, 63, 47, 
    63, 48, 47, 
    48, 63, 50, 
    54, 64, 55, 
    64, 65, 55, 
    65, 63, 55, 
    53, 66, 54, 
    66, 64, 54, 
    52, 67, 53, 
    67, 66, 53, 
    51, 68, 52, 
    68, 67, 52, 
    49, 69, 51, 
    69, 68, 51, 
    50, 70, 49, 
    70, 69, 49, 
    63, 71, 50, 
    71, 70, 50, 
    65, 71, 63, 
    64, 25, 65, 
    64, 73, 25, 
    64, 66, 73, 
    67, 24, 66, 
    24, 72, 66, 
    72, 73, 66, 
    68, 75, 67, 
    75, 24, 67, 
    69, 77, 68, 
    77, 75, 68, 
    70, 78, 69, 
    78, 77, 69, 
    62, 78, 30, 
    78, 76, 30, 
    76, 32, 30, 
    32, 76, 60, 
    61, 77, 62, 
    77, 78, 62, 
    56, 75, 61, 
    75, 77, 61, 
    57, 24, 56, 
    24, 75, 56, 
    58, 72, 57, 
    72, 24, 57, 
    59, 25, 58, 
    25, 73, 58, 
    73, 72, 58, 
    60, 74, 59, 
    74, 25, 59, 
    25, 74, 65, 
    65, 74, 71, 
    70, 76, 78, 
    71, 79, 70, 
    79, 76, 70, 
    79, 60, 76, 
    74, 60, 71, 
    60, 79, 71 };
  int vertex_label_save[VERTICES_SAVE] = {
    3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 
    4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
    3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0 };
  double vertex_coordinate_save[3*VERTICES_SAVE] = { 
  1.0,       0.2,         0.0, 
  1.0,       0.141421,    0.141421, 
  1.0,       0.0,         0.2, 
  1.0,      -0.141421,    0.141421, 
  1.0,      -0.2,         0.0, 
  1.0,      -0.141421,   -0.141421, 
  1.0,       0.0,        -0.2, 
  1.0,       0.141421,   -0.141421, 
  1.0,       0.066163,   -0.0302872, 
  1.0,      -0.0615154,  -0.0610739, 
  1.0,      -0.0306985,   0.0668017, 
  0.0,       0.2,         0.0, 
  0.0,       0.141421,   -0.141421, 
  0.0,       0.0,        -0.2, 
  0.0,      -0.141421,   -0.141421, 
  0.0,      -0.2,         0.0, 
  0.0,      -0.141421,    0.141421, 
  0.0,       0.0,         0.2, 
  0.0,       0.141421,    0.141421, 
  0.0,       0.0686748,   0.0255359, 
  0.0,       0.0,        -0.0865993, 
  0.0,      -0.0686749,   0.0255359, 
  0.8816,    0.185522,   -0.0747102, 
  0.642415,  0.187806,   -0.0687668, 
  0.627606, -0.0696445,  -0.187482, 
  0.876431,  0.0811908,  -0.182779, 
  0.881613,  0.186118,    0.0732131, 
  0.872048, -0.0699008,  -0.187387, 
  0.878318,  0.0844232,   0.181308, 
  0.845861, -0.0716063,   0.186742, 
  0.866503, -0.182493,   -0.0818307, 
  0.859402, -0.186751,    0.0715813, 
  0.131355,  0.18477,     0.0765501, 
  0.13317,   0.077694,    0.184292, 
  0.130862,  0.185301,   -0.0752567, 
  0.135181, -0.0749468,   0.185426, 
  0.130839,  0.0781729,  -0.18409, 
  0.131856, -0.0754694,  -0.185214, 
  0.135683, -0.184121,    0.0780993, 
  0.134207, -0.184959,   -0.0760928, 
  0.261923,  0.199982,    0.00264585, 
  0.263928,  0.144161,    0.138627, 
  0.268645,  0.00535339,  0.199928, 
  0.272346, -0.137646,    0.145098, 
  0.26108,   0.144683,   -0.138082, 
  0.260772,  0.00498797, -0.199938, 
  0.264253, -0.139152,   -0.143655, 
  0.270288, -0.199962,    0.00389323, 
  0.408181, -0.0730357,   0.186187, 
  0.411818, -0.184374,    0.0774991, 
  0.397539,  0.080738,    0.182979, 
  0.39192,   0.185619,    0.0744699, 
  0.392192,  0.184438,   -0.0773479, 
  0.389194,  0.0770141,  -0.184577, 
  0.38786,  -0.0747817,  -0.185493, 
  0.762413,  0.199986,   -0.0023425, 
  0.762987,  0.151152,   -0.13097, 
  0.741526,  0.0187858,  -0.199116, 
  0.746899, -0.128364,   -0.153371, 
  0.720076, -0.19917,    -0.0182053, 
  0.7628,    0.152219,    0.129728, 
  0.763882,  0.0434475,   0.195224, 
  0.399903, -0.1841,     -0.0781489, 
  0.506331, -0.00579066, -0.199916, 
  0.514514, -0.133894,   -0.148568, 
  0.526121,  0.135152,   -0.147424, 
  0.517967,  0.199953,   -0.0043215, 
  0.520585,  0.147847,    0.13469, 
  0.533956,  0.0124181,   0.199614, 
  0.558316, -0.136902,    0.145801, 
  0.549126, -0.199624,   -0.0122659, 
  0.657307,  0.117735,   -0.161674, 
  0.611189,  0.041829,   -0.195577, 
  0.631917, -0.164669,   -0.113508, 
  0.641444,  0.187001,    0.0709267, 
  0.720251, -0.155557,    0.125706, 
  0.647345,  0.0932963,   0.176906, 
  0.677484, -0.0430068,   0.195321, 
  0.635293, -0.188734,    0.0661777, 
  0.888023, -0.00868364, -0.00818647, 
  0.112146,  0.0,        -0.0118425, 
  0.676228,  0.0124197,  -0.0856487, 
  0.638436, -0.0639898,   0.0525795, 
  0.452586, -0.0410297,  -0.0704842, 
  0.762004, -0.0188614,   0.0693717, 
  0.463368,  0.0649048,   0.0262133, 
  0.473921, -0.0356443,   0.0388516, 
  0.557002,  0.0123705,  -0.0932599, 
  0.290986, -0.0200898,   0.00857934, 
  0.7038,    0.0856777,   0.0182744, 
  0.576134,  0.0436218,   0.0828782, 
  0.215187,  0.080855,   -0.0314946 };

  r8vec_copy ( 3 * vertices, vertex_coordinate_save, vertex_coordinate );

  i4vec_copy ( vertices, vertex_label_save, vertex_label );

  i4vec_copy ( 3 * triangles, triangle_vertex_save, triangle_vertex );

  i4vec_copy ( triangles, triangle_label_save, triangle_label );

  i4vec_copy ( 4 * tetrahedrons, tetrahedron_vertex_save, tetrahedron_vertex );

  for ( i = 0; i < tetrahedrons; i++ )
  {
    tetrahedron_label[i] = 1;
  }
  return;
# undef TETRAHEDRONS_SAVE
# undef TRIANGLES_SAVE
# undef VERTICES_SAVE
}
/******************************************************************************/

void cyl248_size ( int *dim, int *vertices, int *edges, int *triangles, 
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

/******************************************************************************/
/*
  Purpose:

    CYL248_SIZE defines the sizes for a 3D tetrahedral mesh.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Output, int *DIM, the spatial dimension, which should be 3.

    Output, int *VERTICES, the number of vertices.

    Output, int *EDGES, the number of edges (may be 0).

    Output, int *TRIANGLES, the number of triangles (may be 0).

    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).

    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  *dim = 3;
  *vertices = 92;
  *edges = 0;
  *triangles = 154;
  *quadrilaterals = 0;
  *tetrahedrons = 248;
  *hexahedrons = 0;

  return;
}
/******************************************************************************/

void data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons, 
  double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
  int edge_label[], int triangle_vertex[], int triangle_label[], 
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], 
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    DATA_PRINT prints the data of an ICE grid dataset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Input, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
    of each vertex.

    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Input, int EDGE_LABEL[EDGES], a label for each edge.

    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for 
    each quadrilateral.

    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int i;
  int j;

  printf ( "\n" );
  printf ( "  Vertices:\n" );
  printf ( "\n" );
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      printf ( "  %f", vertex_coordinate[i+j*3] );
    }
    printf ( "  (%d)\n", vertex_label[j] );
  }

  if ( 0 < edges )
  {
    printf ( "\n" );
    printf ( "  Edges:\n" );
    printf ( "\n" );
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        printf ( "  %d", edge_vertex[i+j*2] );
    }
    printf ( "  (%d)\n", edge_label[j] );
    }
  }

  if ( 0 < triangles )
  {
    printf ( "\n" );
    printf ( "  Triangles:\n" );
    printf ( "\n" );
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        printf ( "  %d", triangle_vertex[i+j*3] );
      }
      printf ( "  (%d)\n", triangle_label[j] );
    }
  }

  if ( 0 < quadrilaterals )
  {
    printf ( "\n" );
    printf ( "  Quadrilaterals:\n" );
    printf ( "\n" );
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %d", quadrilateral_vertex[i+j*4] );
      }
      printf ( "  (%d)\n", quadrilateral_label[j] );
    }
  }

  if ( 0 < tetrahedrons )
  {
    printf ( "\n" );
    printf ( "  Tetrahedrons:\n" );
    printf ( "\n" );
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %d", tetrahedron_vertex[i+j*4] );
      }
      printf ( "  (%d)\n", tetrahedron_label[j] );
    }
  }

  if ( 0 < hexahedrons )
  {
    printf ( "\n" );
    printf ( "  Hexahedrons:\n" );
    printf ( "\n" );
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        printf ( "  %d", hexahedron_vertex[i+j*8] );
      }
      printf ( "  (%d)\n", hexahedron_label[j] );
    }
  }
  return;
}
/******************************************************************************/

void data_read ( char *filename, int dim, int vertices, int edges, 
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
  int edge_label[], int triangle_vertex[], int triangle_label[], 
  int quadrilateral_vertex[], int quadrilateral_label[], 
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[], 
  int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    DATA_READ reads ICE data from a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be read.
    Ordinarily, the name should include the extension ".nc".

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
    of each vertex.

    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Output, int EDGE_LABEL[EDGES], a label for each edge.

    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for 
    each quadrilateral.

    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int i;
  int ncid;
  int dimid;
  int mode;
  int varid;
/*
  Open the file.
*/
  mode = NC_NOCLOBBER;
  nc_open ( filename, mode, &ncid );
/*
  Vertices.
*/
  nc_inq_varid ( ncid, "Vertex_Coordinate", &varid );
  nc_get_var_double ( ncid, varid, vertex_coordinate );

  nc_inq_varid ( ncid, "Vertex_Label", &varid );
  nc_get_var_int ( ncid, varid, vertex_label );
/*
  Edges.
*/
  if ( 0 < edges )
  {
    nc_inq_varid ( ncid, "Edge_Vertex", &varid );
    nc_get_var_int ( ncid, varid, edge_vertex );

    nc_inq_varid ( ncid, "Edge_Label", &varid );
    nc_get_var_int ( ncid, varid, edge_label );
  }
/*
  Triangles.
*/
  if ( 0 < triangles )
  {
    nc_inq_varid ( ncid, "Triangle_Vertex", &varid );
    nc_get_var_int ( ncid, varid, triangle_vertex );

    nc_inq_varid ( ncid, "Triangle_Label", &varid );
    nc_get_var_int ( ncid, varid, triangle_label );
  }
/*
  Quadrilaterals.
*/
  if ( 0 < quadrilaterals )
  {
    nc_inq_varid ( ncid, "Quadrilateral_Vertex", &varid );
    nc_get_var_int ( ncid, varid, quadrilateral_vertex );

    nc_inq_varid ( ncid, "Quadrilateral_Label", &varid );
    nc_get_var_int ( ncid, varid, quadrilateral_label );
  }
/*
  Tetrahedrons.
*/
  if ( 0 < tetrahedrons )
  {
    nc_inq_varid ( ncid, "Tetrahedron_Vertex", &varid );
    nc_get_var_int ( ncid, varid, tetrahedron_vertex );

    nc_inq_varid ( ncid, "Tetrahedron_Label", &varid );
    nc_get_var_int ( ncid, varid, tetrahedron_label );
  }
/*
  Hexahedrons.
*/
  if ( 0 < hexahedrons )
  {
    nc_inq_varid ( ncid, "Hexahedron_Vertex", &varid );
    nc_get_var_int ( ncid, varid, hexahedron_vertex );

    nc_inq_varid ( ncid, "Hexahedron_Label", &varid );
    nc_get_var_int ( ncid, varid, hexahedron_label );
  }
/*
  Close the file.
*/
  nc_close ( ncid );

  return;
}
/******************************************************************************/

void hexahexa_2x2x2_data ( int dim, int vertices, int edges, int triangles, 
  int quadrilaterals, int tetrahedrons, int hexahedrons, 
  double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
  int edge_label[], int triangle_vertex[], int triangle_label[], 
  int quadrilateral_vertex[], int quadrilateral_label[], 
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[], 
  int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    HEXAHEXA_2X2X2_DATA defines the data for a 3D hexahedral mesh.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
    of each vertex.

    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Output, int EDGE_LABEL[EDGES], a label for each edge.

    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for 
    each quadrilateral.

    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for 
    each tetrahedron.

    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
# define HEXAHEDRONS_SAVE 8
# define QUADRILATERALS_SAVE 24
# define VERTICES_SAVE 27

  int hexahedron_label_save[HEXAHEDRONS_SAVE] = {
    1, 1, 1, 1, 1, 1, 1, 1 };
  int hexahedron_vertex_save[8*HEXAHEDRONS_SAVE] = { 
      1,  2,  5,  4, 10, 11, 14, 13,
      2,  3,  6,  5, 11, 12, 15, 14, 
      4,  5,  8,  7, 13, 14, 17, 16,
      5,  6,  9,  8, 14, 15, 18, 17, 
     10, 11, 14, 13, 19, 20, 23, 22, 
     11, 12, 15, 14, 20, 21, 24, 23, 
     13, 14, 17, 16, 22, 23, 26, 25, 
     14, 15, 18, 17, 23, 24, 27, 26 };
  int quadrilateral_label_save[QUADRILATERALS_SAVE] = {
     1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 
     3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 
     6, 6, 6, 6 };
  int quadrilateral_vertex_save[4*QUADRILATERALS_SAVE] = {
     1,  4,  5,  2, 
     2,  5,  6,  3,
     4,  7,  8,  5,
     5,  8,  9,  6,
     1,  2, 11, 10,
     2,  3, 12, 11,
    10, 11, 20, 19, 
    11, 12, 21, 20, 
     3,  6, 15, 12, 
     6,  9, 18, 15, 
    12, 15, 24, 21, 
    15, 18, 27, 24, 
     7, 16, 17,  8, 
     8, 17, 18,  9, 
    16, 25, 26, 17, 
    17, 26, 27, 18, 
     1, 10, 13,  4, 
     4, 13, 16,  7, 
    10, 19, 22, 13, 
    13, 22, 25, 16, 
    19, 20, 23, 22, 
    20, 21, 24, 23, 
    22, 23, 26, 25, 
    23, 24, 27, 26 };
  int vertex_label_save[VERTICES_SAVE] = {
    5, 2, 3, 5, 1, 3, 5, 4, 4, 5, 
    2, 3, 5, 0, 3, 5, 4, 4, 6, 6, 
    6, 6, 6, 6, 6, 6, 6 };
  double vertex_coordinate_save[3*VERTICES_SAVE] = { 
    0.0, 0.0, 0.0, 
    0.5, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 0.5, 0.0, 
    0.5, 0.5, 0.0, 
    1.0, 0.5, 0.0, 
    0.0, 1.0, 0.0, 
    0.5, 1.0, 0.0, 
    1.0, 1.0, 0.0, 
    0.0, 0.0, 0.5, 
    0.5, 0.0, 0.5, 
    1.0, 0.0, 0.5, 
    0.0, 0.5, 0.5, 
    0.5, 0.5, 0.5, 
    1.0, 0.5, 0.5, 
    0.0, 1.0, 0.5, 
    0.5, 1.0, 0.5, 
    1.0, 1.0, 0.5, 
    0.0, 0.0, 1.0, 
    0.5, 0.0, 1.0, 
    1.0, 0.0, 1.0, 
    0.0, 0.5, 1.0, 
    0.5, 0.5, 1.0, 
    1.0, 0.5, 1.0, 
    0.0, 1.0, 1.0, 
    0.5, 1.0, 1.0,
    1.0, 1.0, 1.0 };

  r8vec_copy ( 3 * vertices, vertex_coordinate_save, vertex_coordinate );
  i4vec_copy ( vertices, vertex_label_save, vertex_label );
  i4vec_copy ( 4 * quadrilaterals, quadrilateral_vertex_save, quadrilateral_vertex );
  i4vec_copy ( quadrilaterals, quadrilateral_label_save, quadrilateral_label );
  i4vec_copy ( 8 * hexahedrons, hexahedron_vertex_save, hexahedron_vertex );
  i4vec_copy ( hexahedrons, hexahedron_label_save, hexahedron_label );

  return;
# undef HEXAHEDRONS_SAVE
# undef QUADRILATERALS_SAVE
# undef VERTICES_SAVE
}
/******************************************************************************/

void hexahexa_2x2x2_size ( int *dim, int *vertices, int *edges, int *triangles, 
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

/******************************************************************************/
/*
  Purpose:

    HEXAHEXA_2X2X2_SIZE defines the sizes for a 3D hexahedral mesh.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Output, int *DIM, the spatial dimension, which should be 3.

    Output, int *VERTICES, the number of vertices.

    Output, int *EDGES, the number of edges (may be 0).

    Output, int *TRIANGLES, the number of triangles (may be 0).

    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).

    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  *dim = 3;
  *vertices = 27;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 24;
  *tetrahedrons = 0;
  *hexahedrons = 8;

  return;
}
/******************************************************************************/

void i4vec_copy ( int n, int a1[], int a2[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_COPY copies an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 April 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, int A1[N], the vector to be copied.

    Input, int A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

void ice_write ( char *filename, int dim, int vertices, int edges, 
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons, 
  double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
  int edge_label[], int triangle_vertex[], int triangle_label[], 
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], 
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 November 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be created.
    Ordinarily, the name should include the extension ".nc".

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Input, real VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
    of each vertex.

    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Input, int EDGE_LABEL[EDGES], a label for each edge.

    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for 
    each quadrilateral.

    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int dim_dimension;
  int dim_edges;
  int dim_eight;
  int dim_four;
  int dim_hexahedrons;
  int dim_quadrilaterals;
  int dim_tetrahedrons;
  int dim_three;
  int dim_triangles;
  int dim_two;
  int dim_vertices;
  int dimids[2];
  int i;
  int ncid;
  int ndims;
  int mode;
  int var_edge_label;
  int var_edge_vertex;
  int var_hexahedron_label;
  int var_hexahedron_vertex;
  int var_quadrilateral_label;
  int var_quadrilateral_vertex;
  int var_tetrahedron_label;
  int var_tetrahedron_vertex;
  int var_triangle_label;
  int var_triangle_vertex;
  int var_vertex_coordinate;
  int var_vertex_label;
/*
  Create the file.  This automatically "opens" it as well.
*/
  mode = NC_CLOBBER;
  nc_create ( filename, mode, &ncid );
/*
  Put NETCDF into "define" mode.
*/
  nc_redef ( ncid );
/*
  Dimension information.

  If a dimension has length 0, it seems to be taken to be the unlimited
  dimension (not what you want) and then if you have two such dimensions,
  you get a ninny complaint that you have tried to define the unlimited dimension
  twice.  The fix requires the programmer not to write anything whose
  dimension is zero.
*/
  nc_def_dim ( ncid, "Dimension", dim, &dim_dimension );

  nc_def_dim ( ncid, "Vertices", vertices, &dim_vertices );

  if ( 0 < edges )
  {
    nc_def_dim ( ncid, "Edges", edges, &dim_edges );
  }

  if ( 0 < triangles )
  {
    nc_def_dim ( ncid, "Triangles", triangles, &dim_triangles );
  }

  if ( 0 < quadrilaterals )
  {
    nc_def_dim ( ncid, "Quadrilaterals", quadrilaterals, &dim_quadrilaterals );
  }

  if ( 0 < tetrahedrons )
  {
    nc_def_dim ( ncid, "Tetrahedrons", tetrahedrons, &dim_tetrahedrons );
  }

  if ( 0 < hexahedrons )
  {
    nc_def_dim ( ncid, "Hexahedrons", hexahedrons, &dim_hexahedrons );
  }

  nc_def_dim ( ncid, "Two", 2, &dim_two );
  nc_def_dim ( ncid, "Three", 3, &dim_three );
  nc_def_dim ( ncid, "Four", 4, &dim_four );
  nc_def_dim ( ncid, "Eight", 8, &dim_eight );
/*
  Define variables.
*/
  ndims = 2;
  dimids[0] = dim_three;
  dimids[1] = dim_vertices;
  nc_def_var ( ncid, "Vertex_Coordinate", NC_DOUBLE, ndims, dimids,
    &var_vertex_coordinate );

  ndims = 1;
  dimids[0] = dim_vertices;
  nc_def_var ( ncid, "Vertex_Label", NC_INT, ndims, dimids,
    &var_vertex_label );

  if ( 0 < edges )
  {
    ndims = 2;
    dimids[0] = dim_two;
    dimids[1] = dim_edges;
    nc_def_var ( ncid, "Edge_Vertex", NC_INT, ndims, dimids,
      &var_edge_vertex );

    ndims = 1;
    dimids[0] = dim_edges;
    nc_def_var ( ncid, "Edge_Label", NC_INT, ndims, dimids,
      &var_edge_label );
  }

  if ( 0 < triangles )
  {
    ndims = 2;
    dimids[0] = dim_three;
    dimids[1] = dim_triangles;
    nc_def_var ( ncid, "Triangle_Vertex", NC_INT, ndims, dimids,
      &var_triangle_vertex );

    ndims = 1;
    dimids[0] = dim_triangles;
    nc_def_var ( ncid, "Triangle_Label", NC_INT, ndims, dimids,
      &var_triangle_label );
  }

  if ( 0 < quadrilaterals )
  {
    ndims = 2;
    dimids[0] = dim_four;
    dimids[1] = dim_quadrilaterals;
    nc_def_var ( ncid, "Quadrilateral_Vertex", NC_INT, ndims, dimids,
      &var_quadrilateral_vertex );

    ndims = 1;
    dimids[0] = dim_quadrilaterals;
    nc_def_var ( ncid, "Quadrilateral_Label", NC_INT, ndims, dimids,
      &var_quadrilateral_label );
  }

  if ( 0 < tetrahedrons )
  {
    ndims = 2;
    dimids[0] = dim_four;
    dimids[1] = dim_tetrahedrons;
    nc_def_var ( ncid, "Tetrahedron_Vertex", NC_INT, ndims, dimids,
      &var_tetrahedron_vertex );

    ndims = 1;
    dimids[0] = dim_tetrahedrons;
    nc_def_var ( ncid, "Tetrahedron_Label", NC_INT, ndims, dimids,
      &var_tetrahedron_label );
  }

  if ( 0 < hexahedrons )
  {
    ndims = 2;
    dimids[0] = dim_eight;
    dimids[1] = dim_hexahedrons;
    nc_def_var ( ncid, "Hexahedron_Vertex", NC_INT, ndims, dimids,
      &var_hexahedron_vertex );

    ndims = 1;
    dimids[0] = dim_hexahedrons;
    nc_def_var ( ncid, "Hexahedron_Label", NC_INT, ndims, dimids,
      &var_hexahedron_label );
  }
/*
  Terminate the definition phase.
*/
  nc_enddef ( ncid );
/*
  Write the data.
*/
  nc_put_var_double ( ncid, var_vertex_coordinate, vertex_coordinate );
  nc_put_var_int ( ncid, var_vertex_label, vertex_label );
  if ( 0 < edges )
  {
    nc_put_var_int ( ncid, var_edge_vertex, edge_vertex );
    nc_put_var_int ( ncid, var_edge_label, edge_label );
  }
  if ( 0 < triangles )
  {
    nc_put_var_int ( ncid, var_triangle_vertex, triangle_vertex );
    nc_put_var_int ( ncid, var_triangle_label, triangle_label );
  }
  if ( 0 < quadrilaterals )
  {
    nc_put_var_int ( ncid, var_quadrilateral_vertex, quadrilateral_vertex );
    nc_put_var_int ( ncid, var_quadrilateral_label, quadrilateral_label );
  }
  if ( 0 < tetrahedrons )
  {
    nc_put_var_int ( ncid, var_tetrahedron_vertex, tetrahedron_vertex );
    nc_put_var_int ( ncid, var_tetrahedron_label, tetrahedron_label );
  }
  if ( 0 < hexahedrons )
  {
    nc_put_var_int ( ncid, var_hexahedron_vertex, hexahedron_vertex );
    nc_put_var_int ( ncid, var_hexahedron_label, hexahedron_label );
  }
/*
  Close the file.
*/
  nc_close ( ncid );

  return;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

void size_print ( int dim, int vertices, int edges, int triangles, 
  int quadrilaterals, int tetrahedrons, int hexahedrons )

/******************************************************************************/
/*
  Purpose:

    SIZE_PRINT prints the sizes of an ICE grid dataset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  printf ( "\n" );
  printf ( "  Number of dimensions = %d\n", dim );
  printf ( "  Number of vertices = %d\n", vertices );
  printf ( "  Number of edges = %d\n", edges );
  printf ( "  Number of triangles = %d\n", triangles );
  printf ( "  Number of quadrilaterals = %d\n", quadrilaterals );
  printf ( "  Number of tetrahedrons = %d\n", tetrahedrons );
  printf ( "  Number of hexahedrons = %d\n", hexahedrons );

  return;
}
/******************************************************************************/

void size_read ( char *filename, int *dim, int *vertices, int *edges, 
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

/******************************************************************************/
/*
  Purpose:

    SIZE_READ reads ICE sizes from a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be read.
    Ordinarily, the name should include the extension ".nc".

    Output, int *DIM, the spatial dimension, which should be 3.

    Output, int *VERTICES, the number of vertices.

    Output, int *EDGES, the number of edges (may be 0).

    Output, int *TRIANGLES, the number of triangles (may be 0).

    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).

    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  int ncid;
  int dimid;
  int mode;
  size_t value;
/*
  Initialize everything to nothing.
*/
  *dim = 0;
  *vertices = 0;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 0;
  *tetrahedrons = 0;
  *hexahedrons = 0;
/*
  Open the file.
*/
  mode = NC_NOCLOBBER;
  nc_open ( filename, mode, &ncid );
/*
  Get the dimension information.
*/
  if ( nc_inq_dimid ( ncid, "Dimension", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *dim = value;
  }
  if ( nc_inq_dimid ( ncid, "Vertices", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *vertices = value;
  }
  if ( nc_inq_dimid ( ncid, "Edges", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *edges = value;
  }
  if ( nc_inq_dimid ( ncid, "Triangles", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *triangles = value;
  }
  if ( nc_inq_dimid ( ncid, "Quadrilaterals", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *quadrilaterals = value;
  }
  if ( nc_inq_dimid ( ncid, "Tetrahedra", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *tetrahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Tetrahedrons", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *tetrahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Hexahedra", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *hexahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Hexahedrons", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *hexahedrons = value;
  }
/*
  Close the file.
*/
  nc_close ( ncid );

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
