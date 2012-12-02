/*

Create Platonic solids (tetrahedron, octahedron, cube, icosahedron,
dodecahedron.)

Greg Turk

-----------------------------------------------------------------------

Copyright (c) 1998 Georgia Institute of Technology.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <ply_io.h>


/* user's vertex and face definitions for a polygonal object */

typedef struct Vertex {
  float x,y,z;
} Vertex;

typedef struct Face {
  unsigned char nverts;    /* number of vertex indices in list */
  int *verts;              /* vertex index list */
} Face;

char *elem_names[] = { /* list of the kinds of elements in the user's object */
  "vertex", "face"
};

PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", Int32, Int32, offsetof(Face,verts),
   1, Uint8, Uint8, offsetof(Face,nverts)},
};


/*** the PLY object ***/

static int nverts,nfaces;
static Vertex vlist[32];
static Face flist[32];
static int nelems = 2;
static char **elist;
static int num_comments = 0;
static char **comments;
static int num_obj_info = 0;
static char **obj_info;
static int file_type;

static float xtrans = 0;
static float ytrans = 0;
static float ztrans = 0;

static float xscale = 1;
static float yscale = 1;
static float zscale = 1;

#define TETRA  1
#define OCTA   2
#define CUBE   3
#define ICOSA  4
#define DODECA 5


/******************************************************************************
Main program.
******************************************************************************/

main(int argc, char *argv[])
{
  int i,j;
  char *s;
  char *progname;
  int solid = TETRA;

  progname = argv[0];

  while (--argc > 0 && (*++argv)[0]=='-') {
    for (s = argv[0]+1; *s; s++)
      switch (*s) {
        case 's':
          xscale = atof (*++argv);
          yscale = atof (*++argv);
          zscale = atof (*++argv);
          argc -= 3;
          break;
        case 't':
          solid = TETRA;
          break;
        case 'o':
          solid = OCTA;
          break;
        case 'c':
          solid = CUBE;
          break;
        case 'i':
          solid = ICOSA;
          break;
        case 'd':
          solid = DODECA;
          break;
        default:
          usage (progname);
          exit (-1);
          break;
      }
  }

  switch (solid) {
    case TETRA:
      tetra();
      break;
    case OCTA:
      octa();
      break;
    case CUBE:
      cube();
      break;
    case ICOSA:
      icosahedron();
      break;
    case DODECA:
      dodecahedron();
      break;
    default:
      exit (-1);
      break;
  }

  write_file();
}


/******************************************************************************
Print out usage information.
******************************************************************************/

usage(char *progname)
{
  fprintf (stderr, "usage: %s [flags] >out.ply\n", progname);
  fprintf (stderr, "         -t { tetrahedron }\n");
  fprintf (stderr, "         -o { octahedron }\n");
  fprintf (stderr, "         -c { cube }\n");
  fprintf (stderr, "         -i { icosahedron }\n");
  fprintf (stderr, "         -d { dodecahedron }\n");
}


/******************************************************************************
Create a tetrahedron.
******************************************************************************/

tetra()
{
  static float v[4][3] = {
    -1, -1, -1,
     1,  1, -1,
     1, -1,  1,
    -1,  1,  1,
  };
  static int f[4][3] = {
    1, 2, 3,
    1, 0, 2,
    3, 2, 0,
    0, 1, 3,
  };
  int i;

  nverts = 4;
  nfaces = 4;

  for (i = 0; i < nverts; i++) {
    vlist[i].x = v[i][0];
    vlist[i].y = v[i][1];
    vlist[i].z = v[i][2];
  }

  for (i = 0; i < nfaces; i++) {
    flist[i].nverts = 3;
    flist[i].verts = (int *) malloc (sizeof(int) * 3);
    flist[i].verts[0] = f[i][0];
    flist[i].verts[1] = f[i][1];
    flist[i].verts[2] = f[i][2];
  }
}


/******************************************************************************
Create an octahedron.
******************************************************************************/

octa()
{
  static float v[6][3] = {
     1,  0,  0,
     0, -1,  0,
    -1,  0,  0,
     0,  1,  0,
     0,  0,  1,
     0,  0, -1,
  };
  static int f[8][3] = {
    4, 0, 1,
    4, 1, 2,
    4, 2, 3,
    4, 3, 0,
    5, 1, 0,
    5, 2, 1,
    5, 3, 2,
    5, 0, 3,
  };
  int i;

  nverts = 6;
  nfaces = 8;

  for (i = 0; i < nverts; i++) {
    vlist[i].x = v[i][0];
    vlist[i].y = v[i][1];
    vlist[i].z = v[i][2];
  }

  for (i = 0; i < nfaces; i++) {
    flist[i].nverts = 3;
    flist[i].verts = (int *) malloc (sizeof(int) * 3);
    flist[i].verts[0] = f[i][0];
    flist[i].verts[1] = f[i][1];
    flist[i].verts[2] = f[i][2];
  }
}


/******************************************************************************
Create a cube.
******************************************************************************/

cube()
{
  static float v[8][3] = {
    -1, -1, -1,
     1, -1, -1,
     1,  1, -1,
    -1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
     1,  1,  1,
    -1,  1,  1,
  };
  static int f[6][4] = {
    0, 1, 2, 3,
    5, 4, 7, 6,
    6, 2, 1, 5,
    3, 7, 4, 0,
    7, 3, 2, 6,
    5, 1, 0, 4,
  };
  int i;

  nverts = 8;
  nfaces = 6;

  for (i = 0; i < nverts; i++) {
    vlist[i].x = v[i][0];
    vlist[i].y = v[i][1];
    vlist[i].z = v[i][2];
  }

  for (i = 0; i < nfaces; i++) {
    flist[i].nverts = 4;
    flist[i].verts = (int *) malloc (sizeof(int) * 4);
    flist[i].verts[0] = f[i][0];
    flist[i].verts[1] = f[i][1];
    flist[i].verts[2] = f[i][2];
    flist[i].verts[3] = f[i][3];
  }
}


/******************************************************************************
Create an icosahedron.
******************************************************************************/

icosahedron()
{
  static float v[12][3] = {
    0, -2,  1,
    1,  0,  2,
    1,  0, -2,
    -1,  0, -2,
    -1,  0,  2,
    -2,  1,  0,
    2,  1,  0,
    2, -1,  0,
    -2, -1,  0,
    0, -2, -1,
    0,  2, -1,
    0,  2,  1,
  };
  static int f[20][3] = {
    6, 2, 1, 
    2, 7, 1, 
    5, 4, 3, 
    8, 3, 4, 
    11, 5, 6, 
    10, 6, 5, 
    2, 10, 9, 
    3, 9, 10, 
    9, 8, 7, 
    0, 7, 8, 
    1, 0, 11, 
    4, 11, 0, 
    10, 2, 6, 
    11, 6, 1, 
    10, 5, 3, 
    11, 4, 5, 
    9, 7, 2, 
    0, 1, 7, 
    8, 9, 3, 
    0, 8, 4, 
  };
  int i;
  float *p;
  float a,b;
  float len;

  a = 0.5 * (1 + sqrt(5));
  b = 1;

  len = sqrt(a*a + b*b);
  a /= len;
  b /= len;

  for (i = 0, p = (float *) v; i < 12 * 3; i++, p++) {
    if (fabs(*p) == 1)
      *p *= a;
    if (fabs(*p) == 2)
      *p *= b / 2.0;
  }

  nverts = 12;
  nfaces = 20;

  for (i = 0; i < nverts; i++) {
    vlist[i].x = v[i][0];
    vlist[i].y = v[i][1];
    vlist[i].z = v[i][2];
  }

  for (i = 0; i < nfaces; i++) {
    flist[i].nverts = 3;
    flist[i].verts = (int *) malloc (sizeof(int) * 3);
    flist[i].verts[0] = f[i][0];
    flist[i].verts[1] = f[i][1];
    flist[i].verts[2] = f[i][2];
  }
}


/******************************************************************************
Create a dodecahedron.
******************************************************************************/

dodecahedron()
{
  static float v[20][3] = {
    -2, -2, 2,
    3, 1, 0,
    3, -1, 0,
    -3, 1, 0,
    -3, -1, 0,
    0, 3, 1,
    0, 3, -1,
    1, 0, -3,
    -1, 0, -3,
    0, -3, -1,
    0, -3, 1,
    1, 0, 3,
    -1, 0, 3,
    2, 2, -2,
    2, 2, 2,
    -2, 2, -2,
    -2, 2, 2,
    2, -2, -2,
    2, -2, 2,
    -2, -2, -2,
  };
  static int f[12][5] = {
    1, 2, 18, 11, 14, 
    1, 13, 7, 17, 2, 
    3, 4, 19, 8, 15, 
    3, 16, 12, 0, 4, 
    3, 15, 6, 5, 16, 
    1, 14, 5, 6, 13, 
    2, 17, 9, 10, 18, 
    4, 0, 10, 9, 19, 
    7, 8, 19, 9, 17, 
    6, 15, 8, 7, 13, 
    5, 14, 11, 12, 16, 
    10, 0, 12, 11, 18, 
  };
  int i;
  float *p;
  float a,b,c;
  float len;

  a = 1 + 0.5 * (1 + sqrt(5));
  b = 1;
  c = 0.5 * (1 + sqrt(5));

  len = sqrt(a*a + b*b);
  a /= len;
  b /= len;
  c /= len;

  for (i = 0, p = (float *) v; i < 20 * 3; i++, p++) {
    if (fabs(*p) == 1)
      *p *= b;
    if (fabs(*p) == 2)
      *p *= c / 2.0;
    if (fabs(*p) == 3)
      *p *= a / 3.0;
  }

  nverts = 20;
  nfaces = 12;

  for (i = 0; i < nverts; i++) {
    vlist[i].x = v[i][0];
    vlist[i].y = v[i][1];
    vlist[i].z = v[i][2];
  }

  for (i = 0; i < nfaces; i++) {
    flist[i].nverts = 5;
    flist[i].verts = (int *) malloc (sizeof(int) * 5);
    flist[i].verts[0] = f[i][0];
    flist[i].verts[1] = f[i][1];
    flist[i].verts[2] = f[i][2];
    flist[i].verts[3] = f[i][3];
    flist[i].verts[4] = f[i][4];
  }
}


/******************************************************************************
Write out the PLY file to standard out.
******************************************************************************/

write_file()
{
  int i;
  PlyFile *ply;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/

  ply = write_ply (stdout, nelems, elem_names, PLY_ASCII);

  /* describe what properties go into the vertex elements */

  describe_element_ply (ply, "vertex", nverts);
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);

  describe_element_ply (ply, "face", nfaces);
  describe_property_ply (ply, &face_props[0]);

  append_comment_ply (ply, "created by platoply");

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, "vertex");
  for (i = 0; i < nverts; i++)
    put_element_ply (ply, (void *) &vlist[i]);

  /* set up and write the face elements */
  put_element_setup_ply (ply, "face");
  for (i = 0; i < nfaces; i++)
    put_element_ply (ply, (void *) &flist[i]);

  close_ply (ply);
  free_ply (ply);
}

