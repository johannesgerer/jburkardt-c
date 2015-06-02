# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem_basis.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM_BASIS_PRB.

  Discussion:

    FEM_BASIS_PRB tests the FEM_BASIS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM_BASIS_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the FEM_BASIS library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM_BASIS_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests FEM_BASIS_1D

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  double lij;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  FEM_BASIS_1D evaluates an arbitrary\n" );
  printf ( "  basis function over an interval.\n" );

  i1 = 2;
  j1 = 1;
  d = i1 + j1;
  x1 = r8_fraction ( i1, d );
  printf ( "\n" );
  printf ( "   I   J        X      L(I,J)(X)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %10g  %14g\n", i1, j1, x1, 1.0 );
  printf ( "\n" );
  for ( i2 = 0; i2 <= d; i2++ )
  {
    j2 = d - i2;
    x2 = r8_fraction ( i2, d );
    lij = fem_basis_1d ( i1, j1, x2 );
    printf ( "  %2d  %2d  %10g  %14g\n", i2, j2, x2, lij );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests FEM_BASIS_2D

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  int k1;
  int k2;
  double lijk;
  double x1;
  double x2;
  double y1;
  double y2;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  FEM_BASIS_2D evaluates an arbitrary triangular\n" );
  printf ( "  basis function.\n" );

  i1 = 1;
  j1 = 0;
  k1 = 2;
  d = i1 + j1 + k1;
  x1 = r8_fraction ( i1, d );
  y1 = r8_fraction ( j1, d );
  printf ( "\n" );
  printf ( "   I   J   K        X           Y      L(I,J,K)(X,Y)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %10g  %10g  %14g\n", i1, j1, k1, x1, y1, 1.0 );
  printf ( "\n" );
  for ( j2 = 0; j2 <= d; j2++ )
  {
    for ( i2 = 0; i2 <= d - j2; i2++ )
    {
      k2 = d - i2 - j2;
      x2 = r8_fraction ( i2, d );
      y2 = r8_fraction ( j2, d );
      lijk = fem_basis_2d ( i1, j1, k1, x2, y2 );
      printf ( "  %2d  %2d  %2d  %10g  %10g  %14g\n", i2, j2, k2, x2, y2, lijk );
    }
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests FEM_BASIS_3D

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  int k1;
  int k2;
  int l1;
  int l2;
  double lijkl;
  double x1;
  double x2;
  double y1;
  double y2;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  FEM_BASIS_3D evaluates an arbitrary tetrahedral\n" );
  printf ( "  basis function.\n" );

  i1 = 1;
  j1 = 0;
  k1 = 2;
  l1 = 1;
  d = i1 + j1 + k1 + l1;
  x1 = r8_fraction ( i1, d );
  y1 = r8_fraction ( j1, d );
  z1 = r8_fraction ( k1, d );
  printf ( "\n" );
  printf ( "   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %14g\n", 
    i1, j1, k1, l1, x1, y1, z1, 1.0 );
  printf ( "\n" );
  for ( k2 = 0; k2 <= d; k2++ )
  {
    for ( j2 = 0; j2 <= d - k2; j2++ )
    {
      for ( i2 = 0; i2 <= d - j2 - k2; i2++ )
      {
        l2 = d - i2 - j2 - k2;
        x2 = r8_fraction ( i2, d );
        y2 = r8_fraction ( j2, d );
        z2 = r8_fraction ( k2, d );
        lijkl = fem_basis_3d ( i1, j1, k1, l1, x2, y2, z2 );
        printf ( "  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %14g\n", 
          i2, j2, k2, l2, x2, y2, z2, lijkl );
      }
    }
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests FEM_BASIS_MD, repeating TEST01.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int i1[2];
  int i2[2];
  double l;
  int m = 1;
  int p1;
  double x1[1];
  double x2[1];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  FEM_BASIS_MD evaluates an arbitrary\n" );
  printf ( "  basis function over an M-dimensional simplex.\n" );

  i1[0] = 2;
  i1[1] = 1;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  printf ( "\n" );
  printf ( "   I   J        X      L(I,J)(X)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %10g  %14g\n", i1[0], i1[1], x1[0], 1.0 );
  printf ( "\n" );
  for ( p1 = 0; p1 <= d; p1++ )
  {
    i2[0] = p1;
    i2[1] = d - i2[0];
    for ( i = 0; i < m; i++ )
    {
      x2[i] = r8_fraction ( i2[i], d );
    }
    l = fem_basis_md ( m, i1, x2 );
    printf ( "  %2d  %2d  %10g  %14g\n", i2[0], i2[1], x2[0], l );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests FEM_BASIS_MD, repeating TEST02.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int i1[3];
  int i2[3];
  double l;
  int m = 2;
  int p1;
  int p2;
  double x1[2];
  double x2[2];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  FEM_BASIS_MD evaluates an arbitrary\n" );
  printf ( "  basis function over an M-dimensional simplex.\n" );

  i1[0] = 1;
  i1[1] = 0;
  i1[2] = 2;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  printf ( "\n" );
  printf ( "   I   J   K        X           Y      L(I,J,K)(X,Y)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %10g  %10g  %14g\n", i1[0], i1[1], i1[2], x1[0], x1[1], 1.0 );
  printf ( "\n" );
  for ( p2 = 0; p2 <= d; p2++ )
  {
    i2[1] = p2;
    for ( p1 = 0; p1 <= d - p2; p1++ )
    {
      i2[0] = p1;
      i2[2] = d - i2[0] - i2[1];
      for ( i = 0; i < m; i++ )
      {
        x2[i] = r8_fraction ( i2[i], d );
      }
      l = fem_basis_md ( m, i1, x2 );
      printf ( "  %2d  %2d  %2d  %10g  %10g  %14g\n", i2[0], i2[1], i2[2], x2[0], x2[1], l );
    }
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests FEM_BASIS_MD, repeating TEST03.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int i1[4];
  int i2[4];
  double l;
  int m = 3;
  int p1;
  int p2;
  int p3;
  double x1[3];
  double x2[3];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  FEM_BASIS_MD evaluates an arbitrary\n" );
  printf ( "  basis function over an M-dimensional simplex.\n" );

  i1[0] = 1;
  i1[1] = 0;
  i1[2] = 2;
  i1[3] = 1;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  printf ( "\n" );
  printf ( "   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %14g\n", 
    i1[0], i1[1], i1[2], i1[3], x1[0], x1[1], x1[2], 1.0 );
  printf ( "\n" );
  for ( p3 = 0; p3 <= d; p3++ )
  {
    i2[2] = p3;
    for ( p2 = 0; p2 <= d - p3; p2++ )
    {
      i2[1] = p2;
      for ( p1 = 0; p1 <= d - p3 - p2; p1++ )
      {
        i2[0] = p1;
        i2[3] = d - i2[0] - i2[1] - i2[2];
        for ( i = 0; i < m; i++ )
        {
          x2[i] = r8_fraction ( i2[i], d );
        }
        l = fem_basis_md ( m, i1, x2 );
        printf ( "  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %14g\n", 
          i2[0], i2[1], i2[2], i2[3], x2[0], x2[1], x2[2], l );
      }
    }
  }
  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests FEM_BASIS_PRISM_TRIANGLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  double b;
  int di;
  int dj;
  int i1[3] = { 2, 0, 0 };
  int i2[3];
  int i_0;
  int i_1;
  int j1[2] = { 1, 1 };
  int j2[2];
  int j_0;
  double xyz1[3];
  double xyz2[3];

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary\n" );
  printf ( "  basis function over a right triangular prism.\n" );
  printf ( "\n" );
  printf ( "  Here, we generate basis functions which can be\n" );
  printf ( "  up to degree 2 in X and Y, and up to degree 2 in Z.\n" );
  printf ( "\n" );
  printf ( "  Choose a node N1, define the basis function associated\n" );
  printf ( "  with that node, and then evaluate it at all other nodes.\n" );

  di = i4vec_sum ( 3, i1 );
  xyz1[0] = r8_fraction ( i1[0], di );
  xyz1[1] = r8_fraction ( i1[1], di );

  dj = i4vec_sum ( 2, j1 );
  xyz1[2] = r8_fraction ( j1[0], dj );

  printf ( "\n" );
  printf ( "  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %10g\n",
    i1[0], i1[1], i1[2], j1[0], j1[1], xyz1[0], xyz1[1], xyz1[2], 1.0 );

  printf ( "\n" );

  for ( i_0 = 0; i_0 <= di; i_0++ )
  {
    i2[0] = i_0;
    xyz2[0] = r8_fraction ( i2[0], di );
    for ( i_1 = 0; i_1 <= di - i2[0]; i_1++ )
    {
      i2[1] = i_1;
      xyz2[1] = r8_fraction ( i2[1], di );
      i2[2] = di - i2[0] - i2[1];
      for( j_0 = 0; j_0 <= dj; j_0++ )
      {
        j2[0] = j_0;
        j2[1] = dj - j2[0];
        xyz2[2] = r8_fraction ( j2[0], dj );

        b = fem_basis_prism_triangle ( i1, j1, xyz2 );

        printf ( "  %2d  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %10g\n",
          i2[0], i2[1], i2[2], j2[0], j2[1], xyz2[0], xyz2[1], xyz2[2], b );
     }
    }
  }
  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests FEM_BASIS_PRISM_TRIANGLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt
*/
{
  double b;
  int di;
  int dj;
  int i1[3] = { 2, 0, 1 };
  int i2[3];
  int i_0;
  int i_1;
  int j1[2] = { 1, 0 };
  int j2[2];
  int j_0;
  double xyz1[3];
  double xyz2[3];

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary\n" );
  printf ( "  basis function over a right triangular prism.\n" );
  printf ( "\n" );
  printf ( "  Here, we generate basis functions which can be\n" );
  printf ( "  up to degree 3 in X and Y, and up to degree 1 in Z.\n" );
  printf ( "\n" );
  printf ( "  Choose a node N1, define the basis function associated\n" );
  printf ( "  with that node, and then evaluate it at all other nodes.\n" );

  di = i4vec_sum ( 3, i1 );
  xyz1[0] = r8_fraction ( i1[0], di );
  xyz1[1] = r8_fraction ( i1[1], di );

  dj = i4vec_sum ( 2, j1 );
  xyz1[2] = r8_fraction ( j1[0], dj );

  printf ( "\n" );
  printf ( "  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)\n" );
  printf ( "\n" );
  printf ( "  %2d  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %10g\n",
    i1[0], i1[1], i1[2], j1[0], j1[1], xyz1[0], xyz1[1], xyz1[2], 1.0 );
  printf ( "\n" );

  for ( i_0 = 0; i_0 <= di; i_0++ )
  {
    i2[0] = i_0;
    xyz2[0] = r8_fraction ( i2[0], di );
    for ( i_1 = 0; i_1 <= di - i2[0]; i_1++ )
    {
      i2[1] = i_1;
      xyz2[1] = r8_fraction ( i2[1], di );
      i2[2] = di - i2[0] - i2[1];
      for( j_0 = 0; j_0 <= dj; j_0++ )
      {
        j2[0] = j_0;
        j2[1] = dj - j2[0];
        xyz2[2] = r8_fraction ( j2[0], dj );

        b = fem_basis_prism_triangle ( i1, j1, xyz2 );

        printf ( "  %2d  %2d  %2d  %2d  %2d  %10g  %10g  %10g  %10g\n",
          i2[0], i2[1], i2[2], j2[0], j2[1], xyz2[0], xyz2[1], xyz2[2], b );
      }
    }
  }
  return;
}
