# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "fem2d_bvp_quadratic.h"

/******************************************************************************/

double *fem2d_bvp_quadratic ( int nx, int ny, double a ( double x, double y ), 
  double c ( double x, double y ), double f ( double x, double y ), 
  double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    FEM2D_BVP_QUADRATIC solves a boundary value problem on a rectangle.

  Discussion:

    The procedure uses the finite element method, with piecewise quadratic 
    basis functions to solve a 2D boundary value problem over a rectangle

    The following differential equation is imposed inside the region:

      - d/dx a(x,y) du/dx - d/dy a(x,y) du/dy + c(x,y) * u(x,y) = f(x,y)

    where a(x,y), c(x,y), and f(x,y) are given functions.

    On the boundary, the solution is constrained to have the value 0.

    The finite element method will use a regular grid of NX nodes in X, and 
    NY nodes in Y.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of X and Y grid values.
    NX and NY must be odd and at least 3.

    Input, double A ( double X, double Y ), evaluates a(x,y);

    Input, double C ( double X, double Y ), evaluates c(x,y);

    Input, double F ( double X, double Y ), evaluates f(x,y);

    Input, double X[NX], Y[NY], the grid coordinates.

    Output, double FEM1D_BVP_QUADRATIC[NX*NY], the finite element coefficients, 
    which are also the value of the computed solution at the mesh points.
*/
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  double *amat;
  double aq;
  double *b;
  int cc;
  double cq;
  int e;
  int ex;
  int ex_num;
  int ey;
  int ey_num;
  double fq;
  int i;
  int ierror;
  int ii;
  int il;
  int il2;
  int il3;
  int j;
  int jj;
  int jl;
  int jl2;
  int jl3;
  int k;
  int mm;
  int mn;
  int n;
  int node[9];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double t;
  double *u;
  double v[9];
  double vx[9];
  double vy[9];
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xq;
  double xx[3];
  double yq;
  double yy[3];

  mn = nx * ny;

  amat = r8mat_zero_new ( mn, mn );
  b = r8vec_zero_new ( mn );

  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ex = 0; ex < ex_num; ex++ )
  {
    w = 2 * ex;
    cc = 2 * ex + 1;
    e = 2 * ex + 2;

    xx[0] = x[w];
    xx[1] = x[cc];
    xx[2] = x[e];

    for ( ey = 0; ey < ey_num; ey++ )
    {
      s = 2 * ey;
      mm = 2 * ey + 1;
      n = 2 * ey + 2;

      yy[0] = y[s];
      yy[1] = y[mm];
      yy[2] = y[n];
/*
  Node indices

  7  8  9   wn cn en
  4  5  6   wm cm em
  1  2  3   ws cs es
*/
      node[0] = ( 2 * ey     ) * nx + ex * 2 + 0;
      node[1] = ( 2 * ey     ) * nx + ex * 2 + 1;
      node[2] = ( 2 * ey     ) * nx + ex * 2 + 2;
      node[3] = ( 2 * ey + 1 ) * nx + ex * 2 + 0;
      node[4] = ( 2 * ey + 1 ) * nx + ex * 2 + 1;
      node[5] = ( 2 * ey + 1 ) * nx + ex * 2 + 2;
      node[6] = ( 2 * ey + 2 ) * nx + ex * 2 + 0;
      node[7] = ( 2 * ey + 2 ) * nx + ex * 2 + 1;
      node[8] = ( 2 * ey + 2 ) * nx + ex * 2 + 2;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xx[0]
             + ( 1.0 + abscissa[qx] ) * xx[2] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * yy[0]   
               + ( 1.0 + abscissa[qy] ) * yy[2] ) 
                 / 2.0;

          wq = weight[qx] * ( xx[2] - xx[0] ) / 2.0 
             * weight[qy] * ( yy[2] - yy[0] ) / 2.0;

          k = 0;

          for ( jl = 0; jl < 3; jl++ )
          {
            for ( il = 0; il < 3; il++ )
            {
              v[k] = 1.0;
              vx[k] = 0.0;
              vy[k] = 0.0;
              for ( il2 = 0; il2 < 3; il2++ )
              {
                if ( il2 != il )
                {
                  v[k] = v[k] * ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                  t = 1.0 / ( xx[il] - xx[il2] );
                  for ( il3 = 0; il3 < 3; il3++ )
                  {
                    if ( il3 != il && il3 != il2 )
                    {
                      t = t * ( xq - xx[il3] ) / ( xx[il] - xx[il3] );
                    }
                  }
                  for ( jl2 = 0; jl2 < 3; jl2++ )
                  {
                    if ( jl2 != jl )
                    {
                      t = t * ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                    }
                  }
                  vx[k] = vx[k] + t;
                }
              }

              for ( jl2 = 0; jl2 < 3; jl2++ )
              {
                if ( jl2 != jl )
                {
                  v[k] = v[k] * ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                  t = 1.0 / ( yy[jl] - yy[jl2] );
                  for ( il2 = 0; il2 < 3; il2++ )
                  {
                    if ( il2 != il )
                    {
                      t = t * ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                    }
                  }
                  for ( jl3 = 0; jl3 < 3; jl3++ )
                  {
                    if ( jl3 != jl && jl3 != jl2 )
                    {
                      t = t * ( yq - yy[jl3] ) / ( yy[jl] - yy[jl3] );
                    }
                  }
                  vy[k] = vy[k] + t;
                }
              }
              k = k + 1;
            }
          }

          aq = a ( xq, yq );
          cq = c ( xq, yq );
          fq = f ( xq, yq );

          for ( i = 0; i < 9; i++ )
          {
            ii = node[i];
            for ( j = 0; j < 9; j++ )
            {
              jj = node[j];
              amat[ii+jj*mn] = amat[ii+jj*mn] + wq * ( 
                  vx[i] * aq * vx[j] 
                + vy[i] * aq * vy[j] 
                + v[i]  * cq * v[j] );
            }
            b[ii] = b[ii] + wq * ( v[i] * fq );
          }
        }
      }
    }
  }
/*
  Where a node is on the boundary, 
  replace the finite element equation by a boundary condition.
*/
  k = 0;
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        for ( jj = 0; jj < mn; jj++ )
        {
          amat[k+jj*mn] = 0.0;
        }
        for ( ii = 0; ii < mn; ii++ )
        {
          amat[ii+k*mn] = 0.0;
        }
        amat[k+k*mn] = 1.0;
        b[k] = 0.0;
      }
      k = k + 1;
    }
  }
/*
  Solve the linear system.
*/
  u = r8mat_solve2 ( mn, amat, b, &ierror );

  free ( amat );
  free ( b );

  return u;
# undef QUAD_NUM
}
/******************************************************************************/

double fem2d_h1s_error_quadratic ( int nx, int ny, double x[], double y[],
  double u[], double exact_ux ( double x, double y ),
  double exact_uy ( double x, double y ) )

/******************************************************************************/
/*
  Purpose:

    FEM2D_H1S_ERROR_QUADRATIC: seminorm error of a finite element solution.

  Discussion:

    The finite element method has been used, over a rectangle,
    involving a grid of NX*NY nodes, with piecewise quadratic elements used 
    for the basis.

    The finite element solution U(x,y) has been computed, and formulas for the
    exact derivatives Vx(x,y) and Vy(x,y) are known.

    This function estimates the H1 seminorm of the error:

      H1S = sqrt ( integral ( x, y )   ( Ux(x,y) - Vx(x,y) )^2 
                                     + ( Uy(x,y) - Vy(x,y) )^2 dx dy )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of X and Y grid values.

    Input, double X[NX], Y[NY], the grid coordinates.

    Input, double U[NX*NY], the finite element coefficients.

    Input, function EXACT_UX(X,Y), EXACT_UY(X,Y) return the 
    value of the derivatives of the exact solution with respect to
    X and Y, respectively, at the point (X,Y).

    Output, double FEM2D_H1S_ERROR_QUADRATIC, the estimated seminorm of 
    the error.
*/
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  double aq;
  int cc;
  double cq;
  int e;
  int ex;
  int ex_num;
  double exq;
  double eyq;
  int ey;
  int ey_num;
  double fq;
  double h1s;
  int i;
  int ierror;
  int ii;
  int il;
  int il2;
  int il3;
  int j;
  int jj;
  int jl;
  int jl2;
  int jl3;
  int k;
  int mm;
  int mn;
  int n;
  int node[9];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double t;
  double uxq;
  double uyq;
  double vx[9];
  double vy[9];
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xq;
  double xx[3];
  double yq;
  double yy[3];

  mn = nx * ny;

  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ex = 0; ex < ex_num; ex++ )
  {
    w = 2 * ex;
    cc = 2 * ex + 1;
    e = 2 * ex + 2;

    xx[0] = x[w];
    xx[1] = x[cc];
    xx[2] = x[e];

    for ( ey = 0; ey < ey_num; ey++ )
    {
      s = 2 * ey;
      mm = 2 * ey + 1;
      n = 2 * ey + 2;

      yy[0] = y[s];
      yy[1] = y[mm];
      yy[2] = y[n];
/*
  Node indices

  7  8  9   wn cn en
  4  5  6   wm cm em
  1  2  3   ws cs es
*/
      node[0] = ( 2 * ey     ) * nx + ex * 2 + 0;
      node[1] = ( 2 * ey     ) * nx + ex * 2 + 1;
      node[2] = ( 2 * ey     ) * nx + ex * 2 + 2;
      node[3] = ( 2 * ey + 1 ) * nx + ex * 2 + 0;
      node[4] = ( 2 * ey + 1 ) * nx + ex * 2 + 1;
      node[5] = ( 2 * ey + 1 ) * nx + ex * 2 + 2;
      node[6] = ( 2 * ey + 2 ) * nx + ex * 2 + 0;
      node[7] = ( 2 * ey + 2 ) * nx + ex * 2 + 1;
      node[8] = ( 2 * ey + 2 ) * nx + ex * 2 + 2;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xx[0]
             + ( 1.0 + abscissa[qx] ) * xx[2] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * yy[0]   
               + ( 1.0 + abscissa[qy] ) * yy[2] ) 
                 / 2.0;

          wq = weight[qx] * ( xx[2] - xx[0] ) / 2.0 
             * weight[qy] * ( yy[2] - yy[0] ) / 2.0;

          uxq = 0.0;
          uyq = 0.0;

          k = 0;

          for ( jl = 0; jl < 3; jl++ )
          {
            for ( il = 0; il < 3; il++ )
            {
              vx[k] = 0.0;
              vy[k] = 0.0;
              for ( il2 = 0; il2 < 3; il2++ )
              {
                if ( il2 != il )
                {
                  t = 1.0 / ( xx[il] - xx[il2] );
                  for ( il3 = 0; il3 < 3; il3++ )
                  {
                    if ( il3 != il && il3 != il2 )
                    {
                      t = t * ( xq - xx[il3] ) / ( xx[il] - xx[il3] );
                    }
                  }
                  for ( jl2 = 0; jl2 < 3; jl2++ )
                  {
                    if ( jl2 != jl )
                    {
                      t = t * ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                    }
                  }
                  vx[k] = vx[k] + t;
                }
              }

              for ( jl2 = 0; jl2 < 3; jl2++ )
              {
                if ( jl2 != jl )
                {
                  t = 1.0 / ( yy[jl] - yy[jl2] );
                  for ( il2 = 0; il2 < 3; il2++ )
                  {
                    if ( il2 != il )
                    {
                      t = t * ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                    }
                  }
                  for ( jl3 = 0; jl3 < 3; jl3++ )
                  {
                    if ( jl3 != jl && jl3 != jl2 )
                    {
                      t = t * ( yq - yy[jl3] ) / ( yy[jl] - yy[jl3] );
                    }
                  }
                  vy[k] = vy[k] + t;
                }
              }
              uxq = uxq + u[node[k]] * vx[k];
              uyq = uyq + u[node[k]] * vy[k];
              k = k + 1;
            }
          }

          exq = exact_ux ( xq, yq );
          eyq = exact_uy ( xq, yq );

          h1s = h1s + wq * ( pow ( uxq - exq, 2 ) + pow ( uyq - eyq, 2 ) );
        }
      }
    }
  }

  h1s = sqrt ( h1s );

  return h1s;
# undef QUAD_NUM
}
/******************************************************************************/

double fem2d_l2_error_quadratic ( int nx, int ny, double x[], double y[],
  double u[], double exact ( double x, double y ) )

/******************************************************************************/
/*
  Purpose:

    FEM2D_L2_ERROR_QUADRATIC: L2 error norm of a finite element solution.

  Discussion:

    The finite element method has been used, over a rectangle,
    involving a grid of NX*NY nodes, with piecewise quadratic elements used 
    for the basis.

    The finite element coefficients have been computed, and a formula for the
    exact solution is known.

    This function estimates E2, the L2 norm of the error:

      E2 = Integral ( X, Y ) ( U(X,Y) - EXACT(X,Y) )^2 dX dY

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of X and Y grid values.

    Input, double X[NX], Y[NY], the grid coordinates.

    Input, double U[NX*NY], the finite element coefficients.

    Input, function EQ = EXACT(X,Y), returns the value of the exact
    solution at the point (X,Y).

    Output, double FEM2D_L2_ERROR_QUADRATIC, the estimated L2 norm of the error.
*/
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  double aq;
  int cc;
  double cq;
  int e;
  double e2;
  double eq;
  int ex;
  int ex_num;
  int ey;
  int ey_num;
  double fq;
  int i;
  int ierror;
  int ii;
  int il;
  int il2;
  int il3;
  int j;
  int jj;
  int jl;
  int jl2;
  int jl3;
  int k;
  int mm;
  int mn;
  int n;
  int node[9];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double t;
  double uq;
  double v[9];
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xq;
  double xx[3];
  double yq;
  double yy[3];

  mn = nx * ny;

  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ex = 0; ex < ex_num; ex++ )
  {
    w = 2 * ex;
    cc = 2 * ex + 1;
    e = 2 * ex + 2;

    xx[0] = x[w];
    xx[1] = x[cc];
    xx[2] = x[e];

    for ( ey = 0; ey < ey_num; ey++ )
    {
      s = 2 * ey;
      mm = 2 * ey + 1;
      n = 2 * ey + 2;

      yy[0] = y[s];
      yy[1] = y[mm];
      yy[2] = y[n];
/*
  Node indices

  7  8  9   wn cn en
  4  5  6   wm cm em
  1  2  3   ws cs es
*/
      node[0] = ( 2 * ey     ) * nx + ex * 2 + 0;
      node[1] = ( 2 * ey     ) * nx + ex * 2 + 1;
      node[2] = ( 2 * ey     ) * nx + ex * 2 + 2;
      node[3] = ( 2 * ey + 1 ) * nx + ex * 2 + 0;
      node[4] = ( 2 * ey + 1 ) * nx + ex * 2 + 1;
      node[5] = ( 2 * ey + 1 ) * nx + ex * 2 + 2;
      node[6] = ( 2 * ey + 2 ) * nx + ex * 2 + 0;
      node[7] = ( 2 * ey + 2 ) * nx + ex * 2 + 1;
      node[8] = ( 2 * ey + 2 ) * nx + ex * 2 + 2;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xx[0]
             + ( 1.0 + abscissa[qx] ) * xx[2] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * yy[0]   
               + ( 1.0 + abscissa[qy] ) * yy[2] ) 
                 / 2.0;

          wq = weight[qx] * ( xx[2] - xx[0] ) / 2.0 
             * weight[qy] * ( yy[2] - yy[0] ) / 2.0;

          uq = 0.0;
          k = 0;

          for ( jl = 0; jl < 3; jl++ )
          {
            for ( il = 0; il < 3; il++ )
            {
              v[k] = 1.0;
              for ( il2 = 0; il2 < 3; il2++ )
              {
                if ( il2 != il )
                {
                  v[k] = v[k] * ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                }
              }

              for ( jl2 = 0; jl2 < 3; jl2++ )
              {
                if ( jl2 != jl )
                {
                  v[k] = v[k] * ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                }
              }
              uq = uq + u[node[k]] * v[k];
              k = k + 1;
            }
          }

          eq = exact ( xq, yq );

          e2 = e2 + wq * pow ( uq - eq, 2 );
        }
      }
    }
  }

  e2 = sqrt ( e2 );

  return e2;
# undef QUAD_NUM
}

