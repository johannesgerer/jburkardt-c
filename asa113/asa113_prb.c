# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa113.h"

int main ( void );
void test01 ( void );
double crswap ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int i, int j, int l, int m, int iswitch );
double crtran ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int i, int m, int l, int iswitch );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA113_PRB.

  Discussion:

    ASA136_PRB tests the ASA113 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA113_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA113 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA113_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tries out the ASA113 routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 November 2010

  Author:

    John Burkardt
*/
{
  double *a;
  int *c;
  double *c_center;
  int *c_size;
  int ci;
  double critvl;
  int i;
  int ifault;
  FILE *input;
  int j;
  int k = 5;
  int m = 100;
  int n = 2;
  int ntrans1;
  int ntrans2;
  double *wss;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );
  c = ( int * ) malloc ( m * sizeof ( int ) );
  c_center = ( double * ) malloc ( k * n * sizeof ( double ) );
  c_size = ( int * ) malloc ( k * sizeof ( int ) );
  wss = ( double * ) malloc ( k * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the ASA113 classification algorithm.\n" );
/*
  Read the data.
*/
  input = fopen ( "points_100.txt", "rt" );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fscanf ( input, "%lf", a+i+j*m );
    }
  }

  fclose ( input );
/*
  Print a few data values.
*/
  printf ( "\n" );
  printf ( "  First 5 data values:\n" );
  printf ( "\n" );

  for ( i = 0; i < 5; i++ )
  {
    printf ( "  %8d", i );
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %14f", a[i+j*m] );
    }
    printf ( "\n" );
  }
/*
  Assign points randomly to classes.
*/
  for ( i = 1; i <= m; i++ )
  {
    c[i-1] = ( i % k ) + 1;
  }
/*
  Define the critical value as the sum of the squares of the distances
  of the points to their cluster center.
*/
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    wss[i-1] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      wss[ci-1] = wss[ci-1] + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  critvl = 0.0;
  for ( i = 1; i <= k; i++ )
  {
    critvl = critvl + wss[i-1];
  }

  printf ( "\n" );
  printf ( "        Initial CRITVL = %f\n", critvl );
/*
  Compute the clusters.
*/
  ntrans1 = -1;
  ntrans2 = -1;

  for ( ; ; )
  {
    trnsfr ( a, c, c_size, m, k, n, &critvl, &ntrans1, &ifault );

    if ( ntrans1 == 0 && ntrans2 == 0 )
    {
      break;
    }

    printf ( "  After TRNSFR, CRITVL = %f\n", critvl );

    swap ( a, c, c_size, m, k, n, &critvl, &ntrans2, &ifault );

    if ( ntrans1 == 0 && ntrans2 == 0 )
    {
      break;
    }

    printf ( "    After SWAP, CRITVL = %f\n", critvl );
  }
/*
  Define the critical value as the sum of the squares of the distances
  of the points to their cluster center.
*/
  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    wss[i-1] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      wss[ci-1] = wss[ci-1] 
                + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }
  printf ( "\n" );
  printf ( "  Cluster  Population  Energy\n" );
  printf ( "\n" );

  for ( i = 1; i <= k; i++ )
  {
  printf ( "  %8d  %8d  %14f\n", i, c_size[i-1], wss[i-1] );
  }
  printf ( "\n" );
  printf ( "     Total  %8d  %14f\n", m, critvl );

  free ( a );
  free ( c );
  free ( c_center );
  free ( c_size );
  free ( wss );

  return;
}
/******************************************************************************/

double crswap ( double a[], int c[], int c_size[], int m, int k, 
  int n, double *critvl, int i1, int i2, int c1, int c2, int iswitch )

/******************************************************************************/
/*
  Purpose:

    CRSWAP determines the effect of swapping two objects.

  Discussion:

    This computation is very inefficient.  It is only set up so that we
    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
    ASA 136.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 February 2008

  Author:

    John Burkardt

  Reference:

    Colin Banfield, LC Bassill,
    Algorithm AS 113:
    A transfer for non-hierarchichal classification,
    Applied Statistics,
    Volume 26, Number 2, 1977, pages 206-210.

  Parameters:

    Input, double A(M,N), the data values.  There are M objects,
    each having spatial dimension N.

    Input, int C(M), the classification of each object.

    Input, int C_SIZE(K), the number of objects in each class.

    Input, int M, the number of objects.

    Input, int K, the number of classes.

    Input, int N, the number of spatial dimensions, or variates, 
    of the objects.

    Input, double *CRITVL, the current value of the criterion.

    Input, int I1, I2, the objects to be swapped.

    Input, int C1, C2, the current classes of objects I1 and I2.

    Input, int ISWITCH:
    1, indicates that I1 and I2 should be temporarily swapped, the 
       change in CRITVL should be computed, and then I1 and I2 restored.
    2, indicates that I1 and I2 will be swapped.

    Output, double CRSWAP, the change to CRITVL that would occur if I1 and
    I2 were swapped.  This is only computed for ISWITCH = 1.
*/
{ 
  double *c_center;
  int ci;
  double critvl_new;
  int i;
  double inc;
  int j;

  if ( iswitch == 2 )
  {
    inc = 0.0;
    return inc;
  }

  c_center = ( double * ) malloc ( k * n * sizeof ( double ) );
/*
  Move object I1 from class C1 to class C2.
  Move object I2 from class C2 to class C1.
*/
  c[i1-1] = c2;
  c[i2-1] = c1;
/*
  Define the critical value as the sum of the squares of the distances
  of the points to their cluster center.
*/
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  critvl_new= 0.0;

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      critvl_new = critvl_new 
                 + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  inc = critvl_new - *critvl;
/*
  Move object I1 from class C2 to class C1.
  Move object I2 from class C1 to class C2.
*/
  c[i1-1] = c1;
  c[i2-1] = c2;

  free ( c_center );

  return inc;
}
/******************************************************************************/

double crtran ( double a[], int c[], int c_size[], int m, int k, int n, 
  double *critvl, int i1, int c1, int c2, int iswitch )

/******************************************************************************/
/*
  Purpose:

    CRTRAN determines the effect of moving an object to another class.

  Discussion:

    This computation is very inefficient.  It is only set up so that we
    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
    ASA 136.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 February 2008

  Author:

    John Burkardt

  Reference:

    Colin Banfield, LC Bassill,
    Algorithm AS 113:
    A transfer for non-hierarchichal classification,
    Applied Statistics,
    Volume 26, Number 2, 1977, pages 206-210.

  Parameters:

    Input, double AL(M,N), the data values.  There are M objects,
    each having spatial dimension N.

    Input, int C(M), the classification of each object.

    Input, int C_SIZE(K), the number of objects in each class.

    Input, int M, the number of objects.

    Input, int K, the number of classes.

    Input, int N, the number of spatial dimensions, or variates, 
    of the objects.

    Input, double *CRITVL, the current value of the criterion.

    Input, int I1, the object to be transferred.

    Input, int C1, C2, the current class of object I1, and the
    class to which it may be transferred.

    Input, int ISWITCH:
    1, indicates that I1 should be temporarily transferred, the 
       change in CRITVL should be computed, and then I1 restored.
    2, indicates that I1 will be permanently transferred.

    Output, double CRTRAN, the change to CRITVL that would occur if I1 were
    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
*/
{
  double *c_center;
  int ci;
  double critvl_new;
  int i;
  double inc;
  int j;

  if ( iswitch == 2 )
  {
    inc = 0.0;
    return inc;
  }

  c_center = ( double * ) malloc ( k * n * sizeof ( double ) );
/*
  Move object I from class C1 to class C2.
*/
  c[i1-1] = c2;
  c_size[c1-1] = c_size[c1-1] - 1;
  c_size[c2-1] = c_size[c2-1] + 1;
/*
  Define the critical value as the sum of the squares of the distances
  of the points to their cluster center.
*/
  for ( i = 1; i <= k; i++ )
  {
    c_size[i-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    c_size[ci-1] = c_size[ci-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c_center[ci-1+(j-1)*k] = c_center[ci-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c_center[i-1+(j-1)*k] = c_center[i-1+(j-1)*k] / ( double ) ( c_size[i-1] );
    }
  }

  critvl_new = 0.0;

  for ( i = 1; i <= m; i++ )
  {
    ci = c[i-1];
    for ( j = 1; j <= n; j++ )
    {
      critvl_new = critvl_new 
                 + pow ( a[i-1+(j-1)*m] - c_center[ci-1+(j-1)*k], 2 );
    }
  }

  inc = critvl_new - *critvl;
/*
  Move object I1 from class C2 to class C1.
*/
  c[i1-1] = c1;
  c_size[c1-1] = c_size[c1-1] + 1;
  c_size[c2-1] = c_size[c2-1] - 1;

  free ( c_center );

  return inc;
}
