# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <sys/wait.h>
# include <time.h>
# include <sys/types.h>

# define DIM 2
# define maxNeibNodes 32

typedef struct {
  int    ord;
  double crd[DIM];
  int    mark;
  int    n_neib_t;
  int    neib_t[maxNeibNodes];
  int    neib_t_o[maxNeibNodes];
  int    n_neib_v;
  int    neib_v[maxNeibNodes];
  int    rep;
} NODE;

typedef struct {
  int    ord;
  int    vtx[3];
  int    mark;
} TRIG;
/*
  int N_NODES, the number of nodes.

  int N_TRIGS, the number of elements.
*/
NODE *nodes;
TRIG *trigs;
int n_nodes;
int n_trigs;

char infile1[255];
char infile2[255];
char outfile1[255];
char outfile2[255];
/*
  Functions:
*/
int main ( int argc, char **argv );
void comp_circum ( double v0[], double v1[], double v2[], double circum[] );
void comp_neib ( void );
void read_mesh ( void );
void write_neib ( void );
void write_voro ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIG_TO_NEIB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 October 2010

  Author:

    Lili Ju

  Parameters:

    Commandline argument, char *INFILE1, the name of the node file.

    Commandline argument, char *INFILE2, the name of the element file.

    Commandline argument, char *OUTFILE1

    Commandline argument, char *OUTFILE2
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIG_TO_NEIB\n" );
  printf ( "  C version\n" );
  printf ( "  Construct neighbor and Voronoi information for a mesh.\n" );

  if ( ( argc != 4 ) && 
       ( argc != 5 ) ) 
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRIG_TO_NEIB - Fatal error!\n" );
    fprintf ( stderr, "  Command options error.\n");
    exit ( 1 );
  }
  strcpy ( infile1, argv[1] );
  strcpy ( infile2, argv[2] );
  strcpy ( outfile1, argv[3] );
/*
  Read the input files.
*/
  read_mesh ( );

  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n_nodes );
  printf ( "  Number of elements = %d\n", n_trigs );

  comp_neib ( );

  write_neib ( );

  if ( 5 <= argc ) 
  {
    strcpy ( outfile2, argv[4] );
    write_voro ( );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIG_TO_NEIB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void comp_circum ( double v0[], double v1[], double v2[], double circum[] )

/******************************************************************************/
/*
  Purpose:

    COMP_CIRCUM computes the circumcenter of a triangle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    Lili Ju

  Parameters:

    Input, double V0[DIM], V1[DIM], V2[DIM], the vertices.

    Output, double CIRCUM[DIM], the circumcenter.
*/
{
  double cnorm;
  double cu[DIM];
  double den;
  double e1[DIM];
  double e1l;
  double e2[DIM];
  double e2l;
  int i;
  double lg;

  for ( i = 0; i < DIM; i++ ) 
  {
    e1[i] = v1[i] - v0[i];
    e2[i] = v2[i] - v0[i];
  }
/*
  Formula for DIM == 3.
*/
  if ( DIM == 3 ) 
  {
    lg =  sqrt ( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] );
/*
  cu = e1 X e2
*/
    cu[0] = e1[1] * e2[2] - e1[2] * e2[1];
    cu[1] = e1[2] * e2[0] - e1[0] * e2[2];
    cu[2] = e1[0] * e2[1] - e1[1] * e2[0];

    cnorm = sqrt ( cu[0] * cu[0] + cu[1] * cu[1] + cu[2] * cu[2] );

    for ( i = 0; i < DIM; i++ ) 
    {
      circum[i] = lg * cu[i] / cnorm;
    }
  }
/*
  Formula for DIM = 2.
*/
  else
  {
    e1l = e1[0] * e1[0] + e1[1] * e1[1];
    e2l = e2[0] * e2[0] + e2[1] * e2[1];
    den = 0.5 / ( e1[0] * e2[1] - e1[1] * e2[0]);
    circum[0] = v0[0] + ( e2[1] * e1l - e1[1] * e2l ) * den;
    circum[1] = v0[1] + ( e1[0] * e2l - e2[0] * e1l ) * den;
  }
  return;
}
/******************************************************************************/

void comp_neib ( void )

/******************************************************************************/
/*
  Purpose:

    COMP_NEIB computes the neighbor information.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    Lili Ju

  Parameters:

*/
{
  int i;
  int j;
  int k;
  int lt;
  int m;
  int ni;
  int nj;
  int np;
  int nt;
  int ntemp;
  int v[3];

  for ( i = 0; i < n_nodes; i++ )
  {
    lt = (nodes+i)->n_neib_t;

    if ((nodes+i)->mark<=0)
    {
      (nodes+i)->neib_t_o[0] = (nodes+i)->neib_t[0];
      (nodes+i)->n_neib_v = lt;
    }
    else
    {
      for ( j = 0; j < lt; j++ ) 
      {
        nt = (nodes+i)->neib_t[j];
        for ( k = 0; k < 3; k++ ) 
        {
          (nodes+(trigs+nt)->vtx[k])->rep = 0;
        }
      }
      for ( j = 0; j < lt; j++ ) 
      {
        nt = (nodes+i)->neib_t[j];
        for ( k = 0; k < 3; k++ ) 
        {
          (nodes+(trigs+nt)->vtx[k])->rep++;
        }
      }
      for ( j = 0; j < lt; j++ )
      {
        nt = (nodes+i)->neib_t[j];
        if ((trigs+nt)->mark>=2) 
        {
          for ( k = 0; k < 3; k++ ) 
          {
            if (i==(trigs+nt)->vtx[k]) 
            {
              ni = (trigs+nt)->vtx[(k+1)%3];
              if ((nodes+ni)->rep==1)
              {
                (nodes+i)->neib_t_o[0] = nt;
              }
              break;
            }
          }
        }
      }
      (nodes+i)->n_neib_v = lt+1;
    }

    if ( lt == 1 ) 
    {
      nt = (nodes+i)->neib_t_o[0];
      for ( k = 0; k < 3; k++ )
      {
        if ( i == (trigs+nt)->vtx[k])
        {
          (nodes+i)->neib_v[j] = (trigs+nt)->vtx[(k+1)%3];
          np = (trigs+nt)->vtx[(k+1)%3];
          nj = (trigs+nt)->vtx[(k+2)%3];
        }
      }
    }

    for ( j = 0; j < lt - 1; j++ ) 
    {
      nt = (nodes+i)->neib_t_o[j];
      for ( k = 0; k < 3; k++ )
      {
        if ( i == (trigs+nt)->vtx[k]) 
        {
                 (nodes+i)->neib_v[j] = (trigs+nt)->vtx[(k+1)%3];
                 np = (trigs+nt)->vtx[(k+2)%3];
        }
      }
      for ( m = 0; m < lt; m++ ) 
      {
        ni = (nodes+i)->neib_t[m];
        if ( ni != nt )
        {
          for ( k = 0; k < 3; k++ )
          {
            if ((trigs+ni)->vtx[k]==np) 
            {
              (nodes+i)->neib_t_o[j+1] = ni;
              nj = (trigs+ni)->vtx[(k+1)%3];
              break;
            }
          }
        }
      }
    }

    if ((nodes+i)->mark<=0) 
    {
      (nodes+i)->neib_v[lt-1] = np;
    }
    else
    {
      (nodes+i)->neib_v[lt-1] = np;
      (nodes+i)->neib_v[lt] = nj;
    }
  }

  for ( i = 0; i < n_nodes; i++ ) 
  {
    lt = (nodes+i)->n_neib_t;
    for ( j = 0; j < lt; j++ ) 
    {
      (nodes+i)->neib_t[j] = (nodes+i)->neib_t_o[j];
    }
  }
  return;
}
/******************************************************************************/

void read_mesh ( void )

/******************************************************************************/
/*
  Purpose:

    READ_MESH reads the input files that define the mesh.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    Lili Ju

  Parameters:

*/
{
  FILE* fp;
  double ftemp;
  int i;
  int j;
  int k;
  int ntemp;
/*
  Read the node file.
*/
  fp = fopen ( infile1, "rt" );

  fscanf ( fp, "%d", &n_nodes );
  fscanf ( fp, "%d %d %d", &ntemp, &ntemp, &ntemp );
  nodes = ( NODE *) calloc ( n_nodes, sizeof ( NODE ) );

  for ( i = 0; i < n_nodes; i++ ) 
  {
    fscanf ( fp, "%d", &ntemp );
    (nodes+i)->ord = ntemp-1;
    for ( j = 0; j < DIM; j++ )
    {
      fscanf ( fp, "%lf", &ftemp );
      (nodes+i)->crd[j] = ftemp;
    }
    fscanf(fp,"%d",&ntemp);
    (nodes+i)->mark = ntemp;
    (nodes+i)->n_neib_t = 0;
    (nodes+i)->n_neib_v = 0;
  }
  fclose ( fp );
/*
  Get information from the element file.
*/
  fp = fopen(infile2,"rt");
  fscanf ( fp, "%d", &n_trigs );
  fscanf ( fp, "%d %d", &ntemp, &ntemp );
  trigs = (TRIG*)calloc(n_trigs,sizeof(TRIG));

  for ( i = 0; i < n_trigs; i++ ) 
  {
    fscanf(fp,"%d",&ntemp);
    (trigs+i)->ord = ntemp-1;
    for ( j = 0; j < 3; j++ ) 
    {
      fscanf(fp,"%d",&k);
      k = k - 1;
      (trigs+i)->vtx[j] = k;
      if ((nodes+k)->mark>0) 
      {
        (trigs+i)->mark += 1;
      }
      ntemp = (nodes+k)->n_neib_t;
      (nodes+k)->neib_t[ntemp] = (trigs+i)->ord;
      (nodes+k)->n_neib_t = ntemp+1;
    }
  }
  fclose(fp);

  return;
}
/******************************************************************************/

void write_neib ( void )

/******************************************************************************/
/*
  Purpose:

    WRITE_NEIB writes the neighbor information to an output file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    Lili Ju

  Parameters:

    None.
*/
{
  FILE* fp;
  int i;
  int j;
  int lt;
 
  fp = fopen ( outfile1, "wt" );

  for ( i = 0; i < n_nodes; i++ ) 
  {
    lt = (nodes+i)->n_neib_v;
    fprintf ( fp, "%10d ", lt );
    for ( j = 0; j < lt; j++ ) 
    {
      fprintf ( fp, "%10d ", (nodes+i)->neib_v[j] + 1 );
    }
    fprintf ( fp, "\n" );
  }
  fclose ( fp );

  return;
}
/******************************************************************************/

void write_voro ( void )

/******************************************************************************/
/*
  Purpose:

    WRITE_VORO writes the Voronoi information to an output file.

  Discussion:

    The first record gives the number of nodes.
    There follows one record for each node.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    Lili Ju

  Parameters:

    None
*/
{
  double circenter[DIM];
  FILE* fp;
  int  i;
  int j;
  int k;
  int lt;
  int lt1;
  int p1;
  int p2;

  fp = fopen ( outfile2, "wt" );

  fprintf ( fp, "%10d %4d\n", n_nodes, DIM );

  for ( i = 0; i < n_nodes; i++ ) 
  {
    lt = (nodes+i)->n_neib_v;

    if ( (nodes+i)->mark <= 0 ) 
    {
      fprintf ( fp, "%4d  ", lt );
      lt1 = lt;
    }
    else
    {
      fprintf ( fp, "%4d  ", lt + 2 );
      lt1 = lt - 1;
    }

    if ( (nodes+i)->mark > 0 ) 
    {
      p2 = (nodes+i)->neib_v[0];
      for ( j = 0; j < DIM; j++ )
      {
        fprintf ( fp, "%20.10f ", ((nodes+i)->crd[j]+(nodes+p2)->crd[j])/2);
      }
    }

    for ( j = 0; j < lt1; j++ )
    {
      p1 = (nodes+i)->neib_v[j];
      p2 = (nodes+i)->neib_v[(j+1)%lt];
      comp_circum ( (nodes+i)->crd, (nodes+p1)->crd,(nodes+p2)->crd, 
        circenter );
      for ( k = 0; k < DIM; k++ ) 
      {
        fprintf ( fp, "%20.10f ", circenter[k] );
      }
    }

    if ( (nodes+i)->mark > 0 ) 
    {
      p1 = (nodes+i)->neib_v[lt-1];
      for ( j = 0; j < DIM; j++ )
      {
        fprintf ( fp, "%20.10f ", ((nodes+i)->crd[j]+(nodes+p1)->crd[j])/2 );
      }
      for ( j = 0; j < DIM; j++ )
      {
        fprintf ( fp, "%20.10f ", (nodes+i)->crd[j] );
      }
    }
    fprintf ( fp, "\n" );
  }
  fclose ( fp );

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
