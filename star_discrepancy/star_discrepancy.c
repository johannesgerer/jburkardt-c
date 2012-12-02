# include <stdlib.h>
# include <stdio.h>
# include <math.h>

typedef struct som sommet;

struct som
{
  int id;
  sommet *pt;
};
static int s;
static int smin;
static int taille;
static int num;
static int den;
static int num2;
static int subtotala = 0;
static int subtotalb = 0;
static int *lexsizealpha;
static int *maxsizealpha;
static int *lexsizebeta;
static int *maxsizebeta;
static double p;
static double borne = 0.0;
static double borneinf = 0.0;
static double **points;
static sommet *superarbre;
static sommet **lexalpha;
static sommet **lexbeta;
static FILE *fichier;
/*
  Declaration of routines.
*/
int main ( int argc, char *argv[] );
void decomposition ( double *alpha, double *beta, int min, double value );
int explore ( sommet *liste, double *pave, int dim );
int fastexplore ( double *pave, int range, int *maxsize, int *lexsize,
  sommet **lex, int *subtotal );
static void fileformat ( void );
void freetree ( sommet *noeud );
void initlex ( void );
void insertlex ( sommet *noeud, int range, int *maxsize, int *lexsize,
  sommet **lex );
double lowbound ( int npoints, double volume, double *pave );
static void memory ( void );
void quicksort ( sommet *liste, int dim, int l, int r );
void readfile ( char *filename );
sommet *subtree ( sommet *liste, int min, int next, int dim );
void supertree ( void );
void traiter ( double *outputalpha, double *outputbeta, int range );
static void usage ( char *nom );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:
 
    MAIN is the main program for the star discrepancy bound computation.

  Modified:

    30 September 2003

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int j;
  double *oalpha;
  double *obeta;
  
  initparameters(argc,argv);
  readfile(argv[4]);

  printf("x={\n");
  for (i=0;i<taille;i++)
    {
      printf(" (");
      for (j=0;j<s;j++)
  	printf(" %f",points[i][j]);
      printf(" )\n");
    }
  printf("}\n\n");
  supertree();
  initlex();
  
  oalpha = (double *) calloc((unsigned) s,sizeof(double));
  obeta  = (double *) calloc((unsigned) s,sizeof(double));
  for (i=0;i<s;i++)
    obeta[i] = 1.0;

  decomposition(oalpha,obeta,0,1.0);

  printf("s=%d, epsilon=%f, n=%d\n",s,p,taille);
  printf("D_n^*(x) between %.10f and %.10f\n",borneinf,borne);

  return 0;
}
/******************************************************************************/

void decomposition ( double *alpha, double *beta, int min, double value )

/******************************************************************************/
/*
  Purpose:
 
    DECOMPOSITION carries out the decomposition of a subinterval.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  double pbetaminp = 1.0;
  double palpha = 1.0;
  double pbeta;
  double delta;
  double *subalpha;
  double *subbeta;
  double *gamma;
  subalpha = (double *) calloc((unsigned) s,sizeof(double));
  subbeta  = (double *) calloc((unsigned) s+1,sizeof(double));
  gamma  = (double *) calloc((unsigned) s+1,sizeof(double));

  for (i=min;i<s;i++)
    pbetaminp *= beta[i];
  pbeta = pbetaminp;
  for (i=0;i<min;i++)
    {
      pbetaminp *= beta[i];
      palpha *= alpha[i];
    }
  pbetaminp -= p;
  delta = pow(pbetaminp/(pbeta*palpha),1.0/(s-min));

  for (i=0;i<min;i++)
    {
      gamma[i] = alpha[i];
      subalpha[i] = gamma[i];
      subbeta[i] = beta[i];
    }
  for (i=min;i<s;i++)
    {
      gamma[i] = delta*beta[i];
      subalpha[i] = alpha[i];
      subbeta[i] = beta[i];
    }
  subbeta[min] = gamma[min];

  value *= delta;
  if (value>p)
    for (i=min;i<s;i++)
      {
	decomposition(subalpha,subbeta,i,value);
	subalpha[i]  = gamma[i];
	subbeta[i]   = beta[i];
	subbeta[i+1] = gamma[i+1];
      }
  else
    for (i=min;i<s;i++)
      {
	traiter(subalpha,subbeta,(i==0)?0:i-1);
	subalpha[i]  = gamma[i];
	subbeta[i]   = beta[i];
	subbeta[i+1] = gamma[i+1];
      }

  traiter(gamma,beta,smin);

  free(gamma);
  free(subalpha);
  free(subbeta);

  return;
}
/******************************************************************************/

int explore ( sommet *liste, double *pave, int dim )

/******************************************************************************/
/*
  Purpose:
 
    EXPLORE ???

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int min = 1;
  int max;
  int next;
  int total;
  if (pave[dim]<=points[liste[1].id][dim])
    return 0;
  if (liste[0].id==1)
    {
      total = 1;
      next = liste[1].id; 
      for (i=dim;i<s;i++)
	if (points[next][i]>=pave[i])
	  {
	    total = 0;
	    break;
	  }
    }
  else
    {
      total = 0;
      max = liste[0].id;
      if (dim==smin)
	{
	  if (pave[dim]>points[liste[max].id][dim])
	    total = max;
	  else
	    while (max>=min)
	      {
		next = (min+max+1)/2;
		if (points[liste[next].id][dim]<pave[dim])
		  {		
		    total += next-min+1;
		    min = next+1;
		  }
		else
		  max = next-1;
	      }
	}
      else
	{
	  while (max>=min)
	    {
	      next = ((1+min)*num2+max*num)/den;
	      if (points[liste[next].id][dim]<pave[dim])
		{
		  if (liste[next].pt==NULL)
		    liste[next].pt = subtree(liste,min,next,dim+1);
		  total += explore(liste[next].pt,pave,dim+1);
		  min = next+1;
		}
	      else
		max = next-1;
	    }
	}
    }
  return total;
}
/******************************************************************************/

int fastexplore ( double *pave, int range, int *maxsize, int *lexsize,
  sommet **lex, int *subtotal )

/******************************************************************************/
/*
  Purpose:
 
    FASTEXPLORE ???

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int j;
  int i;
  int min;
  int max;
  int next;
  int start;
  int size = lexsize[range];
  int right;
  int total = 0;
  double seuil = pave[range];
  sommet refnoeud;
  sommet *noeud;
  if (range==smin)
    {
      for (i=size-1;i>=0;i--)
	{
	  refnoeud = lex[range][i];
	  noeud = refnoeud.pt;
	  min = refnoeud.id;
	  max = noeud[0].id;
	  if (min>max)
	    {
	      lexsize[range]--;
	      lex[range][i] = lex[range][lexsize[range]];
	      *subtotal += min-1;
	    }
	  else
	    {
	      total += min-1;
	      right = 1;
	      while (max>=min)
		{
		  next = (min+max+1)/2;
		  if (points[noeud[next].id][range]<seuil)
		    {		
		      total += next-min+1;
		      min = next+1;
		      if (right==1)
			lex[range][i].id = min;
		    }
		  else
		    {
		      right = 0;
		      max = next-1;
		    }
		}
	    }
	}
      total += *subtotal;
    }
  else
    {
      *subtotal = 0;
      lexsize[range+1] = 0;
      for (i=0;i<size;i++)
	{
	  refnoeud = lex[range][i];
	  noeud = refnoeud.pt;
	  start = refnoeud.id;
	  min = 1;
	  max = noeud[0].id;
	  while (min!=start)
	    {
	      next = ((1+min)*num2+max*num)/den;
	      insertlex(noeud[next].pt,range+1,maxsize,lexsize,lex);
	      total += explore(noeud[next].pt,pave,range+1);
	      min = next+1;
	    }
	  right = 1;
	  while (max>=min)
	    {
	      next = ((1+min)*num2+max*num)/den;
	      if (points[noeud[next].id][range]<seuil)
		{
		  if (noeud[next].pt==NULL)
		    noeud[next].pt = subtree(noeud,min,next,range+1);
		  insertlex(noeud[next].pt,range+1,maxsize,lexsize,lex);
		  total += explore(noeud[next].pt,pave,range+1);
		  min = next+1;
		  if (right==1)
		    {
		      if (range==0)
			if (lex==lexalpha)
			  for (j=lex[range][i].id;j<next;j++)
			    if (noeud[j].pt!=NULL)
			      freetree(noeud[j].pt);
		      lex[range][i].id = min;
		    }
		}
	      else
		{
		  right = 0;
		  max = next-1;
		}
	    }
	}
    }
  return total;
}
/******************************************************************************/

static void fileformat ( void )

/******************************************************************************/
/*
  Purpose:

    FILEFORMAT reports the legal input file formats.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  fprintf(stderr,"Two FILE FORMATS are possible to give the point set x={x(1),...,x(n)} in I^s:\n\n");
  fprintf(stderr,"1) As real numbers, x(i)=(x(i,1),...,x(i,s))\n");
  fprintf(stderr,"   The file must contain the n+1 following lines:\n\n");
  fprintf(stderr,"s n reals\n");
  fprintf(stderr,"x(1,1) x(1,2) ... x(1,s)\n");
  fprintf(stderr,"...\n");
  fprintf(stderr,"x(n,1) x(n,2) ... x(n,s)\n\n");
  
  fprintf(stderr,"2) As fractions, x(i)=(num(i,1)/den(i),...,num(i,s)/den(i))\n");
  fprintf(stderr,"   The file must contain the n+1 following lines:\n\n");
  fprintf(stderr,"s n fractions\n");
  fprintf(stderr,"num(1,1) num(1,2) ... num(1,s) den(1)\n");
  fprintf(stderr,"...\n");
  fprintf(stderr,"num(n,1) num(n,2) ... num(n,s) den(n)\n");
  exit(EXIT_FAILURE);
}
/******************************************************************************/

void freetree ( sommet *noeud )

/******************************************************************************/
/*
  Purpose:
 
    FREETREE frees storage associated with a tree.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int max = noeud[0].id;
  for (i=1;i<=max;i++)
    if (noeud[i].pt!=NULL)
      freetree(noeud[i].pt);
  free (noeud);

  return;
}
/******************************************************************************/

void initlex ( void )

/******************************************************************************/
/*
  Purpose:
 
    INITLEX initializes the lexicon.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  maxsizealpha = (int *) calloc((unsigned) s,sizeof(int));
  for (i=0;i<s;i++)
    maxsizealpha[i] = 1;
  lexsizealpha = (int *) calloc((unsigned) s,sizeof(int));
  lexsizealpha[0] = 1;
  lexalpha = (sommet **) calloc((unsigned) s,sizeof(sommet *));
  for (i=0;i<s;i++)
    lexalpha[i] = (sommet *) calloc((unsigned) maxsizealpha[i],sizeof(sommet));
  lexalpha[0][0].id = 1;
  lexalpha[0][0].pt = superarbre;

  maxsizebeta = (int *) calloc((unsigned) s,sizeof(int));
  for (i=0;i<s;i++)
    maxsizebeta[i] = 1;
  lexsizebeta = (int *) calloc((unsigned) s,sizeof(int));
  lexsizebeta[0] = 1;
  lexbeta = (sommet **) calloc((unsigned) s,sizeof(sommet *));
  for (i=0;i<s;i++)
    lexbeta[i] = (sommet *) calloc((unsigned) maxsizebeta[i],sizeof(sommet));
  lexbeta[0][0].id = 1;
  lexbeta[0][0].pt = superarbre;

  return;
}
/******************************************************************************/

void initparameters ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:
 
    INITPARAMETERS sets program parameters based on user input and defaults.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  if ((argc!=5) && (argc!=7))
    usage(argv[0]);
  s = atoi(argv[1]);
  smin = s-1;
  p = atof(argv[2]);
  taille = atoi(argv[3]);
  if ((s<2) || (p<=0) || (p>=1) || (taille<1))
    usage(argv[0]);
  if (argc==7)
    {
      num = atoi(argv[5]);
      den = atoi(argv[6]);
      if ((num<1) || (den<=num))
	usage(argv[0]);
    }
  else
    {
      num = 1;
      den = 2;
    }
  num2 = den-num;

  return;
}
/******************************************************************************/

void insertlex ( sommet *noeud, int range, int *maxsize, int *lexsize,
  sommet **lex )

/******************************************************************************/
/*
  Purpose:
 
    INSERTLEX inserts an item into the lexicon.

  Modified:

    30 September 2003

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i = lexsize[range];
  if (i==maxsize[range])
    {
      maxsize[range] *= 2;
      lex[range] = realloc(lex[range],(unsigned)maxsize[range]*sizeof(sommet));
      if (lex[range]==NULL)
	memory();
    }
  lex[range][i].pt = noeud;
  lex[range][i].id = 1;
  lexsize[range] = ++i;

  return;
}
/******************************************************************************/

double lowbound ( int npoints, double volume, double *pave )

/******************************************************************************/
/*
  Purpose:
 
    LOWBOUND computes the lower bound.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int j;
  double tmp;
  if (fabs(volume-((double) npoints/taille))>borneinf)
    {
      if (((double) npoints/taille)>volume)
	{
	  volume = 1.0;
	  for (j=0;j<s;j++)
	    {
	      tmp = 0.0;
	      for (i=0;i<taille;i++)
		if ((points[i][j]>tmp) && (points[i][j]<=pave[j]))
		  tmp = points[i][j];
	      volume *= tmp;
	    }
	}
      else
	{
	  volume = 1.0;
	  for (j=0;j<s;j++)
	    {
	      tmp = 1.0;
	      for (i=0;i<taille;i++)
		if ((points[i][j]<tmp) && (points[i][j]>=pave[j]))
		  tmp = points[i][j];
	      volume *= tmp;
	    }
	}
      return fabs(volume-((double) npoints/taille));
    }
  else
    return borneinf;
}
/******************************************************************************/

static void memory ( void )

/******************************************************************************/
/*
  Purpose:
 
    MEMORY prints a message and terminates on memory allocation errors.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  fprintf(stderr,"Memory allocation problem\n");
  exit(EXIT_FAILURE);
}
/******************************************************************************/

void quicksort ( sommet *liste, int dim, int l, int r )

/******************************************************************************/
/*
  Purpose:
 
    QUICKSORT uses Quicksort to sort an array.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i = l;
  int j = r+1;
  int tmp;
  double pivot = points[liste[l].id][dim];
  while (i<j)
    {
      do
	i++;
      while ((i<r) && (points[liste[i].id][dim]<pivot));
      do
	j--;
      while (points[liste[j].id][dim]>pivot);
      if (j>i)
	{
	  tmp = liste[i].id;
	  liste[i].id = liste[j].id;
	  liste[j].id = tmp;
	}
    }
  tmp = liste[l].id;
  liste[l].id = liste[j].id;
  liste[j].id = tmp;

  if (l<j-1)
    quicksort(liste,dim,l,j-1);
  if (j+1<r)
    quicksort(liste,dim,j+1,r);

  return;
}
/******************************************************************************/

void readfile ( char *filename )

/******************************************************************************/
/*
  Purpose:
 
    READFILE reads the user's input data file.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int j;
  long *tmp;
  char fracreal[256];
  fichier=fopen(filename,"r");
  if (fichier==NULL)
    {
      fprintf(stderr,"PROBLEM: file \"%s\" doesn't exist\n\n",filename);
      fileformat();
    }
  fscanf(fichier,"%d",&i);
  if (i!=s)
    {
      fprintf(stderr,"PROBLEM: this file contains points in dimension %d, not %d\n\n",i,s);
      fileformat();
    }
  fscanf(fichier,"%d",&i);
  if (i<taille)
    {
      fprintf(stderr,"PROBLEM: this file contains only %d (<%d) points\n\n",i,taille);
      fileformat();
    }
  fgets(fracreal,sizeof fracreal,fichier);
  if ((strstr(fracreal,"fractions")==NULL) && (strstr(fracreal,"reals")==NULL))
    fileformat();

  points = (double **) calloc((unsigned) taille,sizeof(double *));
  if (points==NULL)
    memory();
  for (i=0;i<taille;i++)
    {
      points[i] = (double *) calloc((unsigned) s,sizeof(double));
      if (points[i]==NULL)
	memory();
    }

  if (strstr(fracreal,"fractions")!=NULL)
    {
      tmp = (long *) calloc((unsigned) s+1,sizeof(long *));
      for (i=0;i<taille;i++)
	{
	  for (j=0;j<=s;j++)
	    if (fscanf(fichier,"%ld",&tmp[j])==EOF)
	      {
		fprintf(stderr,"PROBLEM: this file contains only %d (<%d) points\n\n",i,taille);
		fileformat();
	      }
	  for (j=0;j<s;j++)
	    {
	      points[i][j] = (double) tmp[j]/tmp[s];
	      if ((points[i][j]<0) || (points[i][j]>1))
		{
		  fprintf(stderr,"PROBLEM: component x(%d,%d) is not in [0,1]\n\n",i+1,j+1);
		  fileformat();
		}
	    }
	}
    }
  else
    {
      for (i=0;i<taille;i++)
	for (j=0;j<s;j++)
	  {
	    if (fscanf(fichier,"%lf",&points[i][j])==EOF)
	      {
		fprintf(stderr,"PROBLEM: this file contains only %d (<%d) points\n\n",i,taille);
		fileformat();
	      }
	    if ((points[i][j]<0) || (points[i][j]>1))
	      {
		fprintf(stderr,"PROBLEM: component x(%d,%d) is not in [0,1]\n\n",i+1,j+1);
		fileformat();
	      }
	  }
    }
  fclose(fichier);

  return;
}
/******************************************************************************/

sommet *subtree ( sommet *liste, int min, int next, int dim )

/******************************************************************************/
/*
  Purpose:
 
    SUBTREE ???

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  int aux;
  int taille;
  sommet *newarbre;
  aux = min-1;
  taille = next-min+1;
  newarbre = (sommet *) calloc((unsigned) taille+1,sizeof(sommet));
  if (newarbre==NULL)
    memory();
  for (i=1;i<=taille;i++)
    newarbre[i].id=liste[i+aux].id;
  newarbre[0].id=taille;
  if (taille>1)
    quicksort(newarbre,dim,1,taille);
  return newarbre;
}
/******************************************************************************/

void supertree ( void )

/******************************************************************************/
/*
  Purpose:
 
    SUPERTREE ???

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  superarbre = (sommet *) calloc((unsigned) taille+1,sizeof(sommet));
  if (superarbre==NULL)
    memory();
  for (i=1;i<=taille;i++)
    superarbre[i].id=i-1;
  superarbre[0].id=taille;
  quicksort(superarbre,0,1,taille);

  return;
}
/******************************************************************************/

void traiter ( double *outputalpha, double *outputbeta, int range )

/******************************************************************************/
/*
  Purpose:
 
    TRAITER ???

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  int i;
  double valpha = 1.0;
  double vbeta = 1.0;
  double newborne;
  int nalpha;
  int nbeta;
  for (i=0;i<s;i++)
    {
      valpha *= outputalpha[i];
      vbeta *= outputbeta[i];
    }
  nalpha = fastexplore(outputalpha,range,maxsizealpha,lexsizealpha,lexalpha,&subtotala);
  nbeta = fastexplore(outputbeta,range,maxsizebeta,lexsizebeta,lexbeta,&subtotalb);
  newborne = ((double) nbeta/taille)-valpha;
  if (newborne>borne)
    borne = newborne;
  newborne = vbeta-((double) nalpha/taille);
  if (newborne>borne)
    borne = newborne;

  borneinf = lowbound(nalpha,valpha,outputalpha);
  borneinf = lowbound(nbeta,vbeta,outputbeta);

  return;
}
/******************************************************************************/

static void usage ( char *nom )

/******************************************************************************/
/*
  Purpose:
 
    USAGE prints a usage message.

  Reference:
 
    Eric Thiemard,
    An Algorithm to Compute Bounds for the Star Discrepancy,
    Journal of Complexity,
    Volume 17, pages 850-880, 2001.
*/
{
  fprintf(stderr,"Usage : %s s epsilon n Filename Num Den\n",nom);
  fprintf(stderr," Dimension s>=2 integer\n");
  fprintf(stderr," Accuracy parameter 0<epsilon<1\n");
  fprintf(stderr," Number of points n>=1 integer\n");
  fprintf(stderr," Filename (the file must contain at least n points)\n");
  fprintf(stderr," (Optional integer balance parameters Num and Den\n");
  fprintf(stderr,"  such that 0<Num/Den<1. Default is Num=1, Den=2)\n\n");
  fileformat();

  return;
}
