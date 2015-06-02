# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "smolpack.h"
/*
  Some global data.
*/
double c[maxdim];
extern int count;
double derf[8][20];
double w[maxdim];
double wp;

int main ( int argc, char *argv[] );
double (*f)(int, double x[]);  
double f1 ( int dim, double x[] );
double f2 ( int dim, double x[] );
double f3 ( int dim, double x[] );
double f4 ( int dim, double x[] );
double f5 ( int dim, double x[] );
double f6 ( int dim, double x[] );
double f7 ( int dim, double x[] );
void init_erf ( void );
double integral ( int fnum, int dim );
double my_erfc ( double x );
void timestamp ( );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    SMOLPACK_INTERACTIVE is an interactive test program for SMOLPACK.

  Modified:

    30 April 2007

  Author:

    Knut Petras

  Reference:

    Knut Petras,
    Smolyak Cubature of Given Polynomial Degree with Few Nodes
    for Increasing Dimension,
    Numerische Mathematik,
    Volume 93, Number 4, February 2003, pages 729-753.

  Parameters:

    Commandline parameter, int FNUM, the index of the test function.
    1 <= FNUM <= 7.

    Commandline parameter, int DIM, the spatial dimension.
    1 <= DIM <= MAXDIM = 40.

    Commandline parameter, int K, the number of integration stages.
    0 <= K <= ?.

    Commandline parameter, int SEED, a seed for the random number generator.
    SEED should not be zero.

    Commandline parameter, int BS, chooses the basic sequence.
    1 = delayed Clenshaw Curtis, (fewer function evaluations),
    2 = standard Clenshaw Curtis. 
*/
{
  int bs;
  int dim;
  double exact;
  int fnum;
  int i;
  int j;
  int k;
  int print_stats;
  int q;
  int qmax = 12;
  double quad;
  int seed;
  double sum;

  printf ( "\n" );
  timestamp ( );

  printf ( "\n" );
  printf ( "SMOLPACK_INTERACTIVE\n" );
  printf ( "  C version\n" );

  if ( argc < 1 )
  {
    printf ( "\n" );
    printf ( "Enter the number of the function to be integrated:\n"); 

    printf ( "1: Oscillatory\n" ); 
    printf ( "2: Product peak\n" ); 
    printf ( "3: Corner peak\n" ); 
    printf ( "4: Gaussian\n" ); 
    printf ( "5: Continuous\n" ); 
    printf ( "6: Discontinuous\n" ); 
    printf ( "7: exp(sum(x[i]))\n" );

    scanf ( "%i", &fnum ); 
  }
  else
  {
    argv++;
    fnum = atoi ( *argv );
  }
    
  if ( fnum == 1 )
  {
    f = f1;
  }
  else if ( fnum == 2 )
  {
    f = f2;
  }
  else if ( fnum == 3 )
  {
    f = f3;
  }
  else if ( fnum == 4 )
  {
    f = f4;
  }
  else if ( fnum == 5 )
  {
    f = f5;
  }
  else if ( fnum == 6 )
  {
    f = f6;
  }
  else if ( fnum == 7 )
  {
    f = f7;
  }
  else
  {
    printf ( "\n" );
    printf ( "SMOLPACK_INTERACTIVE - Fatal error!\n" );
    printf ( "  Illegal value of FNUM = %d\n", fnum );
    printf ( "  1 <= FNUM <= 7 is required.\n" );
    exit ( 1 );
  }

  if ( argc < 2 )
  {
    printf ( "\n" );
    printf ( "Enter DIM, the spatial dimension.\n"); 
    scanf ( "%i", &dim ); 
  }
  else
  {
    argv++;
    dim = atoi ( *argv );
  }

  if ( dim < 0 )
  {
    printf ( "\n" );
    printf ( "SMOLPACK_INTERACTIVE - Fatal error!\n" );
    printf ( "  Illegal value of DIM = %d\n", dim );
    printf ( "  DIM must be 1 or greater.\n" );
    exit ( 1 );
  }

  if ( argc < 3 )
  {
    printf ( "\n" );
    printf ( "Enter K, the number of stages.\n"); 
    scanf ( "%i", &k ); 
  }
  else
  {
    argv++;
    k = atoi ( *argv );
  }

  if ( k < 0 )
  {
    printf ( "\n" );
    printf ( "SMOLPACK_INTERACTIVE - Fatal error!\n" );
    printf ( "  Illegal value of K = %d\n", k );
    printf ( "  K must be 0 or greater.\n" );
    exit ( 1 );
  }

  if ( argc < 4 )
  {
    printf ( "\n" );
    printf ( "Enter a seed for the random number generator:\n"); 
    scanf ( "%i", &seed ); 
  }
  else
  {
    argv++;
    seed = atoi ( *argv );
  }

  if ( argc < 5 )
  {
    printf ( "\n" );
    printf ( "Choose the basic sequence:\n"); 
    printf ( " 1 : delayed Clenshaw-Curtis rule (fewer evaluations)\n"); 
    printf ( " 2 : standard Clenshaw-Curtis rule\n"); 
    scanf ( "%i", &bs ); 
  }
  else
  {
    argv++;
    bs = atoi ( *argv );
  }

  if ( bs != 1 && bs != 2 )
  {
    printf ( "\n" );
    printf ( "SMOLPACK_INTERACTIVE - Fatal error!\n" );
    printf ( "  Illegal value of BS = %d\n", bs );
    printf ( "  BS = 1 or BS = 2 are the legal choices.\n" );
    exit ( 1 );
  }

  printf ( "\n" );
  printf ( "  Test function number FNUM = %d\n", fnum ); 
  printf ( "  Spatial dimension DIM =     %d\n", dim ); 
  printf ( "  Number of states K =        %d\n", k ); 
  printf ( "  Random number SEED =        %d\n", seed ); 
  printf ( "  Basic sequence choice BS =  %d\n", bs ); 
/*
  Set the C vector to a random value.
  This is a parameter vector for the integrand functions.
*/
  srand ( seed );      

  sum = 0.0;

  for ( i = 0; i < dim; i++ )
  {
    c[i] = ( ( double ) rand ( ) ) / ( double ) RAND_MAX; 
    sum = sum + c[i];
  }

  for ( i = 0; i < dim; i++ )
  {
    c[i] = c[i] / sum * 9.0; 
  }
/*
  Set the W vector to a random value.
*/
  for ( i = 0; i <= dim; i++ ) 
  {
    w[i] = ( ( double ) rand ( ) ) / ( double ) RAND_MAX; 
  }  
/*
  Call the chosen Smolyak algorithm.
  Why does this loop run from DIM+K to DIM+K?
  If we ran it from Q = DIM, would we get the nested results maybe???
*/
  print_stats = 1;

  for ( q = dim + k; q <= dim + k; q++ )
  {
    if ( bs == 1 )
    {
      quad = int_smolyak ( dim, q, f, print_stats );
    }
    else if ( bs == 2 )
    {
      quad = cc_int_smolyak ( dim, q, f, print_stats );
    }
    exact = integral ( fnum, dim );
    printf ( "\n" );
    printf ( "  Exact:    %15.10e\n", exact );      

    printf ( "  Estimate: %15.10e \n",  quad ); 
    printf ( "  Error:    %15.10e \n",  fabs ( quad - exact ) ); 
  }

  printf ( "\n" );
  printf ( "SMOLPACK_INTERACTIVE\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;    
} 
/******************************************************************************/

double f1 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F1 evaluates Genz's "oscillatory" test integrand.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F1, the value of the integrand function at X.
*/
{
  int i;
  double pi = 3.141592653589793;
  double value;

  count++;
  value = 2.0 * pi * w[0];
  for ( i = 0; i < dim; i++ )
  {
    value = value + c[i] * x[i];
  }
  value = cos ( value );

  return value;
}
/******************************************************************************/

double f2 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F2 evaluates Genz's "product peak" test integrand.

  Discussion:

    This version of the product peak function uses 1/c[i] where
    normally we would write c[i].  Since c[i] is more or less a
    random parameter, this doesn't really matter, as long as c[i] 
    is nonzero!

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F2, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;
  value = 1.0;
    
  for ( i = 0; i < dim; i++ ) 
  {
    value = value * ( pow ( 1.0 / c[i], 2 ) + pow ( x[i] - w[i], 2 ) );
  }
  value = 1.0 / value;

  return value;
}
/******************************************************************************/

double f3 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F3 evaluates Genz's "corner peak" test integrand.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F3, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;
  value = 1.0;
  for ( i = 0; i < dim; i++ )
  {
    value = value + c[i] * x[i];
  }
  value = 1.0 / pow ( value, dim + 1 );

  return value;
}
/******************************************************************************/

double f4 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F4 evaluates Genz's "Gaussian" test integrand.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F4, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;
  value = 0.0;
  for ( i = 0; i < dim; i++ )  
  {
    value = value - pow ( c[i] * ( x[i] - w[i] ), 2 );
  }
  value = exp ( value );

  return value;
}
/******************************************************************************/

double f5 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F5 evaluates Genz's "Continuous" test integrand.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F5, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;
  value = 0.0;

  for ( i = 0; i < dim; i++ )
  {
    value = value - c[i] * fabs ( x[i] - w[i] );
  }

  value = exp ( value );

  return value;
}
/******************************************************************************/

double f6 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F6 evaluates Genz's "Discontinuous" test integrand.

  Modified:

    30 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F6, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;

  value = 0.0;

  if ( dim == 1 )
  {
    if ( x[0] <= w[0] )
    {
      for ( i = 0; i < dim; i++ )
      {
        value = value + c[i] * x[i];
      }
      value = exp ( value );
    }
  }
  else
  {
    if ( x[0] <= w[0] && x[1] <= w[1] )
    {
      for ( i = 0; i < dim; i++ )
      {
        value = value + c[i] * x[i];
      }
      value = exp ( value );
    }
  }

  return value;
}
/******************************************************************************/

double f7 ( int dim, double x[] )

/******************************************************************************/
/*
  Purpose:

    F7 evaluates the exp(sum(x(i))) test integrand.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int DIM, the spatial dimension.

    Input, double X[DIM], the argument.

    Output, double F7, the value of the integrand function at X.
*/
{
  int i;
  double value;

  count++;
  value = 0.0;
  for ( i = 0; i < dim; i++ )
  {
    value = value + x[i];
  }
  value = exp ( value );

  return value;
}
/******************************************************************************/

void init_erf ( void )

/******************************************************************************/
/*
  Purpose:

    INIT_ERF initializes data for MY_ERFC.

  Modified:

    26 April 2007

  Author:

    Knut Petras
*/
{
  double h;
  double hilf;
  int i;
  int j;
  double pi = 3.141592653589793;

  wp = sqrt ( pi );
  hilf = 2.0 / wp;

  derf[0][0] = 1.0;
  derf[1][0] = 1.5729920705028513065877936491E-01;
  derf[2][0] = 4.6777349810472658379307436327E-03;
  derf[3][0] = 2.2090496998585441372776129582E-05;
  derf[4][0] = 1.5417257900280018852159673486E-08;
  derf[5][0] = 1.5374597944280348501883434853E-12;
  derf[6][0] = 2.1519736712498913116593350399E-17;
  derf[7][0] = 4.1838256077794143986140102238E-23;

  for ( i = 0; i <= 7; i++ )
  {
    derf[i][1] = - exp ( - pow ( ( double ) i, 2 ) ) * hilf;
    derf[i][2] = - ( double ) i * derf[i][1]; 

    for ( j = 3; j <= 19; j++ )
    {
      derf[i][j] = -2.0 * ( ( double ) i * derf[i][j-1] 
                 + ( double ) ( j - 2 ) * derf[i][j-2] 
                 / ( double ) ( j - 1 ) ) / ( double ) j; 
    }
  }
  return;
}
/******************************************************************************/

double integral ( int fnum, int dim )

/******************************************************************************/
/*
  Purpose:

    INTEGRAL returns the exact integral of a test function.

  Modified:

    30 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, int FNUM, the index of the test function.

    Input, int DIM, the spatial dimension.

    Output, double INTEGRAL, the value of the integral of the test function 
    over the DIM-dimensional unit hypercube.
*/
{
  int a;
  double arg;
  double bot;
  double c_prod;
  int d;
  int i;
  int *ivec;
  int ivec_sum;
  double pi = 3.141592653589793;
  double prod;
  int rank;
  double sum;
  double total;
  double value;
/*
  #1: Oscillatory.
*/
  if ( fnum == 1 )
  {
    arg = 0.0;
    for ( i = 0; i < dim; i++ )
    {
      arg = arg + c[i];
    }

    prod = 1.0;
    for ( i = 0; i < dim; i++ )
    {
      prod = prod * sin ( 0.5 * c[i] ) / c[i];
    }

    value = pow ( 2.0, dim ) * cos ( 2.0 * pi * w[0] + 0.5 * arg ) * prod;
  }
/*
  #2: Product Peak.
*/
  else if ( fnum == 2 )
  {
    value = 1.0;
    for ( i = 0; i < dim; i++ )
    {
      value = value * c[i] * ( atan ( c[i] * ( 1.0 - w[i] ) ) 
                             + atan ( c[i] *         w[i]   ) );
    } 
  }
/*
  #3: Corner Peak
*/
  else if ( fnum == 3 )
  {
    ivec = ( int * ) malloc ( dim * sizeof ( int ) );

    total = 0.0;
    rank = 0;

    for ( ; ; )
    {
      tuple_next ( 0, 1, dim, &rank, ivec );

      if ( rank == 0 )
      {
        break;
      }

      ivec_sum = 0;
      for ( i = 0; i < dim; i++ )
      {
        ivec_sum = ivec_sum + ivec[i];
      }

      bot = 1.0;
      for ( i = 0; i < dim; i++ )
      {
        if ( ivec[i] == 1 )
        {
          bot = bot + c[i];
        }
      }
      total = total + pow ( -1.0, ivec_sum ) / bot;
    }

    a = 1;
    for ( i = 1; i <= dim; i++ )
    {
      a = a * i;
    }

    c_prod = 1.0;
    for ( i = 0; i < dim; i++ )
    {
      c_prod = c_prod * c[i];
    }
    value = total / ( ( double ) a * c_prod );

    free ( ivec );
  }
/*
  #4: Gaussian.
*/
  else if ( fnum == 4 )
  {
    init_erf ( ); 
    value = 1.0;
    for ( i = 0; i < dim; i++ )
    {
      value = value * sqrt ( pi ) / ( 2.0 * c[i] ) 
            * ( my_erfc ( -c[i] *         w[i] )
              - my_erfc (  c[i] * ( 1.0 - w[i] ) ) );
    }
  }
/*
  #5: Continuous function.
*/
  else if ( fnum == 5 )
  {
    value = 1.0;
    for ( i = 0; i < dim; i++ )
    {
      value = value / c[i] 
        * ( 2.0 - exp ( -c[i] * w[i] ) - exp ( c[i] * ( w[i] - 1.0 ) ) );
    }
  }
/*
  #6  Discontinuous function.
*/
  else if ( fnum == 6 )
  {
    value = 1.0;
 
    if ( dim < 2 )
    {
      for ( i = 0; i < dim; i++ )
      {
        value = value * ( exp ( c[i] * w[i] ) - 1.0 ) / c[i];
      }
    }
    else
    {
      for ( i = 0; i <= 1; i++ )
      {
        value = value * ( exp ( c[i] * w[i] ) - 1.0 ) / c[i];
      }
      for ( i = 2; i < dim; i++ )
      {
        value = value * ( exp ( c[i] ) - 1.0 ) / c[i];
      }
    }
  }
/*
  #7: exp(sum(x(i)))
*/
  else if ( fnum == 7 )
  {
    value = pow ( exp ( 1.0 ) - 1.0, dim );
  }
/*
  Unexpected call!
*/
  else
  {
    value = 0.0;
    printf ( "\n" );
    printf ( "INTEGRAL - Fatal error!\n" );
    printf ( "  Input function index FNUM must be between 1 and 7.\n" );
    printf ( "  This value was FNUM = %d\n", fnum );
    exit ( 1 );
  }

  return value;
}
/******************************************************************************/

double my_erfc ( double x )

/******************************************************************************/
/*
  Purpose:

    MY_ERFC evaluates the complementary error function.

  Modified:

    26 April 2007

  Author:

    Knut Petras

  Parameters:

    Input, double X, the argument.

    Output, double MY_ERFC, the value of the complementary error function at X.
*/
{
  double h;
  double hilf;
  int i;
  int j;
  double prod[60];
  double sm;
  double smme;
  double value;

  if ( x < 0 )
  {
    value = 2.0 - my_erfc ( - x ); 
  }
  else if ( x < 7.5 )
  {
    i = ( int ) ( x + 0.499999 ); 
    h = x - ( double ) i;
    value = h * derf[i][19] + derf[i][18]; 
    for ( j = 17; 0 <= j; j-- )
    {
      value = h * value + derf[i][j];    
    }
  }
  else
  {
    hilf = 0.5 / pow ( x, 2 );
    prod[0] = 1.0;
    i = 2 + ( int ) ( 190.0 / x );
    for ( j = 1; j <= i; j++ )
    {
      prod[j] = -prod[j-1] * ( 2 * j - 1 ) * hilf; 
    }
    smme = 0.0;
    for ( j = i; 0 <= j; j-- )
    {
      smme = smme + prod[j]; 
    }
    value = ( smme * exp ( - pow ( x, 2 ) ) / ( wp * x ) );
  }
  return value;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

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
/******************************************************************************/

void tuple_next ( int m1, int m2, int n, int *rank, int x[] )

/******************************************************************************/
/*
  Purpose:

    TUPLE_NEXT computes the next element of a tuple space.

  Discussion:

    The elements are N vectors.  Each entry is constrained to lie
    between M1 and M2.  The elements are produced one at a time.
    The first element is
      (M1,M1,...,M1),
    the second element is
      (M1,M1,...,M1+1),
    and the last element is
      (M2,M2,...,M2)
    Intermediate elements are produced in lexicographic order.

  Example:

    N = 2, M1 = 1, M2 = 3

    INPUT        OUTPUT
    -------      -------
    Rank  X      Rank   X
    ----  ---    -----  ---
    0     * *    1      1 1
    1     1 1    2      1 2
    2     1 2    3      1 3
    3     1 3    4      2 1
    4     2 1    5      2 2
    5     2 2    6      2 3
    6     2 3    7      3 1
    7     3 1    8      3 2
    8     3 2    9      3 3
    9     3 3    0      0 0

  Modified:

    29 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, int M1, M2, the minimum and maximum entries.

    Input, int N, the number of components.

    Input/output, int *RANK, counts the elements.
    On first call, set RANK to 0.  Thereafter, the output value of RANK
    will indicate the order of the element returned.  When there are no
    more elements, RANK will be returned as 0.

    Input/output, int X[N], on input the previous tuple.
    On output, the next tuple.
*/
{
  int i;
  int j;

  if ( m2 < m1 )
  {
    *rank = 0;
    return;
  }

  if ( *rank <= 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = m1;
    }
    *rank = 1;
  }
  else
  {
    *rank = *rank + 1;
    i = n - 1;

    for ( ; ; )
    {

      if ( x[i] < m2 )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = m1;

      if ( i == 0 )
      {
        *rank = 0;
        for ( j = 0; j < n; j++ )
        {
          x[j] = m1;
        }
        break;
      }
      i = i - 1;
    }
  }
  return;
}
