# include <complex.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "uniform.h"

/******************************************************************************/

void bvec_print ( int n, int bvec[], char *title )

/******************************************************************************/
/*
  Purpose:

    BVEC_PRINT prints a BVEC, with an optional title.

  Discussion:

    A BVEC is a vector of binary digits representing an integer.  

    BVEC[0] is 0 for positive values and 1 for negative values, which
    are stored in 2's complement form.

    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
    so that printing the digits in order gives the binary form of the number.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 December 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int BVEC[N], the vector to be printed.

    Input, char *TITLE, a title to be printed first.
    TITLE may be blank.
*/
{
  int i;
  int ihi;
  int ilo;

  if ( 0 < strlen ( title ) )
  {
    printf ( "\n" );
    printf ( "%s\n", title );
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 70 )
  {
    ihi = i4_min ( ilo + 70 - 1, n - 1 );
    printf ( "  " );
    for ( i = ilo; i <= ihi; i++ )
    {
      printf ( "%d", bvec[i] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

int *bvec_uniform_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    BVEC_UNIFORM_NEW returns a random binary vector.

  Discussion:

    A BVEC is a vector of binary digits representing an integer.  

    BVEC[0] is 0 for positive values and 1 for negative values, which
    are stored in 2's complement form.

    For positive values, BVEC[N-1] contains the units digit, BVEC[N-2]
    the coefficient of 2, BVEC[N-3] the coefficient of 4 and so on,
    so that printing the digits in order gives the binary form of the number.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the length of the vectors.

    Input/output, int *SEED, a seed for the random number generator.

    Output, int BVEC_UNIFORM_NEW[N], the randomly selected vector.
*/
{
  const int i4_huge      = 2147483647;
  const int i4_huge_half = 1073741823;
  int i;
  int k;
  int *bvec;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "BVEC_UNIFORM_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    return NULL;
  }

  bvec = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }
    if ( i4_huge_half < *seed )
    {
      bvec[i] = 0;
    }
    else
    {
      bvec[i] = 1;
    }
  }

  return bvec;
}
/******************************************************************************/

float complex c4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    C4_UNIFORM_01 returns a unit pseudorandom C4.

  Discussion:

    The angle should be uniformly distributed between 0 and 2 * PI,
    the square root of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C4_UNIFORM_01, a pseudorandom complex value.
*/
{
  const int i4_huge = 2147483647;
  int k;
  float r;
  const float r4_pi = 3.1415926E+00;
  float theta;
  float complex value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C4_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

  value = r * ( cos ( theta ) + I * sin ( theta ) );

  return value;
}
/******************************************************************************/

void c4mat_print ( int m, int n, float complex a[], char *title )

/******************************************************************************/
/*
  Purpose:

    C4MAT_PRINT prints a C4MAT.

  Discussion:

    A C4MAT is a matrix of single precision complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input, float complex A[M*N], the matrix.

    Input, char *TITLE, a title.
*/
{
  c4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void c4mat_print_some ( int m, int n, float complex a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    C4MAT_PRINT_SOME prints some of a C4MAT.

  Discussion:

    A C4MAT is a matrix of float complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input, float complex A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
  float complex c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );
    printf ( "  Col: " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      printf ( "          %10d", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = 1;
    if ( i2lo < ilo )
    {
      i2lo = ilo;
    }
    i2hi = m;
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }
    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) INCX entries in row I, that lie in the current strip.
*/
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        printf ( "  %8g  %8g", creal ( c ), cimag ( c ) );
      }
      printf ( "\n" );
    }
  }
  return;
}
/******************************************************************************/

void c4mat_uniform_01 ( int m, int n, int *seed, float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C[M*N], the pseudorandom complex matrix.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  const float r4_pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C4MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta ) + I * sin ( theta ) );
    }
  }

  return;
}
/******************************************************************************/

float complex *c4mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C4MAT_UNIFORM_01_NEW returns a unit pseudorandom C4MAT.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C4MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
    complex matrix.
*/
{
  float complex *c;
  int i;
  const int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  const float r4_pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C4MAT_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = ( float complex * ) malloc ( m * n * sizeof ( float complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta ) + I * sin ( theta ) );
    }
  }

  return c;
}
/******************************************************************************/

void c4vec_print ( int n, float complex a[], char *title )

/******************************************************************************/
/*
  Purpose:

    C4VEC_PRINT prints a C4VEC.

  Discussion:

    A C4VEC is a vector of float complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, float complex A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal( a[i] ), cimag ( a[i] ) );
  }

  return;
}
/******************************************************************************/

void c4vec_uniform_01 ( int n, int *seed, float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of values to compute.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C[N], the pseudorandom 
    complex vector.
*/
{
  int i;
  const int i4_huge = 2147483647;
  float r;
  int k;
  const float r4_pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C4VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * ( cos ( theta ) + I * sin ( theta ) );
  }

  return;
}
/******************************************************************************/

float complex *c4vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C4VEC_UNIFORM_01_NEW returns a unit pseudorandom C4VEC.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of values to compute.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C4VEC_UNIFORM_01_NEW[N], the pseudorandom 
    complex vector.
*/
{
  float complex *c;
  int i;
  const int i4_huge = 2147483647;
  float r;
  int k;
  const float r4_pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C4VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = ( float complex * ) malloc ( n * sizeof ( float complex ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * ( cos ( theta ) + I * sin ( theta ) );
  }

  return c;
}
/******************************************************************************/

double complex c8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    C8_UNIFORM_01 returns a unit pseudorandom C8.

  Discussion:

    The angle should be uniformly distributed between 0 and 2 * PI,
    the square root of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C8_UNIFORM_01, a pseudorandom complex value.
*/
{
  const int i4_huge = 2147483647;
  int k;
  double r;
  const double r8_pi = 3.141592653589793;
  double theta;
  double complex value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = sqrt ( ( ( double ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  theta = 2.0 * r8_pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = r * ( cos ( theta ) + I * sin ( theta ) );

  return value;
}
/******************************************************************************/

void c8mat_print ( int m, int n, double complex a[], char *title )

/******************************************************************************/
/*
  Purpose:

    C8MAT_PRINT prints a C8MAT.

  Discussion:

    A C8MAT is a matrix of double precision complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input, double complex A[M*N], the matrix.

    Input, char *TITLE, a title.
*/
{
  c8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void c8mat_print_some ( int m, int n, double complex a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    C8MAT_PRINT_SOME prints some of a C8MAT.

  Discussion:

    A C8MAT is a matrix of double precision complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input, double complex A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
  double complex c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );
    printf ( "  Col: " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      printf ( "          %10d", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = 1;
    if ( i2lo < ilo )
    {
      i2lo = ilo;
    }
    i2hi = m;
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }
    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) INCX entries in row I, that lie in the current strip.
*/
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        printf ( "  %8g  %8g", creal ( c ), cimag ( c ) );
      }
      printf ( "\n" );
    }
  }
  return;
}
/******************************************************************************/

void c8mat_uniform_01 ( int m, int n, int *seed, double complex c[] )

/******************************************************************************/
/*
  Purpose:

    C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C[M*N], the pseudorandom 
    complex matrix.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * r8_pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta )+ I * sin ( theta ) );
    }
  }

  return;
}
/******************************************************************************/

double complex *c8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C8MAT_UNIFORM_01_NEW returns a unit pseudorandom C8MAT.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C8MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
    complex matrix.
*/
{
  double complex *c;
  int i;
  const int i4_huge = 2147483647;
  int j;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8MAT_UNIFORM_01_NEW- Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = ( double complex * ) malloc ( m * n * sizeof ( double complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * r8_pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta )+ I * sin ( theta ) );
    }
  }

  return c;
}
/******************************************************************************/

void c8vec_print ( int n, double complex a[], char *title )

/******************************************************************************/
/*
  Purpose:

    C8VEC_PRINT prints a C8VEC.

  Discussion:

    A C8VEC is a vector of double complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double complex A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal( a[i] ), cimag ( a[i] ) );
  }

  return;
}
/******************************************************************************/

void c8vec_uniform_01 ( int n, int *seed, double complex c[] )

/******************************************************************************/
/*
  Purpose:

    C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of values to compute.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C[N], the pseudorandom vector.
*/
{
  int i;
  const int i4_huge = 2147483647;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * r8_pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * ( cos ( theta ) + I * sin ( theta ) );
  }

  return;
}
/******************************************************************************/

double complex *c8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C8VEC_UNIFORM_01_NEW returns a unit pseudorandom C8VEC.

  Discussion:

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of values to compute.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C8VEC_UNIFORM_01_NEW[N], the pseudorandom vector.
*/
{
  double complex *c;
  int i;
  const int i4_huge = 2147483647;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = ( double complex * ) malloc ( n * sizeof ( double complex ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * r8_pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * ( cos ( theta ) + I * sin ( theta ) );
  }

  return c;
}
/*****************************************************************************/

char ch_uniform_ab ( char a, char b, int *seed )

/******************************************************************************/
/*
  Purpose:

    CH_UNIFORM_AB returns a scaled CH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, char A, B, the minimum and maximum acceptable characters.

    Input/output, int *SEED, a seed for the random number generator.

    Output, char CH_UNIFORM, the randomly chosen character.
*/
{
  char c;
  float r;

  r = r4_uniform_01 ( seed );

  c = a + ( char ) ( r * ( float ) ( b + 1 - a ) );

  return c;
}
/******************************************************************************/

int congruence ( int a, int b, int c, int *error )

/******************************************************************************/
/*
  Purpose:

    CONGRUENCE solves a congruence of the form A * X = C ( mod B ).

  Discussion:

    A, B and C are given integers.  The equation is solvable if and only
    if the greatest common divisor of A and B also divides C.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 November 2004

  Author:

    John Burkardt

  Reference:

    Eric Weisstein, editor,
    CRC Concise Encylopedia of Mathematics,
    CRC Press, 1998, page 446.

  Parameters:

    Input, int A, B, C, the coefficients of the Diophantine equation.

    Output, int *ERROR, error flag, is 1 if an error occurred..

    Output, int CONGRUENCE, the solution of the Diophantine equation.
    X will be between 0 and B-1.
*/
{
# define N_MAX 100

  int a_copy;
  int a_mag;
  int a_sign;
  int b_copy;
  int b_mag;
  int b_sign;
  int c_copy;
  int g;
  int k;
  int n;
  int q[N_MAX];
  int swap;
  int x;
  int y;
  int z;
/*
  Defaults for output parameters.
*/
  *error = 0;
  x = 0;
  y = 0;
/*
  Special cases.
*/
  if ( a == 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b == 0 && c != 0 )
  {
    *error = 1;
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c != 0 )
  {
    x = 0;
    if ( ( c % b ) != 0 )
    {
      *error = 2;
    }
    return x;
  }
  else if ( a != 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a != 0 && b == 0 && c != 0 )
  {
    x = c / a;
    if ( ( c % a ) != 0 )
    {
      *error = 3;
      return x;
    }
    return x;
  }
  else if ( a != 0 && b != 0 && c == 0 )
  {
/*  g = i4_gcd ( a, b ); */
/*  x = b / g; */
    x = 0;
    return x;
  }
/*
  Now handle the "general" case: A, B and C are nonzero.

  Step 1: Compute the GCD of A and B, which must also divide C.
*/
  g = i4_gcd ( a, b );

  if ( ( c % g ) != 0 )
  {
    *error = 4;
    return x;
  }

  a_copy = a / g;
  b_copy = b / g;
  c_copy = c / g;
/*
  Step 2: Split A and B into sign and magnitude.
*/
  a_mag = abs ( a_copy );
  a_sign = i4_sign ( a_copy );
  b_mag = abs ( b_copy );
  b_sign = i4_sign ( b_copy );
/*
  Another special case, A_MAG = 1 or B_MAG = 1.
*/
  if ( a_mag == 1 )
  {
    x = a_sign * c_copy;
    return x;
  }
  else if ( b_mag == 1 )
  {
    x = 0;
    return x;
  }
/*
  Step 3: Produce the Euclidean remainder sequence.
*/
  if ( b_mag <= a_mag )
  {
    swap = 0;
    q[0] = a_mag;
    q[1] = b_mag;
  }
  else
  {
    swap = 1;
    q[0] = b_mag;
    q[1] = a_mag;
  }

  n = 3;

  for ( ; ; )
  {
    q[n-1] = ( q[n-3] % q[n-2] );

    if ( q[n-1] == 1 )
    {
      break;
    }

    n = n + 1;

    if ( N_MAX < n )
    {
      *error = 1;
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "CONGRUENCE - Fatal error!\n" );
      fprintf ( stderr, "  Exceeded number of iterations.\n" );
      exit ( 1 );
    }
  }
/*
  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
*/
  y = 0;
  for ( k = n; 2 <= k; k-- )
  {
    x = y;
    y = ( 1 - x * q[k-2] ) / q[k-1];
  }
/*
  Step 5: Undo the swapping.
*/
  if ( swap == 1 )
  {
    z = x;
    x = y;
    y = z;
  }
/*
  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
*/
  x = x * a_sign;
/*
  Step 7: Multiply by C, so that X * A + Y * B = C.
*/
  x = x * c_copy;
/*
  Step 8: Now force 0 <= X < B.
*/
  x = x % b;
/*
  Step 9: Force positivity.
*/
  if ( x < 0 )
  {
    x = x + b;
  }

  return x;
# undef N_MAX
}
/******************************************************************************/

char digit_to_ch ( int i )

/******************************************************************************/
/*
  Purpose:

    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.

  Example:

     I     C
   -----  ---
     0    '0'
     1    '1'
   ...    ...
     9    '9'  
    10    '*'
   -83    '*'

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, int I, the digit, which should be between 0 and 9.

    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
*/
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
/******************************************************************************/

int get_seed ( )

/******************************************************************************/
/*
  Purpose:

    GET_SEED returns a random seed for the random number generator.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2004

  Author:

    John Burkardt

  Parameters:

    Output, int GET_SEED, a random seed value.
*/
{
  time_t clock;
  const int i4_huge = 2147483647;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
/*
  If the internal seed is 0, generate a value based on the time.
*/
  clock = time ( &tloc );
  lt = localtime ( &clock );
/*
  Hours is 1, 2, ..., 12.
*/
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
/*
  Move Hours to 0, 1, ..., 11
*/
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
/*
  We want values in [1,43200], not [0,43199].
*/
  seed = seed + 1;
/*
  Remap SEED from [1,43200] to [1,HUGE].
*/
  seed = ( int ) 
    ( ( ( double ) seed )
    * ( ( double ) i4_huge ) / ( 60.0 * 60.0 * 12.0 ) );
/*
  Never use a seed of 0.
*/
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
}
/******************************************************************************/

int i4_gcd ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_GCD finds the greatest common divisor of I and J.

  Discussion:

    Only the absolute values of I and J are considered, so that the 
    result is always nonnegative.

    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).

    If I and J have no common factor, I4_GCD is returned as 1.

    Otherwise, using the Euclidean algorithm, I4_GCD is the
    largest common factor of I and J.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, two numbers whose greatest common divisor
    is desired.

    Output, int I4_GCD, the greatest common divisor of I and J.
*/
{
  int ip;
  int iq;
  int ir;
/*
  Return immediately if either I or J is zero.
*/
  if ( i == 0 )
  {
    return i4_max ( 1, abs ( j ) );
  }
  else if ( j == 0 )
  {
    return i4_max ( 1, abs ( i ) );
  }
/*
  Set IP to the larger of I and J, IQ to the smaller.
  This way, we can alter IP and IQ as we go.
*/
  ip = i4_max ( abs ( i ), abs ( j ) );
  iq = i4_min ( abs ( i ), abs ( j ) );
/*
  Carry out the Euclidean algorithm.
*/
  for ( ; ; )
  {
    ir = ip % iq;

    if ( ir == 0 )
    {
      break;
    }

    ip = iq;
    iq = ir;
  }

  return iq;
}
/******************************************************************************/

int i4_huge ( )

/******************************************************************************/
/*
  Purpose:

    I4_HUGE returns a "huge" I4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Output, int I4_HUGE, a "huge" integer.
*/
{
  return 2147483647;
}
/******************************************************************************/

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the whole part of the logarithm base 10 of an I4.

  Discussion:

    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
    (except for I = 0).

    The number of decimal digits in I is I4_LOG_10(I) + 1.

  Example:

        I    I4_LOG_10(I)

        0     0
        1     0
        2     0

        9     0
       10     1
       11     1

       99     1
      100     2
      101     2

      999     2
     1000     3
     1001     3

     9999     3
    10000     4
    10001     4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, int I, the integer.

    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
*/
{
  int ten_pow;
  int value;

  i = abs ( i );

  ten_pow = 10;
  value = 0;

  while ( ten_pow <= i )
  {
    ten_pow = ten_pow * 10;
    value = value + 1;
  }

  return value;
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

int i4_seed_advance ( int seed )

/******************************************************************************/
/*
  Purpose:

    I4_SEED_ADVANCE "advances" the seed.

  Discussion:

    This routine implements one step of the recursion

      SEED = ( 16807 * SEED ) mod ( 2^31 - 1 )

    This version of the routine does not check whether the input value of
    SEED is zero.  If the input value is zero, the output value will be zero.

    If we repeatedly use the output of SEED_ADVANCE as the next input, 
    and we start with SEED = 12345, then the first few iterates are:

         Input      Output
          SEED        SEED

         12345   207482415
     207482415  1790989824
    1790989824  2035175616
    2035175616    77048696
      77048696    24794531

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2013

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int SEED, the seed value.

    Output, int I4_SEED_ADVANCE, the "next" seed.
*/
{
  const int i4_huge = 2147483647;
  int k;
  int seed_new;

  seed_new = seed;

  if ( seed_new < 0 )
  {
    seed_new = seed_new + i4_huge;
  }

  k = seed_new / 127773;

  seed_new = 16807 * ( seed_new - k * 127773 ) - k * 2836;

  if ( seed_new < 0 )
  {
    seed_new = seed_new + i4_huge;
  }
  return seed_new;
}
/******************************************************************************/

int i4_sign ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_SIGN returns the sign of an I4.

  Discussion:

    The sign of 0 and all positive integers is taken to be +1.
    The sign of all negative integers is -1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int I, the integer whose sign is desired.

    Output, int I4_SIGN, the sign of I.
*/
{
  int value;

  if ( i < 0 ) 
  {
    value = - 1;
  }
  else
  {
    value = 1;
  }
  return value;
}
/******************************************************************************/

void i4_swap ( int *i, int *j )

/******************************************************************************/
/*
  Purpose:

    I4_SWAP switches two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 January 2002

  Author:

    John Burkardt

  Parameters:

    Input/output, int *I, *J.  On output, the values of I and
    J have been interchanged.
*/
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}
/******************************************************************************/

char *i4_to_s ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_S converts an I4 to a string.

  Example:

    INTVAL  S

         1  1
        -1  -1
         0  0
      1952  1952
    123456  123456
   1234567  1234567

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, int I, an integer to be converted.

    Output, char *I4_TO_S, the representation of the integer.
*/
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;

  length = i4_log_10 ( i );

  ten_power = ( int ) pow ( ( double ) 10, ( double ) length );

  if ( i < 0 )
  {
    length = length + 1;
  }
/*
  Add one position for the trailing null.
*/
  length = length + 1;

  s = ( char * ) malloc ( length * sizeof ( char ) );

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
/*
  Now take care of the sign.
*/
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
/*
  Find the leading digit of I, strip it off, and stick it into the string.
*/
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
/*
  Tack on the trailing NULL.
*/
  s[j] = '\0';
  j = j + 1;

  return s;
}
/******************************************************************************/

int i4_uniform_0i ( int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_0I is a portable random integer generator.

  Formula:

    SEED = SEED * (7^5) mod (2^31 - 1)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int *SEED, the integer "seed" used to generate
    the output value.  SEED should not be 0.  On output, SEED
    has been updated.

    Output, int I4_UNIFORM_0I, a uniform random value between
    1 and 2^31-1.

  Local parameters:

    IA = 7^5
    IB = 2^15
    IB16 = 2^16
    IP = 2^31-1
*/
{
  const int i4_huge = 2147483647;
  int ia = 16807;
  int ib15 = 32768;
  int ib16 = 65536;
  int ip = 2147483647;
  int iprhi;
  int ixhi;
  int k;
  int leftlo;
  int loxa;
  int value;
/*
  Don't let SEED be 0.
*/
  if ( *seed == 0 )
  {
    *seed = i4_huge;
  }
/*
  Get the 15 high order bits of SEED2.
*/
  ixhi = *seed / ib16;
/*
  Get the 16 low bits of SEED and form the low product.
*/
  loxa = ( *seed - ixhi * ib16 ) * ia;
/*
  Get the 15 high order bits of the low product.
*/
  leftlo = loxa / ib16;
/*
  Form the 31 highest bits of the full product.
*/
  iprhi = ixhi * ia + leftlo;
/*
  Get overflow past the 31st bit of full product.
*/
  k = iprhi / ib15;
/*
  Assemble all the parts and presubtract IP.  The parentheses are
  essential.
*/
  value = ( ( ( loxa - leftlo * ib16 ) - ip ) 
            + ( iprhi - k * ib15 ) * ib16 ) + k;
/*
  Add IP back in if necessary.
*/
  if ( value < 0 )
  {
    value = value + i4_huge;
  }
  *seed = value;

  return value;
}
/******************************************************************************/

int i4_uniform_ab ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM_AB, a number between A and B.
*/
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( a ) - 0.5 ) 
    +         r   * ( ( float ) ( b ) + 0.5 );
/*
  Round R to the nearest integer.
*/
  value = round ( r );
/*
  Guarantee that A <= VALUE <= B.
*/
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
/******************************************************************************/

void i4mat_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT prints an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT_SOME prints some of an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:" );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %6d", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to INCX) entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %6d", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void i4mat_uniform_ab ( int m, int n, int a, int b, int *seed, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_UNIFORM_AB returns a scaled pseudorandom I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, int A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, int X[M*N], a matrix of pseudorandom values.
*/
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
      r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
        +         r   * ( ( float ) b + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
      value = round ( r );
/*
  Guarantee A <= VALUE <= B.
*/
      if ( value < a )
      {
        value = a;
      }
      if ( b < value )
      {
        value = b;
      }

      x[i+j*m] = value;
    }
  }

  return;
}
/******************************************************************************/

int *i4mat_uniform_ab_new ( int m, int n, int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4MAT_UNIFORM_AB_NEW returns a scaled pseudorandom I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, int A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, int I4MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  float r;
  int value;
  int *x;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  x = ( int * ) malloc ( m * n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
      r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
        +         r   * ( ( float ) b + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
      value = round ( r );
/*
  Guarantee A <= VALUE <= B.
*/
      if ( value < a )
      {
        value = a;
      }
      if ( b < value )
      {
        value = b;
      }

      x[i+j*m] = value;
    }
  }

  return x;
}
/******************************************************************************/

int i4vec_max ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MAX returns the value of the maximum element in an I4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int IVEC_MAX, the value of the maximum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;  
}
/******************************************************************************/

float i4vec_mean ( int n, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MEAN returns the mean of an I4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 May 1999

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int X[N], the vector whose mean is desired.

    Output, float I4VEC_MEAN, the mean, or average, of the vector entries.
*/
{
  int i;
  float mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + ( float ) x[i];
  }

  mean = mean / ( float ) n;

  return mean;
}
/******************************************************************************/

int i4vec_min ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MIN returns the minimum element in an I4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int A[N], the array to be checked.

    Output, int I4VEC_MIN, the value of the minimum element.  This
    is set to 0 if N <= 0.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value; 
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

void i4vec_uniform_ab ( int n, int a, int b, int *seed, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.

  Discussion:

    The pseudorandom numbers should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, integer N, the dimension of the vector.

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int X[N], a vector of random values between A and B.
*/
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
      +         r   * ( ( float ) b + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
    value = round ( r );
/*
  Guarantee A <= VALUE <= B.
*/
    if ( value < a )
    {
      value = a;
    }
    if ( b < value )
    {
      value = b;
    }

    x[i] = value;
  }

  return;
}
/******************************************************************************/

int *i4vec_uniform_ab_new ( int n, int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.

  Discussion:

    The pseudorandom numbers should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, integer N, the dimension of the vector.

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4VEC_UNIFORM_AB_NEW[N], a vector of random values between A and B.
*/
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;
  int *x;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  x = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
      +         r   * ( ( float ) b + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
    value = round ( r );
/*
  Guarantee A <= VALUE <= B.
*/
    if ( value < a )
    {
      value = a;
    }
    if ( b < value )
    {
      value = b;
    }

    x[i] = value;
  }

  return x;
}
/******************************************************************************/

float i4vec_variance ( int n, int x[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_VARIANCE returns the variance of an I4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 May 1999

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int X[N], the vector whose variance is desired.

    Output, float I4VEC_VARIANCE, the variance of the vector entries.
*/
{
  int i;
  float mean;
  float variance;

  mean = i4vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( ( float ) x[i] - mean ) * ( ( float ) x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( float ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
/******************************************************************************/

int l4_uniform ( int *seed )

/******************************************************************************/
/*
  Purpose:

    L4_UNIFORM returns a pseudorandom L4.

  Discussion:

    An L4 is a LOGICAL value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value, which should
    NOT be 0.  On output, SEED has been updated.

    Output, int L4_UNIFORM, a pseudorandom logical value.
*/
{
  const int i4_huge      = 2147483647;
  const int i4_huge_half = 1073741823;
  int  k;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L4_UNIFORM - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 ) 
  {
    *seed = *seed + i4_huge;
  }
  value = ( i4_huge_half < *seed );

  return value;
}
/******************************************************************************/

void l4mat_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    L4MAT_PRINT prints an L4MAT.

  Discussion:

    An L4MAT is an array of L4 values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the matrix.

    Input, char *TITLE, a title.
*/
{
  l4mat_print_some ( m, n, a, 0, 0, m - 1, n - 1, title );

  return;
}
/******************************************************************************/

void l4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    L4MAT_PRINT_SOME prints some of an L4MAT.

  Discussion:

    An L4MAT is an array of L4 values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, int A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 35;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );

  for ( j2lo = i4_max ( jlo, 0 ); j2lo <= i4_min ( jhi, n - 1 ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    if ( n - 1 < j2hi )
    {
      j2hi = n - 1;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );

    if ( 100 <= j2hi )
    {
      printf ( "      " );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( " %1d", j / 100 );
      }
      printf ( "\n" );
    }

    if ( 10 <= j2hi )
    {
      printf ( "      " );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( " %1d", ( ( j / 10 ) % 10 ) );
      }
      printf ( "\n" );
    }

    printf ( "  Col " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( " %1d", ( j % 10 ) );
    }
    printf ( "\n" );

    printf ( "  Row\n" );
    printf ( "\n" );

    i2lo = 0;
    if ( i2lo < ilo )
    {
      i2lo = ilo;
    }
    i2hi = m - 1;
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
      printf ( "%5d:", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( " %1d", a[i+j*m] );
      }
      printf ( "\n" );
    }
  }
  return;
}
/******************************************************************************/

int *l4mat_uniform_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    L4MAT_UNIFORM_NEW returns a pseudorandom L4MAT.

  Discussion:

    An LMAT is a two dimensional array of LOGICAL values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the order of the matrix.

    Input/output, int *SEED, the "seed" value, which should
    NOT be 0.  On output, SEED has been updated.

    Output, int L4MAT_UNIFORM_NEW[M*N], a pseudorandom logical matrix.
*/
{
  const int i4_huge      = 2147483647;
  const int i4_huge_half = 1073741823;
  int i;
  int j;
  int k;
  int *l4mat;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L4MAT_UNIFORM_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  l4mat = ( int * ) malloc ( m * n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      l4mat[i+j*m] = ( i4_huge_half < *seed );
    }
  }

  return l4mat;
}
/******************************************************************************/

void l4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    L4VEC_PRINT prints an L4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the (logical) vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  for ( i = 0; i < n; i++ ) 
  {
    if ( a[i] == 0 )
    {
      printf ( "  %8d: F\n", i );
    }
    else
    {
      printf ( "  %8d: T\n", i );
    }
  }

  return;
}
/******************************************************************************/

int *l4vec_uniform_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    L4VEC_UNIFORM_NEW returns a pseudorandom L4VEC.

  Discussion:

    An L4VEC is a vector of LOGICAL values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre LEcuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the order of the vector.

    Input/output, int *SEED, the "seed" value, which should
    NOT be 0.  On output, SEED has been updated.

    Output, int L4VEC_UNIFORM_NEW[N], a pseudorandom logical vector.
/*/
{
  const int i4_huge      = 2147483647;
  const int i4_huge_half = 1073741823;
  int i;
  int k;
  int *l4vec;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L4VEC_UNIFORM_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    return NULL;
  }

  l4vec = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }
    l4vec[i] = ( i4_huge_half < *seed );
  }
  return l4vec;
}
/******************************************************************************/

void lcrg_anbn ( int a, int b, int c, int n, int *an, int *bn )

/******************************************************************************/
/*
  Purpose:

    LCRG_ANBN computes the "N-th power" of a linear congruential generator.

  Discussion:

    We are considering a linear congruential random number generator.
    The LCRG takes as input an integer value called SEED, and returns
    an updated value of SEED,

      SEED(out) = ( a * SEED(in) + b ) mod c.

    and an associated pseudorandom real value

      U = SEED(out) / c.

    In most cases, a user is content to call the LCRG repeatedly, with
    the updating of SEED being taken care of automatically.

    The purpose of this routine is to determine the values of AN and BN
    that describe the LCRG that is equivalent to N applications of the
    original LCRG.

    One use for such a facility would be to do random number computations
    in parallel.  If each of N processors is to compute many random values,
    you can guarantee that they work with distinct random values
    by starting with a single value of SEED, using the original LCRG to generate
    the first N-1 "iterates" of SEED, so that you now have N "seed" values,
    and from now on, applying the N-th power of the LCRG to the seeds.

    If the K-th processor starts from the K-th seed, it will essentially
    be computing every N-th entry of the original random number sequence,
    offset by K.  Thus the individual processors will be using a random
    number stream as good as the original one, and without repeating, and
    without having to communicate.

    To evaluate the N-th value of SEED directly, we start by ignoring
    the modular arithmetic, and working out the sequence of calculations
    as follows:

      SEED(0)   =     SEED.
      SEED(1)   = a * SEED      + b
      SEED(2)   = a * SEED(1)   + b = a^2 * SEED           + a * b + b
      SEED(3)   = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
      ...
      SEED(N-1) = a * SEED(N-2) + b

      SEED(N) = a * SEED(N-1) + b = a^N * SEED
                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b

    or, using the geometric series,

      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
              = AN * SEED + BN

    Thus, from any SEED, we can determine the result of N applications of the
    original LCRG directly if we can solve

      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,

    and evaluate:

      AN = a^N

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Barry Wilkinson, Michael Allen,
    Parallel Programming:
    Techniques and Applications Using Networked Workstations and Parallel Computers,
    Prentice Hall,
    ISBN: 0-13-140563-2,
    LC: QA76.642.W54.

  Parameters:

    Input, int A, the multiplier for the LCRG.

    Input, int  B, the added value for the LCRG.

    Input, int  C, the base for the modular arithmetic.
    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
    required that 0 < C.

    Input, int N, the "index", or number of times that the
    LCRG is to be applied.  It is required that 0 <= N.

    Output, int *AN, *BN, the multiplier and added value for
    the LCRG that represent N applications of the original LCRG.
*/
{
  int am1;
  int anm1tb;
  int ierror;

  if ( n < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LCRG_ANBN - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of N = %d\n", n );
    exit ( 1 );
  }

  if ( c <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LCRG_ANBN - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of C = %d\n", c );
    exit ( 1 );
  }

  if ( n == 0 )
  {
    *an = 1;
    *bn = 0;
  }
  else if ( n == 1 )
  {
    *an = a;
    *bn = b;
  }
  else
  {
/*
  Compute A^N.
*/
    *an = power_mod ( a, n, c );
/*
  Solve 
    ( a - 1 ) * BN = ( a^N - 1 ) mod B
  for BN.
*/
    am1 = a - 1;
    anm1tb = ( *an - 1 ) * b;

    *bn = congruence ( am1, c, anm1tb, &ierror );

    if ( ierror )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "LCRG_ANBN - Fatal error!\n" );
      fprintf ( stderr, "  An error occurred in the CONGRUENCE routine.\n" );
      exit ( 1 );
    }
  }

  return;
}
/******************************************************************************/

int lcrg_evaluate ( int a, int b, int c, int x )

/******************************************************************************/
/*
  Purpose:

    LCRG_EVALUATE evaluates an LCRG, y = ( A * x + B ) mod C.

  Discussion:

    This routine cannot be recommended for production use.  Because we want
    to do modular arithmetic, but the base is not a power of 2, we need to
    use "double precision" integers to keep accuracy.

    If we knew the base C, we could try to avoid overflow while not changing
    precision.

    If the base C was a power of 2, we could rely on the usual properties of
    integer arithmetic on computers, in which overflow bits, which are always
    ignored, don ot actually matter.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int A, the multiplier for the LCRG.

    Input, int B, the added value for the LCRG.

    Input, int C, the base for the modular arithmetic.
    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
    required that 0 < C.

    Input, int X, the value to be processed.

    Output, int LCRG_EVALUATE, the processed value.
*/
{
  long long int a8;
  long long int b8;
  long long int c8;
  long long int x8;
  int y;
  long long int y8;
/*
  To avoid roundoff issues, we need to go to "double precision" integers.
  (Not available on all planets.)
*/
  a8 = ( long long int ) a;
  b8 = ( long long int ) b;
  c8 = ( long long int ) c;
  x8 = ( long long int ) x;

  y8 = ( a8 * x8 + b8 ) % c8;

  y = ( int ) ( y8 );

  if ( y < 0 )
  {
    y = y + c;
  }

  return y;
}
/******************************************************************************/

int lcrg_seed ( int a, int b, int c, int n, int seed )

/******************************************************************************/
/*
  Purpose:

    LCRG_SEED computes the N-th seed of a linear congruential generator.

  Discussion:

    We are considering a linear congruential random number generator.
    The LCRG takes as input an integer value called SEED, and returns
    an updated value of SEED,

      SEED(out) = a * SEED(in) + b, mod c.

    and an associated pseudorandom real value

      U = SEED(out) / c.

    In most cases, a user is content to call the LCRG repeatedly, with
    the updating of SEED being taken care of automatically.

    The purpose of this routine is to determine the value of SEED that
    would be output after N successive applications of the LCRG.  This
    allows the user to know, in advance, what the 1000-th value of
    SEED would be, for instance.  Obviously, one way to do this is to
    apply the LCRG formula 1,000 times.  However, it is possible to
    do this in a more direct and efficient way.

    One use for such a facility would be to do random number computations
    in parallel.  If each processor is to compute 1,000 values, you can
    guarantee that they work with distinct random values by starting the
    first processor with SEED, the second with the value of SEED after
    1,000 applications of the LCRG, and so on.

    To evaluate the N-th value of SEED directly, we start by ignoring
    the modular arithmetic, and working out the sequence of calculations
    as follows:

      SEED(0) =     SEED.
      SEED(1) = a * SEED      + b
      SEED(2) = a * SEED(1)   + b = a^2 * SEED + a * b + b
      SEED(3) = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
      ...
      SEED(N) = a * SEED(N-1) + b = a^N * SEED
                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b

    or, using the geometric series,

      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b

    Therefore, we can determine SEED(N) directly if we can solve

      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,

    and evaluated:

      AN = a^N

    Using the formula:

      SEED(N) = AN * SEED + BN, mod c

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, int A, the multiplier for the LCRG.

    Input, int B, the added value for the LCRG.

    Input, int C, the base for the modular arithmetic.  For 32 bit
    arithmetic, this is often 2^31 - 1, or 2147483647.  It is required
    that 0 < C.

    Input, int N, the "index", or number of times that the LCRG
    is to be applied.  It is required that 0 <= N.

    Input, int SEED, the starting value of SEED.  It is customary
    that 0 < SEED.

    Output, int LCRG_SEED, the value of SEED that would be output
    if the LCRG were applied to the starting value N times.
*/
{
  int an;
  int bn;
  int ierror;
  int value;
  long int value2;

  if ( n < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LCRG_SEED - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of N = %d\n",  n );
    exit ( 1 );
  }

  if ( c <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LCRG_SEED - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of C = %d\n", c );
    exit ( 1 );
  }

  if ( n == 0 )
  {
    value = seed % c;
    if ( value < 0 )
    {
      value = value + c;
    }
    return value;
  }
/*
  Get A^N.
*/
  an = power_mod ( a, n, c );
/*
  Solve ( a - 1 ) * BN = ( a^N - 1 ) for BN.

  The LCRG I have been investigating uses B = 0, so this code
  has not been properly tested yet.
*/
  bn = congruence ( a - 1, c, ( an - 1 ) * b, &ierror );

  if ( ierror != 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LCRG_SEED - Fatal error!\n" );
    fprintf ( stderr, "  An error occurred in the CONGRUENCE routine.\n" );
    fprintf ( stderr, "  The error code was IERROR = %d\n", ierror );
    exit ( 1 );
  }
/*
  Set the new SEED.
*/
  value2 = ( long int ) ( an ) * ( long int ) ( seed ) + ( long int ) ( bn );

  value2 = value2 % ( long int ) ( c );
/*
  Guarantee that the value is positive.
*/
  if ( value2 < 0 )
  {
    value2 = value2 + ( long int ) ( c );
  }

  value = ( int ) ( value2 );

  return value;
}
/******************************************************************************/

int power_mod ( int a, int n, int m )

/******************************************************************************/
/*
  Purpose:

    POWER_MOD computes mod ( A^N, M ).

  Discussion:

    Some programming tricks are used to speed up the computation, and to
    allow computations in which A**N is much too large to store in a
    real word.

    First, for efficiency, the power A**N is computed by determining
    the binary expansion of N, then computing A, A^2, A^4, and so on
    by repeated squaring, and multiplying only those factors that
    contribute to A**N.

    Secondly, the intermediate products are immediately "mod'ed", which
    keeps them small.

    For instance, to compute mod ( A^13, 11 ), we essentially compute

       13 = 1 + 4 + 8

       A**13 = A * A^4 * A^8

       mod ( A**13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).

    Fermat's little theorem says that if P is prime, and A is not divisible
    by P, then ( A^(P-1) - 1 ) is divisible by P.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, int A, the base of the expression to be tested.
    A should be nonnegative.

    Input, int N, the power to which the base is raised.
    N should be nonnegative.

    Input, int M, the divisor against which the expression is tested.
    M should be positive.

    Output, int POWER_MOD, the remainder when A**N is divided by M.
*/
{
  long long int a_square2;
  int d;
  long long int m2;
  int x;
  long long int x2;

  if ( a < 0 )
  {
    return -1;
  }

  if ( m <= 0 )
  {
    return -1;
  }

  if ( n < 0 )
  {
    return -1;
  }
/*
  A_SQUARE contains the successive squares of A.
*/
  a_square2 = ( long long int ) a;
  m2 = ( long long int ) m;
  x2 = ( long long int ) 1;

  while ( 0 < n )
  {
    d = n % 2;

    if ( d == 1 )
    {
      x2 = ( x2 * a_square2 ) % m2;
    }

    a_square2 = ( a_square2 * a_square2 ) % m2;
    n = ( n - d ) / 2;
  }
/*
  Ensure that 0 <= X.
*/
  while ( x2 < 0 )
  {
    x2 = x2 + m2;
  }

  x = ( int ) x2;

  return x;
}
/******************************************************************************/

int r4_nint ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_NINT returns the nearest integer to an R4.

  Example:

        X         R4_NINT

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

    Input, float X, the value.

    Output, int R4_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

float r4_uniform_ab ( float a, float b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_AB returns a scaled pseudorandom R4.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 April 2011

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, float A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float R4_UNIFORM_AB, a number strictly between A and B.
*/
{
  const int i4_huge = 2147483647;
  int k;
  float value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  value = ( float ) ( *seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
/******************************************************************************/

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a unit pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r4_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R4_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  const int i4_huge = 2147483647;
  int k;
  float value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  value = ( float ) ( *seed ) * 4.656612875E-10;

  return value;
}
/******************************************************************************/

void r4mat_print ( int m, int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_PRINT prints an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R4's, which
    may be stored as a vector in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, float A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_PRINT_SOME prints some of an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R4's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, float A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r4mat_uniform_01 ( int m, int n, int *seed, float r[] )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, float R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

float *r4mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_01_NEW returns a unit pseudorandom R4MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, float R4MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( float * ) malloc ( m * n * sizeof ( float ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r4mat_uniform_ab ( int m, int n, float b, float c, int *seed, float r[] )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_AB returns a scaled pseudorandom R4MAT.

  Discussion:

    This routine implements the recursion

      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
      u = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float B, C, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, float R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

float *r4mat_uniform_ab_new ( int m, int n, float b, float c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_AB_NEW returns a scaled pseudorandom R4MAT.

  Discussion:

    This routine implements the recursion

      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
      u = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float B, C, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, float R4MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( float * ) malloc ( m * n * sizeof ( float ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r4vec_print ( int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R4VEC_PRINT prints an R4VEC.

  Discussion:

    An R4VEC is a vector of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, float A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r4vec_uniform_ab ( int n, float b, float c, int *seed, float r[] )

/******************************************************************************/
/*
  Purpose:

    R4VEC_UNIFORM_AB returns a scaled pseudorandom R4VEC.

  Discussion:

    This routine implements the recursion

      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
      u = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, float B, C, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4VEC_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

float *r4vec_uniform_ab_new ( int n, float b, float c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R4VEC.

  Discussion:

    This routine implements the recursion

      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
      u = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 April 2008

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, float B, C, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R4VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( float * ) malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( float ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r4vec_uniform_01 ( int n, int *seed, float r[] )

/******************************************************************************/
/*
  Purpose:

    R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( float ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

float *r4vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4VEC_UNIFORM_01_NEW returns a unit pseudorandom R4VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, float R4VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( float * ) malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( float ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

int r8_nint ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_NINT returns the nearest I4 to an R8.

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

    26 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value.

    Output, int R8_NINT, the nearest integer to X.
*/
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
/******************************************************************************/

double r8_uniform_ab ( double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_AB returns a scaled pseudorandom R8.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 April 2011

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, double A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double R8_UNIFORM_AB, a number strictly between A and B.
*/
{
  const int i4_huge = 2147483647;
  int k;
  double value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  value = ( double ) ( *seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  const int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

double *r8col_uniform_abvec_new ( int m, int n, double a[], double b[], 
  int *seed )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIFORM_ABVEC_NEW fills an R8COL with scaled pseudorandom numbers.

  Discussion:

    An R8COL is an array of R8 values, regarded as a set of column vectors.

    The user specifies a minimum and maximum value for each row.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 December 2014

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M], B[M], the upper and lower limits.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8COL_UNIFORM_ABVEC_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = a[i] 
        + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14g", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_uniform_ab_new ( int m, int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB_NEW returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R8MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

double *r8row_uniform_abvec_new ( int m, int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8ROW_UNIFORM_ABVEC_NEW fills an R8ROW with scaled pseudorandom numbers.

  Discussion:

    An R8ROW is an array of R8 values, regarded as a set of row vectors.

    The user specifies a minimum and maximum value for each column.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[N], B[N], the upper and lower limits.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8ROW_UNIFORM_ABVEC_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = a[j] 
        + ( b[j] - a[j] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

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

double *r8vec_normal_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.

    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.

    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.

  Local parameters:

    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.

    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.

    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.

    Local, float Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
# define R8_PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  I'd like to allow the user to reset the internal data.
  But this won't work properly if we have a saved value Y.
  I'm making a crock option that allows the user to signal
  explicitly that any internal memory should be flushed,
  by passing in a negative value for N.
*/
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  Use up the old value, if we have it.
*/
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
/*
  Maybe we don't need any more values.
*/
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * R8_PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * R8_PI * r[1] );

    saved = 1;

    made = made + 2;

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    free ( r );
  }

  return x;
# undef R8_PI
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec_uniform_01 ( int n, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r8vec_uniform_ab ( int n, double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_ab_new ( int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r8vec_uniform_abvec ( int n, double a[], double b[], int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.

  Discussion:

    Dimension I ranges from A[I] to B[I].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], B[N], the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_ABVEC - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_abvec_new ( int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_ABVEC_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Dimension I ranges from A[I] to B[I].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], B[N], the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_ABVEC_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_ABVEC_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *r8vec_uniform_unit_new ( int m, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_UNIT_NEW generates a random unit vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_UNIT_NEW[M], a random direction 
//    vector, with unit norm.
//
{
  double *a;
  int i;
  double norm;
//
//  Take M random samples from the normal distribution.
//
  a = r8vec_normal_01_new ( m, seed );
//
//  Compute the norm.
//
  norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );
//
//  Normalize.
//
  for ( i = 0; i < m; i++ )
  {
    a[i] = a[i] / norm;
  }

  return a;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    May 31 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 October 2003

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
