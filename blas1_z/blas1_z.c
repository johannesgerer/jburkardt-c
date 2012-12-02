# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <complex.h>

# include "blas1_z.h"

/******************************************************************************/

double dzasum ( int n, _Complex double x[], int incx )

/******************************************************************************/
/*
  Purpose:

    DZASUM takes the sum of the absolute values of a vector.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, _Complex double X[], the vector.

    Input, int INCX, the increment between successive entries of X.

    Output, double DZASUM, the sum of the absolute values.
*/
{
  int i;
  int ix;
  double value;

  value = 0.0;

  if ( n <= 0 || incx <= 0 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + r8_abs ( creal ( x[i] ) ) 
                    + r8_abs ( cimagf ( x[i] ) );
    }
  }
  else
  {
    ix = 0;
    for ( i = 0; i < n; i++ )
    {
      value = value + r8_abs ( creal ( x[ix] ) ) 
                    + r8_abs ( cimagf ( x[ix] ) );
      ix = ix + incx;
    }
  }
  return value;
}
/******************************************************************************/

double dznrm2 ( int n, _Complex double x[], int incx )

/******************************************************************************/
/*
  Purpose:

    DZNRM2 returns the euclidean norm of a vector.

  Discussion:

    This routine uses double precision complex arithmetic.

    DZNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
            = sqrt ( dot_product ( x(1:n), x(1:n) ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, _Complex double X[], the vector.

    Input, int INCX, the increment between successive entries of X.

    Output, double DZNRM2, the norm of the vector.
*/
{
  int i;
  int ix;
  double scale;
  double ssq;
  double temp;
  double value;

  if ( n < 1 || incx < 1 )
  {
    value  = 0.0;
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( creal ( x[ix] ) != 0.0 )
      {
        temp = r8_abs ( creal ( x[ix] ) );
        if ( scale < temp )
        {
          ssq = 1.0 + ssq * pow ( scale / temp, 2 );
          scale = temp;
        }
        else
        {
          ssq = ssq + pow ( temp / scale, 2 );
        }
      }

      if ( cimagf ( x[ix] ) != 0.0 )
      {
        temp = r8_abs ( cimagf ( x[ix] ) );
        if ( scale < temp )
        {
          ssq = 1.0 + ssq * pow ( scale / temp, 2 );
          scale = temp;
        }
        else
        {
          ssq = ssq + pow ( temp / scale, 2 );
        }
      }
      ix = ix + incx;
    }
    value  = scale * sqrt ( ssq );
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

int izamax ( int n, _Complex double x[], int incx )

/******************************************************************************/
/*
  Purpose:

    IZAMAX indexes the vector element of maximum absolute value.

  Discussion:

    This routine uses double precision complex arithmetic.

    WARNING: This index is a 1-based index, not a 0-based index!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, _Complex double X[], the vector.

    Input, int INCX, the increment between successive entries of X.

    Output, int IZAMAX, the index of the element of maximum
    absolute value.
*/
{
  int i;
  int ix;
  double smax;
  int value;

  value = 0;

  if ( n < 1 || incx  <=  0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx != 1 )
  {
    ix = 0;
    smax = zabs1 ( x[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( smax < zabs1 ( x[ix] ) )
      {
        value = i + 1;
        smax = zabs1 ( x[ix] );
      }
      ix = ix + incx;
    }
  }
  else
  {
    smax = zabs1 ( x[0] );
    for ( i = 1; i < n; i++ )
    {
      if ( smax < zabs1 ( x[i] ) )
      {
        value = i + 1;
        smax = zabs1 ( x[i] );
      }
    }
  }

  return value;
}
/******************************************************************************/

int lsame ( char ca, char cb )

/******************************************************************************/
/*
  Purpose:

    LSAME returns TRUE if CA is the same letter as CB regardless of case.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, char CA, CB, the characters to compare.

    Output, int LSAME, is 1 if the characters are equal,
    disregarding case.
*/
{
  if ( ca == cb )
  {
    return 1;
  }

  if ( 'A' <= ca && ca <= 'Z' )
  {
    if ( ca - 'A' == cb - 'a' )
    {
      return 1;
    }    
  }
  else if ( 'a' <= ca && ca <= 'z' )
  {
    if ( ca - 'a' == cb - 'A' )
    {
      return 1;
    }
  }

  return 0;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of a number.

  Discussion:

    This routine uses double precision arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two numbers.

  Discussion:

    This routine uses double precision arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of a number.

  Discussion:

    This routine uses double precision arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
}
/******************************************************************************/

void xerbla ( char *srname, int info )

/******************************************************************************/
/*
  Purpose:

    XERBLA is an error handler for the LAPACK routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, char *SRNAME, the name of the routine
    which called XERBLA.

    Input, int INFO, the position of the invalid parameter in
    the parameter list of the calling routine.
*/
{
  printf ( "\n" );
  printf ( "XERBLA - Fatal error!\n" );
  printf ( "  On entry to routine %s\n", srname );
  printf ( "  input parameter number %d had an illegal value.\n", info );
  exit ( 1 );
}
/******************************************************************************/

double zabs1 ( _Complex double z )

/******************************************************************************/
/*
  Purpose:

    ZABS1 returns the L1 norm of a number.

  Discussion:

    This routine uses double precision complex arithmetic.

    The L1 norm of a complex number is the sum of the absolute values
    of the real and imaginary components.

    ZABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, _Complex double Z, the number whose norm is desired.

    Output, double ZABS1, the L1 norm of Z.
*/
{
  double value;

  value = r8_abs ( creal ( z ) ) + r8_abs ( cimagf ( z ) );

  return value;
}
/******************************************************************************/

double zabs2 ( _Complex double z )

/******************************************************************************/
/*
  Purpose:

    ZABS2 returns the L2 norm of a number.

  Discussion:

    This routine uses double precision complex arithmetic.

    The L2 norm of a complex number is the square root of the sum 
    of the squares of the real and imaginary components.

    ZABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, _Complex double Z, the number whose norm is desired.

    Output, float ZABS2, the L2 norm of Z.
*/
{
  double value;

  value = sqrt ( pow ( creal ( z ), 2 ) 
               + pow ( cimag ( z ), 2 ) );

  return value;
}
/******************************************************************************/

void zaxpy ( int n, _Complex double ca, _Complex double cx[], 
  int incx, _Complex double cy[], int incy )

/******************************************************************************/
/*
  Purpose:

    ZAXPY computes a constant times a vector plus a vector.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of elements in CX and CY.

    Input, _Complex double CA, the multiplier of CX.

    Input, _Complex double CX[], the first vector.

    Input, int INCX, the increment between successive entries of CX.

    Input/output, _Complex double CY[], the second vector.
    On output, CY(*) has been replaced by CY(*) + CA * CX(*).

    Input, int INCY, the increment between successive entries of CY.
*/
{
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( zabs1 ( ca ) == 0.0 ) 
  {
    return;
  }

  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      cy[iy] = cy[iy] + ca * cx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cy[i] = cy[i] + ca * cx[i];
    }

  }

  return;
}
/******************************************************************************/

void zcopy ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy )

/******************************************************************************/
/*
  Purpose:

    ZCOPY copies a vector X to a vector Y.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of elements in CX and CY.

    Input, _Complex double CX[], the first vector.

    Input, int INCX, the increment between successive entries of CX.

    Output, _Complex double CY[], the second vector.

    Input, int INCY, the increment between successive entries of CY.
*/
{
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      cy[iy] = cx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cy[i] = cx[i];
    }
  }
  return;
}
/******************************************************************************/

_Complex double zdotc ( int n, _Complex double cx[], int incx, 
  _Complex double cy[], int incy )

/******************************************************************************/
/*
  Purpose:

    ZDOTC forms the conjugated dot product of two vectors.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, _Complex double CX[], the first vector.

    Input, int INCX, the increment between successive entries in CX.

    Input, _Complex double CY[], the second vector.

    Input, int INCY, the increment between successive entries in CY.

    Output, _Complex double ZDOTC, the conjugated dot product of
    the corresponding entries of CX and CY.
*/
{
  int i;
  int ix;
  int iy;
  _Complex double value;

  value = 0.0;

  if ( n <= 0 )
  {
    return value;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + ( ~cx[i] ) * cy[i];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      value = value + ( ~cx[ix] ) * cy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return value;
}
/******************************************************************************/

_Complex double zdotu ( int n, _Complex double cx[], int incx, 
  _Complex double cy[], int incy )

/******************************************************************************/
/*
  Purpose:

    ZDOTU forms the unconjugated dot product of two vectors.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, _Complex double CX[], the first vector.

    Input, int INCX, the increment between successive entries in CX.

    Input, _Complex double CY[], the second vector.

    Input, int INCY, the increment between successive entries in CY.

    Output, _Complex double ZDOTU, the unconjugated dot product of
    the corresponding entries of CX and CY.
*/
{
  int i;
  int ix;
  int iy;
  _Complex double value;

  value = 0.0;

  if ( n <= 0 )
  {
    return value;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + cx[i] * cy[i];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      value = value + cx[ix] * cy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return value;
}
/******************************************************************************/

void zdrot ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy, double c, double s )

/******************************************************************************/
/*
  Purpose:

    ZDROT applies a plane rotation.

  Discussion:

    This routine uses double precision complex arithmetic.

    The cosine C and sine S are real and the vectors CX and CY are complex.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, _Complex double CX[], one of the vectors to be rotated.

    Input, int INCX, the increment between successive entries of CX.

    Input/output, _Complex double CY[], one of the vectors to be rotated.

    Input, int INCY, the increment between successive elements of CY.

    Input, double C, S, parameters (presumably the cosine and sine of
    some angle) that define a plane rotation.
*/
{
  _Complex double ctemp;
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      ctemp = c * cx[i] + s * cy[i];
      cy[i] = c * cy[i] - s * cx[i];
      cx[i] = ctemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      ctemp  = c * cx[ix] + s * cy[iy];
      cy[iy] = c * cy[iy] - s * cx[ix];
      cx[ix] = ctemp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return;
}
/******************************************************************************/

void zdscal ( int n, double sa, _Complex double cx[], int incx )

/******************************************************************************/
/*
  Purpose:

    ZDSCAL scales a vector by a constant.

  Discussion:

    This routine uses double precision complex arithmetic.

    The scaling constant is double precision real.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double SA, the multiplier.

    Input/output, _Complex double CX[], the vector to be scaled.

    Input, int INCX, the increment between successive entries of
    the vector CX.
*/
{
  int i;

  if ( n <= 0 || incx <= 0 )
  {
    return;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i] = sa * cx[i];
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i*incx] = sa * cx[i*incx];
    }
  }
  return;
}
/******************************************************************************/

double zmach ( int job )

/******************************************************************************/
/*
  Purpose:

    ZMACH computes machine parameters for double complex arithmetic.

  Discussion:

    Assume the computer has

      B = base of arithmetic;
      T = number of base B digits;
      L = smallest possible exponent;
      U = largest possible exponent;

    then

      EPS = B**(1-T)
      TINY = 100.0 * B**(-L+T)
      HUGE = 0.01 * B**(U-T)

    If complex division is done by

      1 / (X+i*Y) = (X-i*Y) / (X**2+Y**2)

    then

      TINY = sqrt ( TINY )
      HUGE = sqrt ( HUGE )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input, int JOB:
    1, EPS is desired;
    2, TINY is desired;
    3, HUGE is desired.

    Output, double ZMACH, the requested value.
*/
{
  double eps;
  double huge;
  double s;
  _Complex double temp1;
  _Complex double temp2;
  _Complex double temp3;
  double tiny;
  double value;

  eps = 1.0;

  for ( ; ; )
  {
    eps = eps / 2.0;
    s = 1.0 + eps;
    if ( s <= 1.0 )
    {
      break;
    }
  }

  eps = 2.0 * eps;

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }
  }

  tiny = ( tiny / eps ) * 100.0;
/*
  Had to insert this manually!
*/
  tiny = sqrt ( tiny );

  if ( 0 )
  {
    temp1 = 1.0; 
    temp2 = tiny;
    temp3 = temp1 / temp2;

    s = creal ( temp3 );

    if ( s != 1.0 / tiny )
    {
      tiny = sqrt ( tiny );
    }
  }

  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/******************************************************************************/

void zrotg ( _Complex double *ca, _Complex double cb, double *c, 
  _Complex double *s )

/******************************************************************************/
/*
  Purpose:

    ZROTG determines a Givens rotation.

  Discussion:

    This routine uses double precision complex arithmetic.

    Given values A and B, this routine computes:

    If A = 0:

      R = B
      C = 0
      S = (1,0).

    If A /= 0:

      ALPHA = A / abs ( A )
      NORM  = sqrt ( ( abs ( A ) )**2 + ( abs ( B ) )**2 )
      R     = ALPHA * NORM
      C     = abs ( A ) / NORM
      S     = ALPHA * conj ( B ) / NORM

    In either case, the computed numbers satisfy the equation:

    (         C    S ) * ( A ) = ( R )
    ( -conj ( S )  C )   ( B ) = ( 0 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for FORTRAN usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, pages 308-323, 1979.

  Parameters:

    Input/output, _Complex double *CA, on input, the value A.  On output,
    the value R.

    Input, _Complex double CB, the value B.

    Output, double *C, the cosine of the Givens rotation.

    Output, _Complex double *S, the sine of the Givens rotation.
*/
{
  _Complex double alpha;
  double norm;
  double scale;

  if ( zabs2 ( *ca ) == 0.0 )
  {
    *c = 0.0;
    *s = 1.0;
    *ca = cb;
  }
  else
  {
    scale = zabs2 ( *ca ) + zabs2 ( cb );
    norm = scale * sqrt ( pow ( zabs2 ( *ca / scale ), 2 )
                        + pow ( zabs2 (  cb / scale ), 2 ) );
    alpha = *ca / zabs2 ( *ca );
    *c = zabs2 ( *ca ) / norm;
    *s = alpha * ( ~cb ) / norm;
    *ca = alpha * norm;
  }

  return;
}
/******************************************************************************/

void zscal ( int n, _Complex double ca, _Complex double cx[], int incx )

/******************************************************************************/
/*
  Purpose:

    ZSCAL scales a vector by a constant.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt.

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, _Complex double CA, the multiplier.

    Input/output, _Complex double CX[], the vector to be scaled.

    Input, int INCX, the increment between successive entries of CX.
*/
{
  int i;

  if ( n <= 0 || incx <= 0 )
  {
    return;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i] = ca * cx[i];
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i*incx] = ca * cx[i*incx];
    }
  }
  return;
}
/******************************************************************************/

_Complex double zsign1 ( _Complex double z1, _Complex double z2 )

/******************************************************************************/
/*
  Purpose:

    ZSIGN1 is a transfer-of-sign function.

  Discussion:

    This routine uses double precision complex arithmetic.

    The L1 norm is used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Parameters:

    Input, _Complex double Z1, Z2, the arguments.

    Output, _Complex double ZSIGN1,  a complex value, with the magnitude of
    Z1, and the argument of Z2.
*/
{
  _Complex double value;

  if ( zabs1 ( z2 ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = zabs1 ( z1 ) * ( z2 / zabs1 ( z2 ) );
  }

  return value;
}
/******************************************************************************/

_Complex double zsign2 ( _Complex double z1, _Complex double z2 )

/******************************************************************************/
/*
  Purpose:

    ZSIGN2 is a transfer-of-sign function.

  Discussion:

    This routine uses double precision complex arithmetic.

    The L2 norm is used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt

  Parameters:

    Input, _Complex double Z1, Z2, the arguments.

    Output, _Complex double ZSIGN2,  a complex value, with the magnitude of
    Z1, and the argument of Z2.
*/
{
  _Complex double value;

  if ( zabs2 ( z2 ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = zabs2 ( z1 ) * ( z2 / zabs2 ( z2 ) );
  }

  return value;
}
/******************************************************************************/

void zswap ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy )

/******************************************************************************/
/*
  Purpose:

    ZSWAP interchanges two vectors.

  Discussion:

    This routine uses double precision complex arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, _Complex double CX[], one of the vectors to swap.

    Input, int INCX, the increment between successive entries of CX.

    Input/output, _Complex double CY[], one of the vectors to swap.

    Input, int INCY, the increment between successive elements of CY.
*/
{
  _Complex double ctemp;
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      ctemp = cx[i];
      cx[i] = cy[i];
      cy[i] = ctemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      ctemp = cx[ix];
      cx[ix] = cy[iy];
      cy[iy] = ctemp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }

  return;
}
