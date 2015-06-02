# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "sandia_rules.h"

/******************************************************************************/

void binary_vector_next ( int n, int bvec[] )

/******************************************************************************/
/*
  Purpose:

    BINARY_VECTOR_NEXT generates the next binary vector.

  Discussion:

    A binary vector is a vector whose entries are 0 or 1.

    The user inputs an initial zero vector to start.  The program returns
    the "next" vector.

    The vectors are produced in the order:

    ( 0, 0, 0, ..., 0 )
    ( 1, 0, 0, ..., 0 ) 
    ( 0, 1, 0, ..., 0 )
    ( 1, 1, 0, ..., 0 )
    ( 0, 0, 1, ..., 0 )
    ( 1, 0, 1, ..., 0 )
               ...
    ( 1, 1, 1, ..., 1)

    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
    we allow wrap around.

  Example:

    N = 3

    Input      Output
    -----      ------
    0 0 0  =>  1 0 0
    1 0 0  =>  0 1 0
    0 1 0  =>  1 1 0
    1 1 0  =>  0 0 1
    0 0 1  =>  1 0 1
    1 0 1  =>  0 1 1
    0 1 1  =>  1 1 1
    1 1 1  =>  0 0 0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 September 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input/output, int BVEC[N], on output, the successor 
    to the input vector.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {  
    if ( bvec[i] == 1 )
    {
      bvec[i] = 0;
    }
    else 
    {
      bvec[i] = 1;
      break;
    }
  }
  return;
}
/******************************************************************************/

void chebyshev1_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE computes a Chebyshev type 1 quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0 / sqrt ( 1 - x^2 ).

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV1_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = pi / ( double ) ( order );
  }
  for ( i = 0; i < order; i++ )
  {
    x[i] = cos ( pi * ( double ) ( 2 * order - 1 - 2 * i ) 
                         / ( double ) ( 2 * order ) );
  }
  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  return;
}
/******************************************************************************/

void chebyshev1_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE_NP computes a Chebyshev type 1 quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0 / sqrt ( 1 - x^2 ).

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  chebyshev1_compute ( order, x, w );

  return;
}
/******************************************************************************/

void chebyshev1_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE_POINTS computes Chebyshev type 1 quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.
*/
{
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV1_COMPUTE_POINTS - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    x[i] =  cos ( pi * ( double ) ( 2 * order - 1 - 2 * i ) 
                          / ( double ) ( 2 * order ) );
  }
  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  return;
}
/******************************************************************************/

void chebyshev1_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE_POINTS_NP computes Chebyshev type 1 quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  chebyshev1_compute_points ( order, x );

  return;
}
/******************************************************************************/

void chebyshev1_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE_WEIGHTS computes Chebyshev type 1 quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.
*/
{
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV1_COMPUTE_WEIGHTS - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = pi / ( double ) ( order );
  }

  return;
}
/******************************************************************************/

void chebyshev1_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_COMPUTE_WEIGHTS_NP: Chebyshev type 1 quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  chebyshev1_compute_weights ( order, w );

  return;
}
/******************************************************************************/

double chebyshev1_integral ( int expon )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.

  Discussion:

    To test a Chebyshev type 1 quadrature rule, we use it to approximate the
    integral of a monomial:

      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx

    This routine is given the value of the exponent, and returns the
    exact value of the integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.

    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
*/
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }
	
    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;	
  }

  return exact;
}
/******************************************************************************/

void chebyshev2_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE computes a Chebyshev type 2 quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = sqrt ( 1 - x^2 ).

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV2_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    angle = pi * ( double ) ( order - i ) / ( double ) ( order + 1 );
    w[i] = pi / ( double ) ( order + 1 ) * pow ( sin ( angle ), 2 );
    x[i] =  cos ( angle );
  }

  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  return;
}
/******************************************************************************/

void chebyshev2_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE_NP computes a Chebyshev type 2 quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = sqrt ( 1 - x^2 ).

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  chebyshev2_compute ( order, x, w );

  return;
}
/******************************************************************************/

void chebyshev2_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE_POINTS computes Chebyshev type 2 quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.
*/
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV2_COMPUTE_POINTS - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    angle = pi * ( double ) ( order - i ) / ( double ) ( order + 1 );
    x[i] =  cos ( angle );
  }

  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  return;
}
/******************************************************************************/

void chebyshev2_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE_POINTS_NP computes Chebyshev type 2 quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  chebyshev2_compute_points ( order, x );

  return;
}
/******************************************************************************/

void chebyshev2_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE_WEIGHTS computes Chebyshev type 2 quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double W[ORDER], the weights.
*/
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CHEBYSHEV2_COMPUTE_WEIGHTS - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    angle = pi * ( double ) ( order - i ) / ( double ) ( order + 1 );
    w[i] = pi / ( double ) ( order + 1 ) * pow ( sin ( angle ), 2 );
  }

  return;
}
/******************************************************************************/

void chebyshev2_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_COMPUTE_WEIGHTS_NP: Chebyshev type 2 quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  chebyshev2_compute_weights ( order, w );

  return;
}
/******************************************************************************/

double chebyshev2_integral ( int expon )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.

  Discussion:

    To test a Chebyshev type 2 quadrature rule, we use it to approximate the
    integral of a monomial:

      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx

    This routine is given the value of the exponent, and returns the
    exact value of the integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.

    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
*/
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }

	bot = bot * ( double ) ( expon + 2 );

    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;
  }
  return exact;
}
/******************************************************************************/

void clenshaw_curtis_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      x[i] =  cos ( ( double ) ( order - 1 - i ) * pi 
                       / ( double ) ( order - 1     ) );
    }
    x[0] = -1.0;
    if ( ( order % 2 ) == 1 )
    {
      x[(order-1)/2] = 0.0;
    }
    x[order-1] = +1.0;

    for ( i = 0; i < order; i++ )
    {
      theta = ( double ) ( i ) * pi / ( double ) ( order - 1 );

      w[i] = 1.0;

      for ( j = 1; j <= ( order - 1 ) / 2; j++ )
      {
        if ( 2 * j == ( order - 1 ) )
        {
          b = 1.0;
        }
        else
        {
          b = 2.0;
        }

        w[i] = w[i] - b *  cos ( 2.0 * ( double ) ( j ) * theta ) 
          / ( double ) ( 4 * j * j - 1 );
      }
    }

    w[0] = w[0] / ( double ) ( order - 1 );
    for ( i = 1; i < order - 1; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( order - 1 );
    }
    w[order-1] = w[order-1] / ( double ) ( order - 1 );
  }

  return;
}
/******************************************************************************/

void clenshaw_curtis_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  clenshaw_curtis_compute ( order, x, w );

  return;
}
/******************************************************************************/

void clenshaw_curtis_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    This rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double X[ORDER], the abscissas.
*/
{
  int index;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CLENSHAW_CURTIS_COMPUTE_POINTS - Fatal error!\n" );
    printf ( "  ORDER < 1.\n" );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( index = 1; index <= order; index++ )
    {
      x[index-1] =  cos ( ( double ) ( order - index ) * pi 
                             / ( double ) ( order - 1     ) );
    }
    x[0] = -1.0;
    if ( ( order % 2 ) == 1 )
    {
      x[(order-1)/2] = 0.0;
    }
    x[order-1] = +1.0;
  }
  return;
}
/******************************************************************************/

void clenshaw_curtis_compute_points_np ( int order, int np, double p[], 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE_POINTS_NP: Clenshaw Curtis quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    This rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  clenshaw_curtis_compute_points ( order, x );

  return;
}
/******************************************************************************/

void clenshaw_curtis_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE_WEIGHTS computes Clenshaw Curtis quadrature weights.

  Discussion:

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double W[ORDER], the weights of the rule.
*/
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "CLENSHAW_CURTIS_COMPUTE_WEIGHTS - Fatal error!\n" );
    printf ( "  ORDER < 1.\n" );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    w[0] = 2.0;
    return;
  }

  for ( i = 1; i <= order; i++ )
  {
    theta = ( double ) ( i - 1 ) * pi / ( double ) ( order - 1 );

    w[i-1] = 1.0;

    for ( j = 1; j <= ( order - 1 ) / 2; j++ )
    {
      if ( 2 * j == ( order - 1 ) )
      {
        b = 1.0;
      }
      else
      {
        b = 2.0;
      }

      w[i-1] = w[i-1] - b *  cos ( 2.0 * ( double ) ( j ) * theta ) 
           / ( double ) ( 4 * j * j - 1 );
    }
  }

  w[0] = w[0] / ( double ) ( order - 1 );
  for ( i = 1; i < order - 1; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( order - 1 );
  }
  w[order-1] = w[order-1] / ( double ) ( order - 1 );

  return;
}
/******************************************************************************/

void clenshaw_curtis_compute_weights_np ( int order, int np, double p[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE_WEIGHTS_NP: Clenshaw Curtis quadrature weights.

  Discussion:

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights of the rule.
*/
{
  clenshaw_curtis_compute_weights ( order, w );

  return;
}
/******************************************************************************/

void comp_next ( int n, int k, int a[], int *more, int *h,  int *t )

/******************************************************************************/
/*
  Purpose:

    COMP_NEXT computes the compositions of the integer N into K parts.

  Discussion:

    A composition of the integer N into K parts is an ordered sequence
    of K nonnegative integers which sum to N.  The compositions (1,2,1)
    and (1,1,2) are considered to be distinct.

    The routine computes one composition on each call until there are no more.
    For instance, one composition of 6 into 3 parts is
    3+2+1, another would be 6+0+0.

    On the first call to this routine, set MORE = FALSE.  The routine
    will compute the first element in the sequence of compositions, and
    return it, as well as setting MORE = TRUE.  If more compositions
    are desired, call again, and again.  Each time, the routine will
    return with a new composition.

    However, when the LAST composition in the sequence is computed 
    and returned, the routine will reset MORE to FALSE, signaling that
    the end of the sequence has been reached.

    This routine originally used a STATICE statement to maintain the
    variables H and T.  I have decided (based on an wasting an
    entire morning trying to track down a problem) that it is safer
    to pass these variables as arguments, even though the user should
    never alter them.  This allows this routine to safely shuffle
    between several ongoing calculations.


    There are 28 compositions of 6 into three parts.  This routine will
    produce those compositions in the following order:

     I         A
     -     ---------
     1     6   0   0
     2     5   1   0
     3     4   2   0
     4     3   3   0
     5     2   4   0
     6     1   5   0
     7     0   6   0
     8     5   0   1
     9     4   1   1
    10     3   2   1
    11     2   3   1
    12     1   4   1
    13     0   5   1
    14     4   0   2
    15     3   1   2
    16     2   2   2
    17     1   3   2
    18     0   4   2
    19     3   0   3
    20     2   1   3
    21     1   2   3
    22     0   3   3
    23     2   0   4
    24     1   1   4
    25     0   2   4
    26     1   0   5
    27     0   1   5
    28     0   0   6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2008

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.

  Parameters:

    Input, int N, the integer whose compositions are desired.

    Input, int K, the number of parts in the composition.

    Input/output, int A[K], the parts of the composition.

    Input/output, int *MORE.
    Set MORE = FALSE on first call.  It will be reset to TRUE on return
    with a new composition.  Each new call returns another composition until
    MORE is set to FALSE when the last composition has been computed
    and returned.

    Input/output, int *H, *T, two internal parameters needed for the
    computation.  The user should allocate space for these in the calling
    program, include them in the calling sequence, but never alter them!
*/
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
/*******************************************************************************/

double cpu_time ( void )

/*******************************************************************************/
/*
  Purpose:
 
    CPU_TIME reports the total CPU time for a program.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

void fejer2_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE computes a Fejer type 2 rule.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  int i;
  int j;
  double p;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "FEJER2_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  for ( i = 0; i < order; i++ )
  {
    x[i] =  cos ( ( double ) ( order - i ) * pi 
                     / ( double ) ( order + 1 ) );
  }
  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  if ( order == 2 )
  {
    w[0] = 1.0;
    w[1] = 1.0;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      theta = ( double ) ( order - i ) * pi 
            / ( double ) ( order + 1 );

      w[i] = 1.0;

      for ( j = 1; j <= ( ( order - 1 ) / 2 ); j++ )
      {
        w[i] = w[i] - 2.0 *  cos ( 2.0 * ( double ) ( j ) * theta ) 
          / ( double ) ( 4 * j * j - 1 );
      }
      p = 2.0 * ( double ) ( ( ( order + 1 ) / 2 ) ) - 1.0;
      w[i] = w[i] -  cos ( ( p + 1.0 ) * theta ) / p;
    }
    for ( i = 0; i < order; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( order + 1 );
    }
  }

  return;
}
/******************************************************************************/

void fejer2_compute_np ( int order, int np, double p[], double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE_NP computes a Fejer type 2 rule.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  fejer2_compute ( order, x, w );

  return;
}
/******************************************************************************/

void fejer2_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE_POINTS computes Fejer type 2 quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.
*/
{
  int index;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "FEJER2_COMPUTE_POINTS - Fatal error!\n" );
    printf ( "  ORDER < 1.\n" );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( index = 1; index <= order; index++ )
    {
      x[index-1] =  cos ( ( double ) ( order + 1 - index ) * pi 
                             / ( double ) ( order + 1 ) );
    }
    if ( ( order % 2 ) == 1 )
    {
      x[(order-1)/2] = 0.0;
    }
  }
  return;
}
/******************************************************************************/

void fejer2_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE_POINTS_NP computes Fejer type 2 quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  fejer2_compute_points ( order, x );

  return;
}
/******************************************************************************/

void fejer2_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE_WEIGHTS computes Fejer type 2 quadrature weights.

  Discussion:

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.

    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.

  Parameters:

    Input, int ORDER, the order.

    Output, double W[ORDER], the weights.
*/
{
  int i;
  int j;
  double p;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "FEJER2_COMPUTE_WEIGHTS - Fatal error!\n" );
    printf ( "  ORDER < 1.\n" );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    w[0] = 2.0;
  }
  else if ( order == 2 )
  {
    w[0] = 1.0;
    w[1] = 1.0;
  }
  else
  {
    for ( i = 1; i <= order; i++ )
    {
      theta = ( double ) ( order + 1 - i ) * pi 
            / ( double ) ( order + 1 );

      w[i-1] = 1.0;

      for ( j = 1; j <= ( ( order - 1 ) / 2 ); j++ )
      {
        w[i-1] = w[i-1] - 2.0 *  cos ( 2.0 * ( double ) ( j ) * theta ) 
          / ( double ) ( 4 * j * j - 1 );
      }
      p = 2.0 * ( double ) ( ( ( order + 1 ) / 2 ) ) - 1.0;
      w[i-1] = w[i-1] -  cos ( ( p + 1.0 ) * theta ) / p;
    }
    for ( i = 0; i < order; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( order + 1 );
    }
  }
  return;
}
/******************************************************************************/

void fejer2_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE_WEIGHTS_NP computes Fejer type 2 quadrature weights.

  Discussion:

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.

    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.

  Parameters:

    Input, int ORDER, the order.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  fejer2_compute_weights ( order, w );

  return;
}
/******************************************************************************/

void gegenbauer_compute ( int order, double alpha, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE computes a Gegenbauer quadrature rule.

  Discussion:

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    Thanks to Janiki Raman for pointing out a problem in an earlier
    version of the code that occurred when ALPHA was -0.5.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double ALPHA, the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double an;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x0;
/*
  Check ORDER.
*/
  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "GEGENBAUER_COMPUTE - Fatal error!\n" );
    printf ( "  1 <= ORDER is required.\n" );
    exit ( 1 );
  }
  c = ( double * ) malloc ( order * sizeof ( double ) );
/*
  Check ALPHA.
*/
  if ( alpha <= -1.0 )
  {
    printf ( "\n" );
    printf ( "GEGENBAUER_COMPUTE - Fatal error!\n" );
    printf ( "  -1.0 < ALPHA is required.\n" );
    exit ( 1 );
  }
/*
  Set the recursion coefficients.
*/
  c[0] = 0.0;
  if ( 2 <= order )
  {
    c[1] = 1.0 / ( 2.0 * alpha + 3.0 );
  }

  for ( i = 3; i <= order; i++ )
  {
    c[i-1] = ( double ) ( i - 1 ) 
          * ( alpha + alpha + ( double ) ( i - 1 ) ) / 
          ( ( alpha + alpha + ( double ) ( 2 * i - 1 ) ) 
          * ( alpha + alpha + ( double ) ( 2 * i - 3 ) ) );
  }

  delta = r8_gamma ( alpha         + 1.0 ) 
        * r8_gamma (         alpha + 1.0 ) 
        / r8_gamma ( alpha + alpha + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + alpha + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 2.44 * an + 1.282 * an * an;

      x0 = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * alpha * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x0 = x0 - r1 * r2 * r3 * ( 1.0 - x0 );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) );

      x0 = x0 - r1 * r2 * r3 * ( x[0] - x0 );
    }
    else if ( i < order - 1 )
    {
      x0 = 3.0 * x[i-2] - 3.0 * x[i-3] + x[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * alpha ) / ( 0.766 + 0.119 * alpha );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * alpha ) / ( 1.67 + 0.28 * alpha );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }

    gegenbauer_root ( &x0, order, alpha, &dp2, &p1, c );

    x[i-1] = x0;
    w[i-1] = cc / ( dp2 * p1 );
  }
/*
  Reverse the order of the values.
*/
  for ( i = 1; i <= order/2; i++ )
  {
    temp       = x[i-1];
    x[i-1]     = x[order-i];
    x[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp       = w[i-1];
    w[i-1]     = w[order-i];
    w[order-i] = temp;
  }

  free ( c );

  return;
}
/******************************************************************************/

void gegenbauer_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE_NP computes a Gegenbauer quadrature rule.

  Discussion:

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    Thanks to Janiki Raman for pointing out a problem in an earlier
    version of the code that occurred when ALPHA was -0.5.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA = the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;

  alpha = p[0];

  gegenbauer_compute ( order, alpha, x, w );

  return;
}
/******************************************************************************/

void gegenbauer_compute_points ( int order, double alpha, double x[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE_POINTS computes Gegenbauer quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double ALPHA, the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  gegenbauer_compute ( order, alpha, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void gegenbauer_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE_POINTS_NP computes Gegenbauer quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA = the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double X[ORDER], the abscissas.
*/
{
  double alpha;

  alpha = p[0];

  gegenbauer_compute_points ( order, alpha, x );

  return;
}
/******************************************************************************/

void gegenbauer_compute_weights ( int order, double alpha, double w[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE_WEIGHTS computes Gegenbauer quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double ALPHA, the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  gegenbauer_compute ( order, alpha, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void gegenbauer_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_COMPUTE_WEIGHTS_NP computes Gegenbauer quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double P[1], contains parameters.
    P[0] = ALPHA = the exponent of (1-X^2).  -1.0 < ALPHA is required.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;

  alpha = p[0];

  gegenbauer_compute_weights ( order, alpha, w );

  return;
}
/******************************************************************************/

double gegenbauer_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_INTEGRAL integrates a monomial with Gegenbauer weight.

  Discussion:

    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x^2)^ALPHA dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.

    Input, double ALPHA, the exponent of (1-X^2) in the weight factor.

    Output, double GEGENBAUER_INTEGRAL, the value of the integral.
*/
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double value;
  double value1;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = r8_gamma ( 1.0 + c ) * 2.0 
    * r8_gamma ( 1.0 + alpha  ) * value1 
    / r8_gamma ( 2.0 + alpha  + c );

  return value;
}
/******************************************************************************/

void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double c[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_RECUR evaluates a Gegenbauer polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Output, double *P2, the value of J(ORDER)(X).

    Output, double *DP2, the value of J'(ORDER)(X).

    Output, double *P1, the value of J(ORDER-1)(X).

    Input, double X, the point at which polynomials are evaluated.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, the exponents of (1-X^2).

    Input, double C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x;
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = x *  ( *p1 ) - c[i-1] * p0;
    *dp2 = x * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
/******************************************************************************/

void gegenbauer_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double c[] )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input/output, double *X, the approximate root, which
    should be improved on output.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, the exponents of (1-X^2).

    Output, double *DP2, the value of J'(ORDER)(X).

    Output, double *P1, the value of J(ORDER-1)(X).

    Input, double C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gegenbauer_recur ( &p2, dp2, p1, *x, order, alpha, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
/******************************************************************************/

void gen_hermite_compute ( int order, double alpha, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE computes a Generalized Hermite quadrature rule.

  Discussion:

    The integral to be approximated has the form:

      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.

    Output, double X[ORDER], W[ORDER], the abscissas and weights of the rule.
*/
{
  double alpha_laguerre;
  double arg;
  int i;
  int order_laguerre;
  double *w_laguerre;
  double *x_laguerre;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "GEN_HERMITE_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  if ( order == 1 )
  {
    arg = ( alpha + 1.0 ) / 2.0;
    x[0] = 0.0;
    w[0] = r8_gamma ( arg );
    return;
  }

  if ( ( order % 2 ) == 0 ) 
  {
    order_laguerre = order / 2;
    alpha_laguerre = ( alpha - 1.0 ) / 2.0;
  }
  else
  {
    order_laguerre = ( order - 1 ) / 2;
    alpha_laguerre = ( alpha + 1.0 ) / 2.0;
  }
  
  w_laguerre = ( double * ) malloc ( order_laguerre * sizeof ( double ) );
  x_laguerre = ( double * ) malloc ( order_laguerre * sizeof ( double ) );

  gen_laguerre_compute ( order_laguerre, alpha_laguerre, x_laguerre, 
    w_laguerre );

  if ( ( order % 2 ) == 0 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i];
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+i] = 0.5 * w_laguerre[i];
    }
  }
  else if ( ( order % 2 ) == 1 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    x[order_laguerre] = 0.0;
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+1+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i] / x_laguerre[order_laguerre-1-i];
    }

    arg = ( alpha + 1.0 ) / 2.0;
    w[order_laguerre] = r8_gamma ( arg );
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre] = w[order_laguerre] - w_laguerre[i] / x_laguerre[i];
    }

    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+1+i] = 0.5 * w_laguerre[i] / x_laguerre[i];
    }
  }
  free ( w_laguerre );
  free ( x_laguerre );

  return;
}
/******************************************************************************/

void gen_hermite_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE_NP computes a Generalized Hermite quadrature rule.

  Discussion:

    The integral to be approximated has the form:

      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.

    Output, double X[ORDER], W[ORDER], the abscissas and weights of the rule.
*/
{
  double alpha;

  alpha = p[0];

  gen_hermite_compute ( order, alpha, x, w );

  return;
}
/******************************************************************************/

void gen_hermite_compute_points ( int order, double alpha, double x[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE_POINTS computes Generalized Hermite quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  gen_hermite_compute ( order, alpha, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void gen_hermite_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE_POINTS_NP: Generalized Hermite quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.

    Output, double X[ORDER], the abscissas.
*/
{
  double alpha;

  alpha = p[0];

  gen_hermite_compute_points ( order, alpha, x );

  return;
}
/******************************************************************************/

void gen_hermite_compute_weights ( int order, double alpha, double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE_WEIGHTS computes Generalized Hermite quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  gen_hermite_compute ( order, alpha, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void gen_hermite_compute_weights_np ( int order, int np, double p[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_COMPUTE_WEIGHTS_NP: Generalized Hermite quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;

  alpha = p[0];

  gen_hermite_compute_weights ( order, alpha, w );

  return;
}
/******************************************************************************/

double gen_hermite_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_INTEGRAL evaluates a monomial Generalized Hermite integral.

  Discussion:

    H(n,alpha) = Integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent of the monomial.
    0 <= EXPON.

    Input, double ALPHA, the exponent of |X| in the weight function.
    -1.0 < ALPHA.

    Output, double GEN_HERMITE_INTEGRAL, the value of the integral.
*/
{
  double a;
  double arg;
  double value;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    a = alpha + ( double ) ( expon );
    if ( a <= - 1.0 )
    {
      value = - r8_huge ( );
    }
    else
    {
      arg = ( a + 1.0 ) / 2.0;
      value = r8_gamma ( arg );
    }
  }
  return value;
}
/******************************************************************************/

void gen_laguerre_compute ( int order, double alpha, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE computes a Generalized Laguerre quadrature rule.

  Discussion:

    In the simplest case, ALPHA is 0, and we are approximating the
    integral from 0 to +oo of exp(-X) * F(X).  When this is so,
    it is easy to modify the rule to approximate the integral from
    A to +oo as well.

    If ALPHA is nonzero, then there is no simple way to extend the
    rule to approximate the integral from A to +oo.  The simplest
    procedures would be to approximate the integral from 0 to A.

    The integration interval is [ A, +oo ) or [ 0, +oo ).

    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.


    If the integral to approximate is:

        Integral ( A <= X < +oo ) exp ( - X ) * F(X) dX
      or
        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX

    then the quadrature rule is:

      exp ( - A ) * Sum ( 1 <= I <= ORDER ) W(I) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )


    If the integral to approximate is:

        Integral ( A <= X < +oo ) F(X) dX
      or
        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX

    then the quadrature rule is:

      exp ( - A ) * Sum ( 1 <= I <= ORDER ) 
        W(I) * exp(A+X(I)) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * exp(X(I)) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the quadrature rule to be computed.
    1 <= ORDER.

    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x0;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "GEN_LAGUERRE_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  b = ( double * ) malloc ( order * sizeof ( double ) );
  c = ( double * ) malloc ( order * sizeof ( double ) );
/*
  Set the recursion coefficients.
*/
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( alpha + ( double ) ( 2 * i + 1 ) );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i ) * ( alpha + ( double ) ( i ) );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = r8_gamma ( alpha + 1.0 ) * prod;

  for ( i = 0; i < order; i++ )
  {
/*
  Compute an estimate for the root.
*/
    if ( i == 0 )
    {
      x0 = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
        ( 1.0 + 2.4 * ( double ) ( order ) + 1.8 * alpha );
    }
    else if ( i == 1 )
    {
      x0 = x0 + ( 15.0 + 6.25 * alpha ) / 
        ( 1.0 + 0.9 * alpha + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      r2 = 1.26 * ( double ) ( i - 1 ) * alpha / 
        ( 1.0 + 3.5 * ( double ) ( i - 1 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      x0 = x0 + ratio * ( x0 - x[i-2] );
    }
/*
  Use iteration to find the root.
*/
    gen_laguerre_root ( &x0, order, alpha, &dp2, &p1, b, c );
/*
  Set the abscissa and weight.
*/
    x[i] = x0;
    w[i] = ( cc / dp2 ) / p1;
  }

  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void gen_laguerre_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE_NP computes a Generalized Laguerre quadrature rule.

  Discussion:

    In the simplest case, ALPHA is 0, and we are approximating the
    integral from 0 to +oo of exp(-X) * F(X).  When this is so,
    it is easy to modify the rule to approximate the integral from
    A to +oo as well.

    If ALPHA is nonzero, then there is no simple way to extend the
    rule to approximate the integral from A to +oo.  The simplest
    procedures would be to approximate the integral from 0 to A.

    The integration interval is [ A, +oo ) or [ 0, +oo ).

    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.


    If the integral to approximate is:

        Integral ( A <= X < +oo ) exp ( - X ) * F(X) dX
      or
        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX

    then the quadrature rule is:

      exp ( - A ) * Sum ( 1 <= I <= ORDER ) W(I) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )


    If the integral to approximate is:

        Integral ( A <= X < +oo ) F(X) dX
      or
        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX

    then the quadrature rule is:

      exp ( - A ) * Sum ( 1 <= I <= ORDER ) 
        W(I) * exp(A+X(I)) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * exp(X(I)) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the quadrature rule to be computed.
    1 <= ORDER.

    Input, double P[1], contains parameters.
    P[0] = ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;

  alpha = p[0];

  gen_laguerre_compute ( order, alpha, x, w );

  return;
}
/******************************************************************************/

void gen_laguerre_compute_points ( int order, double alpha, double x[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE_POINTS: Generalized Laguerre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  gen_laguerre_compute ( order, alpha, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void gen_laguerre_compute_points_np ( int order, int np, double p[], 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE_POINTS_NP: Generalized Laguerre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double X[ORDER], the abscissas.
*/
{
  double alpha;

  alpha = p[0];

  gen_laguerre_compute_points ( order, alpha, x );

  return;
}
/******************************************************************************/

void gen_laguerre_compute_weights ( int order, double alpha, double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE_WEIGHTS: Generalized Laguerre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  gen_laguerre_compute ( order, alpha, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void gen_laguerre_compute_weights_np ( int order, int np, double p[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_COMPUTE_WEIGHTS_NP: Generalized Laguerre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;

  alpha = p[0];

  gen_laguerre_compute_weights ( order, alpha, w );

  return;
}
/******************************************************************************/

double gen_laguerre_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_INTEGRAL evaluates a monomial Generalized Laguerre integral.

  Discussion:

    L(n,alpha) = Integral ( 0 <= x < +oo ) x^n * x^alpha exp(-x) dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent of the monomial.
    0 <= EXPON.

    Input, double ALPHA, the exponent of X in the weight function.
    -1.0 < ALPHA.

    Output, double GEN_LAGUERRE_INTEGRAL, the value of the integral.
*/
{
  double arg;
  double value;

  arg = alpha + ( double ) ( expon + 1.0 );
  value = r8_gamma ( arg );

  return value;
}
/******************************************************************************/

void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_RECUR evaluates a Generalized Laguerre polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Output, double *P2, the value of L(ORDER)(X).

    Output, double *DP2, the value of L'(ORDER)(X).

    Output, double *P1, the value of L(ORDER-1)(X).

    Input, double X, the point at which polynomials are evaluated.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, the exponent of the X factor in the
    integrand.

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - alpha - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
/******************************************************************************/

void gen_laguerre_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_ROOT improves a root of a Generalized Laguerre polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input/output, double *X, the approximate root, which
    should be improved on output.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, the exponent of the X factor.

    Output, double *DP2, the value of L'(ORDER)(X).

    Output, double *P1, the value of L(ORDER-1)(X).

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gen_laguerre_recur ( &p2, dp2, p1, *x, order, alpha, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void hermite_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE computes a Hermite quadrature rule.

  Discussion:

    The abscissas are the zeros of the N-th order Hermite polynomial.

    The integration interval is ( -oo, +oo ).

    The weight function is w(x) = exp ( - x * x ).

    The integral to approximate:

      Integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double cc;
  double dp2;
  int i;
  double p1;
  double s;
  double temp;
  double x0;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "HERMITE_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  cc = 1.7724538509 * r8_gamma ( ( double ) ( order ) ) 
    / pow ( 2.0, order - 1 );

  s = pow ( 2.0 * ( double ) ( order ) + 1.0, 1.0 / 6.0 );

  for ( i = 0; i < ( order + 1 ) / 2; i++ )
  {
    if ( i == 0 )
    {
      x0 = s * s * s - 1.85575 / s;
    }
    else if ( i == 1 )
    {
      x0 = x0 - 1.14 * pow ( ( double ) ( order ), 0.426 ) / x0;
    }
    else if ( i == 2 )
    {
      x0 = 1.86 * x0 - 0.86 * x[0];
    }
    else if ( i == 3 )
    {
      x0 = 1.91 * x0 - 0.91 * x[1];
    }
    else
    {
      x0 = 2.0 * x0 - x[i-2];
    }

    hermite_root ( &x0, order, &dp2, &p1 );

    x[i] = x0;
    w[i] = ( cc / dp2 ) / p1;

    x[order-i-1] = -x0;
    w[order-i-1] = w[i];
  }
/*
  Reverse the order of the abscissas.
*/
  for ( i = 1; i <= order/2; i++ )
  {
    temp       = x[i-1];
    x[i-1]     = x[order-i];
    x[order-i] = temp;
  }

  if ( ( order % 2 ) == 1 )
  {
    x[(order-1)/2] = 0.0;
  }

  return;
}
/******************************************************************************/

void hermite_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE_NP computes a Hermite quadrature rule.

  Discussion:

    The abscissas are the zeros of the N-th order Hermite polynomial.

    The integration interval is ( -oo, +oo ).

    The weight function is w(x) = exp ( - x * x ).

    The integral to approximate:

      Integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  hermite_compute ( order, x, w );

  return;
}
/******************************************************************************/

void hermite_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE_POINTS computes Hermite quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  hermite_compute ( order, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void hermite_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE_POINTS_NP computes Hermite quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  hermite_compute_points ( order, x );

  return;
}
/******************************************************************************/

void hermite_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE_WEIGHTS computes Hermite quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  hermite_compute ( order, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void hermite_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_COMPUTE_WEIGHTS_NP computes Hermite quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  hermite_compute_weights ( order, w );

  return;
}
/******************************************************************************/

void hermite_genz_keister_lookup ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_GENZ_KEISTER_LOOKUP looks up a Genz-Keister Hermite rule.

  Discussion:

    The integral:

      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx

    The quadrature rule:

      sum ( 1 <= i <= n ) w(i) * f ( x(i) )

    A nested family of rules for the Hermite integration problem
    was produced by Genz and Keister.  The structure of the nested
    family was denoted by 1+2+6+10+16, that is, it comprised rules 
    of successive orders O = 1, 3, 9, 19, and 35.

    The precisions of these rules are P = 1, 5, 15, 29, and 51.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 June 2010

  Author:

    John Burkardt

  Reference:

    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309

    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    N must be 1, 3, 9, 19, or 35.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  hermite_genz_keister_lookup_points ( n, x );
  hermite_genz_keister_lookup_weights ( n, w );

  return;
}
/******************************************************************************/

void hermite_genz_keister_lookup_points ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_GENZ_KEISTER_LOOKUP_POINTS looks up Genz-Keister Hermite abscissas.

  Discussion:

    The integral:

      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx

    The quadrature rule:

      sum ( 1 <= i <= n ) w(i) * f ( x(i) )

    A nested family of rules for the Hermite integration problem
    was produced by Genz and Keister.  The structure of the nested
    family was denoted by 1+2+6+10+16, that is, it comprised rules 
    of successive orders O = 1, 3, 9, 19, and 35.

    The precisions of these rules are P = 1, 5, 15, 29, and 51.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 June 2010

  Author:

    John Burkardt

  Reference:

    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309

    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    N must be 1, 3, 9, 19, or 35.

    Output, double X[N], the abscissas.
*/
{
  if ( n == 1 )
  {
    x[ 0] =   0.0000000000000000E+00;
  }
  else if ( n == 3 )
  {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;
  }
  else if ( n == 9 )
  {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;
  }
  else if ( n == 19 )
  {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;
  }
  else if ( n == 35 )
  {
    x[ 0] =  -6.3759392709822356E+00;
    x[ 1] =  -5.6432578578857449E+00;
    x[ 2] =  -5.0360899444730940E+00;
    x[ 3] =  -4.4995993983103881E+00;
    x[ 4] =  -4.0292201405043713E+00;
    x[ 5] =  -3.6677742159463378E+00;
    x[ 6] =  -3.3491639537131945E+00;
    x[ 7] =  -2.9592107790638380E+00;
    x[ 8] =  -2.5705583765842968E+00;
    x[ 9] =  -2.2665132620567876E+00;
    x[10] =  -2.0232301911005157E+00;
    x[11] =  -1.8357079751751868E+00;
    x[12] =  -1.5794121348467671E+00;
    x[13] =  -1.2247448713915889E+00;
    x[14] =  -8.7004089535290285E-01;
    x[15] =  -5.2403354748695763E-01;
    x[16] =  -1.7606414208200893E-01;
    x[17] =   0.0000000000000000E+00;
    x[18] =   1.7606414208200893E-01;
    x[19] =   5.2403354748695763E-01;
    x[20] =   8.7004089535290285E-01;
    x[21] =   1.2247448713915889E+00;
    x[22] =   1.5794121348467671E+00;
    x[23] =   1.8357079751751868E+00;
    x[24] =   2.0232301911005157E+00;
    x[25] =   2.2665132620567876E+00;
    x[26] =   2.5705583765842968E+00;
    x[27] =   2.9592107790638380E+00;
    x[28] =   3.3491639537131945E+00;
    x[29] =   3.6677742159463378E+00;
    x[30] =   4.0292201405043713E+00;
    x[31] =   4.4995993983103881E+00;
    x[32] =   5.0360899444730940E+00;
    x[33] =   5.6432578578857449E+00;
    x[34] =   6.3759392709822356E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "HERMITE_GENZ_KEISTER_LOOKUP_POINTS - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of N.\n" );
    fprintf ( stderr, "  N must be 1, 3, 9, 19, or 35.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void hermite_genz_keister_lookup_points_np ( int n, int np, double p[], 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_GENZ_KEISTER_LOOKUP_POINTS_NP looks up Genz-Keister Hermite abscissas.

  Discussion:

    The integral:

      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx

    The quadrature rule:

      sum ( 1 <= i <= n ) w(i) * f ( x(i) )

    A nested family of rules for the Hermite integration problem
    was produced by Genz and Keister.  The structure of the nested
    family was denoted by 1+2+6+10+16, that is, it comprised rules 
    of successive orders O = 1, 3, 9, 19, and 35.

    The precisions of these rules are P = 1, 5, 15, 29, and 51.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 June 2010

  Author:

    John Burkardt

  Reference:

    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309

    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    N must be 1, 3, 9, 19, or 35.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[N], the abscissas.
*/
{
  hermite_genz_keister_lookup_points ( n, x );

  return;
}
/******************************************************************************/

void hermite_genz_keister_lookup_weights ( int n, double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up Genz-Keister Hermite weights.

  Discussion:

    The integral:

      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx

    The quadrature rule:

      sum ( 1 <= i <= n ) w(i) * f ( x(i) )

    A nested family of rules for the Hermite integration problem
    was produced by Genz and Keister.  The structure of the nested
    family was denoted by 1+2+6+10+16, that is, it comprised rules 
    of successive orders O = 1, 3, 9, 19, and 35.

    The precisions of these rules are P = 1, 5, 15, 29, and 51.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 June 2010

  Author:

    John Burkardt

  Reference:

    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309

    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    N must be 1, 3, 9, 19, or 35.

    Output, double W[N], the weights.
*/
{
  if ( n == 1 )
  {
    w[ 0] =   1.7724538509055159E+00;
  }
  else if ( n == 3 )
  {
    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if ( n == 9 )
  {
    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if ( n == 19 )
  {
    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if ( n == 35 )
  {
    w[ 0] =   1.8684014894510604E-18;
    w[ 1] =   9.6599466278563243E-15;
    w[ 2] =   5.4896836948499462E-12;
    w[ 3] =   8.1553721816916897E-10;
    w[ 4] =   3.7920222392319532E-08;
    w[ 5] =   4.3737818040926989E-07;
    w[ 6] =   4.8462799737020461E-06;
    w[ 7] =   6.3328620805617891E-05;
    w[ 8] =   4.8785399304443770E-04;
    w[ 9] =   1.4515580425155904E-03;
    w[10] =   4.0967527720344047E-03;
    w[11] =   5.5928828911469180E-03;
    w[12] =   2.7780508908535097E-02;
    w[13] =   8.0245518147390893E-02;
    w[14] =   1.6371221555735804E-01;
    w[15] =   2.6244871488784277E-01;
    w[16] =   3.3988595585585218E-01;
    w[17] =   9.1262675363737921E-04;
    w[18] =   3.3988595585585218E-01;
    w[19] =   2.6244871488784277E-01;
    w[20] =   1.6371221555735804E-01;
    w[21] =   8.0245518147390893E-02;
    w[22] =   2.7780508908535097E-02;
    w[23] =   5.5928828911469180E-03;
    w[24] =   4.0967527720344047E-03;
    w[25] =   1.4515580425155904E-03;
    w[26] =   4.8785399304443770E-04;
    w[27] =   6.3328620805617891E-05;
    w[28] =   4.8462799737020461E-06;
    w[29] =   4.3737818040926989E-07;
    w[30] =   3.7920222392319532E-08;
    w[31] =   8.1553721816916897E-10;
    w[32] =   5.4896836948499462E-12;
    w[33] =   9.6599466278563243E-15;
    w[34] =   1.8684014894510604E-18;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of N.\n" );
    fprintf ( stderr, "  N must be 1, 3, 9, 19, or 35.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void hermite_genz_keister_lookup_weights_np ( int n, int np, double p[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS_NP looks up Genz-Keister Hermite weights.

  Discussion:

    The integral:

      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx

    The quadrature rule:

      sum ( 1 <= i <= n ) w(i) * f ( x(i) )

    A nested family of rules for the Hermite integration problem
    was produced by Genz and Keister.  The structure of the nested
    family was denoted by 1+2+6+10+16, that is, it comprised rules 
    of successive orders O = 1, 3, 9, 19, and 35.

    The precisions of these rules are P = 1, 5, 15, 29, and 51.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 June 2010

  Author:

    John Burkardt

  Reference:

    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309

    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    N must be 1, 3, 9, 19, or 35.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[N], the weights.
*/
{
  hermite_genz_keister_lookup_weights ( n, w );

  return;
}
/******************************************************************************/

double hermite_integral ( int n )

/******************************************************************************/
/*
  Purpose:

    HERMITE_INTEGRAL evaluates a monomial Hermite integral.

  Discussion:

    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x^2) dx

    H(n) is 0 for n odd.

    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the integral.  
    0 <= N.

    Output, double VALUE, the value of the integral.
*/
{
  double pi = 3.141592653589793;
  double value;

  if ( n < 0 )
  {
    value = - r8_huge ( );
  }
  else if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) 
      / pow ( 2.0, n / 2 );
  }

  return value;
}
/******************************************************************************/

void hermite_recur ( double *p2, double *dp2, double *p1, double x, int order )

/******************************************************************************/
/*
  Purpose:

    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Output, double *P2, the value of H(ORDER)(X).

    Output, double *DP2, the value of H'(ORDER)(X).

    Output, double *P1, the value of H(ORDER-1)(X).

    Input, double X, the point at which polynomials are evaluated.

    Input, int ORDER, the order of the polynomial to be computed.
*/
{
  int i;
  double dq0;
  double dq1;
  double dq2;
  double q0;
  double q1;
  double q2;

  q1 = 1.0;
  dq1 = 0.0;

  q2 = x;
  dq2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    q0 = q1;
    dq0 = dq1;

    q1 = q2;
    dq1 = dq2;

    q2  = x * q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * q0;
    dq2 = x * dq1 + q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * dq0;
  }

  *p2 = q2;
  *dp2 = dq2;
  *p1 = q1;

  return;
}
/******************************************************************************/

void hermite_root ( double *x, int order, double *dp2, double *p1 )

/******************************************************************************/
/*
  Purpose:

    HERMITE_ROOT improves an approximate root of a Hermite polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input/output, double *X, the approximate root, which
    should be improved on output.

    Input, int ORDER, the order of the Hermite polynomial.

    Output, double *DP2, the value of H'(ORDER)(X).

    Output, double *P1, the value of H(ORDER-1)(X).
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    hermite_recur ( &p2, dp2, p1, *x, order );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }

  return;
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

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_WRITE writes an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the table data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "I4MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
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

int i4vec_product ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRODUCT multiplies the entries of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Example:

    Input:

      A = ( 1, 2, 3, 4 )

    Output:

      I4VEC_PRODUCT = 24

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int A[N], the vector

    Output, int I4VEC_PRODUCT, the product of the entries of A.
*/
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
/******************************************************************************/

int i4vec_sum ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_SUM sums the entries of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Example:

    Input:

      A = ( 1, 2, 3, 4 )

    Output:

      I4VEC_SUM = 10

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int A[N], the vector to be summed.

    Output, int I4VEC_SUM, the sum of the entries of A.
*/
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
/******************************************************************************/

void jacobi_compute ( int order, double alpha, double beta, double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE computes a Jacobi quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    Thanks to Xu Xiang of Fudan University for pointing out that
    an earlier implementation of this routine was incorrect!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.  For simple Legendre quadrature,
    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double an;
  double *b;
  double bn;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x0;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "JACOBI_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  b = ( double * ) malloc ( order * sizeof ( double ) );
  c = ( double * ) malloc ( order * sizeof ( double ) );
/*
  Check ALPHA and BETA.
*/
  if ( alpha <= -1.0 )
  {
    printf ( "\n" );
    printf ( "JACOBI_COMPUTE - Fatal error!\n" );
    printf ( "  -1.0 < ALPHA is required.\n" );
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    printf ( "\n" );
    printf ( "JACOBI_COMPUTE - Fatal error!\n" );
    printf ( "  -1.0 < BETA is required.\n" );
    exit ( 1 );
  }
/*
  Set the recursion coefficients.
*/
  for ( i = 1; i <= order; i++ )
  {
    if ( alpha + beta == 0.0 || beta - alpha == 0.0 )
    {
      b[i-1] = 0.0;
    }
    else
    {
      b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
             ( ( alpha + beta + ( double ) ( 2 * i ) ) 
             * ( alpha + beta + ( double ) ( 2 * i - 2 ) ) );
    }

    if ( i == 1 )
    {
      c[i-1] = 0.0;
    }
    else
    {
      c[i-1] = 4.0 * ( double ) ( i - 1 ) 
         * ( alpha + ( double ) ( i - 1 ) ) 
          * ( beta + ( double ) ( i - 1 ) ) 
            * ( alpha + beta + ( double ) ( i - 1 ) ) / 
            ( ( alpha + beta + ( double ) ( 2 * i - 1 ) ) 
            * pow ( alpha + beta + ( double ) ( 2 * i - 2 ), 2 ) 
            * ( alpha + beta + ( double ) ( 2 * i - 3 ) ) );
    }
  }

  delta = r8_gamma ( alpha        + 1.0 ) 
        * r8_gamma (         beta + 1.0 ) 
        / r8_gamma ( alpha + beta + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + beta + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );
      bn = beta / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 1.48 * an + 0.96 * bn 
        + 0.452 * an * an + 0.83 * an * bn;

      x0 = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * beta * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x0 = x0 - r1 * r2 * r3 * ( 1.0 - x0 );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * beta / 
        ( ( 6.28 + beta ) * ( double ) ( order * order ) );

      x0 = x0 - r1 * r2 * r3 * ( x[0] - x0 );
    }
    else if ( i < order - 1 )
    {
      x0 = 3.0 * x[i-2] - 3.0 * x[i-3] + x[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
    }

    jacobi_root ( &x0, order, alpha, beta, &dp2, &p1, b, c );

    x[i-1] = x0;
    w[i-1] = cc / ( dp2 * p1 );
  }
/*
  Reverse the order of the values.
*/
  for ( i = 1; i <= order/2; i++ )
  {
    temp       = x[i-1];
    x[i-1]     = x[order-i];
    x[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp       = w[i-1];
    w[i-1]     = w[order-i];
    w[order-i] = temp;
  }

  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void jacobi_compute_np ( int order, int np, double p[], double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE_NP computes a Jacobi quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    Thanks to Xu Xiang of Fudan University for pointing out that
    an earlier implementation of this routine was incorrect!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameter values.
    P[0] = ALPHA, the exponent of (1-X)
    P[1] = BETA,  the exponent of (1+X).
    -1.0 < ALPHA and -1.0 < BETA are required.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;
  double beta;

  alpha = p[0];
  beta = p[1];

  jacobi_compute ( order, alpha, beta, x, w );

  return;
}
/******************************************************************************/

void jacobi_compute_points ( int order, double alpha, double beta, 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE_POINTS computes Jacobi quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  jacobi_compute ( order, alpha, beta, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void jacobi_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE_POINTS_NP computes Jacobi quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameter values.
    P[0] = ALPHA, the exponent of (1-X)
    P[1] = BETA,  the exponent of (1+X).
    -1.0 < ALPHA and -1.0 < BETA are required.

    Output, double X[ORDER], the abscissas.
*/
{
  double alpha;
  double beta;

  alpha = p[0];
  beta = p[1];

  jacobi_compute_points ( order, alpha, beta, x );

  return;
}
/******************************************************************************/

void jacobi_compute_weights ( int order, double alpha, double beta, 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE_WEIGHTS computes Jacobi quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  jacobi_compute ( order, alpha, beta, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void jacobi_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_COMPUTE_WEIGHTS_NP computes Jacobi quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameter values.
    P[0] = ALPHA, the exponent of (1-X)
    P[1] = BETA,  the exponent of (1+X).
    -1.0 < ALPHA and -1.0 < BETA are required.

    Output, double W[ORDER], the weights.
*/
{
  double alpha;
  double beta;

  alpha = p[0];
  beta = p[1];

  jacobi_compute_weights ( order, alpha, beta, w );

  return;
}
/******************************************************************************/

double jacobi_integral ( int expon, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    JACOBI_INTEGRAL integrates a monomial with Jacobi weight.

  Discussion:

    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.

    Input, double ALPHA, the exponent of (1-X) in the weight factor.

    Input, double BETA, the exponent of (1+X) in the weight factor.

    Output, double JACOBI_INTEGRAL, the value of the integral.
*/
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;
  double value2;

  c = ( double ) ( expon );

  if ( ( expon % 2 ) == 0 )
  {
    s = +1.0;
  }
  else
  {
    s = -1.0;
  }

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + beta + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  arg1 = - beta;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value2 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = r8_gamma ( 1.0 + c ) * ( 
      s * r8_gamma ( 1.0 + beta  ) * value1 
    / r8_gamma ( 2.0 + beta  + c ) 
    +     r8_gamma ( 1.0 + alpha ) * value2 
    / r8_gamma ( 2.0 + alpha + c ) );

  return value;
}
/******************************************************************************/

void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_RECUR evaluates a Jacobi polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Output, double *P2, the value of J(ORDER)(X).

    Output, double *DP2, the value of J'(ORDER)(X).

    Output, double *P1, the value of J(ORDER-1)(X).

    Input, double X, the point at which polynomials are evaluated.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0 );
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( x - b[i-1] ) *  ( *p1 ) - c[i-1] * p0;
    *dp2 = ( x - b[i-1] ) * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
/******************************************************************************/

void jacobi_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_ROOT improves an approximate root of a Jacobi polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input/output, double *X, the approximate root, which
    should be improved on output.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.

    Output, double *DP2, the value of J'(ORDER)(X).

    Output, double *P1, the value of J(ORDER-1)(X).

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    jacobi_recur ( &p2, dp2, p1, *x, order, alpha, beta, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
/******************************************************************************/

void laguerre_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE computes a Laguerre quadrature rule.

  Discussion:

    The integration interval is [ 0, +oo ).

    The weight function is w(x) = exp ( -x );.

    If the integral to approximate is:

        Integral ( 0 <= X < +oo ) exp ( - X ) * F(X) dX

    then the quadrature rule is:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    If the integral to approximate is:

        Integral ( A <= X < +oo ) F(X) dX

    then the quadrature rule is:

      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double x0;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "LAGUERRE_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  b = ( double * ) malloc ( order * sizeof ( double ) );
  c = ( double * ) malloc ( order * sizeof ( double ) );
/*
  Set the recursion coefficients.
*/
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( double ) ( 2 * i + 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i * i );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = prod;

  for ( i = 0; i < order; i++ )
  {
/*
  Compute an estimate for the root.
*/
    if ( i == 0 )
    {
      x0 =  3.0 / ( 1.0 + 2.4 * ( double ) ( order ) );
    }
    else if ( i == 1 )
    {
      x0 = x0 + 15.0 / ( 1.0 + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      x0 = x0 + r1 * ( x0 - x[i-2] );
    }
/*
  Use iteration to find the root.
*/
    laguerre_root ( &x0, order, &dp2, &p1, b, c );
/*
  Set the abscissa and weight.
*/
    x[i] = x0;
    w[i] = ( cc / dp2 ) / p1;
  }

  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void laguerre_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE_NP computes a Laguerre quadrature rule.

  Discussion:

    The integration interval is [ 0, +oo ).

    The weight function is w(x) = exp ( -x );.

    If the integral to approximate is:

        Integral ( 0 <= X < +oo ) exp ( - X ) * F(X) dX

    then the quadrature rule is:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

    If the integral to approximate is:

        Integral ( A <= X < +oo ) F(X) dX

    then the quadrature rule is:

      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  laguerre_compute ( order, x, w );

  return;
}
/******************************************************************************/

void laguerre_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE_POINTS computes Laguerre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  laguerre_compute ( order, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void laguerre_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE_POINTS_NP computes Laguerre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  laguerre_compute_points ( order, x );

  return;
}
/******************************************************************************/

void laguerre_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE_WEIGHTS computes Laguerre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  laguerre_compute ( order, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void laguerre_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_COMPUTE_WEIGHTS_NP computes Laguerre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  laguerre_compute_weights ( order, w );

  return;
}
/******************************************************************************/

double laguerre_integral ( int expon )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.

  Discussion:

    The integral being computed is

      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.
    0 <= EXPON.

    Output, double EXACT, the value of the integral.
*/
{
  double exact;

  exact = r8_factorial ( expon );

  return exact;
}
/******************************************************************************/

void laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_RECUR evaluates a Laguerre polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Output, double *P2, the value of L(ORDER)(X).

    Output, double *DP2, the value of L'(ORDER)(X).

    Output, double *P1, the value of L(ORDER-1)(X).

    Input, double X, the point at which polynomials are evaluated.

    Input, int ORDER, the order of the polynomial to be computed.

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
/******************************************************************************/

void laguerre_root ( double *x, int order, double *dp2, double *p1, 
  double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_ROOT improves a root of a Laguerre polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.

  Reference:

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input/output, double *X, the approximate root, which
    should be improved on output.

    Input, int ORDER, the order of the polynomial to be computed.

    Output, double *DP2, the value of L'(ORDER)(X).

    Output, double *P1, the value of L(ORDER-1)(X).

    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    laguerre_recur ( &p2, dp2, p1, *x, order, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void legendre_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE computes a Legendre quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas of the rule.

    Output, double W[ORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    printf ( "\n" );
    printf ( "LEGENDRE_COMPUTE - Fatal error!\n" );
    printf ( "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 =  cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
/*
  Initial approximation H:
*/
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
/*
  Refine H using one step of Newton's method:
*/
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    x[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    w[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    x[0] = 0.0;
  }
/*
  Shift the data up.
*/
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    x[iback-1] = x[iback-ncopy-1];
    w[iback-1] = w[iback-ncopy-1];
  }
/*
  Reflect values for the negative abscissas.
*/
  for ( i = 1; i <= order - nmove; i++ )
  {
    x[i-1] = - x[order-i];
    w[i-1] = w[order-i];
  }

  return;
}
/******************************************************************************/

void legendre_compute_np ( int order, int np, double p[], double x[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_NP computes a Legendre quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas of the rule.

    Output, double W[ORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
  legendre_compute ( order, x, w );

  return;
}
/******************************************************************************/

void legendre_compute_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_POINTS computes Legendre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double X[ORDER], the abscissas.
*/
{
  double *w;

  w = ( double * ) malloc ( order * sizeof ( double ) );

  legendre_compute ( order, x, w );

  free ( w );

  return;
}
/******************************************************************************/

void legendre_compute_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_POINTS_NP computes Legendre quadrature points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  legendre_compute_points ( order, x );

  return;
}
/******************************************************************************/

void legendre_compute_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_WEIGHTS computes Legendre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Output, double W[ORDER], the weights.
*/
{
  double *x;

  x = ( double * ) malloc ( order * sizeof ( double ) );

  legendre_compute ( order, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void legendre_compute_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_COMPUTE_WEIGHTS_NP computes Legendre quadrature weights.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  legendre_compute_weights ( order, w );

  return;
}
/******************************************************************************/

double legendre_integral ( int expon )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.

  Discussion:

    To test a Legendre quadrature rule, we use it to approximate the
    integral of a monomial:

      integral ( -1 <= x <= +1 ) x^n dx
/
    This routine is given the value of the exponent, and returns the
    exact value of the integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int EXPON, the exponent.

    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
*/
{
  double exact;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    exact = 2.0 / ( double ) ( expon + 1 );
  }
  else
  {
    exact = 0.0;
  }

  return exact;
}
/******************************************************************************/

void level_growth_to_order ( int dim_num, int level[], int rule[], int growth[],
  int order[] )

/******************************************************************************/
/*
  Purpose:

    LEVEL_GROWTH_TO_ORDER: convert Level and Growth to Order.

  Discussion:

    This function is given level, rule, and growth information
    for each dimension of a quadrature rule, and determines the
    corresponding order of the rule in each dimension.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL[DIM_NUM], the 1D levels.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
     2, "F2",  Fejer Type 2, Open Fully Nested.
     3, "GP",  Gauss Patterson, Open Fully Nested.
     4, "GL",  Gauss Legendre, Open Weakly Nested.
     5, "GH",  Gauss Hermite, Open Weakly Nested.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
     7, "LG",  Gauss Laguerre, Open Non Nested.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
     9, "GJ",  Gauss Jacobi, Open Non Nested.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
    11, "HGK", Hermite Genz-Keister, Open Fully Nested.

    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
    0, "DF", default growth associated with this quadrature rule;
    1, "SL", slow linear, L+1;
    2  "SO", slow linear odd, O=1+2((L+1)/2)
    3, "ML", moderate linear, 2L+1;
    4, "SE", slow exponential;
    5, "ME", moderate exponential;
    6, "FE", full exponential.

    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
*/
{
  int dim;
  int l;
  int o;
  static int o_hgk[5] = { 1, 3, 9, 19, 35 };
  int p;
  static int p_hgk[5] = { 1, 5, 15, 29, 51 };
/*
  Check the input.
*/
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
      fprintf ( stderr, "  Negative value of LEVEL[DIM]!\n" );
      fprintf ( stderr, "  LEVEL[%d] = %d\n", dim, level[dim] );
      exit ( 1 );
    }

    if ( rule[dim] < 1 || 11 < rule[dim] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
      fprintf ( stderr, "  Illegal value of RULE[DIM]!\n" );
      fprintf ( stderr, "  RULE[%d] = %d\n", dim, rule[dim] );
      exit ( 1 );
    }

    if ( growth[dim] < 0 || 6 < growth[dim] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
      fprintf ( stderr, "  Illegal value of GROWTH[DIM]!\n" );
      fprintf ( stderr, "  GROWTH[%d] = %d\n", dim, growth[dim] );
      exit ( 1 );
    }
  }
/*
  Compute the order vector.
*/
  for ( dim = 0; dim < dim_num; dim++ )
  {
/*
  CC
  Default is Moderate Exponential Growth.
*/
    if ( rule[dim] == 1 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        if ( level[dim] == 0 )
        {
          o = 1;
        }
        else
        {
          o = 2;
          while ( o < 2 * level[dim] + 1 )
          {
            o = 2 * ( o - 1 ) + 1;
          }
        }
      }
      else if ( growth[dim] == 5 || growth[dim] == 0 )
      {
        if ( level[dim] == 0 )
        {
          o = 1;
        }
        else
        {
          o = 2;
          while ( o < 4 * level[dim] + 1 )
          {
            o = 2 * ( o - 1 ) + 1;
          }
        }
      }
      else if ( growth[dim] == 6 )
      {
        if ( level[dim] == 0 )
        {
          o = 1;
        }
        else
        {
          o = i4_power ( 2, level[dim] ) + 1;
        }
      }
    }
/*
  F2
  Default is Moderate Exponential Growth.
*/
    else if ( rule[dim] == 2 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( o < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 || growth[dim] == 0 )
      {
        o = 1;
        while ( o < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GP
  Default is Moderate Exponential Growth.
*/
    else if ( rule[dim] == 3 )
    {
      if ( growth[dim] == 1 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 1 for rule 3 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 2 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 2 for rule 3 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 3 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 3 for rule 3 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 4 )
      {
        if ( level[dim] == 0 )
        {
          o = 1;
        }
        else
        {
          p = 5;
          o = 3;
          while ( p < 2 * level[dim] + 1 )
          {
            p = 2 * p + 1;
            o = 2 * o + 1;
          }
        }
      }
      else if ( growth[dim] == 5 || growth[dim] == 0 )
      {
        if ( level[dim] == 0 )
        {
          o = 1;
        }
        else
        {
          p = 5;
          o = 3;
          while ( p < 4 * level[dim] + 1 )
          {
            p = 2 * p + 1;
            o = 2 * o + 1;
          }
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GL
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 4 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GH
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 5 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GGH
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 6 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  LG
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 7 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GLG
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 8 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GJ
  Default is Moderate Linear Growth.
*/
    else if ( rule[dim] == 9 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  GW
  Default is Moderate Linear Growth.
  We'll assume the rule is of OPEN type and that it
  has a precision typical of Gauss rules.
*/
    else if ( rule[dim] == 10 )
    {
      if ( growth[dim] == 1 )
      {
        o = level[dim] + 1;
      }
      else if ( growth[dim] == 2 )
      {
        o = 2 * ( ( level[dim] + 1 ) / 2 ) + 1;
      }
      else if ( growth[dim] == 3 || growth[dim] == 0 )
      {
        o = 2 * level[dim] + 1;
      }
      else if ( growth[dim] == 4 )
      {
        o = 1;
        while ( 2 * o - 1 < 2 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 5 )
      {
        o = 1;
        while ( 2 * o - 1 < 4 * level[dim] + 1 )
        {
          o = 2 * o + 1;
        }
      }
      else if ( growth[dim] == 6 )
      {
        o = i4_power ( 2, level[dim] + 1 ) - 1;
      }
    }
/*
  HGK
  Default is Moderate Exponential Growth.
*/
    else if ( rule[dim] == 11 )
    {
      if ( growth[dim] == 1 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 1 for rule 11 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 2 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 2 for rule 11 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 3 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "LEVEL_GROWTH_TO_ORDER - Fatal error!\n" );
        fprintf ( stderr, "  Growth rate 3 for rule 11 not available!\n" );
        exit ( 1 );
      }
      else if ( growth[dim] == 4 )
      {
        l = 0;
        p = p_hgk[l];
        o = o_hgk[l];
        while ( p < 2 * level[dim] + 1 && l < 4 )
        {
          l = l + 1;
          p = p_hgk[l];
          o = o_hgk[l];
        }
      }
      else if ( growth[dim] == 5 || growth[dim] == 0 )
      {
        l = 0;
        p = p_hgk[l];
        o = o_hgk[l];
        while ( p < 4 * level[dim] + 1 && l < 4 )
        {
          l = l + 1;
          p = p_hgk[l];
          o = o_hgk[l];
        }
      }
      else if ( growth[dim] == 6 )
      {
        l = level[dim];
        l = i4_max ( l, 0 );
        l = i4_min ( l, 4 );
        o = o_hgk[l];
      }
    }
    order[dim] = o;
  }
  return;
}
/******************************************************************************/

void level_to_order_default ( int dim_num, int level[], int rule[], 
  int order[] )

/******************************************************************************/
/*
  Purpose:

    LEVEL_TO_ORDER_DEFAULT: default growth.

  Discussion:

    The user must preallocate space for the output array ORDER.

    This function uses:

    * exponential growth rates for fully nested quadrature rules, 
      ( "CC", "F2", "GP");

    * linear growth rates for other rules.
      ( "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW" ).

    * slow exponential growth alternative for fully nested rules:
      ("CCS", "F2S", "GPS").

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL[DIM_NUM], the 1D levels.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis Slow, Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Open Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Open Fully Nested rule.

    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
*/
{
  int dim;
  int o;
  int p;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_DEFAULT - Fatal error!\n" );
      printf ( "  Negative value of LEVEL[DIM]!\n" );
      printf ( "  LEVEL[%d] = %d\n", dim, level[dim] );
      exit ( 1 );
    }
    else if ( rule[dim] == 1 )
    {
      if ( level[dim] == 0 )
      {
        order[dim] = 1;
      }
      else
      {
        order[dim] = i4_power ( 2, level[dim] ) + 1;
      }
    }
    else if ( rule[dim] == 2 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 3 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 4 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 5 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 6 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 7 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 8 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 9 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 10 )
    {
      order[dim] = 2 * level[dim] + 1;
    }
    else if ( rule[dim] == 11 )
    {
      if ( level[dim] == 0 )
      {
        o = 1;
      }
      else
      {
        o = 2;
        while ( o < 2 * level[dim] + 1 )
        {
          o = 2 * ( o - 1 ) + 1;
        }
      }
      order[dim] = o;
    }
    else if ( rule[dim] == 12 )
    {
      o = 1;
      while ( o < 2 * level[dim] + 1 )
      {
        o = 2 * o + 1;
      }
      order[dim] = o;
    }
/*
  Here, we assume that the precision requested is 2 * LEVEL(DIM) + 1,
  but because this is a Gauss-Patterson rule, the order is not equal to
  the precision.  So we first compute the level that guarantees the
  precision, then the order of that level.
*/
    else if ( rule[dim] == 13 )
    {
      if ( level[dim] == 0 )
      {
        order[dim] = 1;
      }
      else
      {
        p = 5;
        o = 3;
        while ( p < 2 * level[dim] + 1 )
        {
          p = 2 * p + 1;
          o = 2 * o + 1;
        }
        order[dim] = o;
      }
    }
    else
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_DEFAULT - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }
  }
  return;
}
/******************************************************************************/

void level_to_order_exponential ( int dim_num, int level[], int rule[], 
  int order[] )

/******************************************************************************/
/*
  Purpose:

    LEVEL_TO_ORDER_EXPONENTIAL: exponential growth.

  Discussion:

    The user must preallocate space for the output array ORDER.

    Closed rules:

      O(0) = 1
      O(L) = 2^L + 1;

      O = 1, 3, 5, 9, 17, 33, ...

    Open rules:

      O(L) = 2^(L+1) - 1;

      O = 1, 3, 7, 15, 31, 63, ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL[DIM_NUM], the 1D levels.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis Slow, Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Open Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Open Fully Nested rule.

    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
*/
{
  int dim;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!\n" );
      printf ( "  Negative value of LEVEL[DIM]!\n" ); 
      printf ( "  LEVEL[%d] = %d\n", dim, level[dim] );
      exit ( 1 );
    }
    else if ( rule[dim] == 1 )
    {
      if ( level[dim] == 0 )
      {
        order[dim] = 1;
      }
      else
      {
        order[dim] = i4_power ( 2, level[dim] ) + 1;
      }
    }
    else if ( rule[dim] == 2 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 3 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 4 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 5 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 6 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 7 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 8 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 9 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 10 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 11 )
    {
      if ( level[dim] == 0 )
      {
        order[dim] = 1;
      }
      else
      {
        order[dim] = i4_power ( 2, level[dim] ) + 1;
      }
    }
    else if ( rule[dim] == 12 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else if ( rule[dim] == 13 )
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1;
    }
    else
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }
  }
  return;
}
/******************************************************************************/

void level_to_order_exponential_slow ( int dim_num, int level[], int rule[], 
  int order[] )

/******************************************************************************/
/*
  Purpose:

    LEVEL_TO_ORDER_EXPONENTIAL_SLOW: slow exponential growth;

  Discussion:

    We seek a sequence of quadrature rules with two opposing constraints:
    * a measured rise in polynomial precision with increasing level;
    * a control on the increase in (new) points per level;

    Essentially, we are trying to keep some of the advantages of nesting,
    while moderating the cost of the explosive growth in order that occurs
    due to the repeated order doubling of nesting.

    We wish the number of points at a given level L to be "about" 2 * L + 1,
    but we also wish the rules to be completely nested.

    One way to do this is to start with a nested family of rules, whose
    order will tend to grow exponentially (doubling from one to the next),
    but simply to REPEAT each rule as many times as possible.  We move to
    the next rule only when the desired precision 2 * L + 1 exceeds the 
    precision of the current rule.

    For both the Clenshaw Curtis and Fejer Type 2 rules, if the order and
    precision are the same if the order is odd.   That is, an 11 point rule 
    will integrate exactly all polynomials up to and including degree 11.

    For Gauss Patterson rules, the relationship between order and precision
    is somewhat more complicated.  For that rule, we take the philosophy
    that at each level L, we wish to choose the rule of smallest order
    so that the precision of 2 * L + 1 is guaranteed.

     L    2*L+1  CC Order    F2 Order    GP Order/Precision

     0        1         1           1        1/1
     1        3         3           3        3/5
     2        5         5           7        3/5
     3        7         9           7        7/11
     4        9         9          15        7/11
     5       11        17          15        7/11
     6       13        17          15       15/23
     7       15        17          15       15/23
     8       17        17          31       15/23
     9       19        33          31       15/23
    10       21        33          31       15/23
    11       23        33          31       15/23
    12       25        33          31       31/47
    13       27        33          31       31/47
    14       29        33          31       31/47
    15       31        33          31       31/47
    16       33        33          63       31/47
    17       35        65          63       31/47
    18       37        65          63       31/47
    19       39        65          63       31/47
    20       41        65          63       31/47

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL[DIM_NUM], the 1D levels.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis Slow, Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Open Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Open Fully Nested rule.

    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
*/
{
  int dim;
  int o;
  int p;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_EXPONENTIAL_SLOW - Fatal error!\n" );
      printf ( "  Negative value of LEVEL[DIM]!\n" );
      printf ( "  LEVEL[%d] = %d\n", dim, level[dim] );
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 || rule[dim] == 11 )
    {
      if ( level[dim] == 0 )
      {
        o = 1;
      }
      else
      {
        o = 2;
        while ( o < 2 * level[dim] + 1 )
        {
          o = 2 * ( o - 1 ) + 1;
        }
      }
    }
    else if ( rule[dim] == 3 || rule[dim] == 13 )
    {
      if ( level[dim] == 0 )
      {
        o = 1;
      }
      else
      {
        p = 5;
        o = 3;
        while ( p < 2 * level[dim] + 1 )
        {
          p = 2 * p + 1;
          o = 2 * o + 1; 
        }
      }
    }
    else
    {
      o = 1;
      while ( o < 2 * level[dim] + 1 )
      {
        o = 2 * o + 1;
      }
    }
    order[dim] = o;
  }

  return;
}
/******************************************************************************/

void level_to_order_linear ( int dim_num, int level[], int rule[], 
  int order[] )

/******************************************************************************/
/*
  Purpose:

    LEVEL_TO_ORDER_LINEAR: linear growth.

  Discussion:

    The user must preallocate space for the output array ORDER.

      O(L) = 2 * L + 1;

      O = 1, 3, 5, 7, 9, ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL[DIM_NUM], the 1D levels.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis Slow, Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Open Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Open Fully Nested rule.

    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
*/
{
  int dim;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      printf ( "\n" );
      printf ( "LEVEL_TO_ORDER_LINEAR - Fatal error!\n" );
      printf ( "  Negative value of LEVEL[DIM]!\n" );
      printf ( "  LEVEL[%d] == %d\n", dim, level[dim] );
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    order[dim] = 2 * level[dim] + 1;
  }

  return;
}
/******************************************************************************/

void nc_compute ( int n, double x_min, double x_max, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    NC_COMPUTE computes a Newton-Cotes quadrature rule.
  
  Discussion:
  
    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
    estimates
  
      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
  
    using N abscissas XTAB(I) and a weight vector WEIGHT(I):
  
      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
  
    For the CLOSED rule, the abscissas include the end points.
    For the OPEN rule, the abscissas do not include the end points.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    17 November 2009
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order.
  
    Input, double X_MIN, X_MAX, the endpoints of the interval.
  
    Input, double XTAB[N], the abscissas.
  
    Output, double WEIGHT[N], the weights.
*/
{
  double *d;
  int i;
  int j;
  int k;
  double yvala;
  double yvalb;

  d = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
/*
  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
  and zero at the other nodes.
*/
    for ( j = 1; j <= n; j++ )
    {
      d[j-1] = 0.0;
    }
    d[i-1] = 1.0;

    for ( j = 2; j <= n; j++ )
    {
      for ( k = j; k <= n; k++ )
      {
        d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( k = 1; k <= n - j; k++ )
      {
        d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];
      }
    }
/*
  Evaluate the antiderivative of the polynomial at the endpoints.
*/
    yvala = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvala = yvala * x_min + d[j] / ( double ) ( j + 1 );
    }
    yvala = yvala * x_min;

    yvalb = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvalb = yvalb * x_max + d[j] / ( double ) ( j + 1 );
    }
    yvalb = yvalb * x_max;

    w[i-1] = yvalb - yvala;
  }

  free ( d );

  return;
}
/******************************************************************************/

void ncc_compute_points ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:
  
    NCC_COMPUTE_POINTS: points of a Newton-Cotes Closed quadrature rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    16 November 2009
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the order.
  
    Output, double X[N], the abscissas.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  if ( n == 1 )
  {
    x[0] = ( x_max + x_min ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * x_min
             + ( double ) (     i     ) * x_max )
             / ( double ) ( n     - 1 );
    }
  }
  return;
}
/******************************************************************************/

void ncc_compute_weights ( int n, double w[] )

/******************************************************************************/
/*
  Purpose:
  
    NCC_COMPUTE_WEIGHTS: weights of a Newton-Cotes Closed quadrature rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    16 November 2009
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the order.
    
    Output, double W[N], the weights.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;
  double *x;

  if ( n == 1 )
  {
    w[0] = x_max - x_min;
  }
  else
  {
    x = ( double * ) malloc ( n * sizeof ( double ) );

    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i - 1 ) * x_min
             + ( double ) (     i     ) * x_max )
             / ( double ) ( n     - 1 );
    }
    nc_compute ( n, x_min, x_max, x, w );

    free ( x );
  }
  return;
}
/******************************************************************************/

void nco_compute_points ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:
  
    NCO_COMPUTE_POINTS: points of a Newton-Cotes Open quadrature rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    17 November 2009
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the order.
  
    Output, double X[N], the abscissas.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * x_min
           + ( double ) (   + i + 1 ) * x_max )
           / ( double ) ( n     + 1 );
  }

  return;
}
/******************************************************************************/

void nco_compute_weights ( int n, double w[] )

/******************************************************************************/
/*
  Purpose:
  
    NCO_COMPUTE_WEIGHTS: weights of a Newton-Cotes Open quadrature rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    17 November 2009
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the order.
    
    Output, double W[N], the weights.
*/
{
  int i;
  double *x;
  double x_max = 1.0;
  double x_min = -1.0;

  x = ( double * ) malloc ( n * sizeof ( double ) );
  
  nco_compute_points ( n, x );

  nc_compute ( n, x_min, x_max, x, w );

  free ( x );

  return;
}
/******************************************************************************/

void patterson_lookup ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    PATTERSON_LOOKUP looks up Patterson quadrature points and weights.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1],

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Reference:

    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int N, the order.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  patterson_lookup_points ( n, x );
  patterson_lookup_weights ( n, w );

  return;
}
/******************************************************************************/

void patterson_lookup_points ( int order, double x[] )

/******************************************************************************/
/*
  Purpose:

    PATTERSON_LOOKUP_POINTS looks up Patterson quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1],

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Reference:

    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int ORDER, the order of the rule.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.

    Output, double X[ORDER], the abscissas.
*/
{
  static double x_001[1] = 
  {
     0.0 
  };
  static double x_003[3] =
  {
    -0.77459666924148337704, 
     0.0, 
     0.77459666924148337704
  };
  static double x_007[7] =
  {
    -0.96049126870802028342, 
    -0.77459666924148337704, 
    -0.43424374934680255800, 
     0.0, 
     0.43424374934680255800, 
     0.77459666924148337704, 
     0.96049126870802028342 
  };
  static double x_015[15] =
  {
    -0.99383196321275502221, 
    -0.96049126870802028342, 
    -0.88845923287225699889, 
    -0.77459666924148337704, 
    -0.62110294673722640294, 
    -0.43424374934680255800, 
    -0.22338668642896688163, 
     0.0, 
     0.22338668642896688163, 
     0.43424374934680255800, 
     0.62110294673722640294, 
     0.77459666924148337704, 
     0.88845923287225699889, 
     0.96049126870802028342, 
     0.99383196321275502221
  };
  static double x_031[31] =
  {
    -0.99909812496766759766, 
    -0.99383196321275502221, 
    -0.98153114955374010687, 
    -0.96049126870802028342, 
    -0.92965485742974005667, 
    -0.88845923287225699889, 
    -0.83672593816886873550, 
    -0.77459666924148337704, 
    -0.70249620649152707861, 
    -0.62110294673722640294, 
    -0.53131974364437562397, 
    -0.43424374934680255800, 
    -0.33113539325797683309, 
    -0.22338668642896688163, 
    -0.11248894313318662575, 
     0.0, 
     0.11248894313318662575, 
     0.22338668642896688163, 
     0.33113539325797683309, 
     0.43424374934680255800, 
     0.53131974364437562397, 
     0.62110294673722640294, 
     0.70249620649152707861, 
     0.77459666924148337704, 
     0.83672593816886873550, 
     0.88845923287225699889, 
     0.92965485742974005667, 
     0.96049126870802028342, 
     0.98153114955374010687, 
     0.99383196321275502221, 
     0.99909812496766759766
  };
  static double x_063[63] =
  {
    -0.99987288812035761194, 
    -0.99909812496766759766, 
    -0.99720625937222195908, 
    -0.99383196321275502221, 
    -0.98868475754742947994, 
    -0.98153114955374010687, 
    -0.97218287474858179658, 
    -0.96049126870802028342, 
    -0.94634285837340290515, 
    -0.92965485742974005667, 
    -0.91037115695700429250, 
    -0.88845923287225699889, 
    -0.86390793819369047715, 
    -0.83672593816886873550, 
    -0.80694053195021761186, 
    -0.77459666924148337704, 
    -0.73975604435269475868, 
    -0.70249620649152707861, 
    -0.66290966002478059546, 
    -0.62110294673722640294, 
    -0.57719571005204581484, 
    -0.53131974364437562397, 
    -0.48361802694584102756, 
    -0.43424374934680255800, 
    -0.38335932419873034692, 
    -0.33113539325797683309, 
    -0.27774982202182431507, 
    -0.22338668642896688163, 
    -0.16823525155220746498, 
    -0.11248894313318662575, 
    -0.056344313046592789972, 
     0.0, 
     0.056344313046592789972, 
     0.11248894313318662575, 
     0.16823525155220746498, 
     0.22338668642896688163, 
     0.27774982202182431507, 
     0.33113539325797683309, 
     0.38335932419873034692, 
     0.43424374934680255800, 
     0.48361802694584102756, 
     0.53131974364437562397, 
     0.57719571005204581484, 
     0.62110294673722640294, 
     0.66290966002478059546, 
     0.70249620649152707861, 
     0.73975604435269475868, 
     0.77459666924148337704, 
     0.80694053195021761186, 
     0.83672593816886873550, 
     0.86390793819369047715, 
     0.88845923287225699889, 
     0.91037115695700429250, 
     0.92965485742974005667, 
     0.94634285837340290515, 
     0.96049126870802028342, 
     0.97218287474858179658, 
     0.98153114955374010687, 
     0.98868475754742947994, 
     0.99383196321275502221, 
     0.99720625937222195908, 
     0.99909812496766759766, 
     0.99987288812035761194
  };
  static double x_127[127] =
  {
    -0.99998243035489159858, 
    -0.99987288812035761194, 
    -0.99959879967191068325, 
    -0.99909812496766759766, 
    -0.99831663531840739253, 
    -0.99720625937222195908, 
    -0.99572410469840718851, 
    -0.99383196321275502221, 
    -0.99149572117810613240, 
    -0.98868475754742947994, 
    -0.98537149959852037111, 
    -0.98153114955374010687, 
    -0.97714151463970571416, 
    -0.97218287474858179658, 
    -0.96663785155841656709, 
    -0.96049126870802028342, 
    -0.95373000642576113641, 
    -0.94634285837340290515, 
    -0.93832039777959288365, 
    -0.92965485742974005667, 
    -0.92034002547001242073, 
    -0.91037115695700429250, 
    -0.89974489977694003664, 
    -0.88845923287225699889, 
    -0.87651341448470526974, 
    -0.86390793819369047715, 
    -0.85064449476835027976, 
    -0.83672593816886873550, 
    -0.82215625436498040737, 
    -0.80694053195021761186, 
    -0.79108493379984836143, 
    -0.77459666924148337704, 
    -0.75748396638051363793, 
    -0.73975604435269475868, 
    -0.72142308537009891548, 
    -0.70249620649152707861, 
    -0.68298743109107922809, 
    -0.66290966002478059546, 
    -0.64227664250975951377, 
    -0.62110294673722640294, 
    -0.59940393024224289297, 
    -0.57719571005204581484, 
    -0.55449513263193254887, 
    -0.53131974364437562397, 
    -0.50768775753371660215, 
    -0.48361802694584102756, 
    -0.45913001198983233287, 
    -0.43424374934680255800, 
    -0.40897982122988867241, 
    -0.38335932419873034692, 
    -0.35740383783153215238, 
    -0.33113539325797683309, 
    -0.30457644155671404334, 
    -0.27774982202182431507, 
    -0.25067873030348317661, 
    -0.22338668642896688163, 
    -0.19589750271110015392, 
    -0.16823525155220746498, 
    -0.14042423315256017459, 
    -0.11248894313318662575, 
    -0.084454040083710883710, 
    -0.056344313046592789972, 
    -0.028184648949745694339, 
     0.0, 
     0.028184648949745694339, 
     0.056344313046592789972, 
     0.084454040083710883710, 
     0.11248894313318662575, 
     0.14042423315256017459, 
     0.16823525155220746498, 
     0.19589750271110015392, 
     0.22338668642896688163, 
     0.25067873030348317661, 
     0.27774982202182431507, 
     0.30457644155671404334, 
     0.33113539325797683309, 
     0.35740383783153215238, 
     0.38335932419873034692, 
     0.40897982122988867241, 
     0.43424374934680255800, 
     0.45913001198983233287, 
     0.48361802694584102756, 
     0.50768775753371660215, 
     0.53131974364437562397, 
     0.55449513263193254887, 
     0.57719571005204581484, 
     0.59940393024224289297, 
     0.62110294673722640294, 
     0.64227664250975951377, 
     0.66290966002478059546, 
     0.68298743109107922809, 
     0.70249620649152707861, 
     0.72142308537009891548, 
     0.73975604435269475868, 
     0.75748396638051363793, 
     0.77459666924148337704, 
     0.79108493379984836143, 
     0.80694053195021761186, 
     0.82215625436498040737, 
     0.83672593816886873550, 
     0.85064449476835027976, 
     0.86390793819369047715, 
     0.87651341448470526974, 
     0.88845923287225699889, 
     0.89974489977694003664, 
     0.91037115695700429250, 
     0.92034002547001242073, 
     0.92965485742974005667, 
     0.93832039777959288365, 
     0.94634285837340290515, 
     0.95373000642576113641, 
     0.96049126870802028342, 
     0.96663785155841656709, 
     0.97218287474858179658, 
     0.97714151463970571416, 
     0.98153114955374010687, 
     0.98537149959852037111, 
     0.98868475754742947994, 
     0.99149572117810613240, 
     0.99383196321275502221, 
     0.99572410469840718851, 
     0.99720625937222195908, 
     0.99831663531840739253, 
     0.99909812496766759766, 
     0.99959879967191068325, 
     0.99987288812035761194, 
     0.99998243035489159858
  };
  static double x_255[255] =
  {
    -0.99999759637974846462, 
    -0.99998243035489159858, 
    -0.99994399620705437576, 
    -0.99987288812035761194, 
    -0.99976049092443204733, 
    -0.99959879967191068325, 
    -0.99938033802502358193, 
    -0.99909812496766759766, 
    -0.99874561446809511470, 
    -0.99831663531840739253, 
    -0.99780535449595727456, 
    -0.99720625937222195908, 
    -0.99651414591489027385, 
    -0.99572410469840718851, 
    -0.99483150280062100052, 
    -0.99383196321275502221, 
    -0.99272134428278861533, 
    -0.99149572117810613240, 
    -0.99015137040077015918, 
    -0.98868475754742947994, 
    -0.98709252795403406719, 
    -0.98537149959852037111, 
    -0.98351865757863272876, 
    -0.98153114955374010687, 
    -0.97940628167086268381, 
    -0.97714151463970571416, 
    -0.97473445975240266776, 
    -0.97218287474858179658, 
    -0.96948465950245923177, 
    -0.96663785155841656709, 
    -0.96364062156981213252, 
    -0.96049126870802028342, 
    -0.95718821610986096274, 
    -0.95373000642576113641, 
    -0.95011529752129487656, 
    -0.94634285837340290515, 
    -0.94241156519108305981, 
    -0.93832039777959288365, 
    -0.93406843615772578800, 
    -0.92965485742974005667, 
    -0.92507893290707565236, 
    -0.92034002547001242073, 
    -0.91543758715576504064, 
    -0.91037115695700429250, 
    -0.90514035881326159519, 
    -0.89974489977694003664, 
    -0.89418456833555902286, 
    -0.88845923287225699889, 
    -0.88256884024734190684, 
    -0.87651341448470526974, 
    -0.87029305554811390585, 
    -0.86390793819369047715, 
    -0.85735831088623215653, 
    -0.85064449476835027976, 
    -0.84376688267270860104, 
    -0.83672593816886873550, 
    -0.82952219463740140018, 
    -0.82215625436498040737, 
    -0.81462878765513741344, 
    -0.80694053195021761186, 
    -0.79909229096084140180, 
    -0.79108493379984836143, 
    -0.78291939411828301639, 
    -0.77459666924148337704, 
    -0.76611781930376009072, 
    -0.75748396638051363793, 
    -0.74869629361693660282, 
    -0.73975604435269475868, 
    -0.73066452124218126133, 
    -0.72142308537009891548, 
    -0.71203315536225203459, 
    -0.70249620649152707861, 
    -0.69281376977911470289, 
    -0.68298743109107922809, 
    -0.67301883023041847920, 
    -0.66290966002478059546, 
    -0.65266166541001749610, 
    -0.64227664250975951377, 
    -0.63175643771119423041, 
    -0.62110294673722640294, 
    -0.61031811371518640016, 
    -0.59940393024224289297, 
    -0.58836243444766254143, 
    -0.57719571005204581484, 
    -0.56590588542365442262, 
    -0.55449513263193254887, 
    -0.54296566649831149049, 
    -0.53131974364437562397, 
    -0.51955966153745702199, 
    -0.50768775753371660215, 
    -0.49570640791876146017, 
    -0.48361802694584102756, 
    -0.47142506587165887693, 
    -0.45913001198983233287, 
    -0.44673538766202847374, 
    -0.43424374934680255800, 
    -0.42165768662616330006, 
    -0.40897982122988867241, 
    -0.39621280605761593918, 
    -0.38335932419873034692, 
    -0.37042208795007823014, 
    -0.35740383783153215238, 
    -0.34430734159943802278, 
    -0.33113539325797683309, 
    -0.31789081206847668318, 
    -0.30457644155671404334, 
    -0.29119514851824668196, 
    -0.27774982202182431507, 
    -0.26424337241092676194, 
    -0.25067873030348317661, 
    -0.23705884558982972721, 
    -0.22338668642896688163, 
    -0.20966523824318119477, 
    -0.19589750271110015392, 
    -0.18208649675925219825, 
    -0.16823525155220746498, 
    -0.15434681148137810869, 
    -0.14042423315256017459, 
    -0.12647058437230196685, 
    -0.11248894313318662575, 
    -0.098482396598119202090, 
    -0.084454040083710883710, 
    -0.070406976042855179063, 
    -0.056344313046592789972, 
    -0.042269164765363603212, 
    -0.028184648949745694339, 
    -0.014093886410782462614, 
    0.0, 
    0.014093886410782462614, 
    0.028184648949745694339, 
    0.042269164765363603212, 
    0.056344313046592789972, 
    0.070406976042855179063, 
    0.084454040083710883710, 
    0.098482396598119202090, 
    0.11248894313318662575, 
    0.12647058437230196685, 
    0.14042423315256017459, 
    0.15434681148137810869, 
    0.16823525155220746498, 
    0.18208649675925219825, 
    0.19589750271110015392, 
    0.20966523824318119477, 
    0.22338668642896688163, 
    0.23705884558982972721, 
    0.25067873030348317661, 
    0.26424337241092676194, 
    0.27774982202182431507, 
    0.29119514851824668196, 
    0.30457644155671404334, 
    0.31789081206847668318, 
    0.33113539325797683309, 
    0.34430734159943802278, 
    0.35740383783153215238, 
    0.37042208795007823014, 
    0.38335932419873034692, 
    0.39621280605761593918, 
    0.40897982122988867241, 
    0.42165768662616330006, 
    0.43424374934680255800, 
    0.44673538766202847374, 
    0.45913001198983233287, 
    0.47142506587165887693, 
    0.48361802694584102756, 
    0.49570640791876146017, 
    0.50768775753371660215, 
    0.51955966153745702199, 
    0.53131974364437562397, 
    0.54296566649831149049, 
    0.55449513263193254887, 
    0.56590588542365442262, 
    0.57719571005204581484, 
    0.58836243444766254143, 
    0.59940393024224289297, 
    0.61031811371518640016, 
    0.62110294673722640294, 
    0.63175643771119423041, 
    0.64227664250975951377, 
    0.65266166541001749610, 
    0.66290966002478059546, 
    0.67301883023041847920, 
    0.68298743109107922809, 
    0.69281376977911470289, 
    0.70249620649152707861, 
    0.71203315536225203459, 
    0.72142308537009891548, 
    0.73066452124218126133, 
    0.73975604435269475868, 
    0.74869629361693660282, 
    0.75748396638051363793, 
    0.76611781930376009072, 
    0.77459666924148337704, 
    0.78291939411828301639, 
    0.79108493379984836143, 
    0.79909229096084140180, 
    0.80694053195021761186, 
    0.81462878765513741344, 
    0.82215625436498040737, 
    0.82952219463740140018, 
    0.83672593816886873550, 
    0.84376688267270860104, 
    0.85064449476835027976, 
    0.85735831088623215653, 
    0.86390793819369047715, 
    0.87029305554811390585, 
    0.87651341448470526974, 
    0.88256884024734190684, 
    0.88845923287225699889, 
    0.89418456833555902286, 
    0.89974489977694003664, 
    0.90514035881326159519, 
    0.91037115695700429250, 
    0.91543758715576504064, 
    0.92034002547001242073, 
    0.92507893290707565236, 
    0.92965485742974005667, 
    0.93406843615772578800, 
    0.93832039777959288365, 
    0.94241156519108305981, 
    0.94634285837340290515, 
    0.95011529752129487656, 
    0.95373000642576113641, 
    0.95718821610986096274, 
    0.96049126870802028342, 
    0.96364062156981213252, 
    0.96663785155841656709, 
    0.96948465950245923177, 
    0.97218287474858179658, 
    0.97473445975240266776, 
    0.97714151463970571416, 
    0.97940628167086268381, 
    0.98153114955374010687, 
    0.98351865757863272876, 
    0.98537149959852037111, 
    0.98709252795403406719, 
    0.98868475754742947994, 
    0.99015137040077015918, 
    0.99149572117810613240, 
    0.99272134428278861533, 
    0.99383196321275502221, 
    0.99483150280062100052, 
    0.99572410469840718851, 
    0.99651414591489027385, 
    0.99720625937222195908, 
    0.99780535449595727456, 
    0.99831663531840739253, 
    0.99874561446809511470, 
    0.99909812496766759766, 
    0.99938033802502358193, 
    0.99959879967191068325, 
    0.99976049092443204733, 
    0.99987288812035761194, 
    0.99994399620705437576, 
    0.99998243035489159858, 
    0.99999759637974846462 
  };

  if ( order == 1 )
  {
    r8vec_copy ( order, x_001, x );
  }
  else if ( order == 3 )
  {
    r8vec_copy ( order, x_003, x );
  }
  else if ( order == 7 )
  {
    r8vec_copy ( order, x_007, x );
  }
  else if ( order == 15 )
  {
    r8vec_copy ( order, x_015, x );
  }
  else if ( order == 31 )
  {
    r8vec_copy ( order, x_031, x );
  }
  else if ( order == 63 )
  {
    r8vec_copy ( order, x_063, x );
  }
  else if ( order == 127 )
  {
    r8vec_copy ( order, x_127, x );
  }
  else if ( order == 255 )
  {
    r8vec_copy ( order, x_255, x );
  }
  else
  {
    printf ( "\n" );
    printf ( "PATTERSON_LOOKUP_POINTS - Fatal error!\n" );
    printf ( "  Unexpected value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void patterson_lookup_points_np ( int order, int np, double p[], double x[] )

/******************************************************************************/
/*
  Purpose:

    PATTERSON_LOOKUP_POINTS_NP looks up Patterson quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    The rule is defined on [-1,1],

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Reference:

    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.

    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.

  Parameters:

    Input, int ORDER, the order of the rule.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double X[ORDER], the abscissas.
*/
{
  patterson_lookup_points ( order, x );

  return;
}
/******************************************************************************/

void patterson_lookup_weights ( int order, double w[] )

/******************************************************************************/
/*
  Purpose:

    PATTERSON_LOOKUP_WEIGHTS looks up Patterson quadrature weights.

  Discussion:

    The allowed orders are 1, 3, 7, 15, 31, 63, 127 and 255.

    The weights are positive, symmetric and should sum to 2.

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 December 2009

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    ORDER must be 1, 3, 7, 15, 31, 63, 127 and 255.

    Output, double W[ORDER], the weights.
*/
{
  static double w_001[1] =
  {
    2.0
  };
  static double w_003[3] = 
  {
    0.555555555555555555556,
    0.888888888888888888889,
    0.555555555555555555556
  };
  static double w_007[7] =
  {
    0.104656226026467265194,
    0.268488089868333440729,
    0.401397414775962222905,
    0.450916538658474142345,
    0.401397414775962222905,
    0.268488089868333440729,
    0.104656226026467265194
  };
  static double w_015[15] =
  {
    0.0170017196299402603390,
    0.0516032829970797396969,
    0.0929271953151245376859,
    0.134415255243784220360,
    0.171511909136391380787,
    0.200628529376989021034,
    0.219156858401587496404,
    0.225510499798206687386,
    0.219156858401587496404,
    0.200628529376989021034,
    0.171511909136391380787,
    0.134415255243784220360,
    0.0929271953151245376859,
    0.0516032829970797396969,
    0.0170017196299402603390
  };
  static double w_031[31] =
  {
    0.00254478079156187441540,
    0.00843456573932110624631,
    0.0164460498543878109338,
    0.0258075980961766535646,
    0.0359571033071293220968,
    0.0464628932617579865414,
    0.0569795094941233574122,
    0.0672077542959907035404,
    0.0768796204990035310427,
    0.0857559200499903511542,
    0.0936271099812644736167,
    0.100314278611795578771,
    0.105669893580234809744,
    0.109578421055924638237,
    0.111956873020953456880,
    0.112755256720768691607,
    0.111956873020953456880,
    0.109578421055924638237,
    0.105669893580234809744,
    0.100314278611795578771,
    0.0936271099812644736167,
    0.0857559200499903511542,
    0.0768796204990035310427,
    0.0672077542959907035404,
    0.0569795094941233574122,
    0.0464628932617579865414,
    0.0359571033071293220968,
    0.0258075980961766535646,
    0.0164460498543878109338,
    0.00843456573932110624631,
    0.00254478079156187441540
  };
  static double w_063[63] =
  {
    0.000363221481845530659694,
    0.00126515655623006801137,
    0.00257904979468568827243,
    0.00421763044155885483908,
    0.00611550682211724633968,
    0.00822300795723592966926,
    0.0104982469096213218983,
    0.0129038001003512656260,
    0.0154067504665594978021,
    0.0179785515681282703329,
    0.0205942339159127111492,
    0.0232314466399102694433,
    0.0258696793272147469108,
    0.0284897547458335486125,
    0.0310735511116879648799,
    0.0336038771482077305417,
    0.0360644327807825726401,
    0.0384398102494555320386,
    0.0407155101169443189339,
    0.0428779600250077344929,
    0.0449145316536321974143,
    0.0468135549906280124026,
    0.0485643304066731987159,
    0.0501571393058995374137,
    0.0515832539520484587768,
    0.0528349467901165198621,
    0.0539054993352660639269,
    0.0547892105279628650322,
    0.0554814043565593639878,
    0.0559784365104763194076,
    0.0562776998312543012726,
    0.0563776283603847173877,
    0.0562776998312543012726,
    0.0559784365104763194076,
    0.0554814043565593639878,
    0.0547892105279628650322,
    0.0539054993352660639269,
    0.0528349467901165198621,
    0.0515832539520484587768,
    0.0501571393058995374137,
    0.0485643304066731987159,
    0.0468135549906280124026,
    0.0449145316536321974143,
    0.0428779600250077344929,
    0.0407155101169443189339,
    0.0384398102494555320386,
    0.0360644327807825726401,
    0.0336038771482077305417,
    0.0310735511116879648799,
    0.0284897547458335486125,
    0.0258696793272147469108,
    0.0232314466399102694433,
    0.0205942339159127111492,
    0.0179785515681282703329,
    0.0154067504665594978021,
    0.0129038001003512656260,
    0.0104982469096213218983,
    0.00822300795723592966926,
    0.00611550682211724633968,
    0.00421763044155885483908,
    0.00257904979468568827243,
    0.00126515655623006801137,
    0.000363221481845530659694
  };
  static double w_127[127] =
  {
    0.0000505360952078625176247,
    0.000180739564445388357820,
    0.000377746646326984660274,
    0.000632607319362633544219,
    0.000938369848542381500794,
    0.00128952408261041739210,
    0.00168114286542146990631,
    0.00210881524572663287933,
    0.00256876494379402037313,
    0.00305775341017553113613,
    0.00357289278351729964938,
    0.00411150397865469304717,
    0.00467105037211432174741,
    0.00524912345480885912513,
    0.00584344987583563950756,
    0.00645190005017573692280,
    0.00707248999543355546805,
    0.00770337523327974184817,
    0.00834283875396815770558,
    0.00898927578406413572328,
    0.00964117772970253669530,
    0.0102971169579563555237,
    0.0109557333878379016480,
    0.0116157233199551347270,
    0.0122758305600827700870,
    0.0129348396636073734547,
    0.0135915710097655467896,
    0.0142448773729167743063,
    0.0148936416648151820348,
    0.0155367755558439824399,
    0.0161732187295777199419,
    0.0168019385741038652709,
    0.0174219301594641737472,
    0.0180322163903912863201,
    0.0186318482561387901863,
    0.0192199051247277660193,
    0.0197954950480974994880,
    0.0203577550584721594669,
    0.0209058514458120238522,
    0.0214389800125038672465,
    0.0219563663053178249393,
    0.0224572658268160987071,
    0.0229409642293877487608,
    0.0234067774953140062013,
    0.0238540521060385400804,
    0.0242821652033365993580,
    0.0246905247444876769091,
    0.0250785696529497687068,
    0.0254457699654647658126,
    0.0257916269760242293884,
    0.0261156733767060976805,
    0.0264174733950582599310,
    0.0266966229274503599062,
    0.0269527496676330319634,
    0.0271855132296247918192,
    0.0273946052639814325161,
    0.0275797495664818730349,
    0.0277407021782796819939,
    0.0278772514766137016085,
    0.0279892182552381597038,
    0.0280764557938172466068,
    0.0281388499156271506363,
    0.0281763190330166021307,
    0.0281888141801923586938,
    0.0281763190330166021307,
    0.0281388499156271506363,
    0.0280764557938172466068,
    0.0279892182552381597038,
    0.0278772514766137016085,
    0.0277407021782796819939,
    0.0275797495664818730349,
    0.0273946052639814325161,
    0.0271855132296247918192,
    0.0269527496676330319634,
    0.0266966229274503599062,
    0.0264174733950582599310,
    0.0261156733767060976805,
    0.0257916269760242293884,
    0.0254457699654647658126,
    0.0250785696529497687068,
    0.0246905247444876769091,
    0.0242821652033365993580,
    0.0238540521060385400804,
    0.0234067774953140062013,
    0.0229409642293877487608,
    0.0224572658268160987071,
    0.0219563663053178249393,
    0.0214389800125038672465,
    0.0209058514458120238522,
    0.0203577550584721594669,
    0.0197954950480974994880,
    0.0192199051247277660193,
    0.0186318482561387901863,
    0.0180322163903912863201,
    0.0174219301594641737472,
    0.0168019385741038652709,
    0.0161732187295777199419,
    0.0155367755558439824399,
    0.0148936416648151820348,
    0.0142448773729167743063,
    0.0135915710097655467896,
    0.0129348396636073734547,
    0.0122758305600827700870,
    0.0116157233199551347270,
    0.0109557333878379016480,
    0.0102971169579563555237,
    0.00964117772970253669530,
    0.00898927578406413572328,
    0.00834283875396815770558,
    0.00770337523327974184817,
    0.00707248999543355546805,
    0.00645190005017573692280,
    0.00584344987583563950756,
    0.00524912345480885912513,
    0.00467105037211432174741,
    0.00411150397865469304717,
    0.00357289278351729964938,
    0.00305775341017553113613,
    0.00256876494379402037313,
    0.00210881524572663287933,
    0.00168114286542146990631,
    0.00128952408261041739210,
    0.000938369848542381500794,
    0.000632607319362633544219,
    0.000377746646326984660274,
    0.000180739564445388357820,
    0.0000505360952078625176247
  };
  static double w_255[255] =
  {
    0.69379364324108267170E-05,
    0.25157870384280661489E-04,
    0.53275293669780613125E-04,
    0.90372734658751149261E-04,
    0.13575491094922871973E-03,
    0.18887326450650491366E-03,
    0.24921240048299729402E-03,
    0.31630366082226447689E-03,
    0.38974528447328229322E-03,
    0.46918492424785040975E-03,
    0.55429531493037471492E-03,
    0.64476204130572477933E-03,
    0.74028280424450333046E-03,
    0.84057143271072246365E-03,
    0.94536151685852538246E-03,
    0.10544076228633167722E-02,
    0.11674841174299594077E-02,
    0.12843824718970101768E-02,
    0.14049079956551446427E-02,
    0.15288767050877655684E-02,
    0.16561127281544526052E-02,
    0.17864463917586498247E-02,
    0.19197129710138724125E-02,
    0.20557519893273465236E-02,
    0.21944069253638388388E-02,
    0.23355251860571608737E-02,
    0.24789582266575679307E-02,
    0.26245617274044295626E-02,
    0.27721957645934509940E-02,
    0.29217249379178197538E-02,
    0.30730184347025783234E-02,
    0.32259500250878684614E-02,
    0.33803979910869203823E-02,
    0.35362449977167777340E-02,
    0.36933779170256508183E-02,
    0.38516876166398709241E-02,
    0.40110687240750233989E-02,
    0.41714193769840788528E-02,
    0.43326409680929828545E-02,
    0.44946378920320678616E-02,
    0.46573172997568547773E-02,
    0.48205888648512683476E-02,
    0.49843645647655386012E-02,
    0.51485584789781777618E-02,
    0.53130866051870565663E-02,
    0.54778666939189508240E-02,
    0.56428181013844441585E-02,
    0.58078616599775673635E-02,
    0.59729195655081658049E-02,
    0.61379152800413850435E-02,
    0.63027734490857587172E-02,
    0.64674198318036867274E-02,
    0.66317812429018878941E-02,
    0.67957855048827733948E-02,
    0.69593614093904229394E-02,
    0.71224386864583871532E-02,
    0.72849479805538070639E-02,
    0.74468208324075910174E-02,
    0.76079896657190565832E-02,
    0.77683877779219912200E-02,
    0.79279493342948491103E-02,
    0.80866093647888599710E-02,
    0.82443037630328680306E-02,
    0.84009692870519326354E-02,
    0.85565435613076896192E-02,
    0.87109650797320868736E-02,
    0.88641732094824942641E-02,
    0.90161081951956431600E-02,
    0.91667111635607884067E-02,
    0.93159241280693950932E-02,
    0.94636899938300652943E-02,
    0.96099525623638830097E-02,
    0.97546565363174114611E-02,
    0.98977475240487497440E-02,
    0.10039172044056840798E-01,
    0.10178877529236079733E-01,
    0.10316812330947621682E-01,
    0.10452925722906011926E-01,
    0.10587167904885197931E-01,
    0.10719490006251933623E-01,
    0.10849844089337314099E-01,
    0.10978183152658912470E-01,
    0.11104461134006926537E-01,
    0.11228632913408049354E-01,
    0.11350654315980596602E-01,
    0.11470482114693874380E-01,
    0.11588074033043952568E-01,
    0.11703388747657003101E-01,
    0.11816385890830235763E-01,
    0.11927026053019270040E-01,
    0.12035270785279562630E-01,
    0.12141082601668299679E-01,
    0.12244424981611985899E-01,
    0.12345262372243838455E-01,
    0.12443560190714035263E-01,
    0.12539284826474884353E-01,
    0.12632403643542078765E-01,
    0.12722884982732382906E-01,
    0.12810698163877361967E-01,
    0.12895813488012114694E-01,
    0.12978202239537399286E-01,
    0.13057836688353048840E-01,
    0.13134690091960152836E-01,
    0.13208736697529129966E-01,
    0.13279951743930530650E-01,
    0.13348311463725179953E-01,
    0.13413793085110098513E-01,
    0.13476374833816515982E-01,
    0.13536035934956213614E-01,
    0.13592756614812395910E-01,
    0.13646518102571291428E-01,
    0.13697302631990716258E-01,
    0.13745093443001896632E-01,
    0.13789874783240936517E-01,
    0.13831631909506428676E-01,
    0.13870351089139840997E-01,
    0.13906019601325461264E-01,
    0.13938625738306850804E-01,
    0.13968158806516938516E-01,
    0.13994609127619079852E-01,
    0.14017968039456608810E-01,
    0.14038227896908623303E-01,
    0.14055382072649964277E-01,
    0.14069424957813575318E-01,
    0.14080351962553661325E-01,
    0.14088159516508301065E-01,
    0.14092845069160408355E-01,
    0.14094407090096179347E-01,
    0.14092845069160408355E-01,
    0.14088159516508301065E-01,
    0.14080351962553661325E-01,
    0.14069424957813575318E-01,
    0.14055382072649964277E-01,
    0.14038227896908623303E-01,
    0.14017968039456608810E-01,
    0.13994609127619079852E-01,
    0.13968158806516938516E-01,
    0.13938625738306850804E-01,
    0.13906019601325461264E-01,
    0.13870351089139840997E-01,
    0.13831631909506428676E-01,
    0.13789874783240936517E-01,
    0.13745093443001896632E-01,
    0.13697302631990716258E-01,
    0.13646518102571291428E-01,
    0.13592756614812395910E-01,
    0.13536035934956213614E-01,
    0.13476374833816515982E-01,
    0.13413793085110098513E-01,
    0.13348311463725179953E-01,
    0.13279951743930530650E-01,
    0.13208736697529129966E-01,
    0.13134690091960152836E-01,
    0.13057836688353048840E-01,
    0.12978202239537399286E-01,
    0.12895813488012114694E-01,
    0.12810698163877361967E-01,
    0.12722884982732382906E-01,
    0.12632403643542078765E-01,
    0.12539284826474884353E-01,
    0.12443560190714035263E-01,
    0.12345262372243838455E-01,
    0.12244424981611985899E-01,
    0.12141082601668299679E-01,
    0.12035270785279562630E-01,
    0.11927026053019270040E-01,
    0.11816385890830235763E-01,
    0.11703388747657003101E-01,
    0.11588074033043952568E-01,
    0.11470482114693874380E-01,
    0.11350654315980596602E-01,
    0.11228632913408049354E-01,
    0.11104461134006926537E-01,
    0.10978183152658912470E-01,
    0.10849844089337314099E-01,
    0.10719490006251933623E-01,
    0.10587167904885197931E-01,
    0.10452925722906011926E-01,
    0.10316812330947621682E-01,
    0.10178877529236079733E-01,
    0.10039172044056840798E-01,
    0.98977475240487497440E-02,
    0.97546565363174114611E-02,
    0.96099525623638830097E-02,
    0.94636899938300652943E-02,
    0.93159241280693950932E-02,
    0.91667111635607884067E-02,
    0.90161081951956431600E-02,
    0.88641732094824942641E-02,
    0.87109650797320868736E-02,
    0.85565435613076896192E-02,
    0.84009692870519326354E-02,
    0.82443037630328680306E-02,
    0.80866093647888599710E-02,
    0.79279493342948491103E-02,
    0.77683877779219912200E-02,
    0.76079896657190565832E-02,
    0.74468208324075910174E-02,
    0.72849479805538070639E-02,
    0.71224386864583871532E-02,
    0.69593614093904229394E-02,
    0.67957855048827733948E-02,
    0.66317812429018878941E-02,
    0.64674198318036867274E-02,
    0.63027734490857587172E-02,
    0.61379152800413850435E-02,
    0.59729195655081658049E-02,
    0.58078616599775673635E-02,
    0.56428181013844441585E-02,
    0.54778666939189508240E-02,
    0.53130866051870565663E-02,
    0.51485584789781777618E-02,
    0.49843645647655386012E-02,
    0.48205888648512683476E-02,
    0.46573172997568547773E-02,
    0.44946378920320678616E-02,
    0.43326409680929828545E-02,
    0.41714193769840788528E-02,
    0.40110687240750233989E-02,
    0.38516876166398709241E-02,
    0.36933779170256508183E-02,
    0.35362449977167777340E-02,
    0.33803979910869203823E-02,
    0.32259500250878684614E-02,
    0.30730184347025783234E-02,
    0.29217249379178197538E-02,
    0.27721957645934509940E-02,
    0.26245617274044295626E-02,
    0.24789582266575679307E-02,
    0.23355251860571608737E-02,
    0.21944069253638388388E-02,
    0.20557519893273465236E-02,
    0.19197129710138724125E-02,
    0.17864463917586498247E-02,
    0.16561127281544526052E-02,
    0.15288767050877655684E-02,
    0.14049079956551446427E-02,
    0.12843824718970101768E-02,
    0.11674841174299594077E-02,
    0.10544076228633167722E-02,
    0.94536151685852538246E-03,
    0.84057143271072246365E-03,
    0.74028280424450333046E-03,
    0.64476204130572477933E-03,
    0.55429531493037471492E-03,
    0.46918492424785040975E-03,
    0.38974528447328229322E-03,
    0.31630366082226447689E-03,
    0.24921240048299729402E-03,
    0.18887326450650491366E-03,
    0.13575491094922871973E-03,
    0.90372734658751149261E-04,
    0.53275293669780613125E-04,
    0.25157870384280661489E-04,
    0.69379364324108267170E-05
  };
  if ( order == 1 )
  {
    r8vec_copy ( order, w_001, w );
  }
  else if ( order == 3 )
  {
    r8vec_copy ( order, w_003, w );
  }
  else if ( order == 7 )
  {
    r8vec_copy ( order, w_007, w );
  }
  else if ( order == 15 )
  {
    r8vec_copy ( order, w_015, w );
  }
  else if ( order == 31 )
  {
    r8vec_copy ( order, w_031, w );
  }
  else if ( order == 63 )
  {
    r8vec_copy ( order, w_063, w );
  }
  else if ( order == 127 )
  {
    r8vec_copy ( order, w_127, w );
  }
  else if ( order == 255 )
  {
    r8vec_copy ( order, w_255, w );
  }
  else
  {
    printf ( "\n" );
    printf ( "PATTERSON_LOOKUP_WEIGHTS - Fatal error!\n" );
    printf ( "  Unexpected value of ORDER = %d.\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void patterson_lookup_weights_np ( int order, int np, double p[], double w[] )

/******************************************************************************/
/*
  Purpose:

    PATTERSON_LOOKUP_WEIGHTS_NP looks up Patterson quadrature weights.

  Discussion:

    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.

    The weights are positive, symmetric and should sum to 2.

    The user must preallocate space for the output array W.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.

  Parameters:

    Input, int ORDER, the order of the rule.
    ORDER must be 1, 3, 7, 15, 31, 63 or 127.

    Input, int NP, the number of parameters.

    Input, double P[NP], parameters which are not needed by this function.

    Output, double W[ORDER], the weights.
*/
{
  patterson_lookup_weights ( order, w );

  return;
}
/******************************************************************************/

int point_radial_tol_unique_count ( int m, int n, double a[], double tol, 
  int *seed )

/******************************************************************************/
/*
  Purpose:

    POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.

  Discussion:

    The input data is an M x N array A, representing the M-dimensional
    coordinates of N points.

    The output is the number of tolerably unique points in the list.

    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
    But that program is guaranteed to use N^2 comparisons.

    It is hoped that this function, on the other hand, will tend
    to use O(N) comparisons after an O(NLog(N)) sort.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, double A[M*N], the array of N columns of data.

    Input, double TOL, a tolerance for equality.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, int POINT_RADIAL_TOL_UNIQUE_COUNT, the number of tolerably
    unique points.
*/
{
  double dist;
  int hi;
  int i;
  int *indx;
  int j;
  int k;
  double *r;
  int *unique;
  int unique_num;
  double *w;
  double w_sum;
  double *z;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
/*
  Assign a base point Z randomly in the convex hull.
*/
  w = r8vec_uniform_01_new ( n, seed );
  w_sum = r8vec_sum ( n, w );
  for ( j = 0; j < n; j++ )
  {
    w[j] = w[j] / w_sum;
  }

  z = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    z[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      z[i] = z[i] + a[i+j*m] * w[j];
    }
  }
/*
  Compute the radial distance R of each point to Z.
*/
  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r[j] = r[j] + pow ( a[i+j*m] - z[i], 2 );
    }
    r[j] = sqrt ( r[j] );
  }
/*
  Implicitly sort the R array.
*/
  indx = r8vec_sort_heap_index_a ( n, r );
/*
  To determine if a point I is tolerably unique, we only have to check
  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
*/
  unique_num = 0;

  unique = ( int * ) malloc ( n * sizeof ( int ) );
  for ( i = 0; i < n; i++ )
  {
    unique[i] = 1;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( unique[indx[i]] )
    {
/*
  Point INDX(I) is unique, in that no earlier point is near it.
*/
      unique_num = unique_num + 1;
/*
  Look for later points which are close to point INDX(I)
  in terms of R.
*/
      hi = i;

      while ( hi < n - 1 )
      {
        if ( r[indx[i]] + tol < r[indx[hi+1]] )
        {
          break;
        }
        hi = hi + 1;
      }
/*
  Points INDX(I+1) through INDX(HI) have an R value close to 
  point INDX(I).  Are they truly close to point INDEX(I)?
*/
      for ( j = i + 1; j <= hi; j++ )
      {
        if ( unique[indx[j]] )
        {
          dist = 0.0;
          for ( k = 0; k < m; k++ )
          {
            dist = dist + pow ( a[k+indx[i]*m] - a[k+indx[j]*m], 2 );
          }
          dist = sqrt ( dist );

          if ( dist <= tol )
          {
            unique[indx[j]] = 0;
          }
        }
      }
    }
  }

  free ( indx );
  free ( r );
  free ( unique );
  free ( w );
  free ( z );

  return unique_num;
}
/******************************************************************************/

int point_radial_tol_unique_index ( int m, int n, double a[], double tol, 
  int *seed, int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.

  Discussion:

    The input data is an M x N array A, representing the M-dimensional
    coordinates of N points.

    The output is:
    * the number of tolerably unique points in the list;
    * the index, in the list of unique items, of the representatives 
      of each point;
    * the index, in A, of the tolerably unique representatives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, double A[M*N], the array of N columns of data.

    Input, double TOL, a tolerance for equality.

    Input/output, int SEED, a seed for the random
    number generator.

    Output, int UNDX[UNIQUE_NUM], the index, in A, of the 
    tolerably unique points.

    Output, int XDNU[N], the index, in UNDX, of the 
    tolerably unique point that "represents" this point.

    Output, int POINT_RADIAL_TOL_UNIQUE_INDEX, the number of tolerably
    unique points.
*/
{
  double dist;
  int hi;
  int i;
  int *indx;
  int j;
  int k;
  double *r;
  int *unique;
  int unique_num;
  double *w;
  double w_sum;
  double *z;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
/*
  Assign a base point Z randomly in the convex hull.
*/
  w = r8vec_uniform_01_new ( n, seed );
  w_sum = r8vec_sum ( n, w );
  for ( j = 0; j < n; j++ )
  {
    w[j] = w[j] / w_sum;
  }

  z = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    z[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      z[i] = z[i] + a[i+j*m] * w[j];
    }
  }
/*
  Compute the radial distance R of each point to Z.
*/
  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r[j] = r[j] + pow ( a[i+j*m] - z[i], 2 );
    }
    r[j] = sqrt ( r[j] );
  }
/*
  Implicitly sort the R array.
*/
  indx = r8vec_sort_heap_index_a ( n, r );
/*
  To determine if a point I is tolerably unique, we only have to check
  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
*/
  unique_num = 0;

  unique = ( int * ) malloc ( n * sizeof ( int ) );
  for ( i = 0; i < n; i++ )
  {
    unique[i] = 1;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( unique[indx[i]] )
    {
/*
  Point INDX(I) is unique, in that no earlier point is near it.
*/
      xdnu[indx[i]] = unique_num;
      undx[unique_num] = indx[i];
      unique_num = unique_num + 1;
/*
  Look for later points which are close to point INDX(I)
  in terms of R.
*/
      hi = i;

      while ( hi < n - 1 )
      {
        if ( r[indx[i]] + tol < r[indx[hi+1]] )
        {
          break;
        }
        hi = hi + 1;
      }
/*
  Points INDX(I+1) through INDX(HI) have an R value close to 
  point INDX(I).  Are they truly close to point INDEX(I)?
*/
      for ( j = i + 1; j <= hi; j++ )
      {
        if ( unique[indx[j]] )
        {
          dist = 0.0;
          for ( k = 0; k < m; k++ )
          {
            dist = dist + pow ( a[k+indx[i]*m] - a[k+indx[j]*m], 2 );
          }
          dist = sqrt ( dist );

          if ( dist <= tol )
          {
            unique[indx[j]] = 0;
            xdnu[indx[j]] = xdnu[indx[i]];
          }
        }
      }
    }
  }

  free ( indx );
  free ( r );
  free ( unique );
  free ( w );
  free ( z );

  return unique_num;
}
/******************************************************************************/

void point_unique_index ( int m, int n, double a[], int unique_num, int undx[], 
  int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    POINT_UNIQUE_INDEX indexes unique points.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    The first step of the algorithm requires the indexed sorting
    of A, which creates arrays INDX and XDNI.  (If all the entries
    of A are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of A in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      X(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     A  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is A(3).
    XDNI(2) = 5 means that A(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
    XDNU(8) = 2 means that A(8) is at unique sorted item(2).

    XU(XDNU(I))) = A(I).
    XU(I)        = A(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  int base = 0;
  double diff;
  int i;
  int *indx;
  int j;
  int k;
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Walk through the implicitly sorted array X.
*/
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r8_max ( diff, r8_abs ( a[k+indx[i]*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  free ( indx );

  return;
}
/******************************************************************************/

void product_mixed_weight ( int dim_num, int order_1d[], int order_nd, 
  int rule[], double alpha[], double beta[], double weight_nd[] )

/******************************************************************************/
/*
  Purpose:

    PRODUCT_MIXED_WEIGHT computes the weights of a mixed product rule.

  Discussion:

    This routine computes the weights for a quadrature rule which is
    a product of 1D rules of varying order and kind.

    The user must preallocate space for the output array WEIGHT_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.

    Input, int ORDER_ND, the order of the product rule.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis Slow, Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Open Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Open Fully Nested rule.

    Input, double ALPHA[DIM_NUM], BETA[DIM_NUM], parameters used for
    Generalized Gauss Hermite, Generalized Gauss Laguerre, 
    and Gauss Jacobi rules.

    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
*/
{
  int dim;
  int i;
  double *weight_1d;

  for ( i = 0; i < order_nd; i++ )
  {
    weight_nd[i] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    weight_1d = ( double * ) malloc ( order_1d[dim] * sizeof ( double ) );

    if ( rule[dim] == 1 )
    {
      clenshaw_curtis_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 2 )
    {
      fejer2_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 3 )
    {
      patterson_lookup_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 4 )
    {
      legendre_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 5 )
    {
      hermite_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 6 )
    {
      gen_hermite_compute_weights ( 
        order_1d[dim], alpha[dim], weight_1d );
    }
    else if ( rule[dim] == 7 )
    {
      laguerre_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 8 )
    {
      gen_laguerre_compute_weights ( 
        order_1d[dim], alpha[dim], weight_1d );
    }
    else if ( rule[dim] == 9 )
    {
      jacobi_compute_weights ( 
        order_1d[dim], alpha[dim], beta[dim], weight_1d );
    }
    else if ( rule[dim] == 10 )
    {
      printf ( "\n" );
      printf ( "PRODUCT_MIXED_WEIGHT - Fatal error!\n" );
      printf ( "  Do not know how to set weights for rule 10.\n" );
      exit ( 1 );
    }
    else if ( rule[dim] == 11 )
    {
      clenshaw_curtis_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 12 )
    {
      fejer2_compute_weights ( 
        order_1d[dim], weight_1d );
    }
    else if ( rule[dim] == 13 )
    {
      patterson_lookup_weights ( 
        order_1d[dim], weight_1d );
    }
    else
    {
      printf ( "\n" );
      printf ( "PRODUCT_MIXED_WEIGHT - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }

    r8vec_direct_product2 ( dim, order_1d[dim], weight_1d, 
      dim_num, order_nd, weight_nd );

    free ( weight_1d );
  }
  return;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

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
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_ceiling ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_CEILING returns the "ceiling" of an R8.

  Discussion:

    The "ceiling" of X is the value of X rounded towards plus infinity.

  Example:

    X        R8_CEILING(X)

   -1.1      -1.0
   -1.0      -1.0
   -0.9       0.0
   -0.1       0.0
    0.0       0.0
    0.1       1.0
    0.9       1.0
    1.0       1.0
    1.1       2.0
    2.9       3.0
    3.0       3.0
    3.14159   4.0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose ceiling is desired.

    Output, double R8_CEILING, the ceiling of X.
*/
{
  double value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1.0;
  }

  return value;
}
/******************************************************************************/

double r8_choose ( int n, int k )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.

  Discussion:

    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in R8 arithmetic.

    The formula used is:

      C(N,K) = N / ( K * (N-K) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Reference:

    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.

  Parameters:

    Input, int N, K, the values of N and K.

    Output, double R8_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_factorial ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL computes the factorial of N.

  Discussion:

    factorial ( N ) = product ( 1 <= I <= N ) I

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.

    Output, double R8_FACTORIAL, the factorial of N.
*/
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
/******************************************************************************/

double r8_factorial2 ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL2 computes the double factorial function.

  Discussion:

    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)

  Example:

     N    Factorial2(N)

     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the double factorial 
    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.

    Output, double R8_FACTORIAL2, the value of Factorial2(N).
*/
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
/******************************************************************************/

double r8_floor ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_FLOOR rounds an R8 "down" (towards -infinity) to the next integer.

  Example:

    X        R8_FLOOR(X)

   -1.1      -2.0
   -1.0      -1.0
   -0.9      -1.0
   -0.1      -1.0
    0.0       0.0
    0.1       0.0
    0.9       0.0
    1.0       1.0
    1.1       1.0
    2.9       2.0
    3.0       3.0
    3.14159   3.0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose floor is desired.

    Output, double R8_FLOOR, the floor of X.
*/
{
  double value;

  value = ( int ) x;

  if ( x < value )
  {
    value = value - 1.0;
  }

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2008

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
/*
  Coefficients for minimax approximation over (12, INF).
*/
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = 0;
  fact = one;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_huge ( )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

double r8_hyper_2f1 ( double a, double b, double c, double x )

/******************************************************************************/
/*
  Purpose:

    R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).

  Discussion:

    A bug was corrected.  A line which read
      c1 = - ( - 1.0, m ) * gc / ( gam * gbm * rm );
    was corrected to read
      c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );
    JVB, 05 July 2009.

    A minor bug was corrected.  The HW variable, used in several places as
    the "old" value of a quantity being iteratively improved, was not
    being initialized.  JVB, 11 February 2008.

    The FORTRAN77 original version of this routine is copyrighted by
    Shanjie Zhang and Jianming Jin.  However, they give permission to
    incorporate this routine into a user program provided that the copyright
    is acknowledged.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    C version by John Burkardt.

  Reference:

    Shanjie Zhang, Jianming Jin,
    Computation of Special Functions,
    Wiley, 1996,
    ISBN: 0-471-11963-6,
    LC: QA351.C45

  Parameters:

    Input, double A, B, C, X, the arguments of the function.
    C must not be equal to a nonpositive integer.
    X < 1.

    Output, double R8_HYPER_2F1, the value of the function.
*/
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  int l0;
  int l1;
  int l2;
  int l3;
  int l4;
  int l5;
  int m;
  int nm;
  double pa;
  double pb;
  double pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 || l1 )
  {
    printf ( "\n" );
    printf ( "R8_HYPER_2F1 - Fatal error!\n" );
    printf ( "  The hypergeometric series is divergent.\n" );
    hf = 0.0;
    return hf;
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = r8_gamma ( c );
    gcab = r8_gamma ( c - a - b );
    gca = r8_gamma ( c - a );
    gcb = r8_gamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && r8_abs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( pi ) * pow ( 2.0, - a );
    g1 = r8_gamma ( c );
    g2 = r8_gamma ( 1.0 + a / 2.0 - b );
    g3 = r8_gamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( r8_abs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( r8_abs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( r8_abs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( r8_abs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( r8_abs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = ( int ) ( c - a - b );
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gam = r8_gamma ( a + m );
      gbm = r8_gamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 ) 
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a ) 
              / ( ( j + k ) * ( a + j + k - 1.0 ) ) 
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) 
            / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gca = r8_gamma ( c - a );
      gcb = r8_gamma ( c - b );
      gcab = r8_gamma ( c - a - b );
      gabc = r8_gamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( r8_abs ( hf - hw ) < r8_abs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( r8_abs ( hf - hw ) <= r8_abs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    printf ( "\n" );
    printf ( "R8_HYPER_2F1 - Warning!\n" );
    printf ( "  A large number of iterations were needed.\n" );
    printf ( "  The accuracy of the results should be checked.\n" );
  }

  return hf;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

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

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

double r8_mop ( int i )

/******************************************************************************/
/*
  Purpose:

    R8_MOP returns the I-th power of -1 as an R8 value.

  Discussion:

    An R8 is an double value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int I, the power of -1.

    Output, double R8_MOP, the I-th power of -1.
*/
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
/******************************************************************************/

double r8_psi ( double xx )

/******************************************************************************/
/*
  Purpose:

    R8_PSI evaluates the function Psi(X).

  Discussion:

    This routine evaluates the logarithmic derivative of the
    Gamma function,

      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
             = d/dX LN ( GAMMA(X) )

    for real X, where either

      - XMAX1 < X < - XMIN, and X is not a negative integer,

    or

      XMIN < X.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.

  Reference:

    William Cody, Anthony Strecok, Henry Thacher,
    Chebyshev Approximations for the Psi Function,
    Mathematics of Computation,
    Volume 27, Number 121, January 1973, pages 123-127.

  Parameters:

    Input, double XX, the argument of the function.

    Output, double R8_PSI, the value of the function.
*/
{
  double aug;
  double den;
  double four = 4.0;
  double fourth = 0.25;
  double half = 0.5;
  int i;
  int n;
  int nq;
  double one = 1.0;
  double p1[9] = { 
   4.5104681245762934160E-03, 
   5.4932855833000385356, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = { 
  -2.7103228277757834192, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821, 
  -1.4479614616899842986, 
  -7.3689600332394549911E-02, 
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = { 
   9.6141654774222358525E+01, 
   2.6287715790581193330E+03, 
   2.9862497022250277920E+04, 
   1.6206566091533671639E+05, 
   4.3487880712768329037E+05, 
   5.4256384537269993733E+05, 
   2.4242185002017985252E+05, 
   6.4155223783576225996E-08 };
  double q2[6] = { 
   4.4992760373789365846E+01, 
   2.0240955312679931159E+02, 
   2.4736979003315290057E+02, 
   1.0742543875702278326E+02, 
   1.7463965060678569906E+01, 
   8.8427520398873480342E-01 };
  double sgn;
  double three = 3.0;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;
  double zero = 0.0;

  x = xx;
  w = r8_abs ( x );
  aug = zero;
/*
  Check for valid arguments, then branch to appropriate algorithm.
*/
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( zero < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < half )
  {
/*
  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
*/
    if ( w <= xsmall )
    {
      aug = - one / x;
    }
/*
  Argument reduction for cotangent.
*/
    else
    {
      if ( x < zero )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = ( int ) ( w * four );
      w = four * ( w - ( double ) ( nq ) * fourth );
/*
  W is now related to the fractional part of 4.0 * X.
  Adjust argument to correspond to values in the first
  quadrant and determine the sign.
*/
      n = nq / 2;

      if ( n + n != nq )
      {
        w = one - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
/*
  Determine the final value for  -pi * cotan(pi*x).
*/
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
/*
  Check for singularity.
*/
        if ( z == zero )
        {
          if ( zero < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
          }
          return value;
        }
        aug = sgn * ( four / tan ( z ) );
      }
      else
      {
        aug = sgn * ( four * tan ( z ) );
      }
    }
    x = one - x;
  }
/*
  0.5 <= X <= 3.0.

  if ( x <= three )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
/*
  3.0 < X.
*/
  if ( x < xlarge )
  {
    w = one / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - half / x + aug;
  }

  value = aug + log ( x );

  return value;
}
/******************************************************************************/

int r8col_compare ( int m, int n, double a[], int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8COL_COMPARE compares two columns in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      M = 3, N = 4, I = 2, J = 4

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

    Output:

      R8COL_COMPARE = -1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Input, int I, J, the columns to be compared.
    I and J must be between 1 and N.

    Output, int R8COL_COMPARE, the results of the comparison:
    -1, column I < column J,
     0, column I = column J,
    +1, column J < column I.
*/
{
  int k;
  int value;
/*
  Check.
*/
  if ( i < 1 || n < i )
  {
    printf ( "\n" );
    printf ( "R8COL_COMPARE - Fatal error!\n" );
    printf ( "  Column index I is out of bounds.\n" );
    printf ( "  I = %d\n", i );
    exit ( 1 );
  }

  if ( j < 1 || n < j )
  {
    printf ( "\n" );
    printf ( "R8COL_COMPARE - Fatal error!\n" );
    printf ( "  Column index J is out of bounds.\n" );
    printf ( "  J = %d\n", j );
    exit ( 1 );
  }

  value = 0;

  if ( i == j )
  {
    return value;
  }

  k = 0;

  while ( k < m )
  {
    if ( a[k+(i-1)*m] < a[k+(j-1)*m] )
    {
      value = -1;
      return value;
    }
    else if ( a[k+(j-1)*m] < a[k+(i-1)*m] )
    {
      value = +1;
      return value;
    }
    k = k + 1;
  }

  return value;
}
/******************************************************************************/

void r8col_sort_heap_a ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_A ascending heapsorts an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    In lexicographic order, the statement "X < Y", applied to two real
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that

      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).

    In other words, the first time they differ, X is smaller.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N].
    On input, the array of N columns of M-vectors.
    On output, the columns of A have been sorted in lexicographic order.
*/
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( n <= 1 )
  {
    return;
  }
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      r8col_swap ( m, n, a, i, j );
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      isgn = r8col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

int *r8col_sort_heap_index_a ( int m, int n, int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      A(*,INDX(*)) is sorted,

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in each column of A.

    Input, int N, the number of columns in A.

    Input, int BASE, the desired indexing for the sort index:
    0 for 0-based indexing, 
    1 for 1-based indexing.

    Input, double A[M*N], the array.

    Output, int R8COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th column of the sorted array is A(*,INDX(I)).
*/
{
  double *column;
  int i;
  int *indx;
  int indxt;
  int ir;
  int isgn;
  int j;
  int k;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  column = ( double * ) malloc ( m * sizeof ( double ) );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
    }
    else
    {
      indxt = indx[ir-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        isgn = r8vec_compare ( m, a+indx[j-1]*m, a+indx[j]*m );

        if ( isgn < 0 )
        {
          j = j + 1;
        }
      }

      isgn = r8vec_compare ( m, column, a+indx[j-1]*m );

      if ( isgn < 0 )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
  free ( column );
/*
  Take care of the base.
*/
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
/******************************************************************************/

int r8col_sorted_unique_count ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

    If the tolerance is large enough, then the concept of uniqueness
    can become ambiguous.  If we have a tolerance of 1.5, then in the
    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
    one unique entry?  That would be because 1 may be regarded as unique,
    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
    be unique and so on.

    This seems wrongheaded.  So I prefer the idea that an item is not
    unique under a tolerance only if it is close to something that IS unique.
    Thus, the unique items are guaranteed to cover the space if we include
    a disk of radius TOL around each one.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], a sorted array, containing
    N columns of data.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_SORTED_UNIQUE_COUNT, the number of unique columns.
*/
{
  double diff;
  int i;
  int j1;
  int j2;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    diff = 0.0;
    for ( i = 0; i < m; i++ )
    {
      diff = r8_max ( diff,  r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
    }
    if ( tol < diff )
    {
      unique_num = unique_num + 1;
      j1 = j2;
    }
  }

  return unique_num;
}
/******************************************************************************/

void r8col_swap ( int m, int n, double a[], int j1, int j2 )

/******************************************************************************/
/*
  Purpose:

    R8COL_SWAP swaps columns J1 and J2 of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      M = 3, N = 4, J1 = 2, J2 = 4

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

    Output:

      A = (
        1.  4.  3.  2.
        5.  8.  7.  6.
        9. 12. 11. 10. )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the M by N array.

    Input, int J1, J2, the columns to be swapped.
    These columns are 1-based.
*/
{
  int i;
  double temp;

  if ( j1 < 1 || n < j1 || j2 < 1 || n < j2 )
  {
    printf ( "\n" );
    printf ( "R8COL_SWAP - Fatal error!\n" );
    printf ( "  J1 or J2 is out of bounds.\n" );
    printf ( "  J1 =   %d\n", j1 );
    printf ( "  J2 =   %d\n", j2 );
    printf ( "  NCOL = %d\n", n );
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    return;
  }

  for ( i = 0; i < m; i++ )
  {
    temp          = a[i+(j1-1)*m];
    a[i+(j1-1)*m] = a[i+(j2-1)*m];
    a[i+(j2-1)*m] = temp;
  }

  return;
}
/******************************************************************************/

void r8col_tol_undex ( int m, int n, double a[], int unique_num, double tol, 
  int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNDEX indexes tolerably unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    The first step of the algorithm requires the indexed sorting
    of A, which creates arrays INDX and XDNI.  (If all the entries
    of A are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of A in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      A(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     A  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is A(3).
    XDNI(2) = 5 means that A(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
    XDNU(8) = 2 means that A(8) is at unique sorted item(2).

    XU(XDNU(I))) = A(I).
    XU(I)        = A(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  int base = 0;
  double diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  int unique;
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = indx[i];
  xdnu[indx[i]] = k;
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        xdnu[indx[i]] = j;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
      xdnu[indx[i]] = k;
    }
  }
  free ( indx );

  return;
}
/******************************************************************************/

int r8col_tol_unique_count ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

    If the tolerance is large enough, then the concept of uniqueness
    can become ambiguous.  If we have a tolerance of 1.5, then in the
    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
    one unique entry?  That would be because 1 may be regarded as unique,
    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
    be unique and so on.

    This seems wrongheaded.  So I prefer the idea that an item is not
    unique under a tolerance only if it is close to something that IS unique.
    Thus, the unique items are guaranteed to cover the space if we include
    a disk of radius TOL around each one.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array of N columns of data.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_TOL_UNIQUE_COUNT, the number of unique columns.
*/
{
  int base = 0;
  double diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  int *undx;
  int unique;

  undx = ( int * ) malloc ( n * sizeof ( int ) );
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = indx[i];
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
    }
  }
  free ( indx );
  free ( undx );

  k = k + 1;

  return k;
}
/******************************************************************************/

void r8col_undex ( int x_dim, int x_num, double x_val[], int x_unique_num, 
  double tol, int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNDEX returns unique sorted indexes for an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of X, in sorted order,
    and a vector XDNU, which identifies, for each entry of X, the index of
    the unique sorted element of X.

    This is all done with index vectors, so that the elements of
    X are never moved.

    The first step of the algorithm requires the indexed sorting
    of X, which creates arrays INDX and XDNI.  (If all the entries
    of X are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of X in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector X could be
    replaced by a compressed vector XU, containing the unique entries
    of X in sorted order, using the formula

      XU(*) = X(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector X, or
    any element of it, by index, as follows:

      X(I) = XU(XDNU(I)).

    We could then replace X by the combination of XU and XDNU.

    Later, when we need the I-th entry of X, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector X, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     X  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).

    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int X_DIM, the dimension of the data values.
    (the number of rows in the R8COL).

    Input, int X_NUM, the number of data values,
    (the number of columns in the R8COL).

    Input, double X_VAL[X_DIM*X_NUM], the data values.

    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
  int base = 0;
  double diff;
  int i;
  int *indx;
  int j;
  int k;
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( x_dim, x_num, base, x_val );
/*
  Walk through the implicitly sorted array X.
*/
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < x_num; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < x_dim; k++ )
    {
      diff = r8_max ( diff,
        r8_abs ( x_val[k+indx[i]*x_dim] - x_val[k+undx[j]*x_dim] ) );
    }
    if ( tol < diff )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  free ( indx );

  return;
}
/******************************************************************************/

int *r8col_unique_index ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIQUE_INDEX indexes the unique occurrence of values in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A, 
    gathered in order, then 

      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".

    Input, double A[M*N], the array.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_UNIQUE_INDEX[N], the unique index.
*/
{
  double diff;
  int i;
  int j1;
  int j2;
  int *unique_index;
  int unique_num;

  unique_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_index[j1] = -1;
  }
  unique_num = 0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( unique_index[j1] == -1 )
    {
      unique_index[j1] = unique_num;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r8_max ( diff, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff <= tol )
        {
          unique_index[j2] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14f", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the table data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16e", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

double *r8vec_chebyshev_new ( int n, double a_first, double a_last )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CHEBYSHEV_NEW creates a vector of Chebyshev spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A_FIRST, A_LAST, the first and last entries.

    Output, double R8VEC_CHEBYSHEV_NEW[N], a vector of Chebyshev spaced data.
*/
{
  double *a;
  double c;
  int i;
  double pi = 3.141592653589793;
  double theta;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      a[i] = ( ( 1.0 - c ) * a_first  
             + ( 1.0 + c ) * a_last ) 
             /   2.0;

    }
  }

  return a;
}
/******************************************************************************/

int r8vec_compare ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COMPARE compares two R8VEC's.

  Discussion:

    The lexicographic ordering is used.

  Example:

    Input:

      A1 = ( 2.0, 6.0, 2.0 )
      A2 = ( 2.0, 8.0, 12.0 )

    Output:

      ISGN = -1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A[N], B[N], the vectors to be compared.

    Output, int R8VEC_COMPARE, the results of the comparison:
    -1, A is lexicographically less than B,
     0, A is equal to B,
    +1, A is lexicographically greater than B.
*/
{
  int isgn;
  int k;

  isgn = 0;

  for ( k = 0; k < n; k++ )
  {
    if ( a[k] < b[k] )
    {
      isgn = -1;
      return isgn;
    }
    else if ( b[k] < a[k] )
    {
      isgn = +1;
      return isgn;
    }
  }
  return isgn;
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

void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    To explain what is going on here, suppose we had to construct
    a multidimensional quadrature rule as the product of K rules
    for 1D quadrature.

    The product rule will be represented as a list of points and weights.

    The J-th item in the product rule will be associated with
      item J1 of 1D rule 1,
      item J2 of 1D rule 2, 
      ..., 
      item JK of 1D rule K.

    In particular, 
      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    and
      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)

    So we can construct the quadrature rule if we can properly
    distribute the information in the 1D quadrature rules.

    This routine carries out that task for the weights W.

    Another way to do this would be to compute, one by one, the
    set of all possible indices (J1,J2,...,JK), and then index
    the appropriate information.  An advantage of the method shown
    here is that you can process the K-th set of information and
    then discard it.

  Example:

    Rule 1: 
      Order = 4
      W(1:4) = ( 2, 3, 5, 7 )

    Rule 2:
      Order = 3
      W(1:3) = ( 11, 13, 17 )

    Rule 3:
      Order = 2
      W(1:2) = ( 19, 23 )

    Product Rule:
      Order = 24
      W(1:24) =
        ( 2 * 11 * 19 )
        ( 3 * 11 * 19 )
        ( 4 * 11 * 19 )
        ( 7 * 11 * 19 )
        ( 2 * 13 * 19 )
        ( 3 * 13 * 19 )
        ( 5 * 13 * 19 )
        ( 7 * 13 * 19 )
        ( 2 * 17 * 19 )
        ( 3 * 17 * 19 )
        ( 5 * 17 * 19 )
        ( 7 * 17 * 19 )
        ( 2 * 11 * 23 )
        ( 3 * 11 * 23 )
        ( 5 * 11 * 23 )
        ( 7 * 11 * 23 )
        ( 2 * 13 * 23 )
        ( 3 * 13 * 23 )
        ( 5 * 13 * 23 )
        ( 7 * 13 * 23 )
        ( 2 * 17 * 23 )
        ( 3 * 17 * 23 )
        ( 5 * 17 * 23 )
        ( 7 * 17 * 23 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.

    Input, int FACTOR_ORDER, the order of the factor.

    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
    factor FACTOR_INDEX.

    Input, int FACTOR_NUM, the number of factors.

    Input, int POINT_NUM, the number of elements in the direct product.

    Input/output, double W[POINT_NUM], the elements of the
    direct product, which are built up gradually.  

  Local Parameters:

    Local, integer START, the first location of a block of values to set.

    Local, integer CONTIG, the number of consecutive values to set.

    Local, integer SKIP, the distance from the current value of START
    to the next location of a block of values to set.

    Local, integer REP, the number of blocks of values to set.
*/
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}
/******************************************************************************/

void r8vec_index_sorted_range ( int n, double r[], int indx[], double r_lo, 
  double r_hi, int *i_lo, int *i_hi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items in the vector.

    Input, double R[N], the index sorted vector.

    Input, int INDX[N], the vector used to sort R.
    The vector R[INDX[*]] is sorted.

    Input, double R_LO, R_HI, the limits of the range.

    Output, int *I_LO, *I_HI, the range of indices
    so that I_LO <= I <= I_HI => R_LO <= R[INDX[I]] <= R_HI.  If no
    values in R lie in the range, then I_HI < I_LO will be returned.
*/
{
  int i1;
  int i2;
  int j1;
  int j2;
/*
  Cases we can handle immediately.
*/
  if ( r[indx[n-1]] < r_lo )
  {
    *i_lo = n;
    *i_hi = n - 1;
    return;
  }

  if ( r_hi < r[indx[0]] )
  {
    *i_lo = 0;
    *i_hi = -1;
    return;
  }
/*
  Are there are least two intervals?
*/
  if ( n == 1 )
  {
    if ( r_lo <= r[indx[0]] && r[indx[0]] <= r_hi )
    {
      *i_lo = 1;
      *i_hi = 1;
    }
    else
    {
      *i_lo = 0;
      *i_hi = -1;
    }
    return;
  }
/*
  Bracket R_LO.
*/
  if ( r_lo <= r[indx[0]] )
  {
    i_lo = 0;
  }
  else
  {
/*
  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
  Does R_LO lie here, or below or above?
*/
    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;
 
    for ( ; ; )
    {
      if ( r_lo < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_lo = i1;
        break;
      }
    }
  }
/*
  Bracket R_HI.
*/
  if ( r[indx[n-1]] <= r_hi )
  {
    *i_hi = n - 1;
  }
  else
  {
    j1 = *i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;
 
    for ( ; ; )
    {
      if ( r_hi < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_hi = i2;
        break;
      }
    }
  }
/*
  We expect to have computed the largest I_LO and smallest I_HI such that
    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
  but what we want is actually
    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
  which we can usually get simply by incrementing I_LO and decrementing I_HI.
*/
  if ( r[indx[*i_lo]] < r_lo )
  {
    *i_lo = *i_lo + 1;
    if ( n - 1 < *i_lo )
    {
      *i_hi = *i_lo - 1;
    }
  }

  if ( r_hi < r[indx[*i_hi]] )
  {
    *i_hi = *i_hi - 1;
    if ( i_hi < 0 )
    {
      *i_lo = *i_hi + 1;
    }
  }

  return;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min_pos ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN_POS returns the minimum positive value of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double A[N], the array.

    Output, double R8VEC_MIN_POS, the smallest positive entry, 
    or R8_HUGE if no entry is positive.
*/
{
  int i;
  double r8_huge = 1.0E+30;
  double value;

  value = r8_huge;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      if ( a[i] < value )
      {
        value = a[i];
      }
    }
  }
  return value;
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
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

int *r8vec_sort_heap_index_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int R8VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
/******************************************************************************/

double r8vec_sum ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SUM returns the sum of an R8VEC.

  Discussion:

    An R8VEC is a double precision vector.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_SUM, the sum of the vector.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
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
  int i4_huge = 2147483647;
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

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

/******************************************************************************/
/*
  Purpose:

    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.

  Discussion:

    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.

    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 February 2004

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the length of the input list.

    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.

    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.

    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
/*
  INDX = 0: This is the first call.
*/
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
/*
  INDX < 0: The user is returning the results of a comparison.
*/
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
/*
  0 < INDX: the user was asked to make an interchange.
*/
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
/******************************************************************************/

void timestamp ( )

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
/******************************************************************************/

void vec_colex_next3 ( int dim_num, int base[], int a[], int *more )

/******************************************************************************/
/*
  Purpose:

    VEC_COLEX_NEXT3 generates vectors in colex order.

  Discussion:

    The vectors are produced in colexical order, starting with

    (1,        1,        ...,1),
    (2,        1,        ...,1),
     ...
    (BASE(1),  1,        ...,1)

    (1,        2,        ...,1)
    (2,        2,        ...,1)
    ...
    (BASE(1),  2,        ...,1)

    (1,        3,        ...,1)
    (2,        3,        ...,1)
    ...
    (BASE(1),  BASE(2), ...,BASE(DIM_NUM)).

  Example:

    DIM_NUM = 2,
    BASE = { 3, 3 }

    1   1
    2   1
    3   1
    1   2
    2   2
    3   2
    1   3
    2   3
    3   3

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
    In dimension I, entries will range from 1 to BASE[I].

    Output, int A[DIM_NUM], the next vector.

    Input/output, int *MORE.  Set this variable 0 before
    the first call.  On return, MORE is TRUE if another vector has
    been computed.  If MORE is returned FALSE, ignore the output 
    vector and stop calling the routine.
*/
{
  int i;

  if ( !( *more ) )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 1;
    }
    *more = 1;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] <= base[i] )
      {
        return;
      }
      a[i] = 1;
    }
    *more = 0;
  }
  return;
}

