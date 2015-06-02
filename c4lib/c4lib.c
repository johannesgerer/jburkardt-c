# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

# include "c4lib.h"
# include "r4lib.h"

/******************************************************************************/

float c4_abs ( float complex c )

/******************************************************************************/
/*
  Purpose:

    C4_ABS returns the absolute value of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the argument.

    Output, float C4_ABS, the function value.
*/
{
  float value;

  value = sqrt ( pow ( creal ( c ), 2 ) 
               + pow ( cimag ( c ), 2 ) );

  return value;
}
/******************************************************************************/

float complex c4_acos ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ACOS evaluates the inverse cosine of a C4.

  Discussion:

    Here we use the relationship:

      C4_ACOS ( Z ) = pi/2 - C4_ASIN ( Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ACOS, the function value.
*/
{
  float complex c2;
  float c2_imag;
  float c2_real;
  float r4_pi_half = 1.5707963;

  c2 = c4_asin ( c1 );

  c2_real = r4_pi_half - creal ( c2 );
  c2_imag =            - cimag ( c2 );

  c2 = c2_real + c2_imag * I;

  return c2;
}
/******************************************************************************/

float complex c4_acosh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ACOSH evaluates the inverse hyperbolic cosine of a C4.

  Discussion:

    Here we use the relationship:

      C4_ACOSH ( Z ) = i * C4_ACOS ( Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ACOSH, the function value.
*/
{
  float complex c2;

  c2 = c4_i ( ) * c4_acos ( c1 );
  
  return c2;
}
/******************************************************************************/

float complex c4_add ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_ADD adds two C4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, C2, the arguments.

    Output, float complex C4_ADD, the sum of C1 and C2.
*/
{
  float complex c3;
  float c3_imag;
  float c3_real;

  c3_real = creal ( c1 ) + creal ( c2 );
  c3_imag = cimag ( c1 ) + cimag ( c2 );

  c3 = c3_real + c3_imag * I;

  return c3;
}
/******************************************************************************/

float c4_arg ( float complex c )

/******************************************************************************/
/*
  Purpose:

    C4_ARG returns the argument of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the complex number.

    Output, float C4_ARG, the function value.
*/
{
  float value;

  if ( cimag ( c ) == 0.0 && creal ( c ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = atan2 ( cimag ( c ), creal ( c ) );
  }
  return value;
}
/******************************************************************************/

float complex c4_asin ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ASIN evaluates the inverse sine of a C4.

  Discussion:

    Here we use the relationship:

      C4_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ASIN, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex ce;

  c2 = c4_i ( );
  c3 = c4_sqrt ( 1.0 - c1 * c1 );
  c4 = c4_log ( c3 + c2 * c1 );
  ce = - c2 * c4;

  return ce;
}
/******************************************************************************/

float complex c4_asinh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ASINH evaluates the inverse hyperbolic sine of a C4.

  Discussion:

    Here we use the relationship:

      C4_ASINH ( Z ) = - i * C4_ASIN ( i * Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ASINH, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;

  c2 = c4_i ( );
  c3 = c2 * c1;
  c4 = c4_asin ( c3 );
  c5 = c2 * c4;
  c6 = - c5;

  return c6;
}
/******************************************************************************/

float complex c4_atan ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ATAN evaluates the inverse tangent of a C4.

  Discussion:

    Here we use the relationship:

      C4_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ATAN, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;
  float complex c7;
  float complex c8;
  float complex c9;
  float complex cx;

  c2 = c4_i ( );
  c3 = c4_one ( );
  c4 = c4_mul ( c2, c1 );
  c5 = c4_sub ( c3, c4 );
  c6 = c4_add ( c3, c4 );
  c7 = c4_div ( c5, c6 );

  c8 = c4_log ( c7 );
  c9 =  c2 * c8;
  cx = c9 / 2.0;

  return cx;
}
/******************************************************************************/

float complex c4_atanh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_ATANH evaluates the inverse hyperbolic tangent of a C4.

  Discussion:

    Here we use the relationship:

      C4_ATANH ( Z ) = - i * C4_ATAN ( i * Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ATANH, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;

  c2 = c4_i ( );

  c3 = c4_mul ( c2, c1 );
  c4 = c4_atan ( c3 );
  c5 = c4_mul ( c2, c4 );
  c6 = c4_neg ( c5 );

  return c6;
}
/******************************************************************************/

float complex c4_conj ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_CONJ conjugates a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_CONJ, the function value.
*/
{
  float complex c2;

  c2 = conj ( c1 );

  return c2;
}
/******************************************************************************/

void c4_copy ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_COPY copies a C4.

  Discussion:

    The order of the arguments may seem unnatural, but it is arranged so
    that the call

      c4_copy ( c1, c2 )

    mimics the assignment

      c1 = c2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, float complex C1, the copy of C2.

    Input, float complex C2, the value to be copied.
*/
{
  c1 = c2;

  return;
}
/******************************************************************************/

float complex c4_cos ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_COS evaluates the cosine of a C4.

  Discussion:

    We use the relationship:

      C4_COS ( C ) = ( C4_EXP ( i * C ) + C4_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_COS, the function value.
*/
{
  float complex c2;

  c2 = ( cexp ( c1 * I ) + cexp ( - c1 * I ) ) / 2.0;

  return c2;
}
/******************************************************************************/

float complex c4_cosh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_COSH evaluates the hyperbolic cosine of a C4.

  Discussion:

    A C4 is a C4_COMPLEX value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_COSH, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_add ( c2, c4 );
  c6 = c5 / 2.0;

  return c6;
}
/******************************************************************************/

float complex c4_cube_root ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_CUBE_ROOT computes the principal cube root of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_CUBE_ROOT, the function value.
*/
{
  float complex c2;
  float r;
  float t;

  c4_to_polar ( c1, &r, &t );

  r = pow ( r, 1.0 / 3.0 );
  t = t / 3.0;

  c2 = polar_to_c4 ( r, t );

  return c2;
}
/******************************************************************************/

float complex c4_div ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_DIV divides two C4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, C2, the arguments.

    Output, float complex C4_DIV, the function value.
*/
{
  float c2_norm;
  float complex c3;
  float c3_imag;
  float c3_real;

  c2_norm = c4_abs ( c2 );

  c3_real = ( creal ( c1 ) * creal ( c2 ) 
            + cimag ( c1 ) * cimag ( c2 ) ) / c2_norm / c2_norm;

  c3_imag = ( cimag ( c1 ) * creal ( c2 ) 
            - creal ( c1 ) * cimag ( c2 ) ) / c2_norm / c2_norm;

  c3 = c3_real + c3_imag * I;

  return c3;
}
/******************************************************************************/

float complex c4_div_r4 ( float complex c1, float r )

/******************************************************************************/
/*
  Purpose:

    C4_DIV_R4 divides a C4 by an R4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the value to be divided.

    Input, float R, the divisor.

    Output, float complex C4_DIV_R4, the function value.
*/
{
  float complex c2;

  c2 = c1 / r;

  return c2;
}
/******************************************************************************/

float complex c4_exp ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_EXP exponentiates a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_EXP, the function value.
*/
{
  float complex c2;

  c2 = cexp ( c1 );

  return c2;
}
/******************************************************************************/

float complex c4_i ( )

/******************************************************************************/
/*
  Purpose:

    C4_I returns the value of I as a C4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, float complex C4_I, the value of I.
*/
{
  float complex c1;

  c1 = I;

  return c1;
}
/******************************************************************************/

float c4_imag ( float complex c )

/******************************************************************************/
/*
  Purpose:

    C4_IMAG returns the imaginary part of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the argument.

    Output, float C4_IMAG, the function value.
*/
{
  float value;

  value = cimag ( c );

  return value;
}
/******************************************************************************/

float complex c4_inv ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_INV inverts a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_INV, the function value;
*/
{
  float complex c2;

  c2 = 1.0 / c1;

  return c2;
}
/******************************************************************************/

int c4_le_l1 ( float complex x, float complex y )

/******************************************************************************/
/*
  Purpose:

    C4_LE_L1 := X <= Y for C4 values, and the L1 norm.

  Discussion:

    A C4 is a float complex value.

    The L1 norm can be defined here as:

      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, Y, the values to be compared.

    Output, int C4_LE_L1, is 1 if X <= Y.
*/
{
  int value;

  if ( r4_abs ( creal ( x ) ) + r4_abs ( cimag ( x ) ) <= 
       r4_abs ( creal ( y ) ) + r4_abs ( cimag ( y ) ) )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int c4_le_l2 ( float complex x, float complex y )

/******************************************************************************/
/*
  Purpose:

    C4_LE_L2 := X <= Y for C4 values, and the L2 norm.

  Discussion:

    A C8 is a float complex value.

    The L2 norm can be defined here as:

      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, Y, the values to be compared.

    Output, int C4_LE_L2, is 1 if X <= Y.
*/
{
  int value;

  if ( pow ( creal ( x ), 2 ) + pow ( cimag ( x ), 2 ) <= 
       pow ( creal ( y ), 2 ) + pow ( cimag ( y ), 2 ) )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int c4_le_li ( float complex x, float complex y )

/******************************************************************************/
/*
  Purpose:

    C4_LE_LI := X <= Y for C4 values, and the L-oo norm.

  Discussion:

    A C4 is a float complex value.

    The L-oo norm can be defined here as:

      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, Y, the values to be compared.

    Output, int C4_LE_LI, is 1 if X <= Y.
*/
{
  int value;

  if ( r4_max ( r4_abs ( creal ( x ) ), r4_abs ( cimag ( x ) ) ) <= 
       r4_max ( r4_abs ( creal ( y ) ), r4_abs ( cimag ( y ) ) ) )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

float complex c4_log ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_LOG evaluates the logarithm of a C4.

  Discussion:

    Here we use the relationship:

      C4_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_ATAN, the function value.
*/
{
  float arg;
  float complex c2;
  float mag;

  arg = c4_arg ( c1 );
  mag = c4_mag ( c1 );

  c2 = log ( mag ) + arg * I;

  return c2;
}
/******************************************************************************/

float c4_mag ( float complex c )

/******************************************************************************/
/*
  Purpose:

    C4_MAG returns the magnitude of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the argument.

    Output, float C4_MAG, the function value.
*/
{
  float value;

  value = sqrt ( pow ( creal ( c ), 2 ) 
               + pow ( cimag ( c ), 2 ) );

  return value;
}
/******************************************************************************/

float complex c4_mul ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_MUL multiplies two C4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, C2, the arguments.

    Output, float complex C4_MUL, the function value.
*/
{
  float complex c3;

  c3 = c1 * c2;

  return c3;
}
/******************************************************************************/

float complex c4_neg ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_NEG negates a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_NEG, the function value.
*/
{
  float complex c2;

  c2 = - c1;

  return c2;
}
/******************************************************************************/

float complex c4_nint ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_NINT returns the nearest complex integer of a C4.

  Discussion:

    A C4 is a float complex value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the value to be NINT'ed.

    Output, float complex C4_NINT, the NINT'ed value.
*/
{
  double r;
  double r_min;
  double x;
  double x_min;
  double xc;
  double y;
  double y_min;
  double yc;
  double complex value;

  xc = creal ( c1 );
  yc = cimag ( c1 );
/*
  Lower left.
*/
  x = r4_floor ( creal ( c1 ) );
  y = r4_floor ( cimag ( c1 ) );
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  r_min = r;
  x_min = x;
  y_min = y;
/*
  Lower right.
*/
  x = r4_floor ( creal ( c1 ) ) + 1.0;
  y = r4_floor ( cimag ( c1 ) );
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }
/*
  Upper right.
*/
  x = r4_floor ( creal ( c1 ) ) + 1.0;
  y = r4_floor ( cimag ( c1 ) ) + 1.0;
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }
/*
  Upper left.
*/
  x = r4_floor ( creal ( c1 ) );
  y = r4_floor ( cimag ( c1 ) ) + 1.0;
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }

  value = x_min + y_min * I;

  return value;
}
/******************************************************************************/

float c4_norm_l1 ( float complex x )

/******************************************************************************/
/*
  Purpose:

    C4_NORM_L1 evaluates the L1 norm of a C4.

  Discussion:

    A C4 is a float complex value.

    Numbers of equal norm lie along diamonds centered at (0,0).

    The L1 norm can be defined here as:

      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, the value whose norm is desired.

    Output, float C4_NORM_L1, the norm of X.
*/
{
  double value;

  value = fabs ( creal ( x ) )
        + fabs ( cimag ( x ) );

  return value;
}
/******************************************************************************/

float c4_norm_l2 ( float complex x )

/******************************************************************************/
/*
  Purpose:

    C4_NORM_L2 evaluates the L2 norm of a C4.

  Discussion:

    A C4 is a float complex value.

    Numbers of equal norm lie on circles centered at (0,0).

    The L2 norm can be defined here as:

      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, the value whose norm is desired.

    Output, float C4_NORM_L2, the 2-norm of X.
*/
{
  double value;

  value = sqrt ( pow ( creal ( x ), 2 )
               + pow ( cimag ( x ), 2 ) );

  return value;
}
/******************************************************************************/

float c4_norm_li ( float complex x )

/******************************************************************************/
/*
  Purpose:

    C4_NORM_LI evaluates the L-infinity norm of a C4.

  Discussion:

    A C4 is a float complex value.

    Numbers of equal norm lie along squares whose centers are at (0,0).

    The L-infinity norm can be defined here as:

      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex X, the value whose norm is desired.

    Output, float C4_NORM_LI, the infinity norm of X.
*/
{
  double value;

  value = r4_max ( fabs ( creal ( x ) ), fabs ( cimag ( x ) ) );

  return value;
}
/******************************************************************************/

float complex c4_normal_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    C4_NORMAL_01 returns a unit pseudonormal C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, float complex C4_NORMAL_01, a unit pseudonormal value.
*/
{
  float r4_pi = 3.1415926;
  float v1;
  float v2;
  float complex value;
  float x_c;
  float x_r;

  v1 = r4_uniform_01 ( seed );
  v2 = r4_uniform_01 ( seed );

  x_r = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * r4_pi * v2 );
  x_c = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * r4_pi * v2 );

  value = cartesian_to_c4 ( x_r, x_c );

  return value;
}
/******************************************************************************/

float complex c4_one ( )

/******************************************************************************/
/*
  Purpose:

    C4_ONE returns the value of 1 as a C4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, float complex C4_ONE, the value of 1.
*/
{
  float complex c1;

  c1 = 1.0;

  return c1;
}
/******************************************************************************/

void c4_print ( float complex a, char *title )

/******************************************************************************/
/*
  Purpose:

    C4_PRINT prints a C4.

  Discussion:

    A C4 is a float complex value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, float complex A, the value to be printed.

    Input, char *TITLE, a title.
*/
{
  printf ( "%s  ( %g, %g )\n", title, creal ( a ), cimag ( a ) );

  return;
}
/******************************************************************************/

float c4_real ( float complex c )

/******************************************************************************/
/*
  Purpose:

    C4_REAL returns the real part of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the complex number.

    Output, float C4_REAL, the function value.
*/
{
  float value;

  value = creal ( c );

  return value;
}
/******************************************************************************/

float complex c4_sin ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_SIN evaluates the sine of a C4.

  Discussion:

    We use the relationship:

      C4_SIN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_SIN, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;
  float complex c7;
  float complex c8;
  float complex c9;
  float complex cx;
  float r;

  c2 = c4_i ( );

  c3 = c4_mul ( c2, c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_neg ( c3 );
  c6 = c4_exp ( c5 );

  c7 = c4_sub ( c4, c6 );

  r = 2.0;
  c8 = c4_div_r4 ( c7, r );
  c9 = c8 * c2;
  cx = - c9;

  return cx;
}
/******************************************************************************/

float complex c4_sinh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_SINH evaluates the hyperbolic sine of a C4.

  Discussion:

    We use the relationship:

      C4_SINH ( C ) = ( C4_EXP ( C ) - C4_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_SINH, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_sub ( c2, c4 );

  c6 = c5 / 2.0;

  return c6;
}
/******************************************************************************/

float complex c4_sqrt ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_SQRT computes a square root of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_SQRT, the function value.
*/
{
  float complex c2;
  float r;
  float t;

  c4_to_polar ( c1, &r, &t );

  r = sqrt ( r );
  t = t / 2.0;

  c2 = polar_to_c4 ( r, t );

  return c2;
}
/******************************************************************************/

float complex c4_sub ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_SUB subtracts two C4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, C2, the arguments.

    Output, float complex C4_SUB, the function value.
*/
{
  float complex c3;

  c3 = c1 - c2;

  return c3;
}
/******************************************************************************/

void c4_swap ( float complex c1, float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4_SWAP swaps two C4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, float complex C1, C2, the arguments.
*/
{
  float complex c3;

  c3 = c1;
  c1 = c2;
  c2 = c3;

  return;
}
/******************************************************************************/

float complex c4_tan ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_TAN evaluates the tangent of a C4.

  Discussion:

    We use the relationship:

      C4_TAN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) 
                         / ( C4_EXP ( I * C ) + C4_EXP ( - i * C ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_TAN, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;
  float complex c7;
  float complex c8;
  float complex c9;
  float complex cx;
  float complex ce;

  c2 = c4_i ( );
  c3 = c4_mul ( c2, c1 );
  c4 = c4_neg ( c3 );
  
  c5 = c4_exp ( c3 );
  c6 = c4_exp ( c4 );

  c7 = c4_sub ( c5, c6 );
  c8 = c4_add ( c5, c6 );

  c9 = c4_div ( c7, c8 );
  cx = c4_mul ( c2, c9 );
  ce = - cx;

  return ce;
}
/******************************************************************************/

float complex c4_tanh ( float complex c1 )

/******************************************************************************/
/*
  Purpose:

    C4_TANH evaluates the hyperbolic tangent of a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_TANH, the function value.
*/
{
  float complex c2;
  float complex c3;
  float complex c4;
  float complex c5;
  float complex c6;
  float complex c7;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_sub ( c2, c4 );
  c6 = c4_add ( c2, c4 );

  c7 = c4_div ( c5, c6 );

  return c7;
}
/******************************************************************************/

void c4_to_cartesian ( float complex c, float *x, float *y )

/******************************************************************************/
/*
  Purpose:

    C4_TO_CARTESIAN converts a C4 to Cartesian form.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the argument.

    Output, float *X, *Y, the Cartesian form.
*/
{
  *x = creal ( c );
  *y = cimag ( c );

  return;
}
/******************************************************************************/

void c4_to_polar ( float complex c, float *r, float *theta )

/******************************************************************************/
/*
  Purpose:

    C4_TO_POLAR converts a C4 to polar form.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C, the argument.

    Output, float *R, *THETA, the polar form.
*/
{
  *r = c4_abs ( c );
  *theta = c4_arg ( c );

  return;
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

    06 October 2010

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
  int i4_huge = 2147483647;
  int k;
  float r;
  float r4_pi = 3.1415926;
  float theta;
  float complex value;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "C4_UNIFORM_01 - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = sqrt ( ( ( float ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  theta = 2.0 * r4_pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

  value = r * cos ( theta ) + r * sin ( theta ) * I;
  
  return value;
}
/******************************************************************************/

float complex c4_zero ( )

/******************************************************************************/
/*
  Purpose:

    C4_ZERO returns the value of 0 as a C4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, float complex C4_ZERO, the value of 0.
*/
{
  float complex c1;

  c1 = 0.0;

  return c1;
}
/******************************************************************************/

void c4mat_add ( int m, int n, float complex alpha, float complex a[],
  float complex beta, float complex b[], float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_ADD combines two C4MAT's with complex scale factors.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float complex ALPHA, the first scale factor.

    Input, float complex A[M*N], the first matrix.

    Input, float complex BETA, the second scale factor.

    Input, float complex B[M*N], the second matrix.

    Output, float complex C[M*N], the result.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
/******************************************************************************/

void c8mat_add_r4 ( int m, int n, float alpha, float complex a[],
  float beta, float complex b[], float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_ADD_R4 combines two C4MAT's with real scale factors.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float ALPHA, the first scale factor.

    Input, float complex A[M*N], the first matrix.

    Input, float BETA, the second scale factor.

    Input, float complex B[M*N], the second matrix.

    Output, float complex C[M*N], the result.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
/******************************************************************************/

void c4mat_copy ( int m, int n, float complex a[], float complex b[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_COPY copies one C4MAT to another.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float complex A[M*N], the matrix to be copied.

    Output, float complex B[M*N], the copy.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+j*m] = a[i+j*m];
    }
  }

  return;
}
/******************************************************************************/

float complex *c4mat_copy_new ( int m, int n, float complex a1[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_COPY_NEW copies one C4MAT to a "new" C4MAT.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float complex A1[M*N], the matrix to be copied.

    Output, float complex C4MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  float complex *a2;
  int i;
  int j;

  a2 = ( float complex * ) malloc ( m * n * sizeof ( float complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return a2;
}
/******************************************************************************/

void c4mat_fss ( int n, float complex a[], int nb, float complex x[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_FSS factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, float complex A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, int NB, the number of right hand sides.

    Input/output, float complex X[N*NB], on input, the right hand sides of the
    linear systems.  On output, the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  float piv;
  float complex t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = c4_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < c4_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = c4_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "C4MAT_FSS - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}
/******************************************************************************/

float complex *c4mat_fss_new ( int n, float complex a[], int nb, 
  float complex b[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_FSS_NEW factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, float complex A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, int NB, the number of right hand sides.

    Input, float complex B[N*NB], the right hand sides of the linear systems.

    Output, float complex C4MAT_FSS_NEW[N*NB], the solutions of the 
    linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  float piv;
  float complex t;
  float complex *x;

  x = ( float complex * ) malloc ( n * nb * sizeof ( float complex ) );

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = c4_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol + 1; i <= n; i++ )
    {
      if ( piv < c4_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = c4_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "C4MAT_FSS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol + 1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
/******************************************************************************/

float complex *c4mat_identity_new ( int n )

/******************************************************************************/
/*
  Purpose:

    C4MAT_IDENTITY_NEW sets a C4MAT to the identity.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Output, float complex C4MAT_IDENTITY_NEW[N*N], the matrix.
*/
{
  float complex *a;
  int i;
  int j;

  a = ( float complex * ) malloc ( n * n * sizeof ( float complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
  return a;
}
/******************************************************************************/

float complex *c4mat_indicator_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    C4MAT_INDICATOR_NEW returns the C4MAT indicator matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, float complex C4MAT_INDICATOR_NEW[M*N], the matrix.
*/
{
  float complex *a;
  int i;
  int j;

  a = ( float complex * ) malloc ( m * n * sizeof ( float complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = i + j * I;
    }
  }
  return a;
}
/******************************************************************************/

void c4mat_minvm ( int n1, int n2, float complex a[], 
  float complex b[], float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_MINVM returns inverse(A) * B for C4MAT's.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, the order of the matrices.

    Input, float complex A[N1*N1], B[N1*N2], the matrices.

    Output, float complex C[N1*N2], the result, C = inverse(A) * B.
*/
{
  float complex *alu;

  alu = c4mat_copy_new ( n1, n1, a );

  c4mat_copy ( n1, n2, b, c );

  c4mat_fss ( n1, alu, n2, c );
 
  free ( alu );

  return;
}
/******************************************************************************/

float complex *c4mat_minvm_new ( int n1, int n2, float complex a[], 
  float complex b[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_MINVM_NEW returns inverse(A) * B for C4MAT's.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, the order of the matrices.

    Input, float complex A[N1*N1], B[N1*N2], the matrices.

    Output, float complex C4MAT_MINVM_NEW[N1*N2], the result, C = inverse(A) * B.
*/
{
  float complex *alu;
  float complex *c;

  alu = c4mat_copy_new ( n1, n1, a );

  c = c4mat_fss_new ( n1, alu, n2, b );
 
  free ( alu );

  return c;
}
/******************************************************************************/

void c4mat_mm ( int n1, int n2, int n3, float complex a[], float complex b[], 
  float complex c[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_MM multiplies two C4MAT's.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, float complex A[N1*N2], float complex B[N2*N3], 
    the matrices to multiply.

    Output, float complex C[N1*N3], the product matrix C = A * B.
*/
{
  float complex *c1;
  int i;
  int j;
  int k;

  c1 = ( float complex * ) malloc ( n1 * n3 * sizeof ( float complex ) );

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c1[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c1[i+j*n1] = c1[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  c4mat_copy ( n1, n3, c1, c );

  free ( c1 );

  return;
}
/******************************************************************************/

float complex *c4mat_mm_new ( int n1, int n2, int n3, float complex a[], 
  float complex b[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_MM_NEW multiplies two C4MAT's.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, float complex A[N1*N2], float complex B[N2*N3], 
    the matrices to multiply.

    Output, float complex C4MAT_MM_NEW[N1*N3], the product matrix C = A * B.
*/
{
  float complex *c;
  int i;
  int j;
  int k;

  c = ( float complex * ) malloc ( n1 * n3 * sizeof ( float complex ) );

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }
  return c;
}
/******************************************************************************/

void c4mat_nint ( int m, int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_NINT rounds the entries of a C4MAT.

  Discussion:

    A C4MAT is an array of float complex values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.

    Input/output, float complex A[M*N], the matrix to be NINT'ed.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = c4_nint ( a[i+j*m] );
    }
  }
  return;
}
/******************************************************************************/

float c4mat_norm_fro ( int m, int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_NORM_FRO returns the Frobenius norm of a C4MAT.

  Discussion:

    A C4MAT is an array of C4 values.

    The Frobenius norm is defined as

      C4MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) |A(I,J)| )

    The matrix Frobenius-norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      c4vec_norm_l2 ( A*x ) <= c4mat_norm_fro ( A ) * c4vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of the matrix.

    Input, float complex A[M*N], the matrix.

    Output, double C4MAT_NORM_FRO, the Frobenius norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( creal ( a[i+j*m] ), 2 )
                    + pow ( cimag ( a[i+j*m] ), 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

float c4mat_norm_l1 ( int m, int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_NORM_L1 returns the matrix L1 norm of a C4MAT.

  Discussion:

    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].

    The matrix L1 norm is defined as:

      C4MAT_NORM_L1 = max ( 1 <= J <= N )
        sum ( 1 <= I <= M ) abs ( A(I,J) ).

    The matrix L1 norm is derived from the vector L1 norm, and
    satisifies:

      c4vec_norm_l1 ( A * x ) <= c4mat_norm_l1 ( A ) * c4vec_norm_l1 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, float complex A[M*N], the matrix whose L1 norm is desired.

    Output, float C4MAT_NORM_L1, the L1 norm of A.
*/
{
  float col_sum;
  int i;
  int j;
  float value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    col_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      col_sum = col_sum + c4_abs ( a[i+j*m] );
    }
    value = r4_max ( value, col_sum );
  }

  return value;
}
/******************************************************************************/

float c4mat_norm_li ( int m, int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_NORM_LI returns the L-infinity norm of a C4MAT.

  Discussion:

    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
    in column-major order.

    The matrix L-oo norm is defined as:

      C4MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).

    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:

      c4vec_norm_li ( A * x ) <= c4mat_norm_li ( A ) * c4vec_norm_li ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of the matrix.

    Input, float complex A[M*N], the matrix.

    Output, float C4MAT_NORM_LI, the L-infinity norm of A.
*/
{
  int i;
  int j;
  float row_sum;
  float value;

  value = 0.0;
  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + c4_abs ( a[i+j*m] );
    }
    value = r4_max ( value, row_sum );
  }
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

void c4mat_scale ( int m, int n, float complex alpha, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_SCALE scales a C4MAT by a complex scale factor.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float complex ALPHA, the scale factor.

    Input/output, float complex A[M*N], the matrix to be scaled.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * alpha;
    }
  }
  return;
}
/******************************************************************************/

void c4mat_scale_r4 ( int m, int n, float alpha, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4MAT_SCALE_R4 scales a C4MAT by a real scale factor

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float ALPHA, the scale factor.

    Input/output, float complex A[M*N], the matrix to be scaled.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * alpha;
    }
  }
  return;
}
/******************************************************************************/

float complex *c4mat_uniform_01 ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.

  Discussion:

    A C4MAT is a matrix of float complex values.

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float complex C[M*N], the pseudorandom complex matrix.
*/
{
  float complex *c;
  int i;
  int j;
  float r;
  int k;
  float pi = 3.141592653589793;
  float theta;

  c = ( float complex * ) malloc ( m * n * sizeof ( float complex ) );

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

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta ) + sin ( theta ) * I );
    }
  }

  return c;
}
/******************************************************************************/

float complex *c4mat_zero_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    C4MAT_ZERO_NEW returns a new zeroed C4MAT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, float complex C4MAT_ZERO[M*N], the new zeroed matrix.
*/
{
  float complex *a;
  int i;
  int j;

  a = ( float complex * ) malloc ( m * n * sizeof ( float complex ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

void c4vec_copy ( int n, float complex a[], float complex b[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_COPY copies a C4VEC to another.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, float complex A[N], the vector to be copied.

    Output, float complex B[N], the copy.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[i];
  }
  return;
}
/******************************************************************************/

float complex *c4vec_copy_new ( int n, float complex a1[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_COPY_NEW copies a C4VEC to a "new" C4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, float complex A1[N], the vector to be copied.

    Output, float complex C4VEC_COPY_NEW[N], the copy of A1.
*/
{
  float complex *a2;
  int i;

  a2 = ( float complex * ) malloc ( n * sizeof ( float complex ) );

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
/******************************************************************************/

void c4vec_nint ( int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_NINT rounds the entries of a C4VEC.

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

    Input/output, float complex A[N], the vector to be nint'ed.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = c4_nint ( a[i] );
  }

  return;
}
/******************************************************************************/

float c4vec_norm_l2 ( int n, float complex a[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_NORM_L2 returns the L2 norm of a C4VEC.

  Discussion:

    The vector L2 norm is defined as:

      value = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, float complex A[N], the vector whose L2 norm is desired.

    Output, float C4VEC_NORM_L2, the L2 norm of A.
*/
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value 
          + creal ( a[i] ) * creal ( a[i] ) 
          + cimag ( a[i] ) * cimag ( a[i] );
  }
  value = sqrt ( value );

  return value;
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

void c4vec_print_part ( int n, float complex a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    C4VEC_PRINT_PART prints "part" of a C4VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, float complex A[N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: (%14f,  %14f)\n", i, creal ( a[i] ), cimag ( a[i] ) );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: (%14f,  %14f)\n", i, creal( a[i] ), cimag ( a[i] ) );
    }
    fprintf ( stdout, "  ........   ..............   ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: (%14f,  %14f)\n", i, creal( a[i] ), cimag ( a[i] ) );
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: (%14f,  %14f)\n", i, creal( a[i] ), cimag ( a[i] ) );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: (%14f,  %14f )  ...more entries...\n", 
      i, creal( a[i] ), cimag ( a[i] ) );
  }

  return;
}
/******************************************************************************/

void c4vec_print_some ( int n, float complex a[], int i_lo, int i_hi, 
  char *title )

/******************************************************************************/
/*
  Purpose:

    C4VEC_PRINT_SOME prints "some" of a C4VEC.

  Discussion:

    A C4VEC is a vector of C4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, float complex A[N], the vector to be printed.

    Input, integer I_LO, I_HI, the first and last indices to print.
    The routine expects 1 <= I_LO <= I_HI <= N.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    fprintf ( stdout, "  %8d: (%14f, %14f)\n", i, creal( a[i] ), cimag ( a[i] ) );
  }

  return;
}
/******************************************************************************/

void c4vec_sort_a_l2 ( int n, float complex x[] )

/******************************************************************************/
/*
  Purpose:

    C4VEC_SORT_A_L2 ascending sorts a float complex array by L2 norm.

  Discussion:

    The L2 norm of A+Bi is sqrt ( A * A + B * B ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input/output, float complex X[N].
    On input, an unsorted array.
    On output, X has been sorted.
*/
{
  int i;
  int indx;
  int isgn;
  int j;
  float normsq_i;
  float normsq_j;
  float complex temp;

  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;

  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );

    if ( 0 < indx )
    {
      temp = x[i-1];
      x[i-1] = x[j-1];
      x[j-1] = temp;
    }
    else if ( indx < 0 )
    {
      normsq_i = pow ( creal ( x[i-1] ), 2 )
               + pow ( cimag ( x[i-1] ), 2 );

      normsq_j = pow ( creal ( x[j-1] ), 2 )
               + pow ( cimag ( x[j-1] ), 2 );

      if ( normsq_i < normsq_j )
      {
        isgn = -1;
      }
      else
      {
        isgn = +1;
      }
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

float complex *c4vec_spiral ( int n, int m, float complex c1, 
  float complex c2 )

/******************************************************************************/
/*
  Purpose:

    C4VEC_SPIRAL returns N points on a spiral between C1 and C2.

  Discussion:

    A C4VEC is a vector of C4's.

    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
    
    Then the polar form of the I-th point C(I) is:

      R(I) = ( ( N - I     ) * R1 
             + (     I - 1 ) * R2 ) 
              / ( N    - 1 )

      T(I) = ( ( N - I     ) * T1 
             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
             / ( N     - 1 )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points on the spiral.

    Input, int M, the number of full circuits the 
    spiral makes.

    Input, float complex C1, C2, the first and last points 
    on the spiral.

    Output, float complex C4VEC_SPIRAL_NEW[N], the points.
*/
{
  float complex *c;
  int i;
  float r1;
  float r2;
  float ri;
  float r8_pi = 3.141592653589793;
  float t1;
  float t2;
  float ti;

  c = ( float complex * ) malloc ( n * sizeof ( float complex ) );

  r1 = c4_abs ( c1 );
  r2 = c4_abs ( c2 );

  t1 = c4_arg ( c1 );
  t2 = c4_arg ( c2 );

  if ( m == 0 )
  {
    if ( t2 < t1 )
    {
      t2 = t2 + 2.0 * r8_pi;
    }
  }
  else if ( 0 < m )
  {
    if ( t2 < t1 )
    {
      t2 = t2 + 2.0 * r8_pi;
    }
    t2 = t2 + ( float ) ( m ) * 2.0 * r8_pi;
  }
  else if ( m < 0 )
  {
    if ( t1 < t2 )
    {
      t2 = t2 - 2.0 * r8_pi;
    }
    t2 = t2 - ( float ) ( m ) * 2.0 * r8_pi;
  }

  for ( i = 0; i < n; i++ )
  {
    ri = ( ( float ) ( n - i - 1 ) * r1
         + ( float ) (     i     ) * r2 )
         / ( float ) ( n     - 1 );

    ti = ( ( float ) ( n - i - 1 ) * t1
         + ( float ) (     i     ) * t2 )
         / ( float ) ( n     - 1 );

    c[i] = polar_to_c4 ( ri, ti );
  }

  return c;
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

    09 March 2014

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

    Output, float complex C4VEC_UNIFORM_01_NEW[N], the pseudorandom vector.
*/
{
  float complex *c;
  int i;
  int i4_huge = 2147483647;
  float r;
  int k;
  float pi = 3.141592653589793;
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

    theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * cos ( theta ) + r * sin ( theta ) * I;
  }

  return c;
}
/******************************************************************************/

float complex *c4vec_unity ( int n )

/******************************************************************************/
/*
  Purpose:

    C4VEC_UNITY returns the N roots of unity in a C4VEC.

  Discussion:

    A C4VEC is a vector of C4's.

    X(1:N) = exp ( 2 * PI * (0:N-1) / N )

    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, float complex C4VEC_UNITY[N], the N roots of unity.
*/
{
  float complex *a;
  int i;
  float pi = 3.141592653589793;
  float theta;

  a = ( float complex * ) malloc ( n * sizeof ( float complex ) );

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( float ) ( 2 * i ) / ( float ) ( n );
    a[i] = cos ( theta ) + sin ( theta ) * I;
  }

  return a;
}
/******************************************************************************/

float complex cartesian_to_c4 ( float x, float y )

/******************************************************************************/
/*
  Purpose:

    CARTESIAN_TO_C4 converts a Cartesian form to a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float X, Y, the Cartesian form.

    Output, float complex CARTESIAN_TO_C4, the complex number.
*/
{
  float complex c;

  c = x + y * I;

  return c;
}
/******************************************************************************/

float complex polar_to_c4 ( float r, float theta )

/******************************************************************************/
/*
  Purpose:

    POLAR_TO_C4 converts a polar form to a C4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float R, THETA, the polar form.

    Output, float complex POLAR_TO_C4, the complex number.
*/
{
  float complex c;

  c = r * cos ( theta ) + r * sin ( theta ) * I;

  return c;
}
/******************************************************************************/

float complex r4_csqrt ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_CSQRT returns the complex square root of an R4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, float X, the number whose square root is desired.

    Output, float R4_CSQRT, the square root of X:
*/
{
  float argument;
  float magnitude;
  float r8_pi = 3.141592653589793;
  float complex value;

  if ( 0.0 < x )
  {
    magnitude = x;
    argument = 0.0;
  }
  else if ( 0.0 == x )
  {
    magnitude = 0.0;
    argument = 0.0;
  }
  else if ( x < 0.0 )
  {
    magnitude = -x;
    argument = r8_pi;
  }

  magnitude = sqrt ( magnitude );
  argument = argument / 2.0;

  value = magnitude * ( cos ( argument ) + sin ( argument ) * I );

  return value;
}
/******************************************************************************/

void r4poly2_root ( float a, float b, float c, float complex *r1,
  float complex *r2 )

/******************************************************************************/
/*
  Purpose:

    R4POLY2_ROOT returns the two roots of a quadratic polynomial.

  Discussion:

    The polynomial has the form:

      A * X^2 + B * X + C = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Parameters:

    Input, float A, B, C, the coefficients of the polynomial.
    A must not be zero.

    Output, float complex *R1, *R2, the roots of the polynomial, which
    might be real and distinct, real and equal, or complex conjugates.
*/
{
  float disc;
  float complex q;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4POLY2_ROOT - Fatal error!\n" );
    fprintf ( stderr, "  The coefficient A is zero.\n" );
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  q = -0.5 * ( b + r4_sign ( b ) * r4_csqrt ( disc ) );
  *r1 = q / a;
  *r2 = c / q;

  return;
}
/******************************************************************************/

void r4poly3_root ( float a, float b, float c, float d,
  float complex *r1, float complex *r2, float complex *r3 )

/******************************************************************************/
/*
  Purpose:

    R4POLY3_ROOT returns the three roots of a cubic polynomial.

  Discussion:

    The polynomial has the form

      A * X^3 + B * X^2 + C * X + D = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Parameters:

    Input, float A, B, C, D, the coefficients of the polynomial.
    A must not be zero.

    Output, float complex *R1, *R2, *R3, the roots of the polynomial, which
    will include at least one real root.
*/
{
  float pi = 3.141592653589793;
  float q;
  float r;
  float s1;
  float s2;
  float temp;
  float theta;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4POLY3_ROOT - Fatal error!\n" );
    fprintf ( stderr, "  A must not be zero.\n" );
    exit ( 1 );
  }

  q = ( pow ( b / a, 2 ) - 3.0 * ( c / a ) ) / 9.0;

  r = ( 2.0 * pow ( b / a, 3 ) - 9.0 * ( b / a ) * ( c / a )
      + 27.0 * ( d / a ) ) / 54.0;

  if ( r * r < q * q * q )
  {
    theta = acos ( r / sqrt ( pow ( q, 3 ) ) );
    *r1 = - 2.0 * sqrt ( q ) * cos (   theta              / 3.0 );
    *r2 = - 2.0 * sqrt ( q ) * cos ( ( theta + 2.0 * pi ) / 3.0 );
    *r3 = - 2.0 * sqrt ( q ) * cos ( ( theta + 4.0 * pi ) / 3.0 );
  }
  else if ( q * q * q <= r * r )
  {
    temp = - r + sqrt ( r * r - q * q * q );
    s1 = r4_sign ( temp ) * pow ( r4_abs ( temp ), 1.0 / 3.0 );

    temp = - r - sqrt ( r * r - q * q * q );
    s2 = r4_sign ( temp ) * pow ( r4_abs ( temp ), 1.0 / 3.0 );

    *r1 = s1 + s2;
    *r2 = -0.5 * ( s1 + s2 ) + I * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
    *r3 = -0.5 * ( s1 + s2 ) - I * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
  }

  *r1 = *r1 - b / ( 3.0 * a );
  *r2 = *r2 - b / ( 3.0 * a );
  *r3 = *r3 - b / ( 3.0 * a );

  return;
}
/******************************************************************************/

void r4poly4_root ( float a, float b, float c, float d, float e,
  float complex *r1, float complex *r2, float complex *r3,
  float complex *r4 )

/******************************************************************************/
/*
  Purpose:

    R4POLY4_ROOT returns the four roots of a quartic polynomial.

  Discussion:

    The polynomial has the form:

      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2014

  Parameters:

    Input, float A, B, C, D, the coefficients of the polynomial.
    A must not be zero.

    Output, float complex *R1, *R2, *R3, *R4, the roots of the polynomial.
*/
{
  float a3;
  float a4;
  float b3;
  float b4;
  float c3;
  float c4;
  float d3;
  float d4;
  float complex p;
  float complex q;
  float complex r;
  float complex zero;

  zero = 0.0;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4POLY4_ROOT - Fatal error!\n" );
    fprintf ( stderr, "  A must not be zero.\n" );
    exit ( 1 );
  }

  a4 = b / a;
  b4 = c / a;
  c4 = d / a;
  d4 = e / a;
/*
  Set the coefficients of the resolvent cubic equation.
*/
  a3 = 1.0;
  b3 = -b4;
  c3 = a4 * c4 - 4.0 * d4;
  d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
/*
  Find the roots of the resolvent cubic.
*/
  r4poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
/*
  Choose one root of the cubic, here R1.

  Set R = sqrt ( 0.25 * A4^2 - B4 + R1 )
*/
  r = c4_sqrt ( 0.25 * a4 * a4 - b4  + (*r1) );

  if ( creal ( r ) != 0.0 || cimag ( r ) != 0.0 )
  {
    p = c4_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        + 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );

    q = c4_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        - 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );
  }
  else
  {
    p = c4_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      + 2.0 * c4_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );

    q = c4_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      - 2.0 * c4_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );
  }
/*
  Set the roots.
*/
  *r1 = -0.25 * a4 + 0.5 * r + 0.5 * p;
  *r2 = -0.25 * a4 + 0.5 * r - 0.5 * p;
  *r3 = -0.25 * a4 - 0.5 * r + 0.5 * q;
  *r4 = -0.25 * a4 - 0.5 * r - 0.5 * q;

  return;
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
