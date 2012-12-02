# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

# include "c4lib.h"

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

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, float complex C1, the argument.

    Output, float complex C4_CONJ, the function value.
*/
{
  float c1_norm;
  float complex c2;
  float c2_imag;
  float c2_real;

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
  float r;

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

float complex c4_i ( void )

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

float complex c4_one ( void )

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

float complex c4_zero ( void )

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

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a unit pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r4_uniform_01 = seed / ( 2**31 - 1 )

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

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int i4_huge = 2147483647;
  int k;
  float r;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "R4_UNIFORM_01 - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( float ) ( *seed ) * 4.656612875E-10;

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
    C version by John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.

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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
