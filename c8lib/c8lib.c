# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

# include "c8lib.h"

/******************************************************************************/

double c8_abs ( double complex c )

/******************************************************************************/
/*
  Purpose:

    C8_ABS returns the absolute value of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the argument.

    Output, double C8_ABS, the function value.
*/
{
  double value;

  value = sqrt ( pow ( creal ( c ), 2 ) 
               + pow ( cimag ( c ), 2 ) );

  return value;
}
/******************************************************************************/

double complex c8_acos ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ACOS evaluates the inverse cosine of a C8.

  Discussion:

    Here we use the relationship:

      C8_ACOS ( Z ) = pi/2 - C8_ASIN ( Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ACOS, the function value.
*/
{
  double complex c2;
  double c2_imag;
  double c2_real;
  double r8_pi_half = 1.57079632679489661923;

  c2 = c8_asin ( c1 );

  c2_real = r8_pi_half - creal ( c2 );
  c2_imag =            - cimag ( c2 );

  c2 = c2_real + c2_imag * I;

  return c2;
}
/******************************************************************************/

double complex c8_acosh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ACOSH evaluates the inverse hyperbolic cosine of a C8.

  Discussion:

    Here we use the relationship:

      C8_ACOSH ( Z ) = i * C8_ACOS ( Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ACOSH, the function value.
*/
{
  double complex c2;

  c2 = c8_i ( ) * c8_acos ( c1 );
  
  return c2;
}
/******************************************************************************/

double complex c8_add ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_ADD adds two C8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, C2, the arguments.

    Output, double complex C8_ADD, the sum of C1 and C2.
*/
{
  double complex c3;
  double c3_imag;
  double c3_real;

  c3_real = creal ( c1 ) + creal ( c2 );
  c3_imag = cimag ( c1 ) + cimag ( c2 );

  c3 = c3_real + c3_imag * I;

  return c3;
}
/******************************************************************************/

double c8_arg ( double complex c )

/******************************************************************************/
/*
  Purpose:

    C8_ARG returns the argument of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the complex number.

    Output, double C8_ARG, the function value.
*/
{
  double value;

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

double complex c8_asin ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ASIN evaluates the inverse sine of a C8.

  Discussion:

    Here we use the relationship:

      C8_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ASIN, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex ce;

  c2 = c8_i ( );
  c3 = c8_sqrt ( 1.0 - c1 * c1 );
  c4 = c8_log ( c3 + c2 * c1 );
  ce = - c2 * c4;

  return ce;
}
/******************************************************************************/

double complex c8_asinh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ASINH evaluates the inverse hyperbolic sine of a C8.

  Discussion:

    Here we use the relationship:

      C8_ASINH ( Z ) = - i * C8_ASIN ( i * Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ASINH, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;

  c2 = c8_i ( );
  c3 = c2 * c1;
  c4 = c8_asin ( c3 );
  c5 = c2 * c4;
  c6 = - c5;

  return c6;
}
/******************************************************************************/

double complex c8_atan ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ATAN evaluates the inverse tangent of a C8.

  Discussion:

    Here we use the relationship:

      C8_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ATAN, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;
  double complex c7;
  double complex c8;
  double complex c9;
  double complex cx;

  c2 = c8_i ( );
  c3 = c8_one ( );
  c4 = c8_mul ( c2, c1 );
  c5 = c8_sub ( c3, c4 );
  c6 = c8_add ( c3, c4 );
  c7 = c8_div ( c5, c6 );

  c8 = c8_log ( c7 );
  c9 = c8_mul ( c2, c8 );
  cx = c9 / 2.0;

  return cx;
}
/******************************************************************************/

double complex c8_atanh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_ATANH evaluates the inverse hyperbolic tangent of a C8.

  Discussion:

    Here we use the relationship:

      C8_ATANH ( Z ) = - i * C8_ATAN ( i * Z ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ATANH, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;

  c2 = c8_i ( );

  c3 = c8_mul ( c2, c1 );
  c4 = c8_atan ( c3 );
  c5 = c8_mul ( c2, c4 );
  c6 = c8_neg ( c5 );

  return c6;
}
/******************************************************************************/

double complex c8_conj ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_CONJ conjugates a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_CONJ, the function value.
*/
{
  double c1_norm;
  double complex c2;
  double c2_imag;
  double c2_real;

  c2 = conj ( c1 );

  return c2;
}
/******************************************************************************/

void c8_copy ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_COPY copies a C8.

  Discussion:

    The order of the arguments may seem unnatural, but it is arranged so
    that the call

      c8_copy ( c1, c2 )

    mimics the assignment

      c1 = c2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, double complex C1, the copy of C2.

    Input, double complex C2, the value to be copied.
*/
{
  c1 = c2;

  return;
}
/******************************************************************************/

double complex c8_cos ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_COS evaluates the cosine of a C8.

  Discussion:

    We use the relationship:

      C8_COS ( C ) = ( C8_EXP ( i * C ) + C8_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_COS, the function value.
*/
{
  double complex c2;
  double r;

  c2 = ( cexp ( c1 * I ) + cexp ( - c1 * I ) ) / 2.0;

  return c2;
}
/******************************************************************************/

double complex c8_cosh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_COSH evaluates the hyperbolic cosine of a C8.

  Discussion:

    A C8 is a C8_COMPLEX value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_COSH, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;

  c2 = c8_exp ( c1 );

  c3 = c8_neg ( c1 );
  c4 = c8_exp ( c3 );

  c5 = c8_add ( c2, c4 );
  c6 = c8_div_r8 ( c5, 2.0 );

  return c6;
}
/******************************************************************************/

double complex c8_cube_root ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_CUBE_ROOT computes the principal cube root of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_CUBE_ROOT, the function value.
*/
{
  double complex c2;
  double r;
  double t;

  c8_to_polar ( c1, &r, &t );

  r = pow ( r, 1.0 / 3.0 );
  t = t / 3.0;

  c2 = polar_to_c8 ( r, t );

  return c2;
}
/******************************************************************************/

double complex c8_div ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_DIV divides two C8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, C2, the arguments.

    Output, double complex C8_DIV, the function value.
*/
{
  double c2_norm;
  double complex c3;
  double c3_imag;
  double c3_real;

  c2_norm = c8_abs ( c2 );

  c3_real = ( creal ( c1 ) * creal ( c2 ) 
            + cimag ( c1 ) * cimag ( c2 ) ) / c2_norm / c2_norm;

  c3_imag = ( cimag ( c1 ) * creal ( c2 ) 
            - creal ( c1 ) * cimag ( c2 ) ) / c2_norm / c2_norm;

  c3 = c3_real + c3_imag * I;

  return c3;
}
/******************************************************************************/

double complex c8_div_r8 ( double complex c1, double r )

/******************************************************************************/
/*
  Purpose:

    C8_DIV_R8 divides a C8 by an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the value to be divided.

    Input, double R, the divisor.

    Output, double complex C8_DIV_R8, the function value.
*/
{
  double complex c2;

  c2 = c1 / r;

  return c2;
}
/******************************************************************************/

double complex c8_exp ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_EXP exponentiates a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_EXP, the function value.
*/
{
  double complex c2;

  c2 = cexp ( c1 );

  return c2;
}
/******************************************************************************/

double complex c8_i ( void )

/******************************************************************************/
/*
  Purpose:

    C8_I returns the value of I as a C8

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, double complex C8_I, the value of I.
*/
{
  double complex c1;

  c1 = I;

  return c1;
}
/******************************************************************************/

double c8_imag ( double complex c )

/******************************************************************************/
/*
  Purpose:

    C8_IMAG returns the imaginary part of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the argument.

    Output, double C8_IMAG, the function value.
*/
{
  double value;

  value = cimag ( c );

  return value;
}
/******************************************************************************/

double complex c8_inv ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_INV inverts a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_INV, the function value;
*/
{
  double complex c2;

  c2 = 1.0 / c1;

  return c2;
}
/******************************************************************************/

double complex c8_log ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_LOG evaluates the logarithm of a C8.

  Discussion:

    Here we use the relationship:

      C8_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_ATAN, the function value.
*/
{
  double arg;
  double complex c2;
  double mag;

  arg = c8_arg ( c1 );
  mag = c8_mag ( c1 );

  c2 = log ( mag ) + arg * I;

  return c2;
}
/******************************************************************************/

double c8_mag ( double complex c )

/******************************************************************************/
/*
  Purpose:

    C8_MAG returns the magnitude of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the argument.

    Output, double C8_MAG, the function value.
*/
{
  double value;

  value = sqrt ( pow ( creal ( c ), 2 ) 
               + pow ( cimag ( c ), 2 ) );

  return value;
}
/******************************************************************************/

double complex c8_mul ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_MUL multiplies two C8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, C2, the arguments.

    Output, double complex C8_MUL, the function value.
*/
{
  double complex c3;

  c3 = c1 * c2;

  return c3;
}
/******************************************************************************/

double complex c8_neg ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_NEG negates a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_NEG, the function value.
*/
{
  double complex c2;

  c2 = - c1;

  return c2;
}
/******************************************************************************/

double complex c8_normal_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    C8_NORMAL_01 returns a unit pseudonormal C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, double complex C8_NORMAL_01, a unit pseudonormal value.
*/
{
  double r8_pi = 3.141592653589793;
  double v1;
  double v2;
  double complex value;
  double x_c;
  double x_r;

  v1 = r8_uniform_01 ( seed );
  v2 = r8_uniform_01 ( seed );

  x_r = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * r8_pi * v2 );
  x_c = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * r8_pi * v2 );

  value = cartesian_to_c8 ( x_r, x_c );

  return value;
}
/******************************************************************************/

double complex c8_one ( void )

/******************************************************************************/
/*
  Purpose:

    C8_ONE returns the value of 1 as a C8

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, double complex C8_ONE, the value of 1.
*/
{
  double complex c1;

  c1 = 1.0;

  return c1;
}
/******************************************************************************/

double c8_real ( double complex c )

/******************************************************************************/
/*
  Purpose:

    C8_REAL returns the real part of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the complex number.

    Output, double C8_REAL, the function value.
*/
{
  double value;

  value = creal ( c );

  return value;
}
/******************************************************************************/

double complex c8_sin ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_SIN evaluates the sine of a C8.

  Discussion:

    We use the relationship:

      C8_SIN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_SIN, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;
  double complex c7;
  double complex c8;
  double complex c9;
  double complex cx;
  double r;

  c2 = c8_i ( );

  c3 = c8_mul ( c2, c1 );
  c4 = c8_exp ( c3 );

  c5 = c8_neg ( c3 );
  c6 = c8_exp ( c5 );

  c7 = c8_sub ( c4, c6 );

  r = 2.0;
  c8 = c8_div_r8 ( c7, r );
  c9 = c8_mul ( c8, c2 );
  cx = c8_neg ( c9 );

  return cx;
}
/******************************************************************************/

double complex c8_sinh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_SINH evaluates the hyperbolic sine of a C8.

  Discussion:

    We use the relationship:

      C8_SINH ( C ) = ( C8_EXP ( C ) - C8_EXP ( - i * C ) ) / 2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_SINH, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;
  double r;

  c2 = c8_exp ( c1 );

  c3 = c8_neg ( c1 );
  c4 = c8_exp ( c3 );

  c5 = c8_sub ( c2, c4 );

  r = 2.0;
  c6 = c8_div_r8 ( c5, r );

  return c6;
}
/******************************************************************************/

double complex c8_sqrt ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_SQRT computes a square root of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_SQRT, the function value.
*/
{
  double complex c2;
  double r;
  double t;

  c8_to_polar ( c1, &r, &t );

  r = sqrt ( r );
  t = t / 2.0;

  c2 = polar_to_c8 ( r, t );

  return c2;
}
/******************************************************************************/

double complex c8_sub ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_SUB subtracts two C8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, C2, the arguments.

    Output, double complex C8_SUB, the function value.
*/
{
  double complex c3;

  c3 = c1 - c2;

  return c3;
}
/******************************************************************************/

void c8_swap ( double complex c1, double complex c2 )

/******************************************************************************/
/*
  Purpose:

    C8_SWAP swaps two C8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, double complex C1, C2, the arguments.
*/
{
  double complex c3;

  c3 = c1;
  c1 = c2;
  c2 = c3;

  return;
}
/******************************************************************************/

double complex c8_tan ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_TAN evaluates the tangent of a C8.

  Discussion:

    We use the relationship:

      C8_TAN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) 
                         / ( C8_EXP ( I * C ) + C8_EXP ( - i * C ) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_TAN, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;
  double complex c7;
  double complex c8;
  double complex c9;
  double complex cx;
  double complex ce;

  c2 = c8_i ( );
  c3 = c8_mul ( c2, c1 );
  c4 = c8_neg ( c3 );
  
  c5 = c8_exp ( c3 );
  c6 = c8_exp ( c4 );

  c7 = c8_sub ( c5, c6 );
  c8 = c8_add ( c5, c6 );

  c9 = c8_div ( c7, c8 );
  cx = c8_mul ( c2, c9 );
  ce = c8_neg ( cx );

  return ce;
}
/******************************************************************************/

double complex c8_tanh ( double complex c1 )

/******************************************************************************/
/*
  Purpose:

    C8_TANH evaluates the hyperbolic tangent of a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C1, the argument.

    Output, double complex C8_TANH, the function value.
*/
{
  double complex c2;
  double complex c3;
  double complex c4;
  double complex c5;
  double complex c6;
  double complex c7;

  c2 = c8_exp ( c1 );

  c3 = c8_neg ( c1 );
  c4 = c8_exp ( c3 );

  c5 = c8_sub ( c2, c4 );
  c6 = c8_add ( c2, c4 );

  c7 = c8_div ( c5, c6 );

  return c7;
}
/******************************************************************************/

void c8_to_cartesian ( double complex c, double *x, double *y )

/******************************************************************************/
/*
  Purpose:

    C8_TO_CARTESIAN converts a C8 to Cartesian form.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the argument.

    Output, double *X, *Y, the Cartesian form.
*/
{
  *x = creal ( c );
  *y = cimag ( c );

  return;
}
/******************************************************************************/

void c8_to_polar ( double complex c, double *r, double *theta )

/******************************************************************************/
/*
  Purpose:

    C8_TO_POLAR converts a C8 to polar form.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double complex C, the argument.

    Output, double *R, *THETA, the polar form.
*/
{
  *r = c8_abs ( c );
  *theta = c8_arg ( c );

  return;
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

    Output, double complex C8_UNIFORM_01, a pseudorandom complex value.
*/
{
  int i4_huge = 2147483647;
  int k;
  double r;
  double r8_pi = 3.141592653589793;
  double theta;
  double complex value;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "C8_UNIFORM_01 - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
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

  value = r * cos ( theta ) + r * sin ( theta ) * I;
  
  return value;
}
/******************************************************************************/

double complex c8_zero ( void )

/******************************************************************************/
/*
  Purpose:

    C8_ZERO returns the value of 0 as a C8

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Output, double complex C8_ZERO, the value of 0.
*/
{
  double complex c1;

  c1 = 0.0;

  return c1;
}
/******************************************************************************/

double complex *c8mat_copy_new ( int m, int n, double complex a1[] )

/******************************************************************************/
/*
  Purpose:

    C8MAT_COPY_NEW copies one C8MAT to a "new" C8MAT.

  Discussion:

    A C8MAT is a doubly dimensioned array of C8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double complex A1[M*N], the matrix to be copied.

    Output, double complex C8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  double complex *a2;
  int i;
  int j;

  a2 = ( double complex * ) malloc ( m * n * sizeof ( double complex ) );

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

double complex *c8mat_identity ( int n )

/******************************************************************************/
/*
  Purpose:

    C8MAT_IDENTITY sets a C8MAT to the identity.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Output, double complex C8MAT_IDENTITY[N*N], the matrix.
*/
{
  double complex *a;
  int i;
  int j;

  a = ( double complex * ) malloc ( n * n * sizeof ( double complex ) );

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

double complex *c8mat_mm_new ( int n1, int n2, int n3, double complex a[], 
  double complex b[] )

/******************************************************************************/
/*
  Purpose:

    C8MAT_MM_NEW multiplies two C8MAT's.

  Discussion:

    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double complex A[N1*N2], double complex B[N2*N3], 
    the matrices to multiply.

    Output, double complex C[N1*N3], the product matrix C = A * B.
*/
{
  double complex *c;
  int i;
  int j;
  int k;

  c = ( double complex * ) malloc ( n1 * n3 * sizeof ( double complex ) );

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

double c8mat_norm_fro ( int m, int n, double complex a[] )

/******************************************************************************/
/*
  Purpose:

    C8MAT_NORM_FRO returns the Frobenius norm of a C8MAT.

  Discussion:

    The Frobenius norm is defined as

      C8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) |A(I,J)| )

    The matrix Frobenius-norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      c8vec_norm_l2 ( A*x ) <= c8mat_norm_fro ( A ) * c8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of the matrix.

    Input, double complex A[M*N], the matrix.

    Output, double C8MAT_NORM_FRO, the Frobenius norm of A.
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
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

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
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

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

double complex *c8mat_uniform_01 ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.

  Discussion:

    A C8MAT is a matrix of double precision complex values.

    The angles should be uniformly distributed between 0 and 2 * PI,
    the square roots of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double complex C[M*N], the pseudorandom complex matrix.
*/
{
  double complex *c;
  int i;
  int j;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  c = ( double complex * ) malloc ( m * n * sizeof ( double complex ) );

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

      r = sqrt ( ( double ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * ( cos ( theta ) + sin ( theta ) * I );
    }
  }

  return c;
}
/******************************************************************************/

double complex *c8mat_zero_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    C8MAT_ZERO_NEW returns a new zeroed C8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double complex C8MAT_ZERO[M*N], the new zeroed matrix.
*/
{
  double complex *a;
  int i;
  int j;

  a = ( double complex * ) malloc ( m * n * sizeof ( double complex ) );

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

double complex *c8vec_copy_new ( int n, double complex a1[] )

/******************************************************************************/
/*
  Purpose:

    C8VEC_COPY_NEW copies a C8VEC to a "new" C8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double complex A1[N], the vector to be copied.

    Output, double complex C8VEC_COPY_NEW[N], the copy of A1.
*/
{
  double complex *a2;
  int i;

  a2 = ( double complex * ) malloc ( n * sizeof ( double complex ) );

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
/******************************************************************************/

void c8vec_print_part ( int n, double complex a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    C8VEC_PRINT_PART prints "part" of a C8VEC.

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

    09 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double complex A[N], the vector to be printed.

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
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal ( a[i] ), cimag ( a[i] ) );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal( a[i] ), cimag ( a[i] ) );
    }
    fprintf ( stdout, "  ........  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal( a[i] ), cimag ( a[i] ) );
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14f  %14f\n", i, creal( a[i] ), cimag ( a[i] ) );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14f  %14f  ...more entries...\n", 
      i, creal( a[i] ), cimag ( a[i] ) );
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

    09 July 2011

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
  int i4_huge = 2147483647;
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "C8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  c = malloc ( n * sizeof ( double complex ) );

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

    theta = 2.0 * pi * ( ( double ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * cos ( theta ) + r * sin ( theta ) * I;
  }

  return c;
}
/******************************************************************************/

double complex *c8vec_unity ( int n )

/******************************************************************************/
/*
  Purpose:

    C8VEC_UNITY returns the N roots of unity in a C8VEC.

  Discussion:

    A C8VEC is a vector of C8's.

    X(1:N) = exp ( 2 * PI * (0:N-1) / N )

    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, double complex C8VEC_UNITY[N], the N roots of unity.
*/
{
  double complex *a;
  int i;
  double pi = 3.141592653589793;
  double theta;

  a = ( double complex * ) malloc ( n * sizeof ( double complex ) );

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( double ) ( 2 * i ) / ( double ) ( n );
    a[i] = cos ( theta ) + sin ( theta ) * I;
  }

  return a;
}
/******************************************************************************/

double complex cartesian_to_c8 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    CARTESIAN_TO_C8 converts a Cartesian form to a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the Cartesian form.

    Output, double complex CARTESIAN_TO_C8, the complex number.
*/
{
  double complex c;

  c = x + y * I;

  return c;
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

double complex polar_to_c8 ( double r, double theta )

/******************************************************************************/
/*
  Purpose:

    POLAR_TO_C8 converts a polar form to a C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double R, THETA, the polar form.

    Output, double complex POLAR_TO_C8, the complex number.
*/
{
  double complex c;

  c = r * cos ( theta ) + r * sin ( theta ) * I;

  return c;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double complex r8_csqrt ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_CSQRT returns the complex square root of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose square root is desired.

    Output, double complex R8_CSQRT, the square root of X:
*/
{
  double argument;
  double magnitude;
  double pi = 3.141592653589793;
  double complex value;

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
    argument = pi;
  }

  magnitude = sqrt ( magnitude );
  argument = argument / 2.0;

  value = magnitude * ( double complex ) ( cos ( argument ), sin ( argument ) );

  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

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
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
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

      seed = 16807 * seed mod ( 2**31 - 1 )
      r8_uniform_01 = seed / ( 2**31 - 1 )

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
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "R8_UNIFORM_01 - Fatal error!\n" );
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
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void r8poly2_root ( double a, double b, double c, double complex *r1,
  double complex *r2 )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_ROOT returns the two roots of a quadratic polynomial.

  Discussion:

    The polynomial has the form:

      A * X^2 + B * X + C = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2005

  Parameters:

    Input, double A, B, C, the coefficients of the polynomial.
    A must not be zero.

    Output, double complex *R1, *R2, the roots of the polynomial, which
    might be real and distinct, real and equal, or complex conjugates.
*/
{
  double disc;
  double complex q;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY2_ROOT - Fatal error!\n" );
    fprintf ( stderr, "  The coefficient A is zero.\n" );
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  q = -0.5 * ( b + r8_sign ( b ) * r8_csqrt ( disc ) );
  *r1 = q / a;
  *r2 = c / q;

  return;
}
/******************************************************************************/

void r8poly3_root ( double a, double b, double c, double d,
  double complex *r1, double complex *r2, double complex *r3 )

/******************************************************************************/
/*
  Purpose:

    R8POLY3_ROOT returns the three roots of a cubic polynomial.

  Discussion:

    The polynomial has the form

      A * X^3 + B * X^2 + C * X + D = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Parameters:

    Input, double A, B, C, D, the coefficients of the polynomial.
    A must not be zero.

    Output, double complex *R1, *R2, *R3, the roots of the polynomial, which
    will include at least one real root.
*/
{
  double complex i;
  double pi = 3.141592653589793;
  double q;
  double r;
  double s1;
  double s2;
  double temp;
  double theta;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY3_ROOT - Fatal error!\n" );
    fprintf ( stderr, "  A must not be zero.\n" );
    exit ( 1 );
  }

  i = ( double complex ) ( 0.0, 1.0 );

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
    s1 = r8_sign ( temp ) * pow ( r8_abs ( temp ), 1.0 / 3.0 );

    temp = - r - sqrt ( r * r - q * q * q );
    s2 = r8_sign ( temp ) * pow ( r8_abs ( temp ), 1.0 / 3.0 );

    *r1 = s1 + s2;
    *r2 = -0.5 * ( s1 + s2 ) + i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
    *r3 = -0.5 * ( s1 + s2 ) - i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
  }

  *r1 = *r1 - b / ( 3.0 * a );
  *r2 = *r2 - b / ( 3.0 * a );
  *r3 = *r3 - b / ( 3.0 * a );

  return;
}
/******************************************************************************/

void r8poly4_root ( double a, double b, double c, double d, double e,
  double complex *r1, double complex *r2, double complex *r3,
  double complex *r4 )

/******************************************************************************/
/*
  Purpose:

    R8POLY4_ROOT returns the four roots of a quartic polynomial.

  Discussion:

    The polynomial has the form:

      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Parameters:

    Input, double A, B, C, D, the coefficients of the polynomial.
    A must not be zero.

    Output, double complex *R1, *R2, *R3, *R4, the roots of the polynomial.
*/
{
  double a3;
  double a4;
  double b3;
  double b4;
  double c3;
  double c4;
  double d3;
  double d4;
  double complex p;
  double complex q;
  double complex r;
  double complex zero;

  zero = 0.0;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY4_ROOT - Fatal error!\n" );
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
  r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
/*
  Choose one root of the cubic, here R1.

  Set R = sqrt ( 0.25 * A4^2 - B4 + R1 )
*/
  r = c8_sqrt ( 0.25 * a4 * a4 - b4  + (*r1) );

  if ( creal ( r ) != 0.0 || cimag ( r ) != 0.0 )
  {
    p = c8_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        + 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );

    q = c8_sqrt ( 0.75 * a4 * a4 - r * r - 2.0 * b4
        - 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );
  }
  else
  {
    p = c8_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      + 2.0 * c8_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );

    q = c8_sqrt ( 0.75 * a4 * a4 - 2.0 * b4
      - 2.0 * c8_sqrt ( (*r1) * (*r1) - 4.0 * d4 ) );
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
