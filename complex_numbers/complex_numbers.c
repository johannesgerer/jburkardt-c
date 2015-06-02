# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    COMPLEX_NUMBERS is a program which demonstrates the use of complex numbers.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "COMPLEX_NUMBERS\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrate complex number usage.\n" );
/*
  Single precision complex numbers: "complex".
*/
  test01 ( );
  test02 ( );
  test03 ( );
/*
  Double precision complex numbers: "double complex".
*/
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "COMPLEX_NUMBERS\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 demonstrate declaration and assignment for complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
/*
  Declare a complex number A.
  Declare a complex vector B.
  Declare a complex array C.
*/
  complex a;
  complex b[3];
  complex c[2][2];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Declare a COMPLEX variable.\n" );
  printf ( "  Assign value with an = statement.\n" );
/*
  Assign values to A, B, and C.
*/
  a = 1.0 + 2.0 * I;

  b[0] = 1.0 + 2.0 * I;
  b[1] = 3.0 + 4.0 * I;
  b[2] = 5.0 + 6.0 * I;

  c[0][0] = 1.0 + 0.1 * I;
  c[0][1] = 1.0 + 0.2 * I;
  c[1][0] = 2.0 + 0.1 * I;
  c[1][1] = 2.0 + 0.2 * I;
/*
  Print them.
*/
  printf ( "\n" );
  printf ( "  Scalar A:\n" );
  printf ( "\n" );

  printf ( "  (%g, %g)\n", creal ( a ), cimag ( a ) );

  printf ( "\n" );
  printf ( "  Vector B:\n" );
  printf ( "\n" );

  for ( i = 0; i < 3; i++ )
  {
    printf ( "  (%g, %g)\n", creal ( b[i] ), cimag ( b[i] ) );
  }

  printf ( "\n" );
  printf ( "  Array C:\n" );
  printf ( "\n" );

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      printf ( "  (%g, %g)", creal ( c[i][j] ), cimag ( c[i][j] ) );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrate declaration and initialization for complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
/*
  Declare and initialize a complex number A.
  Declare and initialize a complex vector B.
  Declare and initialize a complex array C.
*/
  complex a = {1.0 + 2.0 * I};
  complex b[3] = {1.0 + 2.0 * I, 3.0 + 4.0 * I, 5.0 + 6.0 * I };
  complex c[2][2] = { { 1.0 + 0.1 * I, 1.0 + 0.2 * I},
                      { 2.0 + 0.1 * I, 2.0 + 0.2 * I} };
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Declare a COMPLEX variable.\n" );
  printf ( "  Initialize value as part of the declaration.\n" );
/*
  Print them.
*/
  printf ( "\n" );
  printf ( "  Scalar A:\n" );
  printf ( "\n" );

  printf ( "  (%g, %g)\n", creal ( a ), cimag ( a ) );

  printf ( "\n" );
  printf ( "  Vector B:\n" );
  printf ( "\n" );

  for ( i = 0; i < 3; i++ )
  {
    printf ( "  (%g, %g)\n", creal ( b[i] ), cimag ( b[i] ) );
  }

  printf ( "\n" );
  printf ( "  Array C:\n" );
  printf ( "\n" );

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      printf ( "  (%g, %g)", creal ( c[i][j] ), cimag ( c[i][j] ) );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03: intrinsic functions for complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 November 2010

  Author:

    John Burkardt
*/
{
  complex a = {1.0 + 2.0 * I};
  int i;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Apply intrinsic functions to COMPLEX variables\n" );
/*
  Print them.
*/
  printf ( "\n" );
/*
  Note that "I" by itself is NOT a complex number, nor is it the
  imaginary unit.  You have to cast it to ( complex ) or ( double complex )
  or multiply it by a float or double before it results in a numerical
  result.
*/
  printf ( "  ( complex ) I =      (%14.6g,%14.6g)\n", ( complex ) I );
  printf ( "  a =                  (%14.6g,%14.6g)\n", a );
  printf ( "  - a =                (%14.6g,%14.6g)\n", - a );
  printf ( "  a + 3 =              (%14.6g,%14.6g)\n", a + 3 );
  printf ( "  a + (0,5) =          (%14.6g,%14.6g)\n", a + 5 * I );
  printf ( "  4 * a =              (%14.6g,%14.6g)\n", 4 * a );
  printf ( "  a / 8 =              (%14.6g,%14.6g)\n", a / 8 );
  printf ( "  a * a =              (%14.6g,%14.6g)\n", a * a );
  printf ( "  cpow ( a, 2 ) =      (%14.6g,%14.6g)\n", cpow ( a, 2 ) );
  printf ( "  cpow ( 2, a ) =      (%14.6g,%14.6g)\n", cpow ( 2, a ) );
  printf ( "  cpow ( a, a ) =      (%14.6g,%14.6g)\n", cpow ( a, a ) );
  printf ( "  1/a =                (%14.6g,%14.6g)\n", 1.0 / a );
  printf ( "\n" );
  printf ( "  cabs(a) =             %14.6g\n",         cabs ( a ) );
  printf ( "  cacos(a) =           (%14.6g,%14.6g)\n", cacos ( a ) );
  printf ( "  cacosh(a) =          (%14.6g,%14.6g)\n", cacosh ( a ) );
  printf ( "  carg(a) =             %14.6g\n",         carg ( a ) );
  printf ( "  casin(a) =           (%14.6g,%14.6g)\n", casin ( a ) );
  printf ( "  casinh(a) =          (%14.6g,%14.6g)\n", casinh ( a ) );
  printf ( "  catan(a) =           (%14.6g,%14.6g)\n", catan ( a ) );
  printf ( "  catanh(a) =          (%14.6g,%14.6g)\n", catanh ( a ) );
  printf ( "  ccos(a) =            (%14.6g,%14.6g)\n", ccos ( a ) );
  printf ( "  ccosh(a) =           (%14.6g,%14.6g)\n", ccosh ( a ) );
  printf ( "  cexp(a) =            (%14.6g,%14.6g)\n", cexp ( a ) );
  printf ( "  cimag(a) =            %14.6g\n",         cimag ( a ) );
  printf ( "  clog(a) =            (%14.6g,%14.6g)\n", clog ( a ) );
  printf ( "  (complex)(1) =       (%14.6g,%14.6g)\n", (complex) ( 1 ) );
  printf ( "  (complex)(4.0) =     (%14.6g,%14.6g)\n", (complex) ( 4.0 ) );
  printf ( "  conj(a) =            (%14.6g,%14.6g)\n", conj ( a ) );
  printf ( "  cproj(a) =           (%14.6g,%14.6g)\n", cproj ( a ) );
  printf ( "  creal(a) =            %14.6g\n",         creal ( a ) );
  printf ( "  csin(a) =            (%14.6g,%14.6g)\n", csin ( a ) );
  printf ( "  csinh(a) =           (%14.6g,%14.6g)\n", csinh ( a ) );
  printf ( "  csqrt(a) =           (%14.6g,%14.6g)\n", csqrt ( a ) );
  printf ( "  ctan(a) =            (%14.6g,%14.6g)\n", ctan ( a ) );
  printf ( "  ctanh(a) =           (%14.6g,%14.6g)\n", ctanh ( a ) );
  printf ( "  (int)(a) =            %10d\n",           ( int ) ( a ) );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 demonstrate declaration and assignment for double complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
/*
  Declare a double complex number A.
  Declare a double complex vector B.
  Declare a double complex array C.
*/
  double complex a;
  double complex b[3];
  double complex c[2][2];
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Declare a DOUBLE COMPLEX variable.\n" );
  printf ( "  Assign value with an = statement.\n" );
/*
  Assign values to A, B, and C.
*/
  a = 1.0 + 2.0 * I;

  b[0] = 1.0 + 2.0 * I;
  b[1] = 3.0 + 4.0 * I;
  b[2] = 5.0 + 6.0 * I;

  c[0][0] = 1.0 + 0.1 * I;
  c[0][1] = 1.0 + 0.2 * I;
  c[1][0] = 2.0 + 0.1 * I;
  c[1][1] = 2.0 + 0.2 * I;
/*
  Print them.
*/
  printf ( "\n" );
  printf ( "  Scalar A:\n" );
  printf ( "\n" );

  printf ( "  (%g, %g)\n", creal ( a ), cimag ( a ) );

  printf ( "\n" );
  printf ( "  Vector B:\n" );
  printf ( "\n" );

  for ( i = 0; i < 3; i++ )
  {
    printf ( "  (%g, %g)\n", creal ( b[i] ), cimag ( b[i] ) );
  }

  printf ( "\n" );
  printf ( "  Array C:\n" );
  printf ( "\n" );

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      printf ( "  (%g, %g)", creal ( c[i][j] ), cimag ( c[i][j] ) );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05: declaration and initialization for double complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
/*
  Declare and initialize a double complex number A.
  Declare and initialize a double complex vector B.
  Declare and initialize a double complex array C.
*/
  double complex a = {1.0 + 2.0 * I};
  double complex b[3] = {1.0 + 2.0 * I, 3.0 + 4.0 * I, 5.0 + 6.0 * I };
  double complex c[2][2] = { { 1.0 + 0.1 * I, 1.0 + 0.2 * I},
                             { 2.0 + 0.1 * I, 2.0 + 0.2 * I} };
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Declare a DOUBLE COMPLEX variable.\n" );
  printf ( "  Initialize value as part of the declaration.\n" );
/*
  Print them.
*/
  printf ( "\n" );
  printf ( "  Scalar A:\n" );
  printf ( "\n" );

  printf ( "  (%g, %g)\n", a );

  printf ( "\n" );
  printf ( "  Vector B:\n" );
  printf ( "\n" );

  for ( i = 0; i < 3; i++ )
  {
    printf ( "  (%g, %g)\n", b[i] );
  }

  printf ( "\n" );
  printf ( "  Array C:\n" );
  printf ( "\n" );

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      printf ( "  (%g, %g)", c[i][j] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06: intrinsic functions for double complex variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 November 2010

  Author:

    John Burkardt
*/
{
  double complex a = {1.0 + 2.0 * I};

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Apply intrinsic functions to DOUBLE COMPLEX variables\n" );
/*
  Print them.
*/
  printf ( "\n" );
/*
  Note that "I" by itself is NOT a complex number, nor is it the
  imaginary unit.  You have to cast it to ( complex ) or ( double complex )
  or multiply it by a float or double before it results in a numerical
  result.
*/
  printf ( "  ( double complex ) I =  (%14.6g,%14.6g)\n", ( double complex ) I );
  printf ( "  a =                     (%14.6g,%14.6g)\n", a );
  printf ( "  - a =                   (%14.6g,%14.6g)\n", - a );
  printf ( "  a + 3 =                 (%14.6g,%14.6g)\n", a + 3 );
  printf ( "  a + (0,5) =             (%14.6g,%14.6g)\n", a + ( 0, 5 ) );
  printf ( "  4 * a =                 (%14.6g,%14.6g)\n", 4 * a );
  printf ( "  a / 8 =                 (%14.6g,%14.6g)\n", a / 8 );
  printf ( "  a * a =                 (%14.6g,%14.6g)\n", a * a );
  printf ( "  cpow ( a, 2 ) =         (%14.6g,%14.6g)\n", cpow ( a, 2 ) );
  printf ( "  cpow ( 2, a ) =         (%14.6g,%14.6g)\n", cpow ( 2, a ) );
  printf ( "  cpow ( a, a ) =         (%14.6g,%14.6g)\n", cpow ( a, a ) );
  printf ( "  1/a =                   (%14.6g,%14.6g)\n", 1.0 / a );
  printf ( "\n" );
  printf ( "  cabs(a) =                %14.6g\n",         cabs ( a ) );
  printf ( "  cacos(a) =              (%14.6g,%14.6g)\n", cacos ( a ) );
  printf ( "  cacosh(a) =             (%14.6g,%14.6g)\n", cacosh ( a ) );
  printf ( "  carg(a) =                %14.6g\n",         carg ( a ) );
  printf ( "  casin(a) =              (%14.6g,%14.6g)\n", casin ( a ) );
  printf ( "  casinh(a) =             (%14.6g,%14.6g)\n", casinh ( a ) );
  printf ( "  catan(a) =              (%14.6g,%14.6g)\n", 
    creal ( catan ( a ) ), cimag ( catan ( a ) ) );
  printf ( "  catanh(a) =             (%14.6g,%14.6g)\n", 
    creal ( catanh ( a ) ), cimag ( catanh ( a ) ) );
  printf ( "  ccos(a) =               (%14.6g,%14.6g)\n", 
    creal ( ccos ( a ) ), cimag ( ccos ( a ) ) );
  printf ( "  ccosh(a) =              (%14.6g,%14.6g)\n", 
    creal ( ccosh ( a ) ), cimag ( ccosh ( a ) ) );
  printf ( "  cexp(a) =               (%14.6g,%14.6g)\n", 
    creal ( cexp ( a ) ), cimag ( cexp ( a ) ) );
  printf ( "  cimag(a) =               %14.6g\n",         cimag ( a ) );
  printf ( "  clog(a) =               (%14.6g,%14.6g)\n", 
    creal ( clog ( a ) ), cimag ( clog ( a ) ) );
  printf ( "  (double complex)(1) =   (%14.6g,%14.6g)\n", 
    creal ( ( double complex ) ( 1 ) ), cimag ( ( double complex ) ( 1 ) ) );
  printf ( "  (double complex)(4.0) = (%14.6g,%14.6g)\n", 
    creal ( ( double complex ) ( 4.0 ) ), cimag ( ( double complex ) ( 4.0 ) ) );
  printf ( "  conj(a) =               (%14.6g,%14.6g)\n", 
    creal ( conj ( a ) ), cimag ( conj ( a ) ) );
  printf ( "  cproj(a) =              (%14.6g,%14.6g)\n", 
    creal ( cproj ( a ) ), cimag ( cproj ( a ) ) );
  printf ( "  creal(a) =               %14.6g\n",         creal ( a ) );
  printf ( "  csin(a) =               (%14.6g,%14.6g)\n", 
    creal ( csin ( a ) ), cimag ( csin ( a ) ) );
  printf ( "  csinh(a) =              (%14.6g,%14.6g)\n", 
    creal ( csinh ( a ) ), cimag ( csinh ( a ) ) );
  printf ( "  csqrt(a) =              (%14.6g,%14.6g)\n", 
    creal ( csqrt ( a ) ), cimag ( csqrt ( a ) ) );
  printf ( "  ctan(a) =               (%14.6g,%14.6g)\n", 
    creal ( ctan ( a ) ), cimag ( ctan ( a ) ) );
  printf ( "  ctanh(a) =              (%14.6g,%14.6g)\n", 
    creal ( ctanh ( a ) ), cimag ( ctanh ( a ) ) );
  printf ( "  (int)(a) =               %10d\n",           ( int ) ( a ) );

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
