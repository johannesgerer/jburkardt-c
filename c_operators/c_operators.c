# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test_amp ( void );
void test_ampamp ( void );
void test_bang ( void );
void test_bar ( void );
void test_barbar ( void );
void test_caret ( void );
void test_lshiftlshift ( void );
void test_minus ( void );
void test_plus ( void );
void test_plusplus ( void );
void test_rshiftrshift ( void );
void test_slash ( void );
void test_star ( void );
void test_twiddle ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for C_OPERATORS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "C_OPERATORS:\n" );
  printf ( "  Demonstrate the operators available in C.\n" );

  test_amp ( );
  test_ampamp ( );
  test_bang ( );
  test_bar ( );
  test_barbar ( );
  test_caret ( );
  test_lshiftlshift ( );
  test_minus ( );
  test_plus ( );
  test_plusplus ( );
  test_rshiftrshift ( );
  test_slash ( );
  test_star ( );
  test_twiddle ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "C_OPERATORS:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_amp ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_AMP demonstrates &.

  Discussion:

    & implements the bitwise AND function.

    The operation of this function is easiest to see when we look at data
    that is unsigned characters, that is, integers between 0 and 255.
    In binary, some examples include

      00000000 =   0
      00000011 =   3
      00000101 =   5
      00001101 =  13
      01010000 =  80
      01010101 =  85
      10000000 = 128
      11110010 = 242
      11111111 = 255
      
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_AMP:\n" );
  printf ( "  Demonstrate &, which implements the bitwise AND operator.\n" );
  printf ( "  This is most useful when the data is unsigned characters,\n" );
  printf ( "  that is, binary integers from 0 to 255.\n" );
  printf ( "\n" );
  printf ( "  Type                A      B    A&B\n" );
  printf ( "\n" );

  a_uchar = 3;
  b_uchar = 5;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 242;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 13;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 13;
  b_uchar = 85;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 85;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 255;
  b_uchar = 80;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 0;
  b_uchar = 128;
  c_uchar = a_uchar & b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_ampamp ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_AMPAMP demonstrates &&.

  Discussion:

    && implements the logical AND function.

    In C, 0 counts as "false", and nonzero value as true.
    
    The bool datatype can also be used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  bool a_bool;
  bool b_bool;
  bool c_bool;
  double a_double;
  double b_double;
  double c_double;
  float a_float;
  float b_float;
  float c_float;
  int a_int;
  int b_int;
  int c_int;

  printf ( "\n" );
  printf ( "TEST_AMPAMP:\n" );
  printf ( "  Demonstrate &&, which implements the logical AND operator.\n" );
  printf ( "  For numeric data, 0 acts as false, nonzero as true.\n" );
  printf ( "\n" );
  printf ( "  Type              A    B    A&&B\n" );
  printf ( "\n" );

  a_bool = true;
  b_bool = true;
  c_bool = a_bool && b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = true;
  b_bool = false;
  c_bool = a_bool && b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = false;
  b_bool = true;
  c_bool = a_bool && b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = false;
  b_bool = false;
  c_bool = a_bool && b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  printf ( "\n" );

  a_double = 1.0;
  b_double = 1.0;
  c_double = a_double && b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 1.0;
  b_double = 0.0;
  c_double = a_double && b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 0.0;
  b_double = 1.0;
  c_double = a_double && b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 0.0;
  b_double = 0.0;
  c_double = a_double && b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  printf ( "\n" );

  a_float = 1.0;
  b_float = 1.0;
  c_float = a_float && b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 1.0;
  b_float = 0.0;
  c_float = a_float && b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 0.0;
  b_float = 1.0;
  c_float = a_float && b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 0.0;
  b_float = 0.0;
  c_float = a_float && b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  printf ( "\n" );

  a_int = 1;
  b_int = 1;
  c_int = a_int && b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 1;
  b_int = 0;
  c_int = a_int && b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 0;
  b_int = 1;
  c_int = a_int && b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 0;
  b_int = 0;
  c_int = a_int && b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  return;
}
/******************************************************************************/

void test_bang ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_BANG demonstrates !.

  Discussion:

    ! implements logical negation.

    In C, 0 counts as "false", and nonzero value as true.
    
    The bool datatype can also be used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  bool a_bool;
  bool b_bool;
  double a_double;
  double b_double;
  float a_float;
  float b_float;
  int a_int;
  int b_int;

  printf ( "\n" );
  printf ( "TEST_BANG:\n" );
  printf ( "  Demonstrate !, which implements the logical negation operator.\n" );
  printf ( "  For numeric data, 0 acts as false, nonzero as true.\n" );
  printf ( "\n" );
  printf ( "  Type              A     !A\n" );
  printf ( "\n" );

  a_bool = true;
  b_bool = !a_bool;
  printf ( "  bool              %d    %d\n", a_bool, b_bool );
  a_bool = false;
  b_bool = !a_bool;
  printf ( "  bool              %d    %d\n", a_bool, b_bool );

  a_double = -1.0;
  b_double = !a_double;
  printf ( "  double            %f    %f\n", a_double, b_double );
  a_double = 0.0;
  b_double = !a_double;
  printf ( "  double            %f    %f\n", a_double, b_double );
  a_double = 1.0;
  b_double = !a_double;
  printf ( "  double            %f    %f\n", a_double, b_double );
  a_double = 3.14159265;
  b_double = !a_double;
  printf ( "  double            %f    %f\n", a_double, b_double );

  a_float = -1.0;
  b_float = !a_float;
  printf ( "  float             %f    %f\n", a_float, b_float );
  a_float = 0.0;
  b_float = !a_float;
  printf ( "  float             %f    %f\n", a_float, b_float );
  a_float = 1.0;
  b_float = !a_float;
  printf ( "  float             %f    %f\n", a_float, b_float );
  a_float = 3.14159265;
  b_float = !a_float;
  printf ( "  float             %f    %f\n", a_float, b_float );

  a_int = -1;
  b_int = !a_int;
  printf ( "  int               %d    %d\n", a_int, b_int );
  a_int = 0;
  b_int = !a_int;
  printf ( "  int               %d    %d\n", a_int, b_int );
  a_int = 1;
  b_int = !a_int;
  printf ( "  int               %d    %d\n", a_int, b_int );
  a_int = 2;
  b_int = !a_int;
  printf ( "  int               %d    %d\n", a_int, b_int );

  return;
}
/******************************************************************************/

void test_bar ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_BAR demonstrates |.

  Discussion:

    | implements the bitwise OR function.

    The operation of this function is easiest to see when we look at data
    that is unsigned characters, that is, integers between 0 and 255.
    In binary, some examples include

      00000000 =   0
      00000011 =   3
      00000101 =   5
      00001101 =  13
      01010000 =  80
      01010101 =  85
      10000000 = 128
      11110010 = 242
      11111111 = 255
      
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_BAR:\n" );
  printf ( "  Demonstrate |, which implements the bitwise OR operator.\n" );
  printf ( "  This is most useful when the data is unsigned characters,\n" );
  printf ( "  that is, binary integers from 0 to 255.\n" );
  printf ( "\n" );
  printf ( "  Type                A      B    A|B\n" );
  printf ( "\n" );

  a_uchar = 3;
  b_uchar = 5;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 242;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 13;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 13;
  b_uchar = 85;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 85;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 255;
  b_uchar = 80;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 0;
  b_uchar = 128;
  c_uchar = a_uchar | b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_barbar ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_BARBAR demonstrates ||.

  Discussion:

    || implements the logical OR function.

    In C, 0 counts as "false", and nonzero value as true.
    
    The bool datatype can also be used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  bool a_bool;
  bool b_bool;
  bool c_bool;
  double a_double;
  double b_double;
  double c_double;
  float a_float;
  float b_float;
  float c_float;
  int a_int;
  int b_int;
  int c_int;

  printf ( "\n" );
  printf ( "TEST_BARBAR:\n" );
  printf ( "  Demonstrate ||, which implements the logical OR operator.\n" );
  printf ( "  For numeric data, 0 acts as false, nonzero as true.\n" );
  printf ( "\n" );
  printf ( "  Type              A    B    A||B\n" );
  printf ( "\n" );

  a_bool = true;
  b_bool = true;
  c_bool = a_bool || b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = true;
  b_bool = false;
  c_bool = a_bool || b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = false;
  b_bool = true;
  c_bool = a_bool || b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  a_bool = false;
  b_bool = false;
  c_bool = a_bool || b_bool;
  printf ( "  bool              %d    %d    %d\n", a_bool, b_bool, c_bool );

  printf ( "\n" );

  a_double = 1.0;
  b_double = 1.0;
  c_double = a_double || b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 1.0;
  b_double = 0.0;
  c_double = a_double || b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 0.0;
  b_double = 1.0;
  c_double = a_double || b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_double = 0.0;
  b_double = 0.0;
  c_double = a_double || b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  printf ( "\n" );

  a_float = 1.0;
  b_float = 1.0;
  c_float = a_float || b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 1.0;
  b_float = 0.0;
  c_float = a_float || b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 0.0;
  b_float = 1.0;
  c_float = a_float || b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_float = 0.0;
  b_float = 0.0;
  c_float = a_float || b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  printf ( "\n" );

  a_int = 1;
  b_int = 1;
  c_int = a_int || b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 1;
  b_int = 0;
  c_int = a_int || b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 0;
  b_int = 1;
  c_int = a_int || b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_int = 0;
  b_int = 0;
  c_int = a_int || b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  return;
}
/******************************************************************************/

void test_caret ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_CARET demonstrates ^.

  Discussion:

    ^ implements the bitwise exclusive OR function.

    The operation of this function is easiest to see when we look at data
    that is unsigned characters, that is, integers between 0 and 255.
    In binary, some examples include

      00000000 =   0
      00000011 =   3
      00000101 =   5
      00000110 =   6
      00001101 =  13
      01010000 =  80
      01010101 =  85
      01011000 =  88
      10000000 = 128
      10100111 = 167
      10101111 = 175
      11110010 = 242
      11111111 = 255
      
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_CARET:\n" );
  printf ( "  Demonstrate ^, which implements the bitwise exclusive OR operator.\n" );
  printf ( "  This is most useful when the data is unsigned characters,\n" );
  printf ( "  that is, binary integers from 0 to 255.\n" );
  printf ( "\n" );
  printf ( "  Type                A      B    A^B\n" );
  printf ( "\n" );

  a_uchar = 3;
  b_uchar = 5;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 242;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 13;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 13;
  b_uchar = 85;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = 85;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 255;
  b_uchar = 80;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 0;
  b_uchar = 128;
  c_uchar = a_uchar ^ b_uchar;
  printf ( "  uchar             %3d    %3d    %3d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_lshiftlshift ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LSHIFTLSHIFT demonstrates <<.

  Discussion:

    << implements left shift.

    A << B has the value A * 2^B.
    
    A and B must be of an integer type.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  int a_int;
  int b_int;
  int c_int;
  int i;

  printf ( "\n" );
  printf ( "TEST_LSHIFTLSHIFT:\n" );
  printf ( "  Demonstrate <<, which implements the left shift:\n" );
  printf ( "  A << B results in A * 2^B.\n" );
  printf ( "  Generally, B must be nonnegative.\n" );
  printf ( "  Moreover, B should not be so large that the result overflows.\n" );
  printf ( "\n" );
  printf ( "  Type              A     B    A<<B\n" );
  printf ( "\n" );

  for ( i = -1; i < 32; i++ )
  {
    a_int = 11;
    b_int = i;
    c_int = a_int << b_int;

    printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );
  }

  return;
}
/******************************************************************************/

void test_minus ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_MINUS demonstrates -.

  Discussion:

    - implements subtraction.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  char a_char;
  char b_char;
  char c_char;

  double a_double;
  double b_double;
  double c_double;

  float a_float;
  float b_float;
  float c_float;

  int a_int;
  int b_int;
  int c_int;

  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_MINUS:\n" );
  printf ( "  Demonstrate -, which carries out subtraction:\n" );
  printf ( "\n" );
  printf ( "  Type              A     B     A-B\n" );
  printf ( "\n" );

  a_char = 'Z';
  b_char = 'b';
  c_char = a_char - b_char;
  printf ( "  char              %c   %c   %d\n", a_char, b_char, c_char );

  a_double = 9.25;
  b_double = 3.12;
  c_double = a_double - b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_float = 9.25;
  b_float = 3.12;
  c_float = a_float - b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_int = 11;
  b_int = 2;
  c_int = a_int - b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_uchar = 'Z';
  b_uchar = 'b';
  c_uchar = a_uchar - b_uchar;
  printf ( "  unsigned char     %c    %c    %d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_plus ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_PLUS demonstrates +.

  Discussion:

    + implements addition.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  char a_char;
  char b_char;
  char c_char;

  double a_double;
  double b_double;
  double c_double;

  float a_float;
  float b_float;
  float c_float;

  int a_int;
  int b_int;
  int c_int;

  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_PLUS:\n" );
  printf ( "  Demonstrate +, which carries out addition:\n" );
  printf ( "\n" );
  printf ( "  Type              A     B     A+B\n" );
  printf ( "\n" );

  a_char = 'a';
  b_char = 'b';
  c_char = a_char + b_char;
  printf ( "  char              %c   %c   %d\n", a_char, b_char, c_char );

  a_double = 1.25;
  b_double = 3.12;
  c_double = a_double + b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_float = 1.25;
  b_float = 3.12;
  c_float = a_float + b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_int = 1;
  b_int = 2;
  c_int = a_int + b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_uchar = 'a';
  b_uchar = 'b';
  c_uchar = a_uchar + b_uchar;
  printf ( "  unsigned char     %c    %c    %d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_plusplus ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_PLUSPLUS demonstrates ++.

  Discussion:

    ++ implements increment by 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  char a_char;
  char b_char;
  char c_char;

  double a_double;
  double b_double;
  double c_double;

  float a_float;
  float b_float;
  float c_float;

  int a_int;
  int b_int;
  int c_int;

  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_PLUSPLUS:\n" );
  printf ( "  Demonstrate ++, which increments by 1:\n" );
  printf ( "\n" );
  printf ( "  We execute the following statements:\n" );
  printf ( "\n" );
  printf ( "  A = value;\n" );
  printf ( "  B = A;\n" );
  printf ( "  C = A++;\n" );
  printf ( "\n" );
  printf ( "  Type              A     B    C\n" );
  printf ( "\n" );

  a_char = 'a';
  b_char = a_char;
  c_char = a_char++;
  printf ( "  char              %c   %c  %c\n", a_char, b_char, c_char );

  a_double = 1.25;
  b_double = a_double;
  c_double = a_double++;
  printf ( "  double            %f    %f  %f\n", a_double, b_double, c_double );

  a_float = 1.25;
  b_float = a_float;
  c_float = a_float++;
  printf ( "  float             %f    %f  %f\n", a_float, b_float, c_float );

  a_int = 1;
  b_int = a_int;
  c_int = a_int++;
  printf ( "  int               %d    %d  %d\n", a_int, b_int, c_int );

  a_uchar = 'a';
  b_uchar = a_uchar;
  c_uchar = a_uchar++;
  printf ( "  unsigned char     %c    %c    %c\n", a_uchar, b_uchar, c_uchar );

  printf ( "\n" );
  printf ( "  We execute the following statements:\n" );
  printf ( "\n" );
  printf ( "  A = value;\n" );
  printf ( "  B = A;\n" );
  printf ( "  C = ++A;\n" );
  printf ( "\n" );
  printf ( "  Type              A     B    C\n" );
  printf ( "\n" );

  a_char = 'a';
  b_char = a_char;
  c_char = ++a_char;
  printf ( "  char              %c   %c  %c\n", a_char, b_char, c_char );

  a_double = 1.25;
  b_double = a_double;
  c_double = ++a_double;
  printf ( "  double            %f    %f  %f\n", a_double, b_double, c_double );

  a_float = 1.25;
  b_float = a_float;
  c_float = ++a_float;
  printf ( "  float             %f    %f  %f\n", a_float, b_float, c_float );

  a_int = 1;
  b_int = a_int;
  c_int = ++a_int;
  printf ( "  int               %d    %d  %d\n", a_int, b_int, c_int );

  a_uchar = 'a';
  b_uchar = a_uchar;
  c_uchar = ++a_uchar;
  printf ( "  unsigned char     %c    %c    %c\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_rshiftrshift ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_RSHIFTRSHIFT demonstrates >>.

  Discussion:

    >> implements right shift.

    A >> B has the value A / 2^B.
    
    A and B must be of an integer type.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  int a_int;
  int b_int;
  int c_int;
  int i;

  printf ( "\n" );
  printf ( "TEST_RSHIFTRSHIFT:\n" );
  printf ( "  Demonstrate >>, which implements the right shift:\n" );
  printf ( "  A >> B results in A / 2^B.\n" );
  printf ( "  Generally, B must be nonnegative.\n" );
  printf ( "  Moreover, B should not be so large that the result underflows.\n" );
  printf ( "\n" );
  printf ( "  When A is negative, the shift might be logical or arithmetic.\n" );
  printf ( "\n" );
  printf ( "  Type              A     B    A>>B\n" );
  printf ( "\n" );

  a_int = 3239;
  for ( i = -1; i < 32; i++ )
  {
    b_int = i;
    c_int = a_int >> b_int;

    printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );
  }

  printf ( "\n" );
  a_int = -3239;
  for ( i = -1; i < 32; i++ )
  {
    b_int = i;
    c_int = a_int >> b_int;

    printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );
  }

  return;
}
/******************************************************************************/

void test_slash ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_SLASH demonstrates /.

  Discussion:

    * implements division.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  char a_char;
  char b_char;
  char c_char;

  double a_double;
  double b_double;
  double c_double;

  float a_float;
  float b_float;
  float c_float;

  int a_int;
  int b_int;
  int c_int;

  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_SLASH:\n" );
  printf ( "  Demonstrate /, which carries out division:\n" );
  printf ( "\n" );
  printf ( "  Type              A     B     A/B\n" );
  printf ( "\n" );

  a_char = 'a';
  b_char = 'b';
  c_char = a_char / b_char;
  printf ( "  char              %c   %c   %d\n", a_char, b_char, c_char );

  a_char = 'a';
  b_char = 'a';
  c_char = a_char / b_char;
  printf ( "  char              %c   %c   %d\n", a_char, b_char, c_char );

  a_double = 20.00;
  b_double = 3.00;
  c_double = a_double / b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_float = 20.00;
  b_float = 3.00;
  c_float = a_float / b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_int = 20;
  b_int = 3;
  c_int = a_int / b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_uchar = 'a';
  b_uchar = 'a';
  c_uchar = a_uchar / b_uchar;
  printf ( "  unsigned char     %c    %c    %d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_star ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_STAR demonstrates *.

  Discussion:

    * implements multiplication.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2012

  Author:

    John Burkardt
*/
{
  char a_char;
  char b_char;
  char c_char;

  double a_double;
  double b_double;
  double c_double;

  float a_float;
  float b_float;
  float c_float;

  int a_int;
  int b_int;
  int c_int;

  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  printf ( "\n" );
  printf ( "TEST_STAR:\n" );
  printf ( "  Demonstrate *, which carries out multiplication:\n" );
  printf ( "\n" );
  printf ( "  Type              A     B     A*B\n" );
  printf ( "\n" );

  a_char = 'a';
  b_char = 'b';
  c_char = a_char * b_char;
  printf ( "  char              %c   %c   %d\n", a_char, b_char, c_char );

  a_double = 1.25;
  b_double = 3.12;
  c_double = a_double * b_double;
  printf ( "  double            %f    %f    %f\n", a_double, b_double, c_double );

  a_float = 1.25;
  b_float = 3.12;
  c_float = a_float * b_float;
  printf ( "  float             %f    %f    %f\n", a_float, b_float, c_float );

  a_int = 10;
  b_int = -7;
  c_int = a_int * b_int;
  printf ( "  int               %d    %d    %d\n", a_int, b_int, c_int );

  a_uchar = 'a';
  b_uchar = 'b';
  c_uchar = a_uchar * b_uchar;
  printf ( "  unsigned char     %c    %c    %d\n", a_uchar, b_uchar, c_uchar );

  return;
}
/******************************************************************************/

void test_twiddle ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_TWIDDLE demonstrates ~.

  Discussion:

    ~ implements the bitwise ones complement function.

    The operation of this function is easiest to see when we look at data
    that is unsigned characters, that is, integers between 0 and 255.
    In binary, some examples include

      00000000 =   0
      00000011 =   3
      00000101 =   5
      00001101 =  13
      01010000 =  80
      01010101 =  85
      10000000 = 128
      11110010 = 242
      11111111 = 255
      
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  unsigned char a_uchar;
  unsigned char b_uchar;
  unsigned char c_uchar;

  char a_char;
  char b_char;
  char c_char;

  printf ( "\n" );
  printf ( "TEST_TWIDDLE:\n" );
  printf ( "  Demonstrate ~, which implements the ones complement operator.\n" );
  printf ( "  This is most useful when the data is unsigned characters,\n" );
  printf ( "  that is, binary integers from 0 to 255.\n" );
  printf ( "\n" );
  printf ( "  Type                A     ~A   A+~A\n" );
  printf ( "\n" );

  a_uchar = 0;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 1;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 2;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 3;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 5;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 13;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 80;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 85;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 128;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 242;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 254;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  a_uchar = 255;
  b_uchar = ~ a_uchar;
  c_uchar = a_uchar + b_uchar;
  printf ( "  uchar             %3d    %3d   %3d\n", a_uchar, b_uchar, c_uchar );

  printf ( "\n" );

  a_char = 0;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 1;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 2;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 3;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 5;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 13;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 80;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 85;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = 127;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = -1;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = -2;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = -3;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

  a_char = -4;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );


  a_char = -128;
  b_char = ~ a_char;
  c_char = a_char + b_char;
  printf ( "  char             %4d   %4d  %4d\n", a_char, b_char, c_char );

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
