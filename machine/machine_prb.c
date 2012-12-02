# include <stdlib.h>
# include <stdio.h>

# include "machine.h"

int main ( void );
void d1mach_prb ( void );
void i1mach_prb ( void );
void r1mach_prb ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MACHINE_PRB.

  Discussion:

    MACHINE_PRB runs the MACHINE tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 April 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "MACHINE_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MACHINE library.\n" );

  d1mach_prb ( );
  i1mach_prb ( );
  r1mach_prb ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MACHINE_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void d1mach_prb ( void )

/******************************************************************************/
/*
  Purpose:

    D1MACH_PRB reports the constants returned by D1MACH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 April 2007

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "D1MACH_PRB\n" );
  printf ( "  D1MACH reports the value of constants associated\n" );
  printf ( "  with real double precision computer arithmetic.\n" );

  printf ( "\n" );
  printf ( "  Assume that double precision numbers are stored\n" );
  printf ( "  with a mantissa of T digits in base B, with an\n" );
  printf ( "  exponent whose value must lie between EMIN and EMAX.\n" );

  printf ( "\n" );
  printf ( "  For input arguments of 1 <= I <= 5,\n" );
  printf ( "  D1MACH will return the following values:\n" );

  printf ( "\n" );
  printf ( "  D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.\n" );
  printf ( "%26.16e\n", d1mach(1) );

  printf ( "\n" );
  printf ( "  D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.\n" );
  printf ( "%26.16e\n", d1mach(2) );

  printf ( "\n" );
  printf ( "  D1MACH(3) = B^(-T), the smallest relative spacing.\n" );
  printf ( "%26.16e\n", d1mach(3) );

  printf ( "\n" );
  printf ( "  D1MACH(4) = B^(1-T), the largest relative spacing.\n" );
  printf ( "%26.16e\n", d1mach(4) );

  printf ( "\n" );
  printf ( "  D1MACH(5) = log10(B).\n" );
  printf ( "%26.16e\n", d1mach(5) );

  return;
}
/******************************************************************************/

void i1mach_prb ( void )

/******************************************************************************/
/*
  Purpose:

    I1MACH_PRB reports the constants returned by I1MACH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 April 2007

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "I1MACH_PRB\n" );
  printf ( "  I1MACH reports the value of constants associated\n" );
  printf ( "  with integer computer arithmetic.\n" );

  printf ( "\n" );
  printf ( "  Numbers associated with input/output units:\n" );

  printf ( "\n" );
  printf ( "  I1MACH(1) = the standard input unit.\n" );
  printf ( "%d\n", i1mach(1) );

  printf ( "\n" );
  printf ( "  I1MACH(2) = the standard output unit.\n" );
  printf ( "%d\n", i1mach(2) );

  printf ( "\n" );
  printf ( "  I1MACH(3) = the standard punch unit.\n" );
  printf ( "%d\n", i1mach(3) );

  printf ( "\n" );
  printf ( "  I1MACH(4) = the standard error message unit.\n" );
  printf ( "%d\n", i1mach(4) );

  printf ( "\n" );
  printf ( "  Numbers associated with words:\n" );

  printf ( "\n" );
  printf ( "  I1MACH(5) = the number of bits per integer.\n" );
  printf ( "%d\n", i1mach(5) );

  printf ( "\n" );
  printf ( "  I1MACH(6) = the number of characters per integer.\n" );
  printf ( "%d\n", i1mach(6) );

  printf ( "\n" );
  printf ( "  Numbers associated with integer values:\n" );

  printf ( "\n" );
  printf ( "  Assume integers are represented in the S digit \n" );
  printf ( "  base A form:\n" );
  printf ( "\n" );
  printf ( "    Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))\n" );
  printf ( "\n" );
  printf ( "  where the digits X satisfy 0 <= X(1:S-1) < A.\n" );

  printf ( "\n" );
  printf ( "  I1MACH(7) = A, the base.\n" );
  printf ( "%d\n", i1mach(7) );

  printf ( "\n" );
  printf ( "  I1MACH(8) = S, the number of base A digits.\n" );
  printf ( "%d\n", i1mach(8) );

  printf ( "\n" );
  printf ( "  I1MACH(9) = A^S-1, the largest integer.\n" );
  printf ( "%d\n", i1mach(9) );

  printf ( "\n" );
  printf ( "  Numbers associated with floating point values:\n" );
  printf ( "\n" );
  printf ( "  Assume floating point numbers are represented \n" );
  printf ( "  in the T digit base B form:\n" );
  printf ( "\n" );
  printf ( "    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B^T) )\n" );
  printf ( "\n" );
  printf ( "  where\n" );
  printf ( "\n" );
  printf ( "    0 <= X(1:T) < B,\n" );
  printf ( "    0 < X(1) (unless the value being represented is 0),\n" );
  printf ( "    EMIN <= E <= EMAX.\n" );

  printf ( "\n" );
  printf ( "  I1MACH(10) = B, the base.\n" );
  printf ( "%d\n", i1mach(10) );

  printf ( "\n" );
  printf ( "  Numbers associated with single precision values:\n" );
  printf ( "\n" );
  printf ( "  I1MACH(11) = T, the number of base B digits.\n" );
  printf ( "%d\n", i1mach(11) );

  printf ( "\n" );
  printf ( "  I1MACH(12) = EMIN, the smallest exponent E.\n" );
  printf ( "%d\n", i1mach(12) );

  printf ( "\n" );
  printf ( "  I1MACH(13) = EMAX, the largest exponent E.\n" );
  printf ( "%d\n", i1mach(13) );

  printf ( "\n" );
  printf ( "  Numbers associated with double precision values:\n" );
  printf ( "\n" );
  printf ( "  I1MACH(14) = T, the number of base B digits.\n" );
  printf ( "%d\n", i1mach(14) );

  printf ( "\n" );
  printf ( "  I1MACH(15) = EMIN, the smallest exponent E.\n" );
  printf ( "%d\n", i1mach(15) );

  printf ( "\n" );
  printf ( "  I1MACH(16) = EMAX, the largest exponent E.\n" );
  printf ( "%d\n", i1mach(16) );

  return;
}
/******************************************************************************/

void r1mach_prb ( void )

/******************************************************************************/
/*
  Purpose:

    R1MACH_PRB reports the constants returned by R1MACH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 April 2007

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "R1MACH_PRB\n" );
  printf ( "  R1MACH reports the value of constants associated\n" );
  printf ( "  with real single precision computer arithmetic.\n" );

  printf ( "\n" );
  printf ( "  Assume that single precision numbers are stored \n" );
  printf ( "  with a mantissa of T digits in base B, with an \n" );
  printf ( "  exponent whose value must lie between EMIN and EMAX.\n" );

  printf ( "\n" );
  printf ( "  For input arguments of 1 <= I <= 5,\n" );
  printf ( "  R1MACH will return the following values:\n" );

  printf ( "\n" );
  printf ( "  R1MACH(1) = B^(EMIN-1), the smallest positive magnitude.\n" );
  printf ( "%26.16e\n", r1mach(1) );

  printf ( "\n" );
  printf ( "  R1MACH(2) = B^EMAX*(1-B**(-T)), the largest magnitude.\n" );
  printf ( "%26.16e\n", r1mach(2) );

  printf ( "\n" );
  printf ( "  R1MACH(3) = B^(-T), the smallest relative spacing.\n" );
  printf ( "%26.16e\n", r1mach(3) );

  printf ( "\n" );
  printf ( "  R1MACH(4) = B^(1-T), the largest relative spacing.\n" );
  printf ( "%26.16e\n", r1mach(4) );

  printf ( "\n" );
  printf ( "  R1MACH(5) = log10(B).\n" );
  printf ( "%26.16e\n", r1mach(5) );

  return;
}
