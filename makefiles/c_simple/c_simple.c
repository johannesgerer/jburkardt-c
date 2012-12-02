# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "c_simple.h"

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for the C_SIMPLE example.

  Modified:

    03 December 2006

  Author:

    John Burkardt

  Parameters:

    Local, float A, B, the endpoints of the interval of integration.

    Local, float F ( float T ), the name of the function to be integrated.

    Local, int INT_NUM, the number of intervals to be used.

    Local, float QUAD, the approximate value of the integral.
*/
{
  float a;
  float b;
  int int_num;
  float quad;
  int test;

  printf ( "\n" );
  printf ( "C_SIMPLE\n" );
  printf ( "  A simple C program to demonstrate\n" );
  printf ( "  the use of makefiles.\n" );

  printf ( "\n" );
  printf ( "  Estimate the integral from 0 to 1000, of\n" );
  printf ( "  F(T) = (4+T/365+1/2 sin(pi*T/91)) * (2+exp(-sin(2*pi*T)))\n" );
  printf ( "  a function which models daily power consumption.\n" );
  printf ( "\n" );
  printf ( "  quad = midpoint ( a, b, f, int_num )\n" );
  printf ( "  estimates the integral using the midpoint rule.\n" );
  printf ( "\n" );
  printf ( "  f ( t )\n" );
  printf ( "  evaluates the integrand.\n" );
  printf ( "\n" );
  printf ( "  Intervals   Estimate\n" );
  printf ( "\n" );

  a = 0.0;
  b = 1000.0;
  int_num = 100;

  for ( test = 1; test <= 3; test++ )
  {
    quad = midpoint ( a, b, &f, int_num );

    printf ( "  %8d  %14e\n", int_num, quad );

    int_num = int_num * 100;
  }
 
  printf ( "\n" );
  printf ( "C_SIMPLE:\n" );
  printf ( "  Normal end of execution.\n" );
 
  return 0;
}
