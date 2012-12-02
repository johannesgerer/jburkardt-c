/*fisher.c - compute the two-tailed probability of correct rejection of the null
  hypothesis with an F-ratio of x, for m degrees of freedom in the numerator and
  n degrees of freedom in the denominator.  In the special case of only two
  populations, this is equivalent to Student's t-test with m=1 and x=t**2.
  Coded by Matthew Belmonte <mkb4@Cornell.edu>, 28 September 1995.  This
  implementation Copyright (c) 1995 by Matthew Belmonte.  Permission for use and
  distribution is hereby granted, subject to the restrictions that this
  copyright notice and reference list be included in its entirety, and that any
  and all changes made to the program be clearly noted in the program text.

  References:

  Egon Dorrer, "Algorithm 322: F-Distribution [S14]", Communications of the
  Association for Computing Machinery 11:2:116-117 (1968).

  J.B.F. Field, "Certification of Algorithm 322 [S14] F-Distribution",
  Communications of the Association for Computing Machinery 12:1:39 (1969).

  Hubert Tolman, "Remark on Algorithm 322 [S14] F-Distribution", Communications
  of the Association for Computing Machinery 14:2:117 (1971).
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms322.h"

/******************************************************************************/

double fisher ( int m, int n, double x )

/******************************************************************************/
{
  int a;
  int b;
  double d;
  int i;
  int j;
  double p;
  double w;
  double y;
  double z;
  double zk;

  a = 2 * ( m / 2 ) - m + 2;
  b = 2 * ( n / 2 ) - n + 2;
  w = ( x * m ) / n;
  z = 1.0 / ( 1.0 + w );

  if ( a == 1) 
  {
    if ( b == 1 )
    {
      p = sqrt ( w );
      y = 0.3183098862;
      d = y * z / p;
      p = 2.0 * y * atan ( p );
    }
    else
    {
      p = sqrt ( w * z );
      d = 0.5 * p * z / w;
    }
  }
  else if ( b == 1 )
  {
    p = sqrt ( z );
    d = 0.5 * z * p;
    p = 1.0 - p;
  }
  else
  {
    d = z * z;
    p = w * z;
  }
  y = 2.0 * w / z;
  if ( a == 1 )
  { 
    for(j = b + 2; j <= n; j = j + 2 )
    {
      d = d * ( 1.0 + 1.0 / ( j - 2 ) ) * z;
      p = p + d * y / ( j - 1 );
    }
  }
  else
  {
    zk = pow ( z, ( double ) ( ( n - 1 ) / 2 ) );
    d = d * ( zk * n ) / b;
    p = p * zk + w * z * ( zk - 1.0 ) / ( z - 1.0 );
  }
  y = w * z;
  z = 2.0 / z;
  b = n - 2;

  for ( i = a + 2; i <= m; i = i + 2 )
  {
    j = i + b;
    d = d * ( y * j ) / ( i - 2 );
    p = p - z * d / j;
  }
  if ( p < 0.0 )
  {
    p = 0.0;
  }
  if ( 1.0 < p ) 
  {
    p = 1.0;
  }
  return p;
}
/******************************************************************************/

double student ( int df, double t )

/******************************************************************************/
{
  return ( fisher ( 1, df, t*t ) );
}
