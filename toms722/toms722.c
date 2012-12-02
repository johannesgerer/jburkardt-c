#include <math.h>

#include "toms722.h"

double copysign_d ( double argVal, double argSign ) {
/*
Purpose:

  COPYSIGN_D returns the first argument with the sign of the second.

Discussion:

  The routine should work properly even if one or both arguments are NaN.

Author:

  W. J. Cody and J. T. Coonen
*/
  double z;

/* 
  Work with local copy of argVal. 
*/
  z = argVal;

  highpart(z) = (highpart(z) & ~SignMask) | (highpart(argSign) & SignMask);

  return z;
}
float copysign_f ( float argVal, float argSign ) {
/*
Purpose:

  COPYSIGN_F returns the first argument with the sign of the second.

Discussion:

  The routine should work properly even if one or both arguments are NaN.

Author:

  W. J. Cody and J. T. Coonen
*/
  float s;
  float v;
/* 
  Work with local copies of arguments. 
*/
  v = argVal;
  s = argSign;

  allof(v) = (allof(v) & ~SignMask) | (allof(s) & SignMask);

  return  v ;
}

int finite_d ( double arg) {
/*
Purpose:

  FINITE_D returns true (1) if its argument is finite, and false (0) otherwise.  
Discussion:

  The routine is independent of float.h.  NaN arguments are filtered with
  ISNAN_D to avoid spurious Invalid Operation exceptions.

Author:

  W. J. Cody and J. T. Coonen
*/
  if ( isnan_d ( arg ) )
    return 0;
  else
    return ( ( ABS ( arg ) < One ) || ( arg * Half != arg ) );
}

int finite_f ( float arg ) {
/*
Purpose:

  FINITE_F returns true (1) if its argument is finite, and false (0) otherwise.  
Discussion:

  The routine is independent of float.h.  NaN arguments are filtered with
  ISNAN_F to avoid spurious Invalid Operation exceptions.

Author:

  W. J. Cody and J. T. Coonen
*/
  if ( isnan_f ( arg ) )
    return 0;
  else
    return ( ABS ( arg ) < One ) || ( arg * Half != arg );
}

int isnan_d ( double arg) {
/*
Purpose:

  ISNAN_D returns true (1) if its argument is NaN, and false (0) otherwise.  

Discussion:

  The routine exploits the IEEE requirement that NaNs compare as unequal 
  to all values, including themselves.

Author:

  W. J. Cody and J. T. Coonen

*/
  return ( arg != arg );
}

int isnan_f ( float arg ) {
/*
Purpose:

  ISNAN_F returns true (1) if its argument is NaN, and false (0) otherwise.  

Discussion:

  The routine exploits the IEEE requirement that NaNs compare as unequal 
  to all values, including themselves.

Author:

  W. J. Cody and J. T. Coonen
*/
  return ( arg != arg );
}

double logb_d ( double arg ) {
 /*
Purpose:

  LOGB_D returns the exponent of a real number.

Discussion:

  For arg positive and finite, 1 <= | SCALB_D(arg,-N) | < 2, where
  N is the long equivalent of LOGB_D(arg).  Special cases
  are as follows:

  1) LOGB_D(NaN) = NaN,
  2) LOGB_D(infinity) = infinity, and
  3) LOGB_D(0) = -infinity  and signals division by zero.
  4) LOGB_D(denormal) = logb(0) on machines that flush to zero.

Author:

  W. J. Cody and J. T. Coonen
*/

        long expAdj;
        double z;

/* 
  Work with local copy of argument with positive sign.
*/
        z = copysign_d(arg,One);

/* 
  logb_d(NaN) is NaN.
*/
        if (isnan_d(z))
           return arg ;
/* 
  logb_d(Infinity) is +Infinity.
*/
        if (!finite_d(z))
           return z ;
/* 
  Renormalize, if necessary, a nonzero value and build
  exponent correction.  Denormals that arise in systems
  that flush underflows to zero may be forced to zero in
  this loop.
*/
        for (expAdj = 0; highpart(z) < DMinNorm && z != 0; expAdj--)
           z += z;

        /* logb_d(0) is -Infinity, with division by zero. */
        if (z == Zero)
           return copysign_d(One/z,-One);

        /* Grab exponent and return. */

        return (double) ((highpart(z) >> DExpShift) - DExpBias + expAdj);

}
float logb_f ( float arg ) {
 /*
Purpose:

  LOGB_F returns the exponent of a real number.

Discussion:

  For arg positive and finite, 1 <= | SCALB_F(arg,-N) | < 2, where
  N is the long equivalent of LOGB_F(arg).  Special cases
  are as follows:

  1) LOGB_F(NaN) = NaN,
  2) LOGB_F(infinity) = infinity, and
  3) LOGB_F(0) = -infinity  and signals division by zero.
  4) LOGB_F(denormal) = logb(0) on machines that flush to zero.

Author:

  W. J. Cody and J. T. Coonen
*/
        long expAdj;
        float z;

/*
  Work with local copy of argument with positive sign. 
*/
        z = copysign_f(arg,ROne);
/* 
  logb_f(NaN) is NaN.
*/
        if (isnan_f(z))
           return arg ;

/*
  logb_f(Infinity) is +Infinity.
*/
        if (!finite_f(z))
           return z ;

/* 
  Renormalize, if necessary, a nonzero value and build exponent
  correction.  Denormals that arise in systems that flush
  underflows to zero may be forced to zero in this loop.
*/
        for (expAdj = 0; allof(z) < RMinNorm && z != 0; expAdj--)
           z += z;

/*
  logb_f(0) is -Infinity, with division by zero.
*/
        if (z == Zero) return copysign_f(ROne/z,-ROne);

/* 
  Grab exponent and return.
*/

        return (float) ((allof(z) >> RExpShift) - RExpBias + expAdj);

}
void machar_d ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
  long int *machep, long int *negep, long int *iexp, long int *minexp,
  long int *maxexp, double *eps, double *epsneg, double *xmin, double *xmax ) {
/*
Purpose:

  MACHAR_D computes machine constants for double floating point arithmetic.

Discussion:

  This routine determines the parameters of the floating-point 
  arithmetic system specified below.  The determination of the first 
  three uses an extension of an algorithm due to Malcolm, 
  incorporating some of the improvements suggested by Gentleman and 
  Marovich.  

  A FORTRAN version of this routine appeared as ACM algorithm 665.

  This routine is a C translation of the FORTRAN code, and appeared
  as part of ACM algorithm 722.

  An earlier version of this program was published in Cody and Waite.

Reference:

  W J Cody,
  ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
    machine parameters,
  ACM Transactions on Mathematical Software,
  Volume 14, Number 4, pages 303-311, 1988.

  W J Cody and W Waite,
  Software Manual for the Elementary Functions,
  Prentice Hall, 1980.

  M Gentleman and S Marovich,
  Communications of the ACM,
  Volume 17, pages 276-277, 1974.

  M. Malcolm,
  Communications of the ACM,
  Volume 15, pages 949-951, 1972.

Author:

  W. J. Cody
  Argonne National Laboratory

Parameters:

  Output, long int * IBETA, the radix for the floating-point representation.

  Output, long int * IT, the number of base IBETA digits in the floating-point
  significand.

  Output, long int * IRND:
  0, if floating-point addition chops.
  1, if floating-point addition rounds, but not in the IEEE style.
  2, if floating-point addition rounds in the IEEE style.
  3, if floating-point addition chops, and there is partial underflow.
  4, if floating-point addition rounds, but not in the IEEE style, and 
    there is partial underflow.
  5, if floating-point addition rounds in the IEEE style, and there is 
    partial underflow.

  Output, long int * NGRD, the number of guard digits for multiplication with
  truncating arithmetic.  It is
  0, if floating-point arithmetic rounds, or if it truncates and only 
    IT base IBETA digits participate in the post-normalization shift of the
    floating-point significand in multiplication;
  1, if floating-point arithmetic truncates and more than IT base IBETA
    digits participate in the post-normalization shift of the floating-point
    significand in multiplication.

  Output, long int * MACHEP, the largest negative integer such that
    1.0 + ( double ) IBETA ^ MACHEP != 1.0, 
  except that MACHEP is bounded below by - ( IT + 3 ).

  Output, long int * NEGEPS, the largest negative integer such that
    1.0 - ( double ) IBETA ) ^ NEGEPS != 1.0, 
  except that NEGEPS is bounded below by - ( IT + 3 ).

  Output, long int * IEXP, the number of bits (decimal places if IBETA = 10)
  reserved for the representation of the exponent (including the bias or
  sign) of a floating-point number.

  Output, long int * MINEXP, the largest in magnitude negative integer such 
  that
    ( double ) IBETA ^ MINEXP 
  is positive and normalized.

  Output, long int * MAXEXP, the smallest positive power of BETA that overflows.
 
  Output, double * EPS, the smallest positive floating-point number such
  that  
    1.0 + EPS != 1.0. 
  in particular, if either IBETA = 2  or IRND = 0, 
    EPS = ( double ) IBETA ^ MACHEP.
  Otherwise,  
    EPS = ( ( double ) IBETA ^ MACHEP ) / 2.

  Output, double * EPSNEG, a small positive floating-point number such that
    1.0 - EPSNEG != 1.0. 
  In particular, if IBETA = 2 or IRND = 0, 
    EPSNEG = ( double ) IBETA ^ NEGEPS.
  Otherwise,  
    EPSNEG = ( double ) IBETA ^ NEGEPS ) / 2.  
  Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
  smallest number that can alter 1.0 by subtraction.

  Output, double * XMIN, the smallest non-vanishing normalized floating-point
  power of the radix:
    XMIN = ( double ) IBETA ^ MINEXP

  Output, float * XMAX, the largest finite floating-point number.  In
  particular,
    XMAX = ( 1.0 - EPSNEG ) * ( double ) IBETA ^ MAXEXP
  On some machines, the computed value of XMAX will be only the second, 
  or perhaps third, largest number, being too small by 1 or 2 units in 
  the last digit of the significand.
*/
  double a;
  double b;
  double beta;
  double betah;
  double betain;
  int i;
  int itmp;
  int iz;
  int j;
  int k;
  int mx;
  int nxres;
  double one;
  double t;
  double tmp;
  double tmp1;
  double tmpa;
  double two;
  double y;
  double z;
  double zero;

  (*irnd) = 1;
  one = (double) (*irnd);
  two = one + one;
  a = two;
  b = a;
  zero = 0.0e0;
/*
  Determine IBETA and BETA ala Malcolm.
*/
  tmp = ( ( a + one ) - a ) - one;

  while ( tmp == zero ) {
    a = a + a;
    tmp = a + one;
    tmp1 = tmp - a;
    tmp = tmp1 - one;
  }

  tmp = a + b;
  itmp = ( int ) ( tmp - a );

  while ( itmp == 0 ) {
    b = b + b;
    tmp = a + b;
    itmp = ( int ) ( tmp - a );
  }

  *ibeta = itmp;
  beta = ( double ) ( *ibeta );
/*
  Determine IRND, IT.
*/
  ( *it ) = 0;
  b = one;
  tmp = ( ( b + one ) - b ) - one;

  while ( tmp == zero ) {
    *it = *it + 1;
    b = b * beta;
    tmp = b + one;
    tmp1 = tmp - b;
    tmp = tmp1 - one;
  }

  *irnd = 0;
  betah = beta / two;
  tmp = a + betah;
  tmp1 = tmp - a;

  if ( tmp1 != zero ) {
    *irnd = 1;
  }

  tmpa = a + beta;
  tmp = tmpa + betah;

  if ( ( *irnd == 0 ) && ( tmp - tmpa != zero ) ) {
    *irnd = 2;
  }
/*
  Determine NEGEP, EPSNEG.
*/
  (*negep) = (*it) + 3;
  betain = one / beta;
  a = one;
 
  for ( i = 1; i <= (*negep); i++ ) {
    a = a * betain;
  }
 
  b = a;
  tmp = ( one - a );
  tmp = tmp - one;

  while ( tmp == zero ) {
    a = a * beta;
    *negep = *negep - 1;
    tmp1 = one - a;
    tmp = tmp1 - one;
  }

  (*negep) = -(*negep);
  (*epsneg) = a;
/*
  Determine MACHEP, EPS.
*/

  (*machep) = -(*it) - 3;
  a = b;
  tmp = one + a;

  while ( tmp - one == zero) {
    a = a * beta;
    *machep = *machep + 1;
    tmp = one + a;
  }

  *eps = a;
/*
  Determine NGRD.
*/
  (*ngrd) = 0;
  tmp = one + *eps;
  tmp = tmp * one;

  if ( ( (*irnd) == 0 ) && ( tmp - one ) != zero ) {
    (*ngrd) = 1;
  }
/*
  Determine IEXP, MINEXP and XMIN.

  Loop to determine largest I such that (1/BETA) ** (2**(I))
  does not underflow.  Exit from loop is signaled by an underflow.
*/

  i = 0;
  k = 1;
  z = betain;
  t = one + *eps;
  nxres = 0;

  for ( ; ; ) {
    y = z;
    z = y * y;
/*
  Check for underflow
*/

    a = z * one;
    tmp = z * t;

    if ( ( a + a == zero ) || ( ABS ( z ) > y ) ) {
      break;
    }

    tmp1 = tmp * betain;

    if ( tmp1 * beta == z ) {
      break;
    }

    i = i + 1;
    k = k + k;
  }
/*
  Determine K such that (1/BETA)**K does not underflow.
  First set  K = 2 ** I.
*/
  (*iexp) = i + 1;
  mx = k + k;

  if ( *ibeta == 10 ) {
/*
  For decimal machines only
*/

    (*iexp) = 2;
    iz = *ibeta;
    while ( k >= iz ) {
      iz = iz * ( *ibeta );
      (*iexp) = (*iexp) + 1;
    }
    mx = iz + iz - 1;
  }
 
/*
  Loop to determine MINEXP, XMIN.
  Exit from loop is signaled by an underflow.
*/
  for ( ; ; ) {
    (*xmin) = y;
    y = y * betain;
    a = y * one;
    tmp = y * t;
    tmp1 = a + a;

    if ( ( tmp1 == zero ) || ( ABS ( y ) >= ( *xmin ) ) ) {
      break;
    }

    k = k + 1;
    tmp1 = tmp * betain;
    tmp1 = tmp1 * beta;

    if ( ( tmp1 == y ) && ( tmp != y ) ) {
      nxres = 3;
      *xmin = y;
      break;
    }

  }

  (*minexp) = -k;

/*
  Determine MAXEXP, XMAX.
*/
  if ( ( mx <= k + k - 3 ) && ( ( *ibeta ) != 10 ) ) {
    mx = mx + mx;
    (*iexp) = (*iexp) + 1;
  }

  (*maxexp) = mx + (*minexp);
/*
  Adjust IRND to reflect partial underflow.
*/
  (*irnd) = (*irnd) + nxres;
/*
  Adjust for IEEE style machines.
*/
  if ( ( *irnd) >= 2 ) {
    (*maxexp) = (*maxexp) - 2;
  }
/*
  Adjust for machines with implicit leading bit in binary
  significand and machines with radix point at extreme
  right of significand.
*/
  i = (*maxexp) + (*minexp);

  if ( ( ( *ibeta ) == 2 ) && ( i == 0 ) ) {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( i > 20 ) {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( a != y ) {
    (*maxexp) = (*maxexp) - 2;
  }

  (*xmax) = one - (*epsneg);
  tmp = (*xmax) * one;

  if ( tmp != (*xmax) ) {
    (*xmax) = one - beta * (*epsneg);
  }

  (*xmax) = (*xmax) / ( beta * beta * beta * (*xmin) );
  i = (*maxexp) + (*minexp) + 3;

  if ( i > 0 ) {
 
    for ( j = 1; j <= i; j++ ) {
      if ( (*ibeta) == 2 ) {
        (*xmax) = (*xmax) + (*xmax);
      }
      if ( (*ibeta) != 2 ) {
        (*xmax) = (*xmax) * beta;
      }
    }

  }

  return;

}
void machar_s ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
  long int *machep, long int *negep, long int *iexp, long int *minexp,
  long int *maxexp, float *eps, float *epsneg, float *xmin, float *xmax ) {
/*
Purpose:

  MACHAR_S computes machine constants for floating point arithmetic.

Discussion:

  This routine determines the parameters of the floating-point 
  arithmetic system specified below.  The determination of the first 
  three uses an extension of an algorithm due to Malcolm, 
  incorporating some of the improvements suggested by Gentleman and 
  Marovich.  

  A FORTRAN version of this routine appeared as ACM algorithm 665.

  This routine is a C translation of the FORTRAN code, and appeared
  as part of ACM algorithm 722.

  An earlier version of this program was published in Cody and Waite.

Reference:

  W J Cody,
  ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
    machine parameters,
  ACM Transactions on Mathematical Software,
  Volume 14, Number 4, pages 303-311, 1988.

  W J Cody and W Waite,
  Software Manual for the Elementary Functions,
  Prentice Hall, 1980.

  M Gentleman and S Marovich,
  Communications of the ACM,
  Volume 17, pages 276-277, 1974.

  M. Malcolm,
  Communications of the ACM,
  Volume 15, pages 949-951, 1972.

Author:

  W. J. Cody
  Argonne National Laboratory

Parameters:

  Output, long int * IBETA, the radix for the floating-point representation.

  Output, long int * IT, the number of base IBETA digits in the floating-point
  significand.

  Output, long int * IRND:
  0, if floating-point addition chops.
  1, if floating-point addition rounds, but not in the IEEE style.
  2, if floating-point addition rounds in the IEEE style.
  3, if floating-point addition chops, and there is partial underflow.
  4, if floating-point addition rounds, but not in the IEEE style, and 
    there is partial underflow.
  5, if floating-point addition rounds in the IEEE style, and there is 
    partial underflow.

  Output, long int * NGRD, the number of guard digits for multiplication with
  truncating arithmetic.  It is
  0, if floating-point arithmetic rounds, or if it truncates and only 
    IT base IBETA digits participate in the post-normalization shift of the
    floating-point significand in multiplication;
  1, if floating-point arithmetic truncates and more than IT base IBETA
    digits participate in the post-normalization shift of the floating-point
    significand in multiplication.

  Output, long int * MACHEP, the largest negative integer such that
    1.0 + ( float ) IBETA ^ MACHEP != 1.0, 
  except that MACHEP is bounded below by - ( IT + 3 ).

  Output, long int * NEGEPS, the largest negative integer such that
    1.0 - ( float ) IBETA ) ^ NEGEPS != 1.0, 
  except that NEGEPS is bounded below by - ( IT + 3 ).

  Output, long int * IEXP, the number of bits (decimal places if IBETA = 10)
  reserved for the representation of the exponent (including the bias or
  sign) of a floating-point number.

  Output, long int * MINEXP, the largest in magnitude negative integer such 
  that
    ( float ) IBETA ^ MINEXP 
  is positive and normalized.

  Output, long int * MAXEXP, the smallest positive power of BETA that overflows.
 
  Output, float * EPS, the smallest positive floating-point number such
  that  
    1.0 + EPS != 1.0. 
  in particular, if either IBETA = 2  or IRND = 0, 
    EPS = ( float ) IBETA ^ MACHEP.
  Otherwise,  
    EPS = ( ( float ) IBETA ^ MACHEP ) / 2.

  Output, float * EPSNEG, a small positive floating-point number such that
    1.0 - EPSNEG != 1.0. 
  In particular, if IBETA = 2 or IRND = 0, 
    EPSNEG = ( float ) IBETA ^ NEGEPS.
  Otherwise,  
    EPSNEG = ( float ) IBETA ^ NEGEPS ) / 2.  
  Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
  smallest number that can alter 1.0 by subtraction.

  Output, float * XMIN, the smallest non-vanishing normalized floating-point
  power of the radix:
    XMIN = ( float ) IBETA ^ MINEXP

  Output, float * XMAX, the largest finite floating-point number.  In
  particular,
    XMAX = ( 1.0 - EPSNEG ) * ( float ) IBETA ^ MAXEXP
  On some machines, the computed value of XMAX will be only the second, 
  or perhaps third, largest number, being too small by 1 or 2 units in 
  the last digit of the significand.
*/
  float a;
  float b;
  float beta;
  float betah;
  float betain;
  int i;
  int itmp;
  int iz;
  int j;
  int k;
  int mx;
  int nxres;
  float one;
  float t;
  float tmp;
  float tmp1;
  float tmpa;
  float two;
  float y;
  float z;
  float zero;

  (*irnd) = 1;
  one = (float) (*irnd);
  two = one + one;
  a = two;
  b = a;
  zero = 0.0e0;
/*
  Determine IBETA and BETA ala Malcolm.
*/
  tmp = ( ( a + one ) - a ) - one;

  while ( tmp == zero ) {
    a = a + a;
    tmp = a + one;
    tmp1 = tmp - a;
    tmp = tmp1 - one;
  }

  tmp = a + b;
  itmp = ( int ) ( tmp - a );

  while ( itmp == 0 ) {
    b = b + b;
    tmp = a + b;
    itmp = ( int ) ( tmp - a );
  }

  *ibeta = itmp;
  beta = ( float ) ( *ibeta );
/*
  Determine IRND, IT.
*/
  ( *it ) = 0;
  b = one;
  tmp = ( ( b + one ) - b ) - one;

  while ( tmp == zero ) {
    *it = *it + 1;
    b = b * beta;
    tmp = b + one;
    tmp1 = tmp - b;
    tmp = tmp1 - one;
  }

  *irnd = 0;
  betah = beta / two;
  tmp = a + betah;
  tmp1 = tmp - a;

  if ( tmp1 != zero ) {
    *irnd = 1;
  }

  tmpa = a + beta;
  tmp = tmpa + betah;

  if ( ( *irnd == 0 ) && ( tmp - tmpa != zero ) ) {
    *irnd = 2;
  }
/*
  Determine NEGEP, EPSNEG.
*/
  (*negep) = (*it) + 3;
  betain = one / beta;
  a = one;
 
  for ( i = 1; i <= (*negep); i++ ) {
    a = a * betain;
  }
 
  b = a;
  tmp = ( one - a );
  tmp = tmp - one;

  while ( tmp == zero ) {
    a = a * beta;
    *negep = *negep - 1;
    tmp1 = one - a;
    tmp = tmp1 - one;
  }

  (*negep) = -(*negep);
  (*epsneg) = a;
/*
  Determine MACHEP, EPS.
*/

  (*machep) = -(*it) - 3;
  a = b;
  tmp = one + a;

  while ( tmp - one == zero) {
    a = a * beta;
    *machep = *machep + 1;
    tmp = one + a;
  }

  *eps = a;
/*
  Determine NGRD.
*/
  (*ngrd) = 0;
  tmp = one + *eps;
  tmp = tmp * one;

  if ( ( (*irnd) == 0 ) && ( tmp - one ) != zero ) {
    (*ngrd) = 1;
  }
/*
  Determine IEXP, MINEXP and XMIN.

  Loop to determine largest I such that (1/BETA) ** (2**(I))
  does not underflow.  Exit from loop is signaled by an underflow.
*/

  i = 0;
  k = 1;
  z = betain;
  t = one + *eps;
  nxres = 0;

  for ( ; ; ) {
    y = z;
    z = y * y;
/*
  Check for underflow
*/

    a = z * one;
    tmp = z * t;

    if ( ( a + a == zero ) || ( ABS ( z ) > y ) ) {
      break;
    }

    tmp1 = tmp * betain;

    if ( tmp1 * beta == z ) {
      break;
    }

    i = i + 1;
    k = k + k;
  }
/*
  Determine K such that (1/BETA)**K does not underflow.
  First set  K = 2 ** I.
*/
  (*iexp) = i + 1;
  mx = k + k;

  if ( *ibeta == 10 ) {
/*
  For decimal machines only
*/

    (*iexp) = 2;
    iz = *ibeta;
    while ( k >= iz ) {
      iz = iz * ( *ibeta );
      (*iexp) = (*iexp) + 1;
    }
    mx = iz + iz - 1;
  }
 
/*
  Loop to determine MINEXP, XMIN.
  Exit from loop is signaled by an underflow.
*/
  for ( ; ; ) {
    (*xmin) = y;
    y = y * betain;
    a = y * one;
    tmp = y * t;
    tmp1 = a + a;

    if ( ( tmp1 == zero ) || ( ABS ( y ) >= ( *xmin ) ) ) {
      break;
    }

    k = k + 1;
    tmp1 = tmp * betain;
    tmp1 = tmp1 * beta;

    if ( ( tmp1 == y ) && ( tmp != y ) ) {
      nxres = 3;
      *xmin = y;
      break;
    }

  }

  (*minexp) = -k;

/*
  Determine MAXEXP, XMAX.
*/
  if ( ( mx <= k + k - 3 ) && ( ( *ibeta ) != 10 ) ) {
    mx = mx + mx;
    (*iexp) = (*iexp) + 1;
  }

  (*maxexp) = mx + (*minexp);
/*
  Adjust IRND to reflect partial underflow.
*/
  (*irnd) = (*irnd) + nxres;
/*
  Adjust for IEEE style machines.
*/
  if ( ( *irnd) >= 2 ) {
    (*maxexp) = (*maxexp) - 2;
  }
/*
  Adjust for machines with implicit leading bit in binary
  significand and machines with radix point at extreme
  right of significand.
*/
  i = (*maxexp) + (*minexp);

  if ( ( ( *ibeta ) == 2 ) && ( i == 0 ) ) {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( i > 20 ) {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( a != y ) {
    (*maxexp) = (*maxexp) - 2;
  }

  (*xmax) = one - (*epsneg);
  tmp = (*xmax) * one;

  if ( tmp != (*xmax) ) {
    (*xmax) = one - beta * (*epsneg);
  }

  (*xmax) = (*xmax) / ( beta * beta * beta * (*xmin) );
  i = (*maxexp) + (*minexp) + 3;

  if ( i > 0 ) {
 
    for ( j = 1; j <= i; j++ ) {
      if ( (*ibeta) == 2 ) {
        (*xmax) = (*xmax) + (*xmax);
      }
      if ( (*ibeta) != 2 ) {
        (*xmax) = (*xmax) * beta;
      }
    }

  }

  return;
}
double nextafter_d ( double argx, double argy ) {
/*
Purpose:

  NEXTAFTER_D returns the next value in a given direction.

  In particular, NEXTAFTER_D(argx,argy) returns the next representable 
  neighbor to argx in the direction toward argy, where the function and
  arguments are double-precision.  The following special cases arise:

  1) if argx != argx or argy != argy, at least one of the arguments
     is a NaN, and one of the input NaNs is returned;
  2) if argx = argy, argx is returned; and
  3) if the result is zero, underflow is signalled.

Author:

  W. J. Cody and J. T. Coonen
*/
  double z, forceException;

/*
  If one argument is NaN, return NaN.
*/
       if (isnan_d(argx) || isnan_d(argy)) return argx + argy;

/*
  If arguments are equal, return argx. 
*/
       if (argx == argy) return argx;

/*
  Work with local copy of argx. 
*/
       z = argx;

/* 
  Initialize special variable used to stimulate side-effects despite 
  a compiler's desire to strip "dead"code.
*/
        forceException = Zero;

/* 
  Return signed smallest non-zero value when argx = zero.
*/
       if (argx == Zero) {
             highpart(z) = HighDMinx;
             lowpart(z) = LowDMinx;
             z = copysign_d(z,argy);
             }
/*
  Otherwise, use integer arithmetic to increment or
  decrement least significant half of z, being careful
  with carries and borrows involving most significant
  half.
*/
          else if (((argx < Zero) && (argx < argy)) ||
                   ((argx > Zero) && (argx > argy))) {
                   --lowpart(z);
                   if (lowpart(z) == -1)
                      --highpart(z);
                   }
                else {
                   ++lowpart(z);
                   if (lowpart(z) == 0)
                      ++highpart(z);
                   }

/*
  Now trigger the underflow signal (with 0 result for
  machines that flush) if the result is denormal, and
  trigger overflow if the result is +/-INF.  If z is
  denormal on a flushing machine, argx must have been
  the tiniest normal number, so just halve argx to force
  the appropriate underflow to zero.  Denormal (or zero)
  z on a graceful underflow machine is correct, but
  requires that the underflow signal be triggered.
*/

       if ((highpart(z) & ~MDMaxDenorm) == DZero)
#if defined(BEF) || defined(LEF)
          z = argx * Half;  
/* 
  Force underflow to 0 with signal. 
*/
#else
/* 
  z is either denormal or zero.  argx is denormal, zero
  or the teeniest normal.  So   z * z + argx * argx
  is guaranteed to underflow.  The subsequent test
  forces the result to a double variable, despite the
  optimizing whims of compilers.
*/
          forceException = z * z + argx * argx;
#endif
/*
  If z is infinite, it's correct, but overflow must be signaled.
 */
       else if (!finite_d(z))
          forceException = argx * argx;

/*
  Now use forceException in a seemingly nontrivial
  expression to be sure the side-effects above are
  generated.  In fact, the test below always fails.
*/
      if ((highpart(forceException) & SignMask) == SignMask)
             z = Zero;
    return  z;
}
float nextafter_f ( float argx, float argy ) { 
/*
Purpose:

  NEXTAFTER_F returns the next value in a given direction.

  In particular, NEXTAFTER_F(argx,argy) returns the next representable 
  neighbor to argx in the direction toward argy, where the function and
  arguments are double-precision.  The following special cases arise:

  1) if argx != argx or argy != argy, at least one of the arguments
     is a NaN, and one of the input NaNs is returned;
  2) if argx = argy, argx is returned; and
  3) if the result is zero, underflow is signalled.

Author:

  W. J. Cody and J. T. Coonen
*/
       float z, forceException;

/* 
  If one argument is NaN, return NaN.
*/
       if (isnan_f(argx) || isnan_f(argy)) return argx + argy ;

/* 
  If arguments are equal, return argx.
*/
       if (argx == argy) return argx ;

/* 
  Work with local copy of argx.
*/
       z = argx;

/* 
  Initialize special variable used to stimulate side-effects despite
  a compiler's desire to strip "dead" code.
*/
           forceException = Zero;

/* 
  Return signed smallest non-zero value when argx = zero.
*/
       if (argx == Zero) {
             allof(z) = RMinx;
             z = copysign_f(z,argy);
             }
/* 
  Otherwise, use integer arithmetic to increment or decrement z.
*/
          else if (argx < Zero) {
             if (argx < argy) 
                   --allof(z);
                else
                   ++allof(z);
             }
          else if (argx < argy)
                   ++allof(z);
                else
                   --allof(z);

/* 
  Now trigger the underflow signal (with 0 result for
  machines that flush) if the result is denormal, and
  trigger overflow if the result is +/-INF.  If z is
  denormal on a flushing machine, argx must have been
  the tiniest normal number, so just halve argx to force
  the appropriate underflow to zero.  Denormal (or zero)
  z on a graceful underflow machine is correct, but
  requires that the underflow signal be triggered.
*/

       if ((allof(z) & ~MRMaxDenorm) == Zero)
#if defined(BEF) || defined(LEF)
          z = argx * RHalf; 
/* 
  Force underflow to 0 with signal.
*/

#else

/* 
  z is either denormal or zero.  argx is denormal, zero
  or the teeniest normal.  So   z * z + argx * argx
  is guaranteed to underflow.  The subsequent test
  forces the result to a float variable, despite the
  optimizing whims of compilers.
*/
          forceException = z * z + argx * argx;
#endif
/* 
  If z is infinite, it's correct, but overflow must be signaled.
*/
       else if (!finite_f(z))
          forceException = argx * argx;

/* 
  Now use forceException in a seemingly nontrivial
  expression to be sure the side-effects above are
  are generated.  In fact, the test below always fails.
*/
          if ((allof(forceException) & SignMask) == SignMask)
                 z = Zero;  
       return  z ;
}
double scalb_d ( double arg, long n ) {
 /*
Purpose:

  SCALB_D returns ARG * 2 ** N.

Discussion:

  The result is computed by adding N to the floating-point exponent of ARG.
  
  This function may generate overflow or underflow.

Author:

  W. J. Cody and J. T. Coonen
*/
  long i;
  long tempExp;
  double z;
  double logbz;
  double scaling;

/* 
  Work with local copy of argument.
*/
  z = arg;

/* 
  Handle special cases of NaN, inf, 0.
*/
  if (isnan_d(z) || !finite_d(z) || z == Zero) return arg;
/*
  Extract exponent using logb, then normalize z, if
  necessary, by doubling.  Watch for infinite return
  from denormals in flushing systems.
*/
  logbz = logb_d(z);
  if (!finite_d(logbz)) return copysign_d(Zero, arg);

  tempExp = logbz;

  for (i=0; i<DMinNormExp-tempExp; i++)
    z += z;

/* 
  Zero out old exponent, and build new one.
*/
  highpart(z) &= ~DExpMask;
  tempExp += DExpBias + n;
/* 
  If new exponent too small, force minimum normalized
  exponent, and denormalize result by multiplying by the
  appropriate power of 1/2, possibly triggering underflow.
  Pin very negative exponents to -60, since that is more
  than sufficient to force denormalization "completely off
  the right-hand side."
*/
  if (tempExp < DMinNormBiasExp) {
    highpart(z) |= DMinNormBiasExp << DExpShift;
    if (tempExp < -60)
      tempExp = -60;
    for (i=tempExp,scaling=One;i<1;i++)
      scaling *= Half;
    z *= scaling;
  } else if (tempExp > DMaxNormBiasExp) {
/* 
  If new exponent too large, trigger overflow. 
*/
    highpart(z) |= DMaxNormBiasExp << DExpShift;
    z += z;  
/* 
  Oveflow according to rounding mode. 
*/
  } else
/* 
  Otherwise, new exponent is okay.  Plant it and continue. 
*/
  highpart(z) |= tempExp << DExpShift;

  return z;
}
float scalb_f ( float arg, long n ) {
 /*
Purpose:

  SCALB_F returns ARG * 2 ** N.

Discussion:

  The result is computed by adding N to the floating-point exponent of ARG.
  
  This function may generate overflow or underflow.

Author:

  W. J. Cody and J. T. Coonen
*/
       long      i, tempExp;
       float     z, logbz, scaling;

/* 
  Work with local copy of argument. 
*/
        z = arg;

/* 
  Handle special cases of NaN, inf, 0. 
*/
        if (isnan_f(z) || !finite_f(z) || z == Zero) return arg;

/* 
  Extract exponent using logb_f, then normalize z, if
  necessary, by doubling.  Watch for infinite return
  from denormals in flushing systems.
*/
        logbz = logb_f(z);
        if (!finite_f(logbz))
           return copysign_f(Zero, arg);

        tempExp = logbz;
        for (i=0; i<RMinNormExp-tempExp; i++)
           z += z;

/* 
  Zero out old exponent, and build new one. 
*/
        allof(z) &= ~RExpMask;
        tempExp += RExpBias + n;

/* 
  If new exponent too small, force minimum normalized
  exponent, and denormalize result by multiplying by the
  appropriate power of 1/2, possibly triggering underflow.
  Pin very negative exponents to -30, since that is more
  than sufficient to force denormalization "completely off
  the right-hand side.
*/
        if (tempExp < RMinNormBiasExp) {
              allof(z) |= RMinNormBiasExp << RExpShift;
              if (tempExp < -30)
                    tempExp = -30;
              for (i=tempExp,scaling=ROne;i<1;i++)
                    scaling *= RHalf;
              z *= scaling;
           }
/* 
  If new exponent too large, trigger overflow.
*/
        else if (tempExp > RMaxNormBiasExp) {
              allof(z) |= RMaxNormBiasExp << RExpShift;
              z += z;    /* Oveflow according to rounding mode. */
           }

/* 
  Otherwise, new exponent is okay.  Plant it and continue.
*/
        else
              allof(z) |= tempExp << RExpShift;

        return z;
}
double test_d ( double argx ) {
/*
  Purpose:

    TEST_D returns its input argument.
*/
  return ( argx );
}

float test_f ( float argx ) {
/*
  Purpose:

    TEST_F returns its input argument.
*/
  return ( argx );
}

