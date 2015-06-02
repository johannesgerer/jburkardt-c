# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "cdflib.h"

int main ( );
void test005 ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CDFLIB_PRB.

  Discussion:

    CDFLIB_PRB tests the CDFLIB library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CDFLIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CDFLIB library.\n" );

  test005 ( );
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CDFLIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests BETA_INC and BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int ierror;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  BETA_INC computes the incomplete Beta ratio.\n" );
  printf ( "  BETA_INC_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         Y         A         B         CDF           CDF\n" );
  printf (
    "                                           (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    y = 1.0 - x;

    beta_inc ( &a, &b, &x, &y, &cdf_compute, &ccdf_compute, &ierror );

    printf ( "  %10g  %10g  %10g  %10g  %14g  %14g\n",
      x, y, a, b, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X         Y         A         B         1-CDF         CCDF\n" );
  printf (
    "                                           (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    y = 1.0 - x;

    beta_inc ( &a, &b, &x, &y, &cdf_compute, &ccdf_compute, &ierror );

    printf ( "  %10g  %10g  %10g  %10g  %14g  %14g\n",
      x, y, a, b, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests CDFBET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double bound;
  double p;
  double q;
  int status;
  int which;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CDFBET computes one missing parameter from the\n" );
  printf ( "  BETA CDF:\n" );
  printf ( "\n" );
  printf ( "   BETA_CDF ( (P,Q), (X,Y), A, B )\n" );
  printf ( "\n" );
  printf ( "      P           Q               X           Y" );
  printf ( "            A           B\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 0.25;
      y = 1.0 - x;
      a = 2.0;
      b = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = -1.0;
      y = -1.0;
      a = 2.0;
      b = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = 0.25;
      y = 1.0 - x;
      a = -1.0;
      b = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = 0.25;
      y = 1.0 - x;
      a = 2.0;
      b = -1.0;
    }

    cdfbet ( &which, &p, &q, &x, &y, &a, &b, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFBET returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      p, q, x, y, a, b );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CDFBIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double ompr;
  double p;
  double pr;
  double q;
  double s;
  int status;
  int which;
  double xn;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  CDFBIN computes one missing parameter from the\n" );
  printf ( "  Binomial CDF:\n" );
  printf ( "\n" );
  printf ( "   BINOMIAL_CDF ( (P,Q), S, XN, (PR,OMPR) )\n" );
  printf ( "\n" );
  printf ( "      P           Q                S          " );
  printf ( "XN         PR         OMPR\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      s = 5.0;
      xn = 8.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 2 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = -1.0;
      xn = 8.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 3 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = 5.0;
      xn = -1.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 4 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = 5.0;
      xn = 8.0;
      pr = -1.0;
      ompr = -1.0;
    }

    cdfbin ( &which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFBIN returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      p, q, s, xn, pr, ompr );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CDFCHI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double df;
  double p;
  double q;
  int status;
  int which;
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  CDFCHI computes one missing parameter from the\n" );
  printf ( "  Chi Square CDF:\n" );
  printf ( "\n" );
  printf ( "   CHI_CDF ( (P,Q), X, DF )\n" );
  printf ( "\n" );
  printf ( "      P           Q                X          DF\n" );
  printf ( "\n" );

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      df = 8.0;
    }
    else if ( which == 2 )
    {
      p = 0.242424;
      q = 1.0 - p;
      x = -1.0;
      df = 8.0;
    }
    else if ( which == 3 )
    {
      p = 0.242424;
      q = 1.0 - p;
      x = 5.0;
      df = -1.0;
    }

    cdfchi ( &which, &p, &q, &x, &df, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFCHI returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g\n", p, q, x, df );
  }
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CDFCHN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double df;
  double p;
  double pnonc;
  double q;
  int status;
  int which;
  double x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CDFCHN computes one missing parameter from the\n" );
  printf ( "  Chi Square CDF:\n" );
  printf ( "\n" );
  printf ( "   CHI_Noncentral_CDF ( (P,Q), X, DF, PNONC )\n" );
  printf ( "\n" );
  printf ( "     P         Q             X        DF     PNONC\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      df = 8.0;
      pnonc = 0.5;
    }
    else if ( which == 2 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = -1.0;
      df = 8.0;
      pnonc = 0.5;
    }
    else if ( which == 3 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = 5.0;
      df = -1.0;
      pnonc = 0.5;
    }
    else if ( which == 4 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = 5.0;
      df = 8.0;
      pnonc = -1.0;
    }

    cdfchn ( &which, &p, &q, &x, &df, &pnonc, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFCHN returned STATUS = %d\n", status );
      continue;
    }

    printf ( "  %10g  %10g  %10g  %10g  %10g\n",
      p, q, x, df, pnonc );
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CDFF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double dfd;
  double dfn;
  double f;
  double p;
  double q;
  int status;
  int which;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CDFF computes one missing parameter from the\n" );
  printf ( "  F CDF:\n" );
  printf ( "\n" );
  printf ( "   F_CDF ( (P,Q), F, DFN, DFD )\n" );
  printf ( "\n" );
  printf ( "     P         Q             F       DFN       DFD\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = -1.0;
      dfn = 8.0;
      dfd = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = 5.0;
      dfn = -1.0;
      dfd = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = -1.0;
    }

    cdff ( &which, &p, &q, &f, &dfn, &dfd, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFF returned STATUS = %d\n", status );
      continue;
    }

    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      p, q, f, dfn, dfd );
  }
  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests CDFFNC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double dfd;
  double dfn;
  double f;
  double p;
  double pnonc;
  double q;
  int status;
  int which;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  CDFFNC computes one missing parameter from the\n" );
  printf ( "  noncentral F CDF:\n" );
  printf ( "\n" );
  printf ( "   F_noncentral_CDF ( (P,Q), F, DFN, DFD, PNONC )\n" );
  printf ( "\n" );
  printf ( "         P         Q         F       DFN       DFD     PNONC\n" );
  printf ( "\n" );

  for ( which = 1; which <= 5; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 2 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = -1.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 3 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = -1.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 4 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = -1.0;
      pnonc = 17.648016;
    }
    else if ( which == 5 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = -1.0;
    }

    cdffnc ( &which, &p, &q, &f, &dfn, &dfd, &pnonc, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFFNC returned STATUS = %d\n", status );
      continue;
    }

    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      p, q, f, dfn, dfd, pnonc );
  }

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests CDFGAM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double p;
  double q;
  double scale;
  double shape;
  int status;
  int which;
  double x;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  CDFGAM computes one missing parameter from the\n" );
  printf ( "  Gamma CDF:\n" );
  printf ( "\n" );
  printf ( "   Gamma_CDF ( (P,Q), X, SHAPE, SCALE )\n" );
  printf ( "\n" );
  printf ( "    P         Q              X     SHAPE     SCALE\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      shape = 8.0;
      scale = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = -1.0;
      shape = 8.0;
      scale = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = 5.0;
      shape = -1.0;
      scale = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = 5.0;
      shape = 8.0;
      scale = -1.0;
    }

    cdfgam ( &which, &p, &q, &x, &shape, &scale, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFGAM returned STATUS = %d\n", status );
      continue;
    }

    printf ( "  %10g  %10g  %10g  %10g  %10g\n",
      p, q, x, shape, scale );
  }

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests CDFNBN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double f;
  double ompr;
  double p;
  double pr;
  double q;
  double s;
  int status;
  int which;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  CDFNBN computes one missing parameter from the\n" );
  printf ( "  Negative_Binomial CDF:\n" );
  printf ( "\n" );
  printf ( "   Negative_BINOMIAL_CDF ( (P,Q), F, S, (PR,OMPR) )\n" );
  printf ( "\n" );
  printf ( "    P         Q               F         S       PR        OMPR\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 3.0;
      s = 5.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 2 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = -1.0;
      s = 5.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 3 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = 3.0;
      s = -1.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 4 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = 3.0;
      s = 5.0;
      pr = -1.0;
      ompr = -1.0;
    }

    cdfnbn ( &which, &p, &q, &f, &s, &pr, &ompr, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFNBN returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      p, q, f, s, pr, ompr );
  }

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests CDFNOR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double mean;
  double p;
  double q;
  double sd;
  int status;
  int which;
  double x;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  CDFNOR computes one missing parameter from the\n" );
  printf ( "  Normal CDF:\n" );
  printf ( "\n" );
  printf ( "   Normal_CDF ( (P,Q), X, MEAN, SD )\n" );
  printf ( "\n" );
  printf ( "    P         Q               X      MEAN       SD\n" );
  printf ( "\n" );

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 3.0;
      mean = 5.0;
      sd = 0.875;
    }
    else if ( which == 2 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = -1.0;
      mean = 5.0;
      sd = 0.875;
    }
    else if ( which == 3 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = 3.0;
      mean = -1.0;
      sd = 0.875;
    }
    else if ( which == 4 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = 3.0;
      mean = 5.0;
      sd = -1.0;
    }

    cdfnor ( &which, &p, &q, &x, &mean, &sd, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFNOR returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g  %10g\n",
      p, q, x, mean, sd );
  }

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests CDFPOI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double p;
  double q;
  double s;
  int status;
  int which;
  double xlam;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  CDFPOI computes one missing parameter from the\n" );
  printf ( "  Poisson CDF:\n" );
  printf ( "\n" );
  printf ( "   POISSON_CDF ( (P,Q), S, XLAM )\n" );
  printf ( "\n" );
  printf ( "     P         Q         S         XLAM\n" );
  printf ( "\n" );

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      s = 3.0;
      xlam = 5.0;
    }
    else if ( which == 2 )
    {
      p = 0.265026;
      q = 1.0 - p;
      s = -1.0;
      xlam = 5.0;
    }
    else if ( which == 3 )
    {
      p = 0.265026;
      q = 1.0 - p;
      s = 3.0;
      xlam = -1.0;
    }

    cdfpoi ( &which, &p, &q, &s, &xlam, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFPOI returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g\n",
      p, q, s, xlam );
  }

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests CDFT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double bound;
  double df;
  double p;
  double q;
  int status;
  double t;
  int which;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  CDFT computes one missing parameter from the\n" );
  printf ( "  T CDF:\n" );
  printf ( "\n" );
  printf ( "   T_CDF ( (P,Q), T, DF )\n" );
  printf ( "\n" );
  printf ( "    P         Q         T         DF\n" );
  printf ( "\n" );

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      t = 3.0;
      df = 5.0;
    }
    else if ( which == 2 )
    {
      p = 0.984950;
      q = 1.0 - p;
      t = -1.0;
      df = 5.0;
    }
    else if ( which == 3 )
    {
      p = 0.984950;
      q = 1.0 - p;
      t = 3.0;
      df = -1.0;
    }

    cdft ( &which, &p, &q, &t, &df, &status, &bound );

    if ( status != 0 )
    {
      printf ( "\n" );
      printf ( "  CDFT returned STATUS = %d\n", status );
      continue;
    }
    printf ( "  %10g  %10g  %10g  %10g\n",
      p, q, t, df );
  }

  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests CUMBET, BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  CUMBET computes the Beta CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  BETA_INC_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         Y         A         B         CDF           CDF\n" );
  printf (
    "                                           (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    y = 1.0 - x;

    cumbet ( &x, &y, &a, &b, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      x, y, a, b, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X         Y         A         B         1-CDF         CCDF\n" );
  printf (
    "                                           (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    y = 1.0 - x;

    cumbet ( &x, &y, &a, &b, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      x, y, a, b, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests CUMBIN, BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double ompr;
  int s;
  double s_double;
  double pr;
  int x;
  double x_double;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  CUMBIN computes the Binomial CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  BINOMIAL_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "   X   S    Pr       CDF           CDF\n" );
  printf ( "                    (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &x, &pr, &s, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ompr = 1.0 - pr;

    s_double = ( double ) s;
    x_double = ( double ) x;

    cumbin ( &s_double, &x_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10d  %10g  %10g  %10g\n",
      s, x, pr, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "   X   S    Pr       1-CDF         CCDF\n" );
  printf ( "                    (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &x, &pr, &s, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    ompr = 1.0 - pr;

    s_double = ( double ) s;
    x_double = ( double ) x;

    cumbin ( &s_double, &x_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10d  %10g  %10g  %10g\n",
      s, x, pr, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests CUMCHI, CHI_SQUARE_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  CUMCHI computes the chi square CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  CHI_SQUARE_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X       DF    CDF           CDF\n" );
  printf ( "                 (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    df_double = ( double ) df;

    cumchi ( &x, &df_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g\n",
      x, df, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X       DF    1-CDF         CCDF\n" );
  printf ( "                 (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumchi ( &x, &df_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g\n",
      x, df, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests CUMCHN, CHI_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  CUMCHN computes the cumulative density\n" );
  printf ( "  function for the noncentral chi-squared\n" );
  printf ( "  distribution.\n" );
  printf ( "  CHI_NONCENTRAL_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    DF    Lambda    X         CDF           CDF\n" );
  printf ( "                             (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_noncentral_cdf_values ( &n_data, &x, &lambda, &df, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    df_double = ( double ) df;

    cumchn ( &x, &df_double, &lambda, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10g  %10g  %10g  %10g\n",
      df, lambda, x, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    DF    Lambda    X         1-CDF         CCDF\n" );
  printf ( "                             (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_noncentral_cdf_values ( &n_data, &x, &lambda, &df, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumchn ( &x, &df_double, &lambda, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10g  %10g  %10g  %10g\n",
      df, lambda, x, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests CUMF, F_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int dfd;
  double dfd_double;
  int dfn;
  double dfn_double;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  CUMF computes the F CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  F_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X      DFN DFD    CDF           CDF\n" );
  printf ( "                     (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &dfn, &dfd, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumf ( &x, &dfn_double, &dfd_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10d  %10d  %10g  %10g\n",
      x, dfn, dfd, cdf_lookup, cdf_compute );

  }

  printf ( "\n" );
  printf ( "    X      DFN DFD    1-CDF         CCDF\n" );
  printf ( "                     (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &dfn, &dfd, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumf ( &x, &dfn_double, &dfd_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10d  %10d  %10g  %10g\n",
      x, dfn, dfd, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests CUMFNC, F_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int dfd;
  double dfd_double;
  int dfn;
  double dfn_double;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  CUMFNC computes the noncentral F CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  F_NONCENTRAL_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X      DFN DFD    LAMBDA    CDF           CDF\n" );
  printf ( "                               (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &dfn, &dfd, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumfnc ( &x, &dfn_double, &dfd_double, &lambda, &cdf_compute,
      &ccdf_compute );

    printf ( "  %10g  %10d  %10d  %10g  %10g  %10g\n",
      x, dfn, dfd, lambda, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X      DFN DFD    LAMBDA    1-CDF         CCDF\n" );
  printf ( "                               (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &dfn, &dfd, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumfnc ( &x, &dfn_double, &dfd_double, &lambda, &cdf_compute,
      &ccdf_compute );

    printf ( "  %10g  %10d  %10d  %10g  %10g  %10g\n",
      x, dfn, dfd, lambda, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests CUMGAM, GAMMA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  CUMGAM computes the Gamma CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  GAMMA_INC_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    A         X         CDF           CDF\n" );
  printf ( "                        (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    cumgam ( &x, &a, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g\n",
      a, x, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    A         X         CDF           CDF\n" );
  printf ( "                        (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    cumgam ( &x, &a, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g  %10g\n",
      a, x, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test19 ( )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests CUMNBN, NEGATIVE_BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int f;
  double f_double;
  int n_data;
  double ompr;
  int s;
  double s_double;
  double pr;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  CUMNBN computes the Negative Binomial CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  NEGATIVE_BINOMIAL_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "   F   S    Pr       CDF           CDF\n" );
  printf ( "                     (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &pr, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ompr = 1.0 - pr;

    f_double = ( double ) f;
    s_double = ( double ) s;

    cumnbn ( &f_double, &s_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10d  %10g  %10g  %10g\n",
      f, s, pr, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "   F   S    Pr       1-CDF         CCDF\n" );
  printf ( "                     (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &pr, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    ompr = 1.0 - pr;

    f_double = ( double ) f;
    s_double = ( double ) s;

    cumnbn ( &f_double, &s_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10d  %10g  %10g  %10g\n",
      f, s, pr, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test20 ( )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests CUMNOR, NORMAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  CUMNOR computes the Normal CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  NORMAL_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         CDF           CDF\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    cumnor ( &x, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g\n",
      x, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X         1-CDF         CCDF\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    cumnor ( &x, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10g  %10g\n",
      x, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test21 ( )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests CUMPOI, POISSON_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double  cdf_lookup;
  double lambda;
  int n_data;
  int x;
  double x_double;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  CUMPOI computes the Poisson CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  POISSON_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "     X    LAMBDA    CDF           CDF\n" );
  printf ( "                   (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    x_double = ( double ) x;

    cumpoi ( &x_double, &lambda, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10g  %10g  %10g\n",
      x, lambda, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "     X    LAMBDA    1-CDF         CCDF\n" );
  printf ( "                   (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    x_double = ( double ) x;
    ccdf_lookup = 1.0 - cdf_lookup;

    cumpoi ( &x_double, &lambda, &cdf_compute, &ccdf_compute );

    printf ( "  %10d  %10g  %10g  %10g\n",
      x, lambda, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests CUMT, STUDENT_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  CUMT computes the Student T CDF\n" );
  printf ( "  and the complementary CDF.\n" );
  printf ( "  STUDENT_CDF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X       DF    CDF           CDF\n" );
  printf ( "                 (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }
    df_double = ( double ) df;

    cumt ( &x, &df_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10d  %10g  %10g\n",
      x, df, cdf_lookup, cdf_compute );
  }

  printf ( "\n" );
  printf ( "    X       DF    1-CDF         CCDF\n" );
  printf ( "                 (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumt ( &x, &df_double, &cdf_compute, &ccdf_compute );

    printf ( "  %10g  %10d  %10g  %10g\n",
      x, df, ccdf_lookup, ccdf_compute );
  }

  return;
}
/******************************************************************************/

void test23 ( )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests BETA, GAMMA_X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double apb;
  double beta1;
  double beta2;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  BETA evaluates the Beta function;\n" );
  printf ( "  GAMMA_X evaluates the Gamma function.\n" );

  a = 2.2;
  b = 3.7;
  apb = a + b;

  beta1 = beta ( a, b );
  beta2 = gamma_x ( &a ) * gamma_x ( &b ) / gamma_x ( &apb );

  printf ( "\n" );
  printf ( "  Argument A =                   %g\n", a );
  printf ( "  Argument B =                   %g\n", b );
  printf ( "  Beta(A,B) =                    %g\n", beta1 );
  printf ( "  (Expected value = 0.0454 )\n" );
  printf ( "\n" );
  printf ( "  Gamma(A)*Gamma(B)/Gamma(A+B) = %g\n", beta2 );

  return;
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests ERROR_F, ERROR_FC, ERF_VALUES..

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2006

  Author:

    John Burkardt
*/
{
  double erf_compute;
  double erf_lookup;
  double erfc_compute;
  double erfc_lookup;
  int ind;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  ERROR_F computes the error function ERF;\n" );
  printf ( "  ERROR_FC the complementary error function ERFC.\n" );
  printf ( "  ERF_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         ERF           ERF\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &erf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    erf_compute = error_f ( &x );

    printf ( "  %10g  %10g  %10g\n",
      x, erf_lookup, erf_compute );
  }

  printf ( "\n" );
  printf ( "    X         ERFC          ERFC\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  ind = 0;
  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &erf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    erfc_lookup = 1.0 - erf_lookup;
    erfc_compute = error_fc ( &ind, &x );

    printf ( "  %10g  %10g  %10g\n",
      x, erfc_lookup, erfc_compute );
  }

  return;
}
/******************************************************************************/

void test25 ( )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests XGAMM, GAMMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double gamma_compute;
  double gamma_lookup;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  XGAMM computes the Gamma function;\n" );
  printf ( "  GAMMA_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         GAMMA         GAMMA\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &gamma_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    gamma_compute = gamma_x ( &x );

    printf ( "  %10g  %10g  %10g\n",
      x, gamma_lookup, gamma_compute );
  }

  return;
}
/******************************************************************************/

void test26 ( )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests GAMMA_INC, GAMMA_INC_INV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  int ierror;
  int ind;
  double p;
  double q;
  int test_num = 10;
  double x;
  double x0;
  double x2;

  a = 3.0;
  ind = 1;
  x0 = 0;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  GAMMA_INC evaluates the incomplete Gamma ratio;\n" );
  printf ( "  GAMMA_INC_INV inverts it.\n" );
  printf ( "\n" );
  printf ( "  Parameters:\n" );
  printf ( "\n" );
  printf ( "    A = %g\n", a );
  printf ( "\n" );
  printf ( "    X             P             Q             Inverse\n" );
  printf ( "\n" );

  for ( i = 0; i <= test_num; i++ )
  {
    x = ( double ) i / ( double ) test_num;

    gamma_inc ( &a, &x, &p, &q, &ind );

    gamma_inc_inv ( &a, &x2, &x0, &p, &q, &ierror );

    printf ( "  %10g  %10g  %10g  %10g\n",
      x, p, q, x2 );
  }

  return;
}
/******************************************************************************/

void test27 ( )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests PSI, PSI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2007

  Author:

    John Burkardt
*/
{
  double psi_compute;
  double psi_lookup;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  PSI computes the Psi function;\n" );
  printf ( "  PSI_VALUES looks up some values.\n" );
  printf ( "\n" );
  printf ( "    X         PSI           PSI\n" );
  printf ( "              (Lookup)      (Computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &psi_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    psi_compute = psi ( &x );

    printf ( "  %10g  %10g  %10g\n",
      x, psi_lookup, psi_compute );
  }

  return;
}
