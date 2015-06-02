# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fn.h"
# include "test_values.h"

int main ( );
void acos_test ( );
void acosh_test ( );
void ai_test ( );
void aid_test ( );
void asin_test ( );
void asinh_test ( );
void atan_test ( );
void atan2_test ( );
void atanh_test ( );
void besi0_test ( );
void besi1_test ( );
void besj0_test ( );
void besj1_test ( );
void besk_test ( );
void besk0_test ( );
void besk1_test ( );
void besy0_test ( );
void besy1_test ( );
void beta_test ( );
void betai_test ( );
void bi_test ( );
void bid_test ( );
void binom_test ( );
void cbrt_test ( );
void chi_test ( );
void chu_test ( );
void ci_test ( );
void cin_test ( );
void cinh_test ( );
void cos_test ( );
void cos_deg_test ( );
void cosh_test ( );
void cot_test ( );
void dawson_test ( );
void e1_test ( );
void ei_test ( );
void erf_test ( );
void erfc_test ( );
void exp_test ( );
void fac_test ( );
void gamma_test ( );
void gamma_inc_test ( );
void gamma_inc_tricomi_test ( );
void int_test ( );
void lbeta_test ( );
void li_test ( );
void lngam_test ( );
void log_test ( );
void log10_test ( );
void poch_test ( );
void psi_test ( );
void rand_test ( );
void shi_test ( );
void si_test ( );
void sin_test ( );
void sin_deg_test ( );
void sinh_test ( );
void spence_test ( );
void sqrt_test ( );
void tan_test ( );
void tanh_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FN_PRB.

  Discussion:

    FN_PRB tests the FN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 November 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FN_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FN library.\n" );

  acos_test ( );
  acosh_test ( );
  ai_test ( );
  aid_test ( );
  asin_test ( );
  asinh_test ( );
  atan_test ( );
  atan2_test ( );
  atanh_test ( );
  besi0_test ( );
  besi1_test ( );
  besj0_test ( );
  besj1_test ( );
  besk_test ( );
  besk0_test ( );
  besk1_test ( );
  besy0_test ( );
  besy1_test ( );
  beta_test ( );
  betai_test ( );
  bi_test ( );
  bid_test ( );
  binom_test ( );
  cbrt_test ( );
  chi_test ( );
  chu_test ( );
  ci_test ( );
  cin_test ( );
  cinh_test ( );
  cos_test ( );
  cos_deg_test ( );
  cosh_test ( );
  cot_test ( );
  dawson_test ( );
  e1_test ( );
  ei_test ( );
  erf_test ( );
  erfc_test ( );
  exp_test ( );
  fac_test ( );
  gamma_test ( );
  gamma_inc_test ( );
  gamma_inc_tricomi_test ( );
  int_test ( );
  lbeta_test ( );
  li_test ( );
  lngam_test ( );
  log_test ( );
  log10_test ( );
  poch_test ( );
  psi_test ( );
  rand_test ( );
  shi_test ( );
  si_test ( );
  sin_test ( );
  sin_deg_test ( );
  sinh_test ( );
  spence_test ( );
  sqrt_test ( );
  tan_test ( );
  tanh_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FN_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void acos_test ( )

/******************************************************************************/
/*
  Purpose:

    ACOS_TEST tests R4_ACOS and R8_ACOS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ACOS_TEST:\n" );
  printf ( "  Test ARCCOS_VALUES, R4_ACOS, R8_ACOS.\n" );
  printf ( "\n" );
  printf ( "             X      ARCCOS(X)\n" );
  printf ( "                   R4_ACOS(X)         Diff\n" );
  printf ( "                   R8_ACOS(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_acos ( ( float ) x );
    fx3 = r8_acos ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void acosh_test ( )

/******************************************************************************/
/*
  Purpose:

    ACOSH_TEST tests R4_ACOSH and R8_ACOSH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ACOSH_TEST:\n" );
  printf ( "  Test ARCCOSH_VALUES, R4_ACOSH, R8_ACOSH\n" );
  printf ( "\n" );
  printf ( "             X      ARCCOSH(X)\n" );
  printf ( "                   R4_ACOSH(X)        Diff\n" );
  printf ( "                   R8_ACOSH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccosh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_acosh ( ( float ) x );
    fx3 = r8_acosh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void ai_test ( )

/******************************************************************************/
/*
  Purpose:

    AI_TEST tests R4_AI and R8_AI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AI_TEST:\n" );
  printf ( "  Test AIRY_AI_VALUES, R4_AI, R8_AI.\n" );
  printf ( "\n" );
  printf ( "             X   AIRY_AI(X)\n" );
  printf ( "                   R4_AI(X)         Diff\n" );
  printf ( "                   R8_AI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ai ( ( float ) x );
    fx3 = r8_ai ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );

  }
  return;
}
/******************************************************************************/

void aid_test ( )

/******************************************************************************/
/*
  Purpose:

    AID_TEST tests R4_AID and R8_AID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AID_TEST:\n" );
  printf ( "  Test AIRY_AI_PRIME_VALUES, R4_AID, R8_AID.\n" );
  printf ( "\n" );
  printf ( "             X   AIRY_AID(X)\n" );
  printf ( "                   R4_AID(X)         Diff\n" );
  printf ( "                   R8_AID(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_prime_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_aid ( ( float ) x );
    fx3 = r8_aid ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );

  }
  return;
}
/******************************************************************************/

void asin_test ( )

/******************************************************************************/
/*
  Purpose:

    ASIN_TEST tests R4_ASIN and R8_ASIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ASIN_TEST:\n" );
  printf ( "  Test ARCSIN_VALUES, R4_ASIN, R8_ASIN.\n" );
  printf ( "\n" );
  printf ( "             X      ARCSIN(X)\n" );
  printf ( "                   R4_ASIN(X)         Diff\n" );
  printf ( "                   R8_ASIN(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_asin ( ( float ) x );
    fx3 = r8_asin ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void asinh_test ( )

/******************************************************************************/
/*
  Purpose:

    ASINH_TEST tests R4_ASINH and R8_ASINH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ASINH_TEST:\n" );
  printf ( "  Test ARCSINH_VALUES, R4_ASINH, R8_ASINH\n" );
  printf ( "\n" );
  printf ( "             X      ARCSINH(X)\n" );
  printf ( "                   R4_ASINH(X)        Diff\n" );
  printf ( "                   R8_ASINH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsinh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_asinh ( ( float ) x );
    fx3 = r8_asinh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void atan_test ( )

/******************************************************************************/
/*
  Purpose:

    ATAN_TEST tests R4_ATAN and R8_ATAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ATAN_TEST:\n" );
  printf ( "  Test ARCTAN_VALUES, R4_ATAN, R8_ATAN.\n" );
  printf ( "\n" );
  printf ( "             X      ARCTAN(X)\n" );
  printf ( "                   R4_ATAN(X)         Diff\n" );
  printf ( "                   R8_ATAN(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atan ( ( float ) x );
    fx3 = r8_atan ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void atan2_test ( )

/******************************************************************************/
/*
  Purpose:

    ATAN2_TEST tests R4_ATAN2 and R8_ATAN2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "ATAN2_TEST:\n" );
  printf ( "  Test ARCTAN2_VALUES, R4_ATAN2, R8_ATAN2.\n" );
  printf ( "\n" );
  printf ( "             X             Y      ARCTAN2(Y,X)\n" );
  printf ( "                                 R4_ATAN2(Y,X)         Diff\n" );
  printf ( "                                 R8_ATAN2(Y,X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan2_values ( &n_data, &x, &y, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atan2 ( ( float ) y, ( float ) x );
    fx3 = r8_atan2 ( y, x );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", x, y, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void atanh_test ( )

/******************************************************************************/
/*
  Purpose:

    ATANH_TEST tests R4_ATANH and R8_ATANH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ATANH_TEST:\n" );
  printf ( "  Test ARCTANH_VALUES, R4_ATANH, R8_ATANH\n" );
  printf ( "\n" );
  printf ( "             X      ARCTANH(X)\n" );
  printf ( "                   R4_ATANH(X)        Diff\n" );
  printf ( "                   R8_ATANH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctanh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atanh ( ( float ) x );
    fx3 = r8_atanh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besi0_test ( )

/******************************************************************************/
/*
  Purpose:

    BESI0_TEST tests R4_BESI0 and R8_BESI0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESI0_TEST:\n" );
  printf ( "  Test BESSEL_I0_VALUES, R4_BESI0, R8_BESI0\n" );
  printf ( "\n" );
  printf ( "             X      BESI0(X)\n" );
  printf ( "                   R4_BESI0(X)        Diff\n" );
  printf ( "                   R8_BESI0(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besi0 ( ( float ) x );
    fx3 = r8_besi0 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besi1_test ( )

/******************************************************************************/
/*
  Purpose:

    BESI1_TEST tests R4_BESI1 and R8_BESI1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESI1_TEST:\n" );
  printf ( "  Test BESSEL_I1_VALUES, R4_BESI1, R8_BESI1\n" );
  printf ( "\n" );
  printf ( "             X      BESI1(X)\n" );
  printf ( "                   R4_BESI1(X)        Diff\n" );
  printf ( "                   R8_BESI1(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besi1 ( ( float ) x );
    fx3 = r8_besi1 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besj0_test ( )

/******************************************************************************/
/*
  Purpose:

    BESJ0_TEST tests R4_BESJ0 and R8_BESJ0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESJ0_TEST:\n" );
  printf ( "  Test BESSEL_J0_VALUES, R4_BESJ0, R8_BESJ0\n" );
  printf ( "\n" );
  printf ( "             X      BESJ0(X)\n" );
  printf ( "                   R4_BESJ0(X)        Diff\n" );
  printf ( "                   R8_BESJ0(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besj0 ( ( float ) x );
    fx3 = r8_besj0 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besj1_test ( )

/******************************************************************************/
/*
  Purpose:

    BESJ1_TEST tests R4_BESJ1 and R8_BESJ1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESJ1_TEST:\n" );
  printf ( "  Test BESSEL_J1_VALUES, R4_BESJ1, R8_BESJ1\n" );
  printf ( "\n" );
  printf ( "             X      BESJ1(X)\n" );
  printf ( "                   R4_BESJ1(X)        Diff\n" );
  printf ( "                   R8_BESJ1(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besj1 ( ( float ) x );
    fx3 = r8_besj1 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besk_test ( )

/******************************************************************************/
/*
  Purpose:

    BESK_TEST tests R4_BESK and R8_BESK.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 November 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double nu;
  double x;

  printf ( "\n" );
  printf ( "BESK_TEST:\n" );
  printf ( "  Test BESSEL_KX_VALUES, R4_BESK, R8_BESK\n" );
  printf ( "\n" );
  printf ( "              NU               X      BESK(X)\n" );
  printf ( "                                   R4_BESK(X)        Diff\n" );
  printf ( "                                   R8_BESK(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_kx_values ( &n_data, &nu, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besk ( ( float ) nu, ( float ) x );
    fx3 = r8_besk ( nu, x );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", nu, x, fx1 );
    printf ( "                                  %14g  %14g\n", 
      fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", 
      fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besk0_test ( )

/******************************************************************************/
/*
  Purpose:

    BESK0_TEST tests R4_BESK0 and R8_BESK0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESK0_TEST:\n" );
  printf ( "  Test BESSEL_K0_VALUES, R4_BESK0, R8_BESK0\n" );
  printf ( "\n" );
  printf ( "             X      BESK0(X)\n" );
  printf ( "                   R4_BESK0(X)        Diff\n" );
  printf ( "                   R8_BESK0(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besk0 ( ( float ) x );
    fx3 = r8_besk0 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besk1_test ( )

/******************************************************************************/
/*
  Purpose:

    BESK1_TEST tests R4_BESK1 and R8_BESK1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESK1_TEST:\n" );
  printf ( "  Test BESSEL_K1_VALUES, R4_BESK1, R8_BESK1\n" );
  printf ( "\n" );
  printf ( "             X      BESK1(X)\n" );
  printf ( "                   R4_BESK1(X)        Diff\n" );
  printf ( "                   R8_BESK1(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k1_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besk1 ( ( float ) x );
    fx3 = r8_besk1 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besy0_test ( )

/******************************************************************************/
/*
  Purpose:

    BESY0_TEST tests R4_BESY0 and R8_BESY0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESY0_TEST:\n" );
  printf ( "  Test BESSEL_Y0_VALUES, R4_BESY0, R8_BESY0\n" );
  printf ( "\n" );
  printf ( "             X      BESY0(X)\n" );
  printf ( "                   R4_BESY0(X)        Diff\n" );
  printf ( "                   R8_BESY0(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besy0 ( ( float ) x );
    fx3 = r8_besy0 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void besy1_test ( )

/******************************************************************************/
/*
  Purpose:

    BESY1_TEST tests R4_BESY1 and R8_BESY1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESY1_TEST:\n" );
  printf ( "  Test BESSEL_Y1_VALUES, R4_BESY1, R8_BESY1\n" );
  printf ( "\n" );
  printf ( "             X      BESY1(X)\n" );
  printf ( "                   R4_BESY1(X)        Diff\n" );
  printf ( "                   R8_BESY1(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y1_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besy1 ( ( float ) x );
    fx3 = r8_besy1 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void beta_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_TEST tests R4_BETA and R8_BETA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;

  printf ( "\n" );
  printf ( "BETA_TEST:\n" );
  printf ( "  Test BETA_VALUES, R4_BETA, R8_BETA.\n" );
  printf ( "\n" );
  printf ( "             A        B        BETA(A,B)\n" );
  printf ( "                           R4_BETA(A,B)       Diff\n" );
  printf ( "                           R8_BETA(A,B)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &a, &b, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_beta ( ( float ) a, ( float ) b );
    fx3 = r8_beta ( a, b );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", a, b, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void betai_test ( )

/******************************************************************************/
/*
  Purpose:

    BETAI_TEST tests R4_BETAI and R8_BETAI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BETAI_TEST:\n" );
  printf ( "  Test BETA_INC_VALUES, R4_BETAI, R8_BETAI.\n" );
  printf ( "\n" );
  printf ( "             X        BETA(A,B,X)\n" );
  printf ( "                   R4_BETAI(A,B,X)       Diff\n" );
  printf ( "                   R8_BETAI(A,B,X)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_betai ( ( float ) x, ( float ) a, ( float ) b );
    fx3 = r8_betai ( x, a, b );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g  %14g\n", a, b, x, fx1 );
    printf ( "                                                    %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                                    %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void bi_test ( )

/******************************************************************************/
/*
  Purpose:

    BI_TEST tests R4_BI and R8_BI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BI_TEST:\n" );
  printf ( "  Test AIRY_BI_VALUES, R4_BI, R8_BI.\n" );
  printf ( "\n" );
  printf ( "             X   AIRY_BI(X)\n" );
  printf ( "                   R4_BI(X)         Diff\n" );
  printf ( "                   R8_BI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_bi ( ( float ) x );
    fx3 = r8_bi ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void bid_test ( )

/******************************************************************************/
/*
  Purpose:

    BID_TEST tests R4_BID and R8_BID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BID_TEST:\n" );
  printf ( "  Test AIRY_BI_PRIME_VALUES, R4_BID, R8_BID.\n" );
  printf ( "\n" );
  printf ( "             X   AIRY_BID(X)\n" );
  printf ( "                   R4_BID(X)         Diff\n" );
  printf ( "                   R8_BID(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_prime_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_bid ( ( float ) x );
    fx3 = r8_bid ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void binom_test ( )

/******************************************************************************/
/*
  Purpose:

    BINOM_TEST tests R4_BINOM and R8_BINOM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double diff;
  int fx1;
  float fx2;
  double fx3;
  int n_data;

  printf ( "\n" );
  printf ( "BINOM_TEST:\n" );
  printf ( "  Test BINOM_VALUES, R4_BINOM, R8_BINOM.\n" );
  printf ( "\n" );
  printf ( "             A    B        BINOM(A,B)\n" );
  printf ( "                        R4_BINOM(A,B)       Diff\n" );
  printf ( "                        R8_BINOM(A,B)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_values ( &n_data, &a, &b, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_binom ( ( float ) a, ( float ) b );
    fx3 = r8_binom ( ( double ) a, ( double ) b );

    printf ( "\n" );
    printf ( "  %14d  %14d  %14d\n", a, b, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );

  }
  return;
}
/******************************************************************************/

void cbrt_test ( )

/******************************************************************************/
/*
  Purpose:

    CBRT_TEST tests R4_CBRT and R8_CBRT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CBRT_TEST:\n" );
  printf ( "  Test CBRT_VALUES, R4_CBRT, R8_CBRT\n" );
  printf ( "\n" );
  printf ( "             X      CBRT(X)\n" );
  printf ( "                   R4_CBRT(X)        Diff\n" );
  printf ( "                   R8_CBRT(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cbrt ( ( float ) x );
    fx3 = r8_cbrt ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void chi_test ( )

/******************************************************************************/
/*
  Purpose:

    CHI_TEST tests R4_CHI and R8_CHI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHI_TEST:\n" );
  printf ( "  Test CHI_VALUES, R4_CHI, R8_CHI.\n" );
  printf ( "\n" );
  printf ( "             X      CHI(X)\n" );
  printf ( "                   R4_CHI(X)         Diff\n" );
  printf ( "                   R8_CHI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_chi ( ( float ) x );
    fx3 = r8_chi ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void chu_test ( )

/******************************************************************************/
/*
  Purpose:

    CHU_TEST tests R4_CHU and R8_CHU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double diff;
  double fx1;
  double fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHU_TEST:\n" );
  printf ( "  Test HYPERGEOMETRIC_U_VALUES, R4_CHU, R8_CHU.\n" );
  printf ( "\n" );
  printf ( "             A               B               X     CHU(A,B,X)\n" );
  printf ( "                                                R4_CHU(A,B,X)" );
  printf ( "         Diff\n" );
  printf ( "                                                R8_CHU(A,B,X)" );
  printf ( "         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_u_values ( &n_data, &a, &b, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_chu ( ( float ) a, ( float ) b, ( float ) x );
    fx3 = r8_chu ( a, b, x );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g  %14g\n", a, b, x, fx1 );
    printf ( "                                                    %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                                    %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }

  return;
}
/******************************************************************************/

void ci_test ( )

/******************************************************************************/
/*
  Purpose:

    CI_TEST tests R4_CI and R8_CI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CI_TEST:\n" );
  printf ( "  Test CI_VALUES, R4_CI, R8_CI.\n" );
  printf ( "\n" );
  printf ( "             X      CI(X)\n" );
  printf ( "                   R4_CI(X)         Diff\n" );
  printf ( "                   R8_CI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    ci_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ci ( ( float ) x );
    fx3 = r8_ci ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cin_test ( )

/******************************************************************************/
/*
  Purpose:

    CIN_TEST tests R4_CIN and R8_CIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CIN_TEST:\n" );
  printf ( "  Test CIN_VALUES, R4_CIN, R8_CIN.\n" );
  printf ( "\n" );
  printf ( "             X      CIN(X)\n" );
  printf ( "                   R4_CIN(X)         Diff\n" );
  printf ( "                   R8_CIN(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cin_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cin ( ( float ) x );
    fx3 = r8_cin ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cinh_test ( )

/******************************************************************************/
/*
  Purpose:

    CINH_TEST tests R4_CINH and R8_CINH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CINH_TEST:\n" );
  printf ( "  Test CINH_VALUES, R4_CINH, R8_CINH.\n" );
  printf ( "\n" );
  printf ( "             X      CINH(X)\n" );
  printf ( "                   R4_CINH(X)         Diff\n" );
  printf ( "                   R8_CINH(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cinh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cinh ( ( float ) x );
    fx3 = r8_cinh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cos_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_TEST tests R4_COS and R8_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COS_TEST:\n" );
  printf ( "  Test COS_VALUES, R4_COS, R8_COS.\n" );
  printf ( "\n" );
  printf ( "             X      COS(X)\n" );
  printf ( "                   R4_COS(X)         Diff\n" );
  printf ( "                   R8_COS(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cos ( ( float ) x );
    fx3 = r8_cos ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cos_deg_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_DEG_TEST tests R4_COS_DEG and R8_COS_DEG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COS_DEG_TEST:\n" );
  printf ( "  Test COS_DEGREE_VALUES, R4_COS_DEG, R8_COS_DEG.\n" );
  printf ( "\n" );
  printf ( "             X      COS_DEG(X)\n" );
  printf ( "                   R4_COS_DEG(X)         Diff\n" );
  printf ( "                   R8_COS_DEG(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_degree_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cos_deg ( ( float ) x );
    fx3 = r8_cos_deg ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cosh_test ( )

/******************************************************************************/
/*
  Purpose:

    COSH_TEST tests R4_COSH and R8_COSH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COSH_TEST:\n" );
  printf ( "  Test COSH_VALUES, R4_COSH, R8_COSH\n" );
  printf ( "\n" );
  printf ( "             X      COSH(X)\n" );
  printf ( "                   R4_COSH(X)        Diff\n" );
  printf ( "                   R8_COSH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cosh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cosh ( ( float ) x );
    fx3 = r8_cosh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void cot_test ( )

/******************************************************************************/
/*
  Purpose:

    COT_TEST tests R4_COT and R8_COT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COT_TEST:\n" );
  printf ( "  Test COT_VALUES, R4_COT, R8_COT.\n" );
  printf ( "\n" );
  printf ( "             X      COT(X)\n" );
  printf ( "                   R4_COT(X)         Diff\n" );
  printf ( "                   R8_COT(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    cot_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cot ( ( float ) x );
    fx3 = r8_cot ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void dawson_test ( )

/******************************************************************************/
/*
  Purpose:

    DAWSON_TEST tests R4_DAWSON and R8_DAWSON.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DAWSON_TEST:\n" );
  printf ( "  Test DAWSON_VALUES, R4_DAWSON, R8_DAWSON.\n" );
  printf ( "\n" );
  printf ( "             X      DAWSON(X)\n" );
  printf ( "                   R4_DAWSON(X)         Diff\n" );
  printf ( "                   R8_DAWSON(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_dawson ( ( float ) x );
    fx3 = r8_dawson ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void e1_test ( )

/******************************************************************************/
/*
  Purpose:

    E1_TEST tests R4_E1 and R8_E1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "E1_TEST:\n" );
  printf ( "  Test E1_VALUES, R4_E1, R8_E1.\n" );
  printf ( "\n" );
  printf ( "             X      E1(X)\n" );
  printf ( "                   R4_E1(X)         Diff\n" );
  printf ( "                   R8_E1(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    e1_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_e1 ( ( float ) x );
    fx3 = r8_e1 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void ei_test ( )

/******************************************************************************/
/*
  Purpose:

    EI_TEST tests R4_EI and R8_EI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EI_TEST:\n" );
  printf ( "  Test EI_VALUES, R4_EI, R8_EI.\n" );
  printf ( "\n" );
  printf ( "             X      EI(X)\n" );
  printf ( "                   R4_EI(X)         Diff\n" );
  printf ( "                   R8_EI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    ei_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ei ( ( float ) x );
    fx3 = r8_ei ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void erf_test ( )

/******************************************************************************/
/*
  Purpose:

    ERF_TEST tests R4_ERF and R8_ERF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ERF_TEST:\n" );
  printf ( "  Test ERF_VALUES, R4_ERF, R8_ERF.\n" );
  printf ( "\n" );
  printf ( "             X      ERF(X)\n" );
  printf ( "                   R4_ERF(X)         Diff\n" );
  printf ( "                   R8_ERF(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_erf ( ( float ) x );
    fx3 = r8_erf ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void erfc_test ( )

/******************************************************************************/
/*
  Purpose:

    ERFC_TEST tests R4_ERFC and R8_ERFC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ERFC_TEST:\n" );
  printf ( "  Test ERFC_VALUES, R4_ERFC, R8_ERFC.\n" );
  printf ( "\n" );
  printf ( "             X      ERFC(X)\n" );
  printf ( "                   R4_ERFC(X)         Diff\n" );
  printf ( "                   R8_ERFC(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_erfc ( ( float ) x );
    fx3 = r8_erfc ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void exp_test ( )

/******************************************************************************/
/*
  Purpose:

    EXP_TEST tests R4_EXP and R8_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EXP_TEST:\n" );
  printf ( "  Test EXP_VALUES, R4_EXP, R8_EXP.\n" );
  printf ( "\n" );
  printf ( "             X      EXP(X)\n" );
  printf ( "                   R4_EXP(X)         Diff\n" );
  printf ( "                   R8_EXP(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_exp ( ( float ) x );
    fx3 = r8_exp ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void fac_test ( )

/******************************************************************************/
/*
  Purpose:

    FAC_TEST tests R4_FAC and R8_FAC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  int fx1;
  float fx2;
  double fx3;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "FAC_TEST:\n" );
  printf ( "  Test FACTORIAL_VALUES, R4_FAC, R8_FAC.\n" );
  printf ( "\n" );
  printf ( "             N      FAC(N)\n" );
  printf ( "                   R4_FAC(N)         Diff\n" );
  printf ( "                   R8_FAC(N)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    factorial_values ( &n_data, &n, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_fac ( n );
    fx3 = r8_fac ( n );

    printf ( "\n" );
    printf ( "  %14d  %14d\n", n, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void gamma_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_TEST tests R4_GAMMA and R8_GAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_TEST:\n" );
  printf ( "  Test GAMMA_VALUES, R4_GAMMA, R8_GAMMA\n" );
  printf ( "\n" );
  printf ( "             X      GAMMA(X)\n" );
  printf ( "                   R4_GAMMA(X)        Diff\n" );
  printf ( "                   R8_GAMMA(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamma ( ( float ) x );
    fx3 = r8_gamma ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void gamma_inc_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_TEST tests R4_GAMIC and R8_GAMIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_TEST:\n" );
  printf ( "  Test GAMMA_INC_VALUES, R4_GAMIC, R8_GAMIC.\n" );
  printf ( "\n" );
  printf ( "             X        GAMIC(A,X)\n" );
  printf ( "                   R4_GAMIC(A,X)       Diff\n" );
  printf ( "                   R8_GAMIC(A,X)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamic ( ( float ) a, ( float ) x );
    fx3 = r8_gamic ( a, x );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", a, x, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void gamma_inc_tricomi_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_TRICOMI_TEST tests R4_GAMIT and R8_GAMIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_TRICOMI_TEST:\n" );
  printf ( "  Test GAMMA_INC_TRICOMI_VALUES, R4_GAMIT, R8_GAMIT.\n" );
  printf ( "\n" );
  printf ( "        A     X        GAMIT(A,X)\n" );
  printf ( "                    R4_GAMIT(A,X)       Diff\n" );
  printf ( "                    R8_GAMIT(A,X)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_tricomi_values ( &n_data, &a, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamit ( ( float ) a, ( float ) x );
    fx3 = r8_gamit ( a, x );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", a, x, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void int_test ( )

/******************************************************************************/
/*
  Purpose:

    INT_TEST tests R4_INT and R8_INT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "INT_TEST:\n" );
  printf ( "  Test INT_VALUES, R4_INT, R8_INT\n" );
  printf ( "\n" );
  printf ( "             X      INT(X)\n" );
  printf ( "                   R4_INT(X)         Diff\n" );
  printf ( "                   R8_INT(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    int_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_int ( ( float ) x );
    fx3 = r8_int ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void lbeta_test ( )

/******************************************************************************/
/*
  Purpose:

    LBETA_TEST tests R4_LBETA and R8_LBETA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;

  printf ( "\n" );
  printf ( "LBETA_TEST:\n" );
  printf ( "  Test BETA_LOG_VALUES, R4_LBETA, R8_LBETA.\n" );
  printf ( "\n" );
  printf ( "             A  B        LBETA(A,B)\n" );
  printf ( "                   R4_LBETA(A,B)       Diff\n" );
  printf ( "                   R8_LBETA(A,B)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_log_values ( &n_data, &a, &b, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_lbeta ( ( float ) a, ( float ) b );
    fx3 = r8_lbeta ( a, b );

    printf ( "\n" );
    printf ( "  %14g  %14g  %14g\n", a, b, fx1 );
    printf ( "                                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void li_test ( )

/******************************************************************************/
/*
  Purpose:

    LI_TEST tests R4_LI and R8_LI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LI_TEST:\n" );
  printf ( "  Test LOGARITHMIC_INTEGRAL_VALUES, R4_LI, R8_LI\n" );
  printf ( "\n" );
  printf ( "             X      LI(X)\n" );
  printf ( "                   R4_LI(X)        Diff\n" );
  printf ( "                   R8_LI(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    logarithmic_integral_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_li ( ( float ) x );
    fx3 = r8_li ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void lngam_test ( )

/******************************************************************************/
/*
  Purpose:

    LNGAM_TEST tests R4_LNGAM and R8_LNGAM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LNGAM_TEST:\n" );
  printf ( "  Test GAMMA_LOG_VALUES, R4_LNGAM, R8_LNGAM\n" );
  printf ( "\n" );
  printf ( "             X        LNGAM(X)\n" );
  printf ( "                   R4_LNGAM(X)        Diff\n" );
  printf ( "                   R8_LNGAM(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_lngam ( ( float ) x );
    fx3 = r8_lngam ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void log_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG_TEST tests R4_LOG and R8_LOG.

  Licensing:

    This code is distributed under the GNU LGPL logcense.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOG_TEST:\n" );
  printf ( "  Test LOG_VALUES, R4_LOG, R8_LOG\n" );
  printf ( "\n" );
  printf ( "             X      LOG(X)\n" );
  printf ( "                   R4_LOG(X)        Diff\n" );
  printf ( "                   R8_LOG(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_log ( ( float ) x );
    fx3 = r8_log ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void log10_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG10_TEST tests R4_LOG10 and R8_LOG10.

  Licensing:

    This code is distributed under the GNU LGPL log10cense.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOG10_TEST:\n" );
  printf ( "  Test LOG10_VALUES, R4_LOG10, R8_LOG10\n" );
  printf ( "\n" );
  printf ( "             X      LOG10(X)\n" );
  printf ( "                   R4_LOG10(X)        Diff\n" );
  printf ( "                   R8_LOG10(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    log10_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_log10 ( ( float ) x );
    fx3 = r8_log10 ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void poch_test ( )

/******************************************************************************/
/*
  Purpose:

    POCH_TEST tests R4_POCH and R8_POCH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "POCH_TEST:\n" );
  printf ( "  Test POCHHAMMER_VALUES, R4_POCH, R8_POCH.\n" );
  printf ( "\n" );
  printf ( "             X        POCH(A,X)\n" );
  printf ( "                   R4_POCH(A,X)       Diff\n" );
  printf ( "                   R8_POCH(A,X)       Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    pochhammer_values ( &n_data, &a, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_poch ( ( float ) a, ( float ) x );
    fx3 = r8_poch ( a, x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void psi_test ( )

/******************************************************************************/
/*
  Purpose:

    PSI_TEST tests R4_PSI and R8_PSI.

  Licensing:

    This code is distributed under the GNU LGPL psicense.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "PSI_TEST:\n" );
  printf ( "  Test PSI_VALUES, R4_PSI, R8_PSI\n" );
  printf ( "\n" );
  printf ( "             X      PSI(X)\n" );
  printf ( "                   R4_PSI(X)        Diff\n" );
  printf ( "                   R8_PSI(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_psi ( ( float ) x );
    fx3 = r8_psi ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void rand_test ( )

/******************************************************************************/
/*
  Purpose:

    RAND_TEST tests R4_RAND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  float average;
  int i;
  int i_value[7] = { 1, 2, 3, 4, 10, 100, 1000 };
  int k;
  float r;
  float r_value[7] = { 
    0.0004127026, 
    0.6750836372, 
    0.1614754200, 
    0.9086198807, 
    0.5527787209, 
    0.3600893021, 
    0.2176990509 };
  float variance;

  printf ( "\n" );
  printf ( "RAND_TEST:\n" );
  printf ( "  Test R4_RAND.\n" );
  printf ( "\n" );
  printf ( "               I       R4_RAND        Expected\n" );
  printf ( "\n" );

  k = 0;

  for ( i = 1; i <= 1000; i++ )
  {
    r = r4_rand ( 0.0 );

    if ( i == i_value[k] )
    {
      printf ( "  %14d  %14g  %14g\n", i, r, r_value[k] );
      k = k + 1;
    }
  }

  average = 0.0;
  for ( i = 1; i <= 1000000; i++ )
  {
    r = r4_rand ( 0.0 );
    average = average + r;
  }
  average = average / 1000000.0;
  printf ( "\n" );
  printf ( "     Average =  %14g  %14g\n", average, 0.5 );

  variance = 0.0;
  for ( i = 1; i <= 1000000; i++ )
  {
    r = r4_rand ( 0.0 );
    variance = variance + pow ( r - average, 2 );
  }
  variance = variance / 1000000.0;
  printf ( "     Variance = %14g  %14g\n", variance, 1.0 / 12.0 );

  return;
}
/******************************************************************************/

void shi_test ( )

/******************************************************************************/
/*
  Purpose:

    SHI_TEST tests R4_SHI and R8_SHI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SHI_TEST:\n" );
  printf ( "  Test SHI_VALUES, R4_SHI, R8_SHI.\n" );
  printf ( "\n" );
  printf ( "             X      SHI(X)\n" );
  printf ( "                   R4_SHI(X)         Diff\n" );
  printf ( "                   R8_SHI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    shi_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_shi ( ( float ) x );
    fx3 = r8_shi ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void si_test ( )

/******************************************************************************/
/*
  Purpose:

    SI_TEST tests R4_SI and R8_SI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SI_TEST:\n" );
  printf ( "  Test SI_VALUES, R4_SI, R8_SI.\n" );
  printf ( "\n" );
  printf ( "             X      SI(X)\n" );
  printf ( "                   R4_SI(X)         Diff\n" );
  printf ( "                   R8_SI(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    si_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_si ( ( float ) x );
    fx3 = r8_si ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void sin_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_TEST tests R4_SIN and R8_SIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SIN_TEST:\n" );
  printf ( "  Test SIN_VALUES, R4_SIN, R8_SIN.\n" );
  printf ( "\n" );
  printf ( "             X      SIN(X)\n" );
  printf ( "                   R4_SIN(X)         Diff\n" );
  printf ( "                   R8_SIN(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sin ( ( float ) x );
    fx3 = r8_sin ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void sin_deg_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_DEG_TEST tests R4_SIN_DEG and R8_SIN_DEG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SIN_DEG_TEST:\n" );
  printf ( "  Test SIN_DEGREE_VALUES, R4_SIN_DEG, R8_SIN_DEG.\n" );
  printf ( "\n" );
  printf ( "             X      SIN_DEG(X)\n" );
  printf ( "                   R4_SIN_DEG(X)         Diff\n" );
  printf ( "                   R8_SIN_DEG(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_degree_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sin_deg ( ( float ) x );
    fx3 = r8_sin_deg ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void sinh_test ( )

/******************************************************************************/
/*
  Purpose:

    SINH_TEST tests R4_SINH and R8_SINH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SINH_TEST:\n" );
  printf ( "  Test SINH_VALUES, R4_SINH, R8_SINH\n" );
  printf ( "\n" );
  printf ( "             X      SINH(X)\n" );
  printf ( "                   R4_SINH(X)        Diff\n" );
  printf ( "                   R8_SINH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    sinh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sinh ( ( float ) x );
    fx3 = r8_sinh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void spence_test ( )

/******************************************************************************/
/*
  Purpose:

    SPENCE_TEST tests R4_SPENCE and R8_SPENCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SPENCE_TEST:\n" );
  printf ( "  Test DILOGARITHM_VALUES, R4_SPENCE, R8_SPENCE\n" );
  printf ( "\n" );
  printf ( "             X      SPENCE(X)\n" );
  printf ( "                   R4_SPENCE(X)        Diff\n" );
  printf ( "                   R8_SPENCE(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    dilogarithm_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_spence ( ( float ) x );
    fx3 = r8_spence ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void sqrt_test ( )

/******************************************************************************/
/*
  Purpose:

    SQRT_TEST tests R4_SQRT and R8_SQRT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SQRT_TEST:\n" );
  printf ( "  Test SQRT_VALUES, R4_SQRT, R8_SQRT\n" );
  printf ( "\n" );
  printf ( "             X      SQRT(X)\n" );
  printf ( "                   R4_SQRT(X)        Diff\n" );
  printf ( "                   R8_SQRT(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    sqrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sqrt ( ( float ) x );
    fx3 = r8_sqrt ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void tan_test ( )

/******************************************************************************/
/*
  Purpose:

    TAN_TEST tests R4_TAN and R8_TAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TAN_TEST:\n" );
  printf ( "  Test TAN_VALUES, R4_TAN, R8_TAN.\n" );
  printf ( "\n" );
  printf ( "             X      TAN(X)\n" );
  printf ( "                   R4_TAN(X)         Diff\n" );
  printf ( "                   R8_TAN(X)         Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_tan ( ( float ) x );
    fx3 = r8_tan ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
/******************************************************************************/

void tanh_test ( )

/******************************************************************************/
/*
  Purpose:

    TANH_TEST tests R4_TANH and R8_TANH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TANH_TEST:\n" );
  printf ( "  Test TANH_VALUES, R4_TANH, R8_TANH\n" );
  printf ( "\n" );
  printf ( "             X      TANH(X)\n" );
  printf ( "                   R4_TANH(X)        Diff\n" );
  printf ( "                   R8_TANH(X)        Diff\n" );

  n_data = 0;

  for ( ; ; )
  {
    tanh_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_tanh ( ( float ) x );
    fx3 = r8_tanh ( x );

    printf ( "\n" );
    printf ( "  %14g  %14g\n", x, fx1 );
    printf ( "                  %14g  %14g\n", fx2, r4_abs ( ( float ) fx1 - fx2 ) );
    printf ( "                  %14g  %14g\n", fx3, r8_abs ( ( float ) fx1 - fx3 ) );
  }
  return;
}
