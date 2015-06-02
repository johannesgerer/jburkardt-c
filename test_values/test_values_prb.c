# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "test_values.h"

int main ( );

void abram0_values_test ( );
void abram1_values_test ( );
void abram2_values_test ( );
void agm_values_test ( );
void airy_ai_values_test ( );
void airy_ai_int_values_test ( );
void airy_ai_prime_values_test ( );
void airy_bi_values_test ( );
void airy_bi_int_values_test ( );
void airy_bi_prime_values_test ( );
void airy_cai_values_test ( );
void airy_cbi_values_test ( );
void airy_gi_values_test ( );
void airy_hi_values_test ( );
void arccos_values_test ( );
void arccosh_values_test ( );
void arcsin_values_test ( );
void arcsinh_values_test ( );
void arctan_values_test ( );
void arctan_int_values_test ( );
void arctan2_values_test ( );
void arctanh_values_test ( );

void bei0_values_test ( );
void bei1_values_test ( );
void bell_values_test ( );
void ber0_values_test ( );
void ber1_values_test ( );
void bernoulli_number_values_test ( );
void bernoulli_poly_values_test ( );
void bernstein_poly_01_values_test ( );
void bessel_i0_values_test ( );
void bessel_i0_int_values_test ( );
void bessel_i0_spherical_values_test ( );
void bessel_i1_values_test ( );
void bessel_i1_spherical_values_test ( );
void bessel_in_values_test ( );
void bessel_ix_values_test ( );
void bessel_j0_values_test ( );
void bessel_j0_int_values_test ( );
void bessel_j0_spherical_values_test ( );
void bessel_j1_values_test ( );
void bessel_j1_spherical_values_test ( );
void bessel_jn_values_test ( );
void bessel_jx_values_test ( );
void bessel_k0_values_test ( );
void bessel_k0_int_values_test ( );
void bessel_k1_values_test ( );
void bessel_kn_values_test ( );
void bessel_kx_values_test ( );
void bessel_y0_values_test ( );
void bessel_y0_int_values_test ( );
void bessel_y0_spherical_values_test ( );
void bessel_y1_values_test ( );
void bessel_y1_spherical_values_test ( );
void bessel_yn_values_test ( );
void bessel_yx_values_test ( );
void beta_cdf_values_test ( );
void beta_inc_values_test ( );
void beta_log_values_test ( );
void beta_noncentral_cdf_values_test ( );
void beta_values_test ( );
void binomial_values_test ( );
void binomial_cdf_values_test ( );
void bivariate_normal_cdf_values_test ( );

void catalan_values_test ( );
void cauchy_cdf_values_test ( );
void cbrt_values_test ( );
void cheby_t_poly_values_test ( );
void cheby_u_poly_values_test ( );
void cheby_v_poly_values_test ( );
void cheby_w_poly_values_test ( );
void chi_values_test ( );
void chi_square_cdf_values_test ( );
void chi_square_noncentral_cdf_values_test ( );
void ci_values_test ( );
void cin_values_test ( );
void cinh_values_test ( );
void clausen_values_test ( );
void clebsch_gordan_values_test ( );
void collatz_count_values_test ( );
void cos_values_test ( );
void cos_degree_values_test ( );
void cos_power_int_values_test ( );
void cosh_values_test ( );
void cot_values_test ( );
void cp_values_test ( );

void dawson_values_test ( );
void debye1_values_test ( );
void debye2_values_test ( );
void debye3_values_test ( );
void debye4_values_test ( );
void dielectric_values_test ( );
void dilogarithm_values_test ( );

void e1_values_test ( );
void ei_values_test ( );
void elliptic_ea_values_test ( );
void elliptic_em_values_test ( );
void elliptic_ka_values_test ( );
void elliptic_km_values_test ( );
void erf_values_test ( );
void erfc_values_test ( );
void euler_number_values_test ( );
void euler_poly_values_test ( );
void exp_values_test ( );
void exp3_int_values_test ( );
void exponential_cdf_values_test ( );
void extreme_values_cdf_values_test ( );

void f_cdf_values_test ( );
void f_noncentral_cdf_values_test ( );
void fresnel_cos_values_test ( );
void fresnel_sin_values_test ( );
void frobenius_number_data_values_test ( );
void frobenius_number_order_values_test ( );
void frobenius_number_order2_values_test ( );

void gamma_values_test ( );
void gamma_cdf_values_test ( );
void gamma_inc_values_test ( );
void gamma_inc_p_values_test ( );
void gamma_inc_q_values_test ( );
void gamma_inc_tricomi_values_test ( );
void gamma_log_values_test ( );
void gegenbauer_poly_values_test ( );
void geometric_cdf_values_test ( );
void goodwin_values_test ( );
void gud_values_test ( );

void hermite_function_values_test ( );
void hermite_poly_phys_values_test ( );
void hermite_poly_prob_values_test ( );
void hyper_1f1_values_test ( );
void hyper_2f1_values_test ( );
void hypergeometric_cdf_values_test ( );
void hypergeometric_pdf_values_test ( );
void hypergeometric_u_values_test ( );

void i0ml0_values_test ( );
void i1ml1_values_test ( );
void i4_factorial_values_test ( );
void i4_factorial2_values_test ( );
void i4_fall_values_test ( );
void i4_rise_values_test ( );
void int_values_test ( );

void jacobi_cn_values_test ( );
void jacobi_dn_values_test ( );
void jacobi_poly_values_test ( );
void jacobi_sn_values_test ( );
void jed_ce_values_test ( );
void jed_mjd_values_test ( );
void jed_rd_values_test ( );
void jed_weekday_values_test ( );

void kei0_values_test ( );
void kei1_values_test ( );
void ker0_values_test ( );
void ker1_values_test ( );

void laguerre_associated_values_test ( );
void laguerre_general_values_test ( );
void laguerre_polynomial_values_test ( );
void lambert_w_values_test ( );
void laplace_cdf_values_test ( );
void legendre_associated_values_test ( );
void legendre_associated_normalized_values_test ( );
void legendre_associated_normalized_sphere_values_test ( );
void legendre_poly_values_test ( );
void legendre_function_q_values_test ( );
void lerch_values_test ( );
void lobachevsky_values_test ( );
void lobatto_polynomial_values_test ( );
void lobatto_polynomial_derivatives_test ( );
void log_values_test ( );
void log_normal_cdf_values_test ( );
void log_series_cdf_values_test ( );
void log10_values_test ( );
void logarithmic_integral_values_test ( );
void logistic_cdf_values_test ( );

void mertens_values_test ( );
void moebius_values_test ( );

void negative_binomial_cdf_values_test ( );
void nine_j_values_test ( );
void normal_cdf_values_test ( );
void normal_01_cdf_values_test ( );

void omega_values_test ( );
void owen_values_test ( );

void partition_count_values_test ( );
void partition_distinct_count_values_test ( );
void phi_values_test ( );
void pi_values_test ( );
void poisson_cdf_values_test ( );
void polylogarithm_values_test ( );
void prandtl_values_test ( );
void prime_values_test ( );
void psat_values_test ( );
void psi_values_test ( );

void r8_factorial_values_test ( );
void r8_factorial_log_values_test ( );
void r8_factorial2_values_test ( );
void r8_fall_values_test ( );
void r8_rise_values_test ( );
void rayleigh_cdf_values_test ( );

void secvir_values_test ( );
void shi_values_test ( );
void si_values_test ( );
void sigma_values_test ( );
void sin_values_test ( );
void sin_degree_values_test ( );
void sin_power_int_values_test ( );
void sinh_values_test ( );
void six_j_values_test ( );
void sound_values_test ( );
void sphere_unit_area_values_test ( );
void sphere_unit_volume_values_test ( );
void spherical_harmonic_values_test ( );
void sqrt_values_test ( );
void stirling1_values_test ( );
void stirling2_values_test ( );
void stromgen_values_test ( );
void struve_h0_values_test ( );
void struve_h1_values_test ( );
void struve_l0_values_test ( );
void struve_l1_values_test ( );
void student_cdf_values_test ( );
void student_noncentral_cdf_values_test ( );
void subfactorial_values_test ( );
void surten_values_test ( );
void synch1_values_test ( );
void synch2_values_test ( );

void tan_values_test ( );
void tanh_values_test ( );
void tau_values_test ( );
void thercon_values_test ( );
void three_j_values_test ( );
void tran02_values_test ( );
void tran03_values_test ( );
void tran04_values_test ( );
void tran05_values_test ( );
void tran06_values_test ( );
void tran07_values_test ( );
void tran08_values_test ( );
void tran09_values_test ( );
void trigamma_values_test ( );
void truncated_normal_ab_cdf_values_test ( );
void truncated_normal_ab_pdf_values_test ( );
void truncated_normal_a_cdf_values_test ( );
void truncated_normal_a_pdf_values_test ( );
void truncated_normal_b_cdf_values_test ( );
void truncated_normal_b_pdf_values_test ( );
void tsat_values_test ( );

void van_der_corput_values_test ( );
void viscosity_values_test ( );
void von_mises_cdf_values_test ( );

void weekday_values_test ( );
void weibull_cdf_values_test ( );

void zeta_values_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_VALUES_PRB.

  Discussion:

    TEST_VALUES_PRB tests the TEST_VALUE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 February 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_VALUES_PRB:\n" );
  printf ( "  C version,\n" );
  printf ( "  Test the TEST_VALUES library.\n" );

  abram0_values_test ( );
  abram1_values_test ( );
  abram2_values_test ( );
  agm_values_test ( );
  airy_ai_values_test ( );
  airy_ai_int_values_test ( );
  airy_ai_prime_values_test ( );
  airy_bi_values_test ( );
  airy_bi_int_values_test ( );
  airy_bi_prime_values_test ( );
  airy_cai_values_test ( );
  airy_cbi_values_test ( );
  airy_gi_values_test ( );
  airy_hi_values_test ( );
  arccos_values_test ( );
  arccosh_values_test ( );
  arcsin_values_test ( );
  arcsinh_values_test ( );
  arctan_values_test ( );
  arctan_int_values_test ( );
  arctan2_values_test ( );
  arctanh_values_test ( );

  bei0_values_test ( );
  bei1_values_test ( );
  bell_values_test ( );
  ber0_values_test ( );
  ber1_values_test ( );
  bernoulli_number_values_test ( );
  bernoulli_poly_values_test ( );
  bernstein_poly_01_values_test ( );
  bessel_i0_values_test ( );
  bessel_i0_int_values_test ( );
  bessel_i0_spherical_values_test ( );
  bessel_i1_values_test ( );
  bessel_i1_spherical_values_test ( );
  bessel_in_values_test ( );
  bessel_ix_values_test ( );
  bessel_j0_values_test ( );
  bessel_j0_int_values_test ( );
  bessel_j0_spherical_values_test ( );
  bessel_j1_values_test ( );
  bessel_j1_spherical_values_test ( );
  bessel_jn_values_test ( );
  bessel_jx_values_test ( );
  bessel_k0_values_test ( );
  bessel_k0_int_values_test ( );
  bessel_k1_values_test ( );
  bessel_kn_values_test ( );
  bessel_kx_values_test ( );
  bessel_y0_values_test ( );
  bessel_y0_int_values_test ( );
  bessel_y0_spherical_values_test ( );
  bessel_y1_values_test ( );
  bessel_y1_spherical_values_test ( );
  bessel_yn_values_test ( );
  bessel_yx_values_test ( );
  beta_cdf_values_test ( );
  beta_inc_values_test ( );
  beta_log_values_test ( );
  beta_noncentral_cdf_values_test ( );
  beta_values_test ( );
  binomial_values_test ( );
  binomial_cdf_values_test ( );
  bivariate_normal_cdf_values_test ( );

  catalan_values_test ( );
  cauchy_cdf_values_test ( );
  cbrt_values_test ( );
  cheby_t_poly_values_test ( );
  cheby_u_poly_values_test ( );
  cheby_v_poly_values_test ( );
  cheby_w_poly_values_test ( );
  chi_values_test ( );
  chi_square_cdf_values_test ( );
  chi_square_noncentral_cdf_values_test ( );
  ci_values_test ( );
  cin_values_test ( );
  cinh_values_test ( );
  clausen_values_test ( );
  clebsch_gordan_values_test ( );
  collatz_count_values_test ( );
  cos_values_test ( );
  cos_degree_values_test ( );
  cos_power_int_values_test ( );
  cosh_values_test ( );
  cot_values_test ( );
  cp_values_test ( );

  dawson_values_test ( );
  debye1_values_test ( );
  debye2_values_test ( );
  debye3_values_test ( );
  debye4_values_test ( );
  dielectric_values_test ( );
  dilogarithm_values_test ( );

  e1_values_test ( );
  ei_values_test ( );
  elliptic_ea_values_test ( );
  elliptic_em_values_test ( );
  elliptic_ka_values_test ( );
  elliptic_km_values_test ( );
  erf_values_test ( );
  erfc_values_test ( );
  euler_number_values_test ( );
  euler_poly_values_test ( );
  exp_values_test ( );
  exp3_int_values_test ( );
  exponential_cdf_values_test ( );
  extreme_values_cdf_values_test ( );

  f_cdf_values_test ( );
  f_noncentral_cdf_values_test ( );
  fresnel_cos_values_test ( );
  fresnel_sin_values_test ( );
  frobenius_number_data_values_test ( );
  frobenius_number_order_values_test ( );
  frobenius_number_order2_values_test ( );

  gamma_values_test ( );
  gamma_cdf_values_test ( );
  gamma_inc_values_test ( );
  gamma_inc_p_values_test ( );
  gamma_inc_q_values_test ( );
  gamma_inc_tricomi_values_test ( );
  gamma_log_values_test ( );
  gegenbauer_poly_values_test ( );
  geometric_cdf_values_test ( );
  goodwin_values_test ( );
  gud_values_test ( );

  hermite_function_values_test ( );
  hermite_poly_phys_values_test ( );
  hermite_poly_prob_values_test ( );
  hyper_1f1_values_test ( );
  hyper_2f1_values_test ( );
  hypergeometric_cdf_values_test ( );
  hypergeometric_pdf_values_test ( );
  hypergeometric_u_values_test ( );

  i0ml0_values_test ( );
  i1ml1_values_test ( );
  i4_factorial_values_test ( );
  i4_factorial2_values_test ( );
  i4_fall_values_test ( );
  i4_rise_values_test ( );
  int_values_test ( );

  jacobi_cn_values_test ( );
  jacobi_dn_values_test ( );
  jacobi_poly_values_test ( );
  jacobi_sn_values_test ( );
  jed_ce_values_test ( );
  jed_mjd_values_test ( );
  jed_rd_values_test ( );
  jed_weekday_values_test ( );

  kei0_values_test ( );
  kei1_values_test ( );
  ker0_values_test ( );
  ker1_values_test ( );

  laguerre_associated_values_test ( );
  laguerre_general_values_test ( );
  laguerre_polynomial_values_test ( );
  lambert_w_values_test ( );
  laplace_cdf_values_test ( );
  legendre_associated_values_test ( );
  legendre_associated_normalized_values_test ( );
  legendre_associated_normalized_sphere_values_test ( );
  legendre_poly_values_test ( );
  legendre_function_q_values_test ( );
  lerch_values_test ( );
  lobachevsky_values_test ( );
  lobatto_polynomial_values_test ( );
  lobatto_polynomial_derivatives_test ( );
  log_values_test ( );
  log_normal_cdf_values_test ( );
  log_series_cdf_values_test ( );
  log10_values_test ( );
  logarithmic_integral_values_test ( );
  logistic_cdf_values_test ( );

  mertens_values_test ( );
  moebius_values_test ( );

  negative_binomial_cdf_values_test ( );
  nine_j_values_test ( );
  normal_cdf_values_test ( );
  normal_01_cdf_values_test ( );

  omega_values_test ( );
  owen_values_test ( );

  partition_count_values_test ( );
  partition_distinct_count_values_test ( );
  phi_values_test ( );
  pi_values_test ( );
  poisson_cdf_values_test ( );
  polylogarithm_values_test ( );
  prandtl_values_test ( );
  prime_values_test ( );
  psat_values_test ( );
  psi_values_test ( );

  r8_factorial_values_test ( );
  r8_factorial_log_values_test ( );
  r8_factorial2_values_test ( );
  r8_fall_values_test ( );
  r8_rise_values_test ( );
  rayleigh_cdf_values_test ( );

  secvir_values_test ( );
  shi_values_test ( );
  si_values_test ( );
  sigma_values_test ( );
  sin_values_test ( );
  sin_degree_values_test ( );
  sin_power_int_values_test ( );
  sinh_values_test ( );
  six_j_values_test ( );
  sound_values_test ( );
  sphere_unit_area_values_test ( );
  sphere_unit_volume_values_test ( );
  spherical_harmonic_values_test ( );
  sqrt_values_test ( );
  stirling1_values_test ( );
  stirling2_values_test ( );
  stromgen_values_test ( );
  struve_h0_values_test ( );
  struve_h1_values_test ( );
  struve_l0_values_test ( );
  struve_l1_values_test ( );
  student_cdf_values_test ( );
  student_noncentral_cdf_values_test ( );
  subfactorial_values_test ( );
  surten_values_test ( );
  synch1_values_test ( );
  synch2_values_test ( );

  tan_values_test ( );
  tanh_values_test ( );
  tau_values_test ( );
  thercon_values_test ( );
  three_j_values_test ( );
  tran02_values_test ( );
  tran03_values_test ( );
  tran04_values_test( );
  tran05_values_test ( );
  tran06_values_test ( );
  tran07_values_test ( );
  tran08_values_test ( );
  tran09_values_test ( );
  trigamma_values_test ( );
  truncated_normal_ab_cdf_values_test ( );
  truncated_normal_ab_pdf_values_test ( );
  truncated_normal_a_cdf_values_test ( );
  truncated_normal_a_pdf_values_test ( );
  truncated_normal_b_cdf_values_test ( );
  truncated_normal_b_pdf_values_test ( );
  tsat_values_test ( );

  van_der_corput_values_test ( );
  viscosity_values_test ( );
  von_mises_cdf_values_test ( );

  weekday_values_test ( );
  weibull_cdf_values_test ( );

  zeta_values_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_VALUES_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void abram0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ABRAM0_VALUES_TEST tests ABRAM0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ABRAM0_VALUES_TEST:\n" );
  printf ( "  ABRAM0_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void abram1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ABRAM1_VALUES_TEST tests ABRAM1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ABRAM1_VALUES_TEST:\n" );
  printf ( "  ABRAM1_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void abram2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ABRAM2_VALUES_TEST tests ABRAM2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ABRAM2_VALUES_TEST:\n" );
  printf ( "  ABRAM2_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM3(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void agm_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AGM_VALUES_TEST tests AGM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;

  printf ( "\n" );
  printf ( "AGM_VALUES_TEST:\n" );
  printf ( "  AGM_VALUES stores values of \n" );
  printf ( "  the arithmetic geometric mean function.\n" );
  printf ( "\n" );
  printf ( "           A          B              AGM(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( &n_data, &a, &b, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %24.16e\n", a, b, fx );
  }
  return;
}
/******************************************************************************/

void airy_ai_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_AI_VALUES_TEST tests AIRY_AI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double ai;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_AI_VALUES_TEST:\n" );
  printf ( "  AIRY_AI_VALUES stores values of \n" );
  printf ( "  the Airy functions Ai(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Ai(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_values ( &n_data, &x, &ai );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, ai );
  }
  return;
}
/******************************************************************************/

void airy_ai_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_AI_INT_VALUES_TEST tests AIRY_AI_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_AI_INT_VALUES_TEST:\n" );
  printf ( "  AIRY_AI_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Airy Ai function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void airy_ai_prime_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_AI_PRIME_VALUES_TEST tests AIRY_AI_PRIME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double aip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_AI_PRIME_VALUES_TEST:\n" );
  printf ( "  AIRY_AI_PRIME_VALUES stores values of \n" );
  printf ( "  the derivative of the Airy function Ai'(X).\n" );
  printf ( "\n" );
  printf ( "                X                    Ai'\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_prime_values ( &n_data, &x, &aip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, aip );
  }
  return;
}
/******************************************************************************/

void airy_bi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_BI_VALUES_TEST tests AIRY_BI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bi;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_BI_VALUES_TEST:\n" );
  printf ( "  AIRY_BI_VALUES stores values of \n" );
  printf ( "  the Airy function Bi.\n" );
  printf ( "\n" );
  printf ( "                X                     Bi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_values ( &n_data, &x, &bi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bi );
  }
  return;
}
/******************************************************************************/

void airy_bi_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_BI_INT_VALUES_TEST tests AIRY_BI_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_BI_INT_VALUES_TEST:\n" );
  printf ( "  AIRY_BI_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Airy Bi function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void airy_bi_prime_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_BI_PRIME_VALUES_TEST tests AIRY_BI_PRIME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_BI_PRIME_VALUES_TEST:\n" );
  printf ( "  AIRY_BI_PRIME_VALUES stores values of \n" );
  printf ( "  the derivative of Airy function Bi'(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Bi'\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_prime_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void airy_cai_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_CAI_VALUES_TEST tests AIRY_CAI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double complex cai;
  int n_data;
  double complex x;

  printf ( "\n" );
  printf ( "AIRY_CAI_VALUES_TEST:\n" );
  printf ( "  AIRY_CAI_VALUES stores values of \n" );
  printf ( "  the Airy functions Ai(X) for complex argument.\n" );
  printf ( "\n" );
  printf ( "                X                     Ai\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_cai_values ( &n_data, &x, &cai );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %14.6f  %14.6f\n", x, cai );
  }
  return;
}
/******************************************************************************/

void airy_cbi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_CBI_VALUES_TEST tests AIRY_CBI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double complex cbi;
  int n_data;
  double complex x;

  printf ( "\n" );
  printf ( "AIRY_CBI_VALUES_TEST:\n" );
  printf ( "  AIRY_CBI_VALUES stores values of \n" );
  printf ( "  the Airy functions Bi(X) for complex argument.\n" );
  printf ( "\n" );
  printf ( "                X                     Bi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_cbi_values ( &n_data, &x, &cbi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %14.6f  %14.6f\n", x, cbi );
  }
  return;
}
/******************************************************************************/

void airy_gi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_GI_VALUES_TEST tests AIRY_GI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_GI_VALUES_TEST:\n" );
  printf ( "  AIRY_GI_VALUES stores values of \n" );
  printf ( "  the modified Airy function Gi(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Gi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_gi_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void airy_hi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    AIRY_HI_VALUES_TEST tests AIRY_HI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "AIRY_HI_VALUES_TEST:\n" );
  printf ( "  AIRY_HI_VALUES stores values of \n" );
  printf ( "  the modified Airy function Hi(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Hi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_hi_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void arccos_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCCOS_VALUES_TEST tests ARCCOS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCCOS_VALUES_TEST:\n" );
  printf ( "  ARCCOS_VALUES stores values of the arc cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void arccosh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCCOSH_VALUES_TEST tests ARCCOSH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCCOSH_VALUES_TEST:\n" );
  printf ( "  ARCCOSH_VALUES stores values of the hyperbolic arc cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccosh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void arcsin_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCSIN_VALUES_TEST tests ARCSIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCSIN_VALUES_TEST:\n" );
  printf ( "  ARCSIN_VALUES stores values of the arc sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void arcsinh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCSINH_VALUES_TEST tests ARCSINH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCSINH_VALUES_TEST:\n" );
  printf ( "  ARCSINH_VALUES stores values of the hyperbolic arc sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsinh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void arctan_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCTAN_VALUES_TEST tests ARCTAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCTAN_VALUES_TEST:\n" );
  printf ( "  ARCTAN_VALUES stores values of the arc tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void arctan_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCTAN_INT_VALUES_TEST tests ARCTAN_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCTAN_INT_VALUES_TEST:\n" );
  printf ( "  ARCTAN_INT_VALUES stores values of \n" );
  printf ( "  the arctangent integral.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void arctan2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCTAN2_VALUES_TEST tests ARCTAN2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 February 2015

  Author:

    John Burkardt
*/
{
  double f;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "ARCTAN2_VALUES_TEST:\n" );
  printf ( "  ARCTAN2_VALUES stores values of the arc tangent function.\n" );
  printf ( "\n" );
  printf ( "        X             Y           F(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan2_values ( &n_data, &x, &y, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12.6f  %12.6f  %24.16e\n", x, y, f );
  }
  return;
}
/******************************************************************************/

void arctanh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ARCTANH_VALUES_TEST tests ARCTANH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ARCTANH_VALUES_TEST:\n" );
  printf ( "  ARCTANH_VALUES stores values of the hyperbolic arc tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctanh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bei0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BEI0_VALUES_TEST tests BEI0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BEI0_VALUES_TEST:\n" );
  printf ( "  BEI0_VALUES stores values of \n" );
  printf ( "  the Kelvin function BEI of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bei0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bei1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BEI1_VALUES_TEST tests BEI1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BEI1_VALUES_TEST:\n" );
  printf ( "  BEI1_VALUES stores values of \n" );
  printf ( "  the Kelvin function BEI of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bei1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bell_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BELL_VALUES_TEST tests BELL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BELL_VALUES_TEST:\n" );
  printf ( "  BELL_VALUES returns values of \n" );
  printf ( "  the Bell numbers.\n" );
  printf ( "\n" );
  printf ( "     N        BELL(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void ber0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BER0_VALUES_TEST tests BER0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BER0_VALUES_TEST:\n" );
  printf ( "  BER0_VALUES stores values of \n" );
  printf ( "  the Kelvin function BER of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ber0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void ber1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BER1_VALUES_TEST tests BER1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BER1_VALUES_TEST:\n" );
  printf ( "  BER1_VALUES stores values of \n" );
  printf ( "  the Kelvin function BER of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ber1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bernoulli_number_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_NUMBER_VALUES_TEST tests BERNOULLI_NUMBER_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BERNOULLI_NUMBER_VALUES_TEST:\n" );
  printf ( "  BERNOULLI_NUMBER_VALUES returns values of \n" );
  printf ( "  the Bernoulli numbers.\n" );
  printf ( "\n" );
  printf ( "     N              B(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12g\n", n, c );
  }
  return;
}
/******************************************************************************/

void bernoulli_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_POLY_VALUES_TEST tests BERNOULLI_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double b;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BERNOULLI_POLY_VALUES_TEST:\n" );
  printf ( "  BERNOULLI_POLY_VALUES returns values of \n" );
  printf ( "  the Bernoulli Polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     X      BERNOULLI(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_poly_values ( &n_data, &n, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, b );
  }
  return;
}
/******************************************************************************/

void bernstein_poly_01_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNSTEIN_POLY_01_VALUES_TEST tests BERNSTEIN_POLY_01_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double b;
  int k;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BERNSTEIN_POLY_01_VALUES_TEST:\n" );
  printf ( "  BERNSTEIN_POLY_01_VALUES returns values of \n" );
  printf ( "  the Bernstein Polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     K       X      BERNSTEIN(N,K)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_01_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %12e\n", n, k, x, b );
  }
  return;
}
/******************************************************************************/

void bessel_i0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_I0_VALUES_TEST tests BESSEL_I0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_I0_VALUES_TEST:\n" );
  printf ( "  BESSEL_I0_VALUES stores values of \n" );
  printf ( "  the Bessel I0 function.\n" );
  printf ( "\n" );
  printf ( "      X         I0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_i0_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_I0_INT_VALUES_TEST tests BESSEL_I0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_I0_INT_VALUES_TEST:\n" );
  printf ( "  BESSEL_I0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel I0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void bessel_i0_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_I0_SPHERICAL_VALUES_TEST tests BESSEL_I0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_I0_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_I0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel i0 function.\n" );
  printf ( "\n" );
  printf ( "      X            i0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_i0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_i1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_I1_VALUES_TEST tests BESSEL_I1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_I1_VALUES_TEST:\n" );
  printf ( "  BESSEL_I1_VALUES stores values of \n" );
  printf ( "  the Bessel I1 function.\n" );
  printf ( "\n" );
  printf ( "      X         I1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_i1_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_I1_SPHERICAL_VALUES_TEST tests BESSEL_I1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_I1_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_I1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel i1 function.\n" );
  printf ( "\n" );
  printf ( "      X            i1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_i1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_in_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_IN_VALUES_TEST tests BESSEL_IN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_IN_VALUES_TEST:\n" );
  printf ( "  BESSEL_IN_VALUES stores values of \n" );
  printf ( "  the Bessel In function.\n" );
  printf ( "\n" );
  printf ( "      N     X         IN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_in_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_ix_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_IX_VALUES_TEST tests BESSEL_IX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_IX_VALUES_TEST:\n" );
  printf ( "  BESSEL_IX_VALUES stores values of \n" );
  printf ( "  the Bessel In function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU      X         IN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_ix_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_j0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J0_VALUES_TEST tests BESSEL_J0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_J0_VALUES_TEST:\n" );
  printf ( "  BESSEL_J0_VALUES stores values of \n" );
  printf ( "  the Bessel J0 function.\n" );
  printf ( "\n" );
  printf ( "      X         J0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_j0_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J0_INT_VALUES_TEST tests BESSEL_J0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_J0_INT_VALUES_TEST:\n" );
  printf ( "  BESSEL_J0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel J0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void bessel_j0_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J0_SPHERICAL_VALUES_TEST tests BESSEL_J0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_J0_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_J0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel j0 function.\n" );
  printf ( "\n" );
  printf ( "      X            j0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_j0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_j1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J1_VALUES_TEST tests BESSEL_J1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_J1_VALUES_TEST:\n" );
  printf ( "  BESSEL_J1_VALUES stores values of \n" );
  printf ( "  the Bessel J1 function.\n" );
  printf ( "\n" );
  printf ( "      X         J1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_j1_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_J1_SPHERICAL_VALUES_TEST tests BESSEL_J1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_J1_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_J1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel j1 function.\n" );
  printf ( "\n" );
  printf ( "      X            j1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_j1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_jn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_JN_VALUES_TEST tests BESSEL_JN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_JN_VALUES_TEST:\n" );
  printf ( "  BESSEL_JN_VALUES stores values of \n" );
  printf ( "  the Bessel Jn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         JN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_jn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_jx_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_JX_VALUES_TEST tests BESSEL_JX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_JX_VALUES_TEST:\n" );
  printf ( "  BESSEL_JX_VALUES stores values of \n" );
  printf ( "  the Bessel Jn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU        X         JN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_jx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_k0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_K0_VALUES_TEST tests BESSEL_K0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_K0_VALUES_TEST:\n" );
  printf ( "  BESSEL_K0_VALUES stores values of \n" );
  printf ( "  the Bessel K0 function.\n" );
  printf ( "\n" );
  printf ( "      X         K0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_k0_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_K0_INT_VALUES_TEST tests BESSEL_K0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_K0_INT_VALUES_TEST:\n" );
  printf ( "  BESSEL_K0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel K0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void bessel_k1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_K1_VALUES_TEST tests BESSEL_K1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_K1_VALUES_TEST:\n" );
  printf ( "  BESSEL_K1_VALUES stores values of \n" );
  printf ( "  the Bessel K1 function.\n" );
  printf ( "\n" );
  printf ( "      X         K1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_kn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_KN_VALUES_TEST tests BESSEL_KN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_KN_VALUES_TEST:\n" );
  printf ( "  BESSEL_KN_VALUES stores values of \n" );
  printf ( "  the Bessel Kn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         KN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_kn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_kx_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_KX_VALUES_TEST tests BESSEL_KX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_KX_VALUES_TEST:\n" );
  printf ( "  BESSEL_KX_VALUES stores values of \n" );
  printf ( "  the Bessel Kn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU      X         KN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_kx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_y0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_Y0_VALUES_TEST tests BESSEL_Y0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_Y0_VALUES_TEST:\n" );
  printf ( "  BESSEL_Y0_VALUES stores values of \n" );
  printf ( "  the Bessel Y0 function.\n" );
  printf ( "\n" );
  printf ( "      X         Y0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_y0_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_Y0_INT_VALUES_TEST tests BESSEL_Y0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_Y0_INT_VALUES_TEST:\n" );
  printf ( "  BESSEL_Y0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel Y0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_y0_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_Y0_SPHERICAL_VALUES_TEST tests BESSEL_Y0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_Y0_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_Y0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel y0 function.\n" );
  printf ( "\n" );
  printf ( "                X                      y0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_y0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_y1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_Y1_VALUES_TEST tests BESSEL_Y1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_Y1_VALUES_TEST:\n" );
  printf ( "  BESSEL_Y1_VALUES stores values of \n" );
  printf ( "  the Bessel Y1 function.\n" );
  printf ( "\n" );
  printf ( "                X                   Y1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_y1_spherical_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_Y1_SPHERICAL_VALUES_TEST tests BESSEL_Y1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_Y1_SPHERICAL_VALUES_TEST:\n" );
  printf ( "  BESSEL_Y1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel y1 function.\n" );
  printf ( "\n" );
  printf ( "                X                      y1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_y1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_yn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_YN_VALUES_TEST tests BESSEL_YN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_YN_VALUES_TEST:\n" );
  printf ( "  BESSEL_YN_VALUES stores values of \n" );
  printf ( "  the Bessel Yn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         YN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_yn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void bessel_yx_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BESSEL_YX_VALUES_TEST tests BESSEL_YX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BESSEL_YX_VALUES_TEST:\n" );
  printf ( "  BESSEL_YX_VALUES stores values of \n" );
  printf ( "  the Bessel Yn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU    X         YN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_yx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void beta_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_CDF_VALUES_TEST tests BETA_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BETA_CDF_VALUES_TEST:\n" );
  printf ( "  BETA_CDF_VALUES stores values of\n" );
  printf ( "  the Beta CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void beta_inc_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_INC_VALUES_TEST tests BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BETA_INC_VALUES_TEST:\n" );
  printf ( "  BETA_INC_VALUES stores values of\n" );
  printf ( "  the incomplete Beta function.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            BETA_INC(A,B)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void beta_log_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_LOG_VALUES_TEST tests BETA_LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "BETA_LOG_VALUES_TEST:\n" );
  printf ( "  BETA_LOG_VALUES stores values of\n" );
  printf ( "  the logarithm of the Beta function.\n" );
  printf ( "\n" );
  printf ( "      X              Y         BETA_LOG(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_log_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", x, y, fxy );
  }
  return;
}
/******************************************************************************/

void beta_noncentral_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_NONCENTRAL_CDF_VALUES_TEST tests BETA_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BETA_NONCENTRAL_CDF_VALUES_TEST:\n" );
  printf ( "  BETA_NONCENTRAL_CDF_VALUES stores values of\n" );
  printf ( "  the noncentral Beta CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B       LAMBDA             X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %12f  %24.16e\n", a, b, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void beta_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BETA_VALUES_TEST tests BETA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "BETA_VALUES_TEST:\n" );
  printf ( "  BETA_VALUES stores values of\n" );
  printf ( "  the Beta function.\n" );
  printf ( "\n" );
  printf ( "      X              Y         BETA(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", x,  y, fxy );
  }
  return;
}
/******************************************************************************/

void binomial_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BINOMIAL_VALUES_TEST tests BINOMIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int n_data;

  printf ( "\n" );
  printf ( "BINOMIAL_VALUES_TEST:\n" );
  printf ( "  BINOMIAL_VALUES returns values of\n" );
  printf ( "  the binomial numbers.\n" );
  printf ( "\n" );
  printf ( "     A     B        C(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_values ( &n_data, &a, &b, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", a, b, c );
  }
  return;
}
/******************************************************************************/

void binomial_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BINOMIAL_CDF_VALUES_TEST tests BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int a;
  double b;
  double fx;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "BINOMIAL_CDF_VALUES_TEST:\n" );
  printf ( "  BINOMIAL_CDF_VALUES returns values of \n" );
  printf ( "  the Binomial Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     A      B        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %4d  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void bivariate_normal_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    BIVARIATE_NORMAL_CDF_VALUES_TEST tests BIVARIATE_NORMAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double r;
  double x;
  double y;

  printf ( "\n" );
  printf ( "BIVARIATE_NORMAL_CDF_VALUES_TEST:\n" );
  printf ( "  BIVARIATE_NORMAL_CDF_VALUES stores values of\n" );
  printf ( "  the bivariate normal CDF.\n" );
  printf ( "\n" );
  printf ( "      X            Y            R            F(R)(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( &n_data, &x, &y, &r, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", x, y, r, fxy );
  }
  return;
}
/******************************************************************************/

void catalan_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CATALAN_VALUES_TEST tests CATALAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "CATALAN_VALUES_TEST:\n" );
  printf ( "  CATALAN_VALUES returns values of \n" );
  printf ( "  the Catalan numbers.\n" );
  printf ( "\n" );
  printf ( "     N        C(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void cauchy_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CAUCHY_CDF_VALUES_TEST tests CAUCHY_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "CAUCHY_CDF_VALUES_TEST:\n" );
  printf ( "  CAUCHY_CDF_VALUES returns values of \n" );
  printf ( "  the Cauchy Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Sigma        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cauchy_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void cbrt_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CBRT_VALUES_TEST tests CBRT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CBRT_VALUES_TEST:\n" );
  printf ( "  CBRT_VALUES stores values of the cube root function.\n" );
  printf ( "\n" );
  printf ( "      X            CBRT(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cheby_t_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_T_POLY_VALUES_TEST tests CHEBY_T_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHEBY_T_POLY_VALUES_TEST:\n" );
  printf ( "  CHEBY_T_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev T polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      T(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void cheby_u_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_U_POLY_VALUES_TEST tests CHEBY_U_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHEBY_U_POLY_VALUES_TEST:\n" );
  printf ( "  CHEBY_U_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev U polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      U(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void cheby_v_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_V_POLY_VALUES_TEST tests CHEBY_V_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHEBY_V_POLY_VALUES_TEST:\n" );
  printf ( "  CHEBY_V_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev V polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      V(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_v_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void cheby_w_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_W_POLY_VALUES_TEST tests CHEBY_W_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHEBY_W_POLY_VALUES_TEST:\n" );
  printf ( "  CHEBY_W_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev W polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      W(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_w_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void chi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHI_VALUES_TEST tests CHI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHI_VALUES_TEST:\n" );
  printf ( "  CHI_VALUES stores values of\n" );
  printf ( "  the Hyperbolic Cosine Integral function CHI(X).\n" );
  printf ( "\n" );
  printf ( "      X            CHI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void chi_square_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHI_SQUARE_CDF_VALUES_TEST tests CHI_SQUARE_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHI_SQUARE_CDF_VALUES_TEST:\n" );
  printf ( "  CHI_SQUARE_CDF_VALUES returns values of \n" );
  printf ( "  the Chi-Squared Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     N       X    CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %12e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void chi_square_noncentral_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST tests CHI_SQUARE_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST:\n" );
  printf ( "  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of\n" );
  printf ( "  the noncentral Chi-Squared Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      X      LAMBDA     DF     CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_noncentral_cdf_values ( &n_data, &df, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %10f  %8f  %4d  %12e\n", x, lambda, df, fx );
  }
  return;
}
/******************************************************************************/

void ci_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CI_VALUES_TEST tests CI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CI_VALUES_TEST:\n" );
  printf ( "  CI_VALUES stores values of\n" );
  printf ( "  the Cosine Integral function CI(X).\n" );
  printf ( "\n" );
  printf ( "      X            CI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ci_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cin_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CIN_VALUES_TEST tests CIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CIN_VALUES_TEST:\n" );
  printf ( "  CIN_VALUES stores values of\n" );
  printf ( "  the Cosine Integral function CIN(X).\n" );
  printf ( "\n" );
  printf ( "      X            CIN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cinh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CINH_VALUES_TEST tests CINH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CINH_VALUES_TEST:\n" );
  printf ( "  CINH_VALUES stores values of\n" );
  printf ( "  the Hyperbolic Cosine Integral function CINH(X).\n" );
  printf ( "\n" );
  printf ( "      X            CINH(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cinh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void clausen_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CLAUSEN_VALUES_TEST tests CLAUSEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "CLAUSEN_VALUES_TEST:\n" );
  printf ( "  CLAUSEN_VALUES stores values of \n" );
  printf ( "  Clausen's integral function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    clausen_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void clebsch_gordan_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CLEBSCH_GORDAN_VALUES_TEST tests CLEBSCH_GORDAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  printf ( "\n" );
  printf ( "CLEBSCH_GORDAN_VALUES_TEST:\n" );
  printf ( "  CLEBSCH_GORDAN_VALUES returns values of\n" );
  printf ( "  the Clebsch Gordan coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      M1      M2      M3        CG\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    clebsch_gordan_values ( &n_data, &j1, &j2, &j3, &m1, &m2, &m3, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, m1, m2, m3, fx );
  }
  return;
}
/******************************************************************************/

void collatz_count_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COLLATZ_COUNT_VALUES_TEST tests COLLATZ_COUNT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 March 2006

  Author:

    John Burkardt
*/
{
  int count;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "COLLATZ_COUNT_VALUES_TEST:\n" );
  printf ( "  COLLATZ_COUNT_VALUES returns values of\n" );
  printf ( "  the length of the Collatz sequence that\n" );
  printf ( "  starts at N.\n" );
  printf ( "\n" );
  printf ( "         N      COLLATZ_COUNT(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %12d\n", n, count );
  }

  return;
}
/******************************************************************************/

void cos_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_VALUES_TEST tests COS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COS_VALUES_TEST:\n" );
  printf ( "   COS_VALUES stores values of the cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cos_degree_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_DEGREE_VALUES_TEST tests COS_DEGREE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COS_DEGREE_VALUES_TEST:\n" );
  printf ( "   COS_DEGREE_VALUES stores values of the cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_degree_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cos_power_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_POWER_INT_VALUES_TEST tests COS_POWER_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "COS_POWER_INT_VALUES_TEST:\n" );
  printf ( "  COS_POWER_INT_VALUES returns values of\n" );
  printf ( "  the integral of the N-th power of the cosine function.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   cos_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %6d  %24.16e\n", a, b, n, fx );
  }
  return;
}
/******************************************************************************/

void cosh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COSH_VALUES_TEST tests COSH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COSH_VALUES_TEST:\n" );
  printf ( "   COSH_VALUES stores values of the hyperbolic cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cosh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cot_values_test ( )

/******************************************************************************/
/*
  Purpose:

    COT_VALUES_TEST tests COT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "COT_VALUES_TEST:\n" );
  printf ( "   COT_VALUES stores values of the cotangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cot_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void cp_values_test ( )

/******************************************************************************/
/*
  Purpose:

    CP_VALUES_TEST tests CP_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double cp;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "CP_VALUES_TEST:\n" );
  printf ( "  CP_VALUES stores values of\n" );
  printf ( "  the specific heat CP\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            CP(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cp_values ( &n_data, &tc, &p, &cp );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, cp );
  }
  return;
}
/******************************************************************************/

void dawson_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DAWSON_VALUES_TEST tests DAWSON_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DAWSON_VALUES_TEST:\n" );
  printf ( "  DAWSON_VALUES stores values of\n" );
  printf ( "  Dawson's integral function.\n" );
  printf ( "\n" );
  printf ( "      X          DAWSON(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void debye1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DEBYE1_VALUES_TEST tests DEBYE1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DEBYE1_VALUES_TEST:\n" );
  printf ( "  DEBYE1_VALUES stores values of \n" );
  printf ( "  the Debye function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void debye2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DEBYE2_VALUES_TEST tests DEBYE2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DEBYE2_VALUES_TEST:\n" );
  printf ( "  DEBYE2_VALUES stores values of \n" );
  printf ( "  the Debye function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void debye3_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DEBYE3_VALUES_TEST tests DEBYE3_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DEBYE3_VALUES_TEST:\n" );
  printf ( "  DEBYE3_VALUES stores values of \n" );
  printf ( "  the Debye function of order 3.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye3_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void debye4_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DEBYE4_VALUES_TEST tests DEBYE4_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DEBYE4_VALUES_TEST:\n" );
  printf ( "  DEBYE4_VALUES stores values of \n" );
  printf ( "  the Debye function of order 4.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye4_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void dielectric_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DIELECTRIC_VALUES_TEST tests DIELECTRIC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double eps;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "DIELECTRIC_VALUES_TEST:\n" );
  printf ( "  DIELECTRIC_VALUES stores values of\n" );
  printf ( "  the dielectric function.\n" );
  printf ( "\n" );
  printf ( "      T           P            EPS(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dielectric_values ( &n_data, &tc, &p, &eps );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, eps );
  }
  return;
}
/******************************************************************************/

void dilogarithm_values_test ( )

/******************************************************************************/
/*
  Purpose:

    DILOGARITHM_VALUES_TEST tests DILOGARITHM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "DILOGARITHM_VALUES_TEST:\n" );
  printf ( "  DILOGARITHM_VALUES stores values of\n" );
  printf ( "  the dilogarithm function.\n" );
  printf ( "\n" );
  printf ( "      X          DILOGARITHM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dilogarithm_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", x, fx );
  }
  return;
}
/******************************************************************************/

void e1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    E1_VALUES_TEST tests E1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "E1_VALUES_TEST:\n" );
  printf ( "  E1_VALUES stores values of\n" );
  printf ( "  the exponential integral function E1(X).\n" );
  printf ( "\n" );
  printf ( "      X          E1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    e1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void ei_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EI_VALUES_TEST tests EI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EI_VALUES_TEST:\n" );
  printf ( "  EI_VALUES stores values of\n" );
  printf ( "  the exponential integral function EI(X).\n" );
  printf ( "\n" );
  printf ( "      X          EI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ei_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void elliptic_ea_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ELLIPTIC_EA_VALUES_TEST tests ELLIPTIC_EA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ELLIPTIC_EA_VALUES_TEST:\n" );
  printf ( "  ELLIPTIC_EA_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the second\n" );
  printf ( "  kind, with parameter angle ALPHA in degrees.\n" );
  printf ( "\n" );
  printf ( "    ALPHA        EA(ALPHA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ea_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void elliptic_em_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ELLIPTIC_EM_VALUES_TEST tests ELLIPTIC_EM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ELLIPTIC_EM_VALUES_TEST:\n" );
  printf ( "  ELLIPTIC_EM_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the second\n" );
  printf ( "  kind, with parameter modulus M.\n" );
  printf ( "\n" );
  printf ( "      M            EM(M)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_em_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void elliptic_ka_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ELLIPTIC_KA_VALUES_TEST tests ELLIPTIC_KA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ELLIIPTIC_KA_VALUES_TEST:\n" );
  printf ( "  ELLIPTIC_KA_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the first\n" );
  printf ( "  kind, with parameter angle ALPHA in degrees.\n" );
  printf ( "\n" );
  printf ( "    ALPHA        KA(ALPHA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ka_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void elliptic_km_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ELLIPTIC_KM_VALUES_TEST tests ELLIPTIC_KM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ELLIPTIC_KM_VALUES_TEST:\n" );
  printf ( "  ELLIPTIC_KM_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the first\n" );
  printf ( "  kind, with parameter modulus M.\n" );
  printf ( "\n" );
  printf ( "      M            KM(M)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_km_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void erf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ERF_VALUES_TEST tests ERF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ERF_VALUES_TEST:\n" );
  printf ( "  ERF_VALUES stores values of\n" );
  printf ( "  the error function ERF(X).\n" );
  printf ( "\n" );
  printf ( "      X          ERF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void erfc_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ERFC_VALUES_TEST tests ERFC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "ERFC_VALUES_TEST:\n" );
  printf ( "  ERFC_VALUES stores values of\n" );
  printf ( "  the complementary error function ERFC(X).\n" );
  printf ( "\n" );
  printf ( "      X          ERFC(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void euler_number_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EULER_NUMBER_VALUES_TEST tests EULER_NUMBER_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "EULER_NUMBER_VALUES_TEST:\n" );
  printf ( "  EULER_NUMBER_VALUES returns values of\n" );
  printf ( "  the Euler numbers.\n" );
  printf ( "\n" );
  printf ( "     N        EULER_NUMBER(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void euler_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EULER_POLY_VALUES_TEST tests EULER_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EULER_POLY_VALUES_TEST:\n" );
  printf ( "  EULER_POLY_VALUES returns values of\n" );
  printf ( "  the Euler numbers.\n" );
  printf ( "\n" );
  printf ( "     N     X       EULER_POLY(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void exp_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EXP_VALUES_TEST tests EXP_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EXP_VALUES_TEST:\n" );
  printf ( "   EXP_VALUES stores values of the exponential function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void exp3_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EXP3_INT_VALUES_TEST tests EXP3_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EXP3_INT_VALUES_TEST:\n" );
  printf ( "  EXP3_INT_VALUES stores values of \n" );
  printf ( "  the exponential integral function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp3_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void exponential_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EXPONENTIAL_CDF_VALUES_TEST tests EXPONENTIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EXPONENTIAL_CDF_VALUES_TEST:\n" );
  printf ( "  EXPONENTIAL_CDF_VALUES stores values of \n" );
  printf ( "  the exponential CDF.\n" );
  printf ( "\n" );
  printf ( "      LAMBDA          X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exponential_cdf_values ( &n_data, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.8f  %24.8f  %24.16e\n", lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void extreme_values_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    EXTREME_VALUES_CDF_VALUES_TEST tests EXTREME_VALUES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "EXTREME_VALUES_CDF_VALUES_TEST:\n" );
  printf ( "  EXTREME_VALUES_CDF_VALUES stores values of \n" );
  printf ( "  the extreme values CDF.\n" );
  printf ( "\n" );
  printf ( "        Alpha    Beta        X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    extreme_values_cdf_values ( &n_data, &alpha, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", alpha, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void f_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    F_CDF_VALUES_TEST tests F_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " F_CDF_VALUES_TEST:\n" );
  printf ( "   F_CDF_VALUES stores values of\n" );
  printf ( "   the F cumulative density function.\n" );
  printf ( "\n" );
  printf ( "     A       B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %12f\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void f_noncentral_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    F_NONCENTRAL_CDF_VALUES_TEST tests F_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "F_NONCENTRAL_CDF_VALUES_TEST:\n" );
  printf ( "  F_NONCENTRAL_CDF_VALUES stores values of\n" );
  printf ( "  the F cumulative density function.\n" );
  printf ( "\n" );
  printf ( "     A       B            LAMBDA    X            CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %8f  %12f  %12e\n", a, b, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void fresnel_cos_values_test ( )

/******************************************************************************/
/*
  Purpose:

    FRESNEL_COS_VALUES_TEST tests FRESNEL_COS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "FRESNEL_COS_VALUES_TEST:\n" );
  printf ( "  FRESNEL_COS_VALUES stores values of\n" );
  printf ( "  the Fresnel cosine integral C(X).\n" );
  printf ( "\n" );
  printf ( "      X           C(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    fresnel_cos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void fresnel_sin_values_test ( )

/******************************************************************************/
/*
  Purpose:

    FRESNEL_SIN_VALUES_TEST tests FRESNEL_SIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "FRESNEL_SIN_VALUES_TEST:\n" );
  printf ( "  FRESNEL_SIN_VALUES stores values of\n" );
  printf ( "  the Fresnel sine integral S(X).\n" );
  printf ( "\n" );
  printf ( "      X           S(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    fresnel_sin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void frobenius_number_order2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    FROBENIUS_NUMBER_ORDER2_VALUES_TEST tests FROBENIUS_NUMBER_ORDER2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  int c1;
  int c2;
  int f;
  int n_data;

  printf ( "\n" );
  printf ( "FROBENIUS_NUMBER_ORDER2_VALUES_TEST:\n" );
  printf ( "  FROBENIUS_NUMBER_ORDER2_VALUES returns values of \n" );
  printf ( "  the Frobenius number of order 2.\n" );
  printf ( "\n" );
  printf ( "         C1        C2          F(C1,C2)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order2_values ( &n_data, &c1, &c2, &f );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %8d  %8d  %8d\n", c1, c2, f );
  }

  return;
}
/******************************************************************************/

void frobenius_number_order_values_test ( )

/******************************************************************************/
/*
  Purpose:

    FROBENIUS_NUMBER_ORDER_VALUES_TEST tests FROBENIUS_NUMBER_ORDER_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2007

  Author:

    John Burkardt
*/
{
  int *c;
  int f;
  int i;
  int n_data;
  int order;

  printf ( "\n" );
  printf ( "FROBENIUS_NUMBER_ORDER_VALUES_TEST:\n" );
  printf ( "  FROBENIUS_NUMBER_ORDER_VALUES returns the order for\n" );
  printf ( "  a Frobenius problem;\n" );

  n_data = 0;
  printf ( "\n" );
  printf ( "       #      Order\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    frobenius_number_order_values ( &n_data, &order );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %4d\n", n_data, order );
  }
  return;
}
/******************************************************************************/

void frobenius_number_data_values_test ( )

/******************************************************************************/
/*
  Purpose:

    FROBENIUS_NUMBER_DATA_VALUES_TEST tests FROBENIUS_NUMBER_DATA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2007

  Author:

    John Burkardt
*/
{
  int *c;
  int f;
  int i;
  int n_data;
  int order;

  printf ( "\n" );
  printf ( "FROBENIUS_NUMBER_DATA_VALUES_TEST:\n" );
  printf ( "  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding\n" );
  printf ( "  coin denominations.\n" );

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order_values ( &n_data, &order );

    if ( n_data == 0 )
    {
      break;
    }

    c = ( int * ) malloc ( order * sizeof ( int ) );

    frobenius_number_data_values ( &n_data, order, c, &f );

    printf ( "\n" );
    printf ( "  Order = %d\n", order );
    for ( i = 0; i < order; i++ )
    {
      printf ( "  %8d", c[i] );
    }
    printf ( "\n" );
    printf ( "  Frobenius number = %d\n", f );

    free ( c );
  }
  return;
}
/******************************************************************************/

void gamma_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_VALUES_TEST tests GAMMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_VALUES_TEST:\n" );
  printf ( "  GAMMA_VALUES stores values of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X            GAMMA(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_CDF_VALUES_TEST tests GAMMA_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_CDF_VALUES_TEST:\n" );
  printf ( "  GAMMA_CDF_VALUES stores values of\n" );
  printf ( "  the Gamma CDF.\n" );
  printf ( "\n" );
  printf ( "      M    Sigma      X            CDF((X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_inc_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_VALUES_TEST tests GAMMA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_VALUES_TEST:\n" );
  printf ( "   GAMMA_INC_VALUES stores values of\n" );
  printf ( "   the incomplete Gamma function.\n" );
  printf ( "\n" );
  printf ( "      A            X            GAMMA_INC(A)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_inc_p_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_P_VALUES_TEST tests GAMMA_INC_P_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_P_VALUES_TEST:\n" );
  printf ( "   GAMMA_INC_P_VALUES stores values of\n" );
  printf ( "   the incomplete Gamma P function.\n" );
  printf ( "\n" );
  printf ( "      A            X            F(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_p_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_inc_q_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_Q_VALUES_TEST tests GAMMA_INC_Q_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_Q_VALUES_TEST:\n" );
  printf ( "   GAMMA_INC_Q_VALUES stores values of\n" );
  printf ( "   the incomplete Gamma Q function.\n" );
  printf ( "\n" );
  printf ( "      A            X            F(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_q_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_inc_tricomi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_INC_TRICOMI_VALUES_TEST tests GAMMA_INC_TRICOMI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_INC_TRICOMI_VALUES_TEST:\n" );
  printf ( "   GAMMA_INC_TRICOMI_VALUES stores values of\n" );
  printf ( "   the incomplete Tricomi Gamma function.\n" );
  printf ( "\n" );
  printf ( "      A            X            F(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_tricomi_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void gamma_log_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GAMMA_LOG_VALUES_TEST tests GAMMA_LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GAMMA_LOG_VALUES_TEST:\n" );
  printf ( "  GAMMA_LOG_VALUES stores values of\n" );
  printf ( "  the logarithm of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X            GAMMA_LOG(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void gegenbauer_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_POLY_VALUES_TEST tests GEGENBAUER_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GEGENBAUER_POLY_VALUES_TEST:\n" );
  printf ( "  GEGENBAUER_POLY_VALUES returns values of\n" );
  printf ( "  the Gegenbauer polynomials.\n" );
  printf ( "\n" );
  printf ( "       N       A       X       G(N,A)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12f  %12f  %24.16e\n", n, a, x, fx );
  }

  return;
}
/******************************************************************************/

void geometric_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GEOMETRIC_CDF_VALUES_TEST tests GEOMETRIC_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double cdf;
  int n_data;
  double p;
  int x;

  printf ( "\n" );
  printf ( "GEOMETRIC_CDF_VALUES_TEST:\n" );
  printf ( "  GEOMETRIC_CDF_VALUES stores values of\n" );
  printf ( "  the Geometric Probability Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      X      P       CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    geometric_cdf_values ( &n_data, &x, &p, &cdf );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12f  %24.16e\n", x, p, cdf );
  }
  return;
}
/******************************************************************************/

void goodwin_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GOODWIN_VALUES_TEST tests GOODWIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GOODWIN_VALUES_TEST:\n" );
  printf ( "  GOODWIN_VALUES stores values of \n" );
  printf ( "  the Goodwin function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    goodwin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void gud_values_test ( )

/******************************************************************************/
/*
  Purpose:

    GUD_VALUES_TEST tests GUD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GUD_VALUES_TEST:\n" );
  printf ( "  GUD_VALUES stores values of\n" );
  printf ( "  the Gudermannian function.\n" );
  printf ( "\n" );
  printf ( "      X            GUD(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void hermite_function_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_FUNCTION_VALUES_TEST tests HERMITE_FUNCTION_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HERMITE_FUNCTION_VALUES_TEST\n" );
  printf ( "  HERMITE_FUNCTION_VALUES stores values of\n" );
  printf ( "  the Hermite function.\n" );
  printf ( "\n" );
  printf ( "     N      X            Hf(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_function_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void hermite_poly_phys_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLY_PHYS_VALUES_TEST tests HERMITE_POLY_PHYS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HERMITE_POLY_PHYS_VALUES_TEST\n" );
  printf ( "  HERMITE_POLY_PHYS_VALUES stores values of\n" );
  printf ( "  the physicist's Hermite polynomials.\n" );
  printf ( "\n" );
  printf ( "     N      X            H(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_phys_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void hermite_poly_prob_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLY_PROB_VALUES_TEST tests HERMITE_POLY_PROB_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HERMITE_POLY_PROB_VALUES_TEST\n" );
  printf ( "  HERMITE_POLY_PROB_VALUES stores values of\n" );
  printf ( "  the probabilist's Hermite polynomials.\n" );
  printf ( "\n" );
  printf ( "     N      X            He(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_prob_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void hyper_1f1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HYPER_1F1_VALUES_TEST tests HYPER_1F1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HYPER_1F1_VALUES_TEST:\n" );
  printf ( "  HYPER_1F1_VALUES stores values of\n" );
  printf ( "  the hypergeometric function 1F1.\n" );
  printf ( "\n" );
  printf ( "      A      B      X   Hyper_1F1(A,B,C,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hyper_1f1_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", a, b,x, fx ); 
  }
  return;
}
/******************************************************************************/

void hyper_2f1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HYPER_2F1_VALUES_TEST tests HYPER_2F1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 September 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HYPER_2F1_VALUES_TEST:\n" );
  printf ( "  HYPER_2F1_VALUES stores values of\n" );
  printf ( "  the hypergeometric function 2F1.\n" );
  printf ( "\n" );
  printf ( "      A      B     C      X   Hyper_2F1(A,B,C,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %8f  %24.16e\n", a, b, c, x, fx ); 
  }
  return;
}
/******************************************************************************/

void hypergeometric_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HYPERGEOMETRIC_CDF_VALUES_TEST tests HYPERGEOMETRIC_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  printf ( "\n" );
  printf ( "HYPERGEOMETRIC_CDF_VALUES_TEST:\n" );
  printf ( "  HYPERGEOMETRIC_CDF_VALUES stores values of\n" );
  printf ( "  the Hypergeometric CDF.\n" );
  printf ( "\n" );
  printf ( "     SAM    SUC   POP     X   HyperCDF(S,S,P)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_cdf_values ( &n_data, &sam, &succ, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %8d  %8d  %24.16e\n", sam, succ, pop, x, fx );
  }
  return;
}
/******************************************************************************/

void hypergeometric_pdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HYPERGEOMETRIC_PDF_VALUES_TEST tests HYPERGEOMETRIC_PDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 January 2008

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  printf ( "\n" );
  printf ( "HYPERGEOMETRIC_PDF_VALUES_TEST:\n" );
  printf ( "  HYPERGEOMETRIC_PDF_VALUES stores values of\n" );
  printf ( "  the Hypergeometric PDF.\n" );
  printf ( "\n" );
  printf ( "     SAM    SUC   POP     X   HyperPDF(S,S,P)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_pdf_values ( &n_data, &sam, &succ, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %8d  %8d  %24.16e\n", sam, succ, pop, x, fx );
  }
  return;
}
/******************************************************************************/

void hypergeometric_u_values_test ( )

/******************************************************************************/
/*
  Purpose:

    HYPERGEOMETRIC_U_VALUES_TEST tests HYPERGEOMETRIC_U_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HYPERGEOMETRIC_U_VALUES_TEST:\n" );
  printf ( "  HYPERGEOMETRIC_U_VALUES stores values of\n" );
  printf ( "  the hypergeometric function 1U.\n" );
  printf ( "\n" );
  printf ( "      A      B      X   HyperU(A,B,C,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_u_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", a, b,x, fx ); 
  }
  return;
}
/******************************************************************************/

void i0ml0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I0ML0_VALUES_TEST tests I0ML0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "I0ML0_VALUES_TEST:\n" );
  printf ( "  I0ML0_VALUES stores values of \n" );
  printf ( "  the I0-L0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i0ml0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void i1ml1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I1ML1_VALUES_TEST tests I1ML1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "I1ML1_VALUES_TEST:\n" );
  printf ( "  I1ML1_VALUES stores values of \n" );
  printf ( "  the I1-L1 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i1ml1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void i4_factorial_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FACTORIAL_VALUES_TEST tests I4_FACTORIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 March 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( " I4_FACTORIAL_VALUES_TEST:\n" );
  printf ( "   I4_FACTORIAL_VALUES return;s values of\n" );
  printf ( "   the factorial function.\n" );
  printf ( "\n" );
  printf ( "      N         Factorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void i4_factorial2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FACTORIAL2_VALUES_TEST tests I4_FACTORIAL2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( " I4_FACTORIAL2_VALUES_TEST:\n" );
  printf ( "   I4_FACTORIAL2_VALUES return;s values of\n" );
  printf ( "   the double factorial function.\n" );
  printf ( "\n" );
  printf ( "      N         DoubleFactorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void i4_fall_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FALL_VALUES_TEST tests I4_FALL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2014

  Author:

    John Burkardt
*/
{
  int fmn;
  int m;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "I4_FALL_VALUES_TEST:\n" );
  printf ( "  I4_FALL_VALUES returns some exact values\n" );
  printf ( "  of the falling factorial function:\n" );
  printf ( "\n" );
  printf ( "     M     N      I4_FALL(M,N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_fall_values ( &n_data, &m, &n, &fmn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", m, n, fmn );
  }
  return;
}
/******************************************************************************/

void i4_rise_values_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_RISE_VALUES_TEST tests I4_RISE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2014

  Author:

    John Burkardt
*/
{
  int fmn;
  int m;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "I4_RISE_VALUES_TEST:\n" );
  printf ( "  I4_RISE_VALUES returns some exact values\n" );
  printf ( "  of the rising factorial function:\n" );
  printf ( "\n" );
  printf ( "     M     N      I4_RISE(M,N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_rise_values ( &n_data, &m, &n, &fmn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", m, n, fmn );
  }
  return;
}
/******************************************************************************/

void int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    INT_VALUES_TEST tests INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "INT_VALUES_TEST:\n" );
  printf ( "  INT_VALUES stores values of the integer part of a real number.\n" );
  printf ( "\n" );
  printf ( "      X            INT(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void jacobi_cn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_CN_VALUES_TEST tests JACOBI_CN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "JACOBI_CN_VALUES_TEST:\n" );
  printf ( "  JACOBI_CN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic CN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       CN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_cn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void jacobi_dn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_DN_VALUES_TEST tests JACOBI_DN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "JACOBI_DN_VALUES_TEST:\n" );
  printf ( "  JACOBI_DN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic DN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       DN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_dn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void jacobi_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_POLY_VALUES_TEST tests JACOBI_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "JACOBI_POLY_VALUES_TEST:\n" );
  printf ( "  JACOBI_POLY_VALUES returns values of\n" );
  printf ( "  the Jacobi polynomial.\n" );
  printf ( "\n" );
  printf ( "       N         A         B      X       J(N,A,B)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_poly_values ( &n_data, &n, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %8f  %12f  %24.16e\n", n, a, b, x, fx );
  }

  return;
}
/******************************************************************************/

void jacobi_sn_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_SN_VALUES_TEST tests JACOBI_SN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "JACOBI_SN_VALUES_TEST:\n" );
  printf ( "  JACOBI_SN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic SN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       SN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_sn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void jed_ce_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JED_CE_VALUES_TEST tests JED_CE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int d;
  double f;
  double jed;
  int n_data;
  int m;
  int y;

  printf ( "\n" );
  printf ( "JED_CE_VALUES_TEST:\n" );
  printf ( "  JED_CE_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  YMDF, the corresponding year, month, day, fraction.\n" );
  printf ( "\n" );
  printf ( "        JED          Y   M   D    F\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_ce_values ( &n_data, &jed, &y, &m, &d, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %6d  %2d  %2d  %6f\n", jed, y, m, d, f );
  }
  return;
}
/******************************************************************************/

void jed_mjd_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JED_MJD_VALUES_TEST tests JED_MJD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  double mjd;

  printf ( "\n" );
  printf ( "JED_MJD_VALUES_TEST:\n" );
  printf ( "  JED_MJD_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  MJD, the corresponding Modified Julian Day count.\n" );
  printf ( "\n" );
  printf ( "   JED      MJD\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_mjd_values ( &n_data, &jed, &mjd );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", jed, mjd );
  }
  return;
}
/******************************************************************************/

void jed_rd_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JED_RD_VALUES_TEST tests JED_RD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  double rd;

  printf ( "\n" );
  printf ( "JED_RD_VALUES_TEST:\n" );
  printf ( "  JED_RD_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  RD, the corresponding Reingold Dershowitz Day count.\n" );
  printf ( "\n" );
  printf ( "   JED      RD\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_rd_values ( &n_data, &jed, &rd );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", jed, rd );
  }
  return;
}
/******************************************************************************/

void jed_weekday_values_test ( )

/******************************************************************************/
/*
  Purpose:

    JED_WEEKDAY_VALUE_TEST tests JED_WEEKDAY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  int weekday;

  printf ( "\n" );
  printf ( "JED_WEEKDAY_VALUES_TEST:\n" );
  printf ( "  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates \n" );
  printf ( "  (JED) and the corresponding weekday\n" );
  printf ( "\n" );
  printf ( "   JED      #  Weekday\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_weekday_values ( &n_data, &jed, &weekday );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %d\n", jed, weekday );
  }
  return;
}
/******************************************************************************/

void kei0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    KEI0_VALUES_TEST tests KEI0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "KEI0_VALUES_TEST:\n" );
  printf ( "  KEI0_VALUES stores values of \n" );
  printf ( "  the Kelvin function KEI of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    kei0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void kei1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    KEI1_VALUES_TEST tests KEI1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "KEI1_VALUES_TEST:\n" );
  printf ( "  KEI1_VALUES stores values of \n" );
  printf ( "  the Kelvin function KEI of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    kei1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void ker0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    KER0_VALUES_TEST tests KER0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "KER0_VALUES_TEST:\n" );
  printf ( "  KER0_VALUES stores values of \n" );
  printf ( "  the Kelvin function KER of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ker0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void ker1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    KER1_VALUES_TEST tests KER1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "KER1_VALUES_TEST:\n" );
  printf ( "  KER1_VALUES stores values of \n" );
  printf ( "  the Kelvin function KER of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ker1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void laguerre_associated_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_ASSOCIATED_VALUES_TEST tests LAGUERRE_ASSOCIATED_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAGUERRE_ASSOCIATED_VALUES_TEST:\n" );
  printf ( "  LAGUERRE_ASSOCIATED_VALUES stores values of\n" );
  printf ( "  the associated Laguerre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             L(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void laguerre_general_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_GENERAL_VALUES_TEST tests LAGUERRE_GENERAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAGUERRE_GENERAL_VALUES_TEST:\n" );
  printf ( "  LAGUERRE_GENERAL_VALUES stores values of\n" );
  printf ( "  the generalized Laguerre function.\n" );
  printf ( "\n" );
  printf ( "     N     A    X             L(N,A)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_general_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %%12.4f  %12f  %24.16e\n", n, a, x, fx );
  }
  return;
}
/******************************************************************************/

void laguerre_polynomial_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_VALUES_TEST tests LAGUERRE_POLYNOMIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_VALUES_TEST:\n" );
  printf ( "  LAGUERRE_POLYNOMIAL_VALUES stores values of \n" );
  printf ( "  the Laguerre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     X            L(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void lambert_w_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LAMBERT_W_VALUES_TEST tests LAMBERT_W_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAMBERT_W_VALUES_TEST:\n" );
  printf ( "  LAMBERT_W_VALUES stores values of \n" );
  printf ( "  the Lambert W function.\n" );
  printf ( "\n" );
  printf ( "                X                     W(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void laplace_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LAPLACE_CDF_VALUES_TEST tests LAPLACE_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAPLACE_CDF_VALUES_TEST:\n" );
  printf ( "  LAPLACE_CDF_VALUES returns values of \n" );
  printf ( "  the Laplace Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Beta         X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laplace_cdf_values ( &n_data, &mu, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void legendre_associated_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ASSOCIATED_VALUES_TEST tests LEGENDRE_ASSOCIATED_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_ASSOCIATED_VALUES_TEST:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_VALUES stores values of\n" );
  printf ( "  the associated Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void legendre_associated_normalized_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES stores values of\n" );
  printf ( "  the normalized associated Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void legendre_associated_normalized_sphere_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES stores values of\n" );
  printf ( "  the associated Legendre polynomials normalized for the unit sphere.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_sphere_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void legendre_poly_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_POLY_VALUES_TEST tests LEGENDRE_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_POLY_VALUES_TEST:\n" );
  printf ( "  LEGENDRE_POLY_VALUES stores values of \n" );
  printf ( "  the Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N    X             P(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void legendre_function_q_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_FUNCTION_Q_VALUES_TEST tests LEGENDRE_FUNCTION_Q_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_FUNCTION_Q_VALUES_TEST:\n" );
  printf ( "  LEGENDRE_FUNCTION_Q_VALUES stores values of\n" );
  printf ( "  the Legendre Q function.\n" );
  printf ( "\n" );
  printf ( "     N    X             Q(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void lerch_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LERCH_VALUES_TEST tests LERCH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  int s;
  double z;

  printf ( "\n" );
  printf ( "LERCH_VALUES_TEST:\n" );
  printf ( "  LERCH_VALUES returns values of\n" );
  printf ( "  the Lerch transcendent function.\n" );
  printf ( "\n" );
  printf ( "      Z      S      A      Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %6d  %12f  %24.16e\n", z, s, a, fx );
  }
  return;
}
/******************************************************************************/

void lobachevsky_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBACHEVSKY_VALUES_TEST tests LOBACHEVSKY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOBACHEVSKY_VALUES_TEST:\n" );
  printf ( "  LOBACHEVSKY_VALUES stores values of \n" );
  printf ( "  the Lobachevsky function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobachevsky_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void lobatto_polynomial_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_VALUES_TEST tests LOBATTO_POLYNOMIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 May 2013

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_VALUES_TEST:\n" );
  printf ( "  LOBATTO_POLYNOMIAL_VALUES stores values of \n" );
  printf ( "  the completed Lobatto polynomials.\n" );
  printf ( "\n" );
  printf ( "     N    X            Lo(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void lobatto_polynomial_derivatives_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_DERIVATIVES_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 November 2014

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_DERIVATIVES_TEST:\n" );
  printf ( "  LOBATTO_POLYNOMIAL_DERIVATIVES stores derivatives of \n" );
  printf ( "  the completed Lobatto polynomials.\n" );
  printf ( "\n" );
  printf ( "     N    X            Lo'(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_derivatives ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void log_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG_VALUES_TEST tests LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOG_VALUES_TEST:\n" );
  printf ( "   LOG_VALUES stores values of the natural logarithm function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void log_normal_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG_NORMAL_CDF_VALUES_TEST tests LOG_NORMAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "LOG_NORMAL_CDF_VALUES_TEST:\n" );
  printf ( "  LOG_NORMAL_CDF_VALUES returns values of \n" );
  printf ( "  the Log Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Sigma        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_normal_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void log_series_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG_SERIES_CDF_VALUES_TEST tests LOG_SERIES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "LOG_SERIES_CDF_VALUES_TEST:\n" );
  printf ( "  LOG_SERIES_CDF_VALUES returns values of \n" );
  printf ( "  the Log Series Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     T      N   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_series_cdf_values ( &n_data, &t, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %6d  %24.16e\n", t, n, fx );
  }
  return;
}
/******************************************************************************/

void log10_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOG10_VALUES_TEST tests LOG10_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOG10_VALUES_TEST:\n" );
  printf ( "   LOG10_VALUES stores values of the natural logarithm function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log10_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void logarithmic_integral_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOGARITHMIC_INTEGRAL_VALUES_TEST tests LOGARITHMIC_INTEGRAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOGARITHMIC_INTEGRAL_VALUES_TEST:\n" );
  printf ( "  LOGARITHMIC_INTEGAL_VALUES stores values of\n" );
  printf ( "  the logarithmic integral function.\n" );
  printf ( "\n" );
  printf ( "      X            LI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    logarithmic_integral_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void logistic_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LOGISTIC_CDF_VALUES_TEST tests LOGISTIC_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LOGISTIC_CDF_VALUES_TEST:\n" );
  printf ( "  LOGISTIC_CDF_VALUES returns values of \n" );
  printf ( "  the Logistic Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Beta         X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    logistic_cdf_values ( &n_data, &mu, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void mertens_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    MERTENS_VALUES_TEST tests MERTENS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 October 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "MERTENS_VALUES_TEST:\n" );
  printf ( "  MERTENS_VALUES returns values of\n" );
  printf ( "  the Mertens function.\n" );
  printf ( "\n" );
  printf ( "     N         MERTENS(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    mertens_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    } 
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void moebius_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    MOEBIUS_VALUES_TEST tests MOEBIUS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "MOEBIUS_VALUES_TEST:\n" );
  printf ( "  MOEBIUS_VALUES returns values of\n" );
  printf ( "  the Moebius function.\n" );
  printf ( "\n" );
  printf ( "     N         MU(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    moebius_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    } 
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void negative_binomial_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    NEGATIVE_BINOMIAL_CDF_VALUES_TEST tests NEGATIVE_BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double cdf;
  int f;
  int n_data;
  double p;
  int s;

  printf ( "\n" );
  printf ( "NEGATIVE_BINOMIAL_CDF_VALUES_TEST:\n" );
  printf ( "  NEGATIVE_BINOMIAL_CDF_VALUES stores values of\n" );
  printf ( "  the Negative Binomial Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     F     S         P         CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &p, &cdf );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %24.16e\n", f, s, p, cdf );
  }
  return;
}
/******************************************************************************/

void nine_j_values_test ( )

/******************************************************************************/
/*
  Purpose:

    NINE_J_VALUES_TEST demonstrates NINE_J_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  double j7;
  double j8;
  double j9;
  int n_data;

  printf ( "\n" );
  printf ( "NINE_J_VALUES_TEST:\n" );
  printf ( "  NINE_J_VALUES returns values of\n" );
  printf ( "  the Wigner 9J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      J4      J5      J6" );
  printf ( "      J7      J8      J9        NINE_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    nine_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, j4, j5, j6, j7, j8, j9, fx );
  }
  return;
}
/******************************************************************************/

void normal_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    NORMAL_CDF_VALUES_TEST tests NORMAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_CDF_VALUES_TEST:\n" );
  printf ( "  NORMAL_CDF_VALUES stores values of\n" );
  printf ( "  the Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "            X                   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void normal_01_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    NORMAL_01_CDF_VALUES_TEST tests NORMAL_01_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "NORMAL_01_CDF_VALUES_TEST:\n" );
  printf ( "  NORMAL_01_CDF_VALUES stores values of\n" );
  printf ( "  the Normal 01 Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "            X                   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void omega_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    OMEGA_VALUES_TEST tests OMEGA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "OMEGA_VALUES_TEST:\n" );
  printf ( "  OMEGA_VALUES returns values of\n" );
  printf ( "  the Omega function.\n" );
  printf ( "\n" );
  printf ( "     N           OMEGA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void owen_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    OWEN_VALUES_TEST tests OWEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double h;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "OWEN_VALUES_TEST\n" );
  printf ( "  OWEN_VALUES stores values of\n" );
  printf ( "  Owen's T function.\n" );
  printf ( "\n" );
  printf ( "          H            A            T\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f %12f  %24.16e\n", h, a, t );
  }
  return;
}
/******************************************************************************/

void partition_count_values_test ( )
 
/******************************************************************************/
/*
  Purpose:

    PARTITION_COUNT_VALUES_TEST tests PARTITION_COUNT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "PARTITION_COUNT_VALUES_TEST:\n" );
  printf ( "  PARTITION_COUNT_VALUES returns values of \n" );
  printf ( "  the integer partition count function.\n" );
  printf ( "\n" );
  printf ( "     N         P(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_count_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void partition_distinct_count_values_test ( )

/******************************************************************************/
/*
  Purpose:

    PARTITION_DISTINCT_COUNT_VALUES_TEST tests PARTITION_DISTINCT_COUNT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "PARTITION_DISTINCT_COUNT_VALUES_TEST:\n" );
  printf ( "  PARTITION_DISTINCT_COUNT_VALUES returns values of \n" );
  printf ( "  the integer distinct partition count function.\n" );
  printf ( "\n" );
  printf ( "     N         Q(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void phi_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    PHI_VALUES_TEST tests PHI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "PHI_VALUES_TEST:\n" );
  printf ( "  PHI_VALUES returns values of\n" );
  printf ( "  the PHI function.\n" );
  printf ( "\n" );
  printf ( "     N         PHI(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void pi_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    PI_VALUES_TEST tests PI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "PI_VALUES_TEST:\n" );
  printf ( "  PI_VALUES returns values of\n" );
  printf ( "  the PI function.\n" );
  printf ( "\n" );
  printf ( "     N         PI(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pi_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void poisson_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    POISSON_CDF_VALUES_TEST tests POISSON_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "POISSON_CDF_VALUES_TEST:\n" );
  printf ( "  POISSON_CDF_VALUES returns values of\n" );
  printf ( "  the Poisson Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      A     X       CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %4d  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void polylogarithm_values_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYLOGARITHM_VALUES_TEST tests POLYLOGARITHM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double z;

  printf ( "\n" );
  printf ( "POLYLOGARITHM_VALUES_TEST:\n" );
  printf ( "  POLYLOGARITHM_VALUES returns values of \n" );
  printf ( "  the polylogarithm function.\n" );
  printf ( "\n" );
  printf ( "     N      Z          Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    polylogarithm_values ( &n_data, &n, &z, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e  %24.16e\n", n, z, fx );
  }
  return;
}
/******************************************************************************/

void prandtl_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    PRANDTL_VALUES_TEST tests PRANDTL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double p;
  double pr;
  double tc;

  printf ( "\n" );
  printf ( "PRANDTL_VALUES_TEST:\n" );
  printf ( "  PRANDTL_VALUES stores values of\n" );
  printf ( "  the Prandtl number of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            Pr(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    prandtl_values ( &n_data, &tc, &p, &pr );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, pr );
  }
  return;
}
/******************************************************************************/

void prime_values_test ( )

/******************************************************************************/
/*
  Purpose:

    PRIME_VALUES_TEST tests PRIME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n;
  int n_data;
  int p;

  printf ( "\n" );
  printf ( "PRIME_VALUES_TEST:\n" );
  printf ( "  PRIME_VALUES returns values of\n" );
  printf ( "  the prime function.\n" );
  printf ( "\n" );
  printf ( "           N          P[N]\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    prime_values ( &n_data, &n, &p );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, p );
  }

  return;
}
/******************************************************************************/

void psat_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    PSAT_VALUES_TEST tests PSAT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double psat;
  double tc;

  printf ( "\n" );
  printf ( "PSAT_VALUES_TEST:\n" );
  printf ( "  PSAT_VALUES stores values of\n" );
  printf ( "  the saturation pressure of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            PSAT(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psat_values ( &n_data, &tc, &psat );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", tc, psat );
  }
  return;
}
/******************************************************************************/

void psi_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    PSI_VALUES_TEST tests PSI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "PSI_VALUES_TEST\n" );
  printf ( "  PSI_VALUES stores values of\n" );
  printf ( "  the PSI function.\n" );
  printf ( "\n" );
  printf ( "      X            PSI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void r8_factorial_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    R8_FACTORIAL_VALUES_TEST tests R8_FACTORIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL_VALUES_TEST:\n" );
  printf ( "  R8_FACTORIAL_VALUES stores values of\n" );
  printf ( "  the factorial function (using double arithmetic).\n" );
  printf ( "\n" );
  printf ( "      N       Factorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %24.16e\n", n, fn );
  }
  return;
}
/******************************************************************************/

void r8_factorial_log_values_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL_LOG_VALUES_TEST tests R8_FACTORIAL_LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL_LOG_VALUES_TEST:\n" );
  printf ( "  R8_FACTORIAL_LOG_VALUES stores values of\n" );
  printf ( "  the logarithm of the factorial function\n" );
  printf ( "  (using real arithmetic).\n" );
  printf ( "\n" );
  printf ( "      N       Log(Factorial(N))\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fn );
  }
  return;
}
/******************************************************************************/

void r8_factorial2_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    R8_FACTORIAL2_VALUES_TEST tests R8_FACTORIAL2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2015

  Author:

    John Burkardt
*/
{
  double f;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL2_VALUES_TEST:\n" );
  printf ( "  R8_FACTORIAL2_VALUES stores values of\n" );
  printf ( "  the double factorial function (using double arithmetic).\n" );
  printf ( "\n" );
  printf ( "      N               F\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( &n_data, &n, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %24.16e\n", n, f );
  }
  return;
}
/******************************************************************************/

void r8_fall_values_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FALL_VALUES_TEST tests R8_FALL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 December 2014

  Author:

    John Burkardt
*/
{
  double f;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_FALL_VALUES_TEST:\n" );
  printf ( "  R8_FALL_VALUES returns some exact values\n" );
  printf ( "  of the falling factorial function:\n" );
  printf ( "\n" );
  printf ( "     X     N      R8_FALL(X,N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_fall_values ( &n_data, &x, &n, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8g  %6d  %12g\n", x, n, f );
  }
  return;
}
/******************************************************************************/

void r8_rise_values_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_RISE_VALUES_TEST tests R8_RISE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 December 2014

  Author:

    John Burkardt
*/
{
  double f;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_RISE_VALUES_TEST:\n" );
  printf ( "  R8_RISE_VALUES returns some exact values\n" );
  printf ( "  of the rising factorial function:\n" );
  printf ( "\n" );
  printf ( "     X     N      R8_RISE(X,N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_rise_values ( &n_data, &x, &n, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8g  %6d  %12g\n", x, n, f );
  }
  return;
}
/******************************************************************************/

void rayleigh_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    RAYLEIGH_CDF_VALUES_TEST tests RAYLEIGH_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "RAYLEIGH_CDF_VALUES_TEST:\n" );
  printf ( "  RAYLEIGH_CDF_VALUES stores values of\n" );
  printf ( "  the Rayleigh CDF.\n" );
  printf ( "\n" );
  printf ( "      SIGMA        X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    rayleigh_cdf_values ( &n_data, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void secvir_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SECVIR_VALUES_TEST tests SECVIR_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double tc;
  double vir;

  printf ( "\n" );
  printf ( "SECVIR_VALUES_TEST:\n" );
  printf ( "  SECVIR_VALUES stores values of\n" );
  printf ( "  the second virial coefficient of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            VIR(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   secvir_values ( &n_data, &tc, &vir );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", tc, vir );
  }
  return;
}
/******************************************************************************/

void shi_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SHI_VALUES_TEST tests SHI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SHI_VALUES_TEST:\n" );
  printf ( "  SHI_VALUES stores values of\n" );
  printf ( "  the hyperbolic sine integral function.\n" );
  printf ( "\n" );
  printf ( "      X            SHI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   shi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void si_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SI_VALUES_TEST tests SI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SI_VALUES_TEST:\n" );
  printf ( "  SI_VALUES stores values of\n" );
  printf ( "  the sine integral function.\n" );
  printf ( "\n" );
  printf ( "      X            SI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   si_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void sigma_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SIGMA_VALUES_TEST tests SIGMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "SIGMA_VALUES_TEST:\n" );
  printf ( "  SIGMA_VALUES returns values of\n" );
  printf ( "  the SIGMA function.\n" );
  printf ( "\n" );
  printf ( "       N         SIGMA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sigma_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void sin_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_VALUES_TEST tests SIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SIN_VALUES_TEST:\n" );
  printf ( "   SIN_VALUES stores values of the sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void sin_degree_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_DEGREE_VALUES_TEST tests SIN_DEGREE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SIN_DEGREE_VALUES_TEST:\n" );
  printf ( "   SIN_DEGREE_VALUES stores values of the sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_degree_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void sin_power_int_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_POWER_INT_VALUES_TEST tests SIN_POWER_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "SIN_POWER_INT_VALUES_TEST:\n" );
  printf ( "  SIN_POWER_INT_VALUES returns values of\n" );
  printf ( "  the integral of the N-th power of the sine function.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %6d  %24.16e\n", a, b, n, fx );
  }
  return;
}
/******************************************************************************/

void sinh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SINH_VALUES_TEST tests SINH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SINH_VALUES_TEST:\n" );
  printf ( "   SINH_VALUES stores values of the hyperbolic sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sinh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void six_j_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SIX_J_VALUES_TEST tests SIX_J_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  int n_data;

  printf ( "\n" );
  printf ( "SIX_J_VALUES_TEST:\n" );
  printf ( "  SIX_J_VALUES returns values of \n" );
  printf ( "  the Wigner 6J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      J4      J5      J6        SIX_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    six_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n", 
    j1, j2, j3, j4, j5, j6, fx );
  }

  return;
}
/******************************************************************************/

void sound_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SOUND_VALUES_TEST tests SOUND_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double c;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "SOUND_VALUES_TEST:\n" );
  printf ( "  SOUND_VALUES stores values of\n" );
  printf ( "  the spead of sound in water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            C(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sound_values ( &n_data, &tc, &p, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, c );
  }
  return;
}
/******************************************************************************/

void sphere_unit_area_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_UNIT_AREA_VALUES_TEST tests SPHERE_UNIT_AREA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "SPHERE_UNIT_AREA_VALUES_TEST:\n" );
  printf ( "  SPHERE_UNIT_AREA_VALUES stores values of\n" );
  printf ( "  the area of the unit sphere in various dimensions.\n" );
  printf ( "\n" );
  printf ( "      N           AREA\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( &n_data, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fx );
  }
  return;
}
/******************************************************************************/

void sphere_unit_volume_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_UNIT_VOLUME_VALUES_TEST tests SPHERE_UNIT_VOLUME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "SPHERE_UNIT_VOLUME_VALUES_TEST:\n" );
  printf ( "  SPHERE_UNIT_VOLUME_VALUES stores values of\n" );
  printf ( "  the volume of the unit sphere in various dimensions.\n" );
  printf ( "\n" );
  printf ( "      N           VOLUME\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( &n_data, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fx );
  }
  return;
}
/******************************************************************************/

void spherical_harmonic_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERICAL_HARMONIC_VALUES_TEST tests SPHERICAL_HARMONIC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{

  int l;
  int m;
  int n_data;
  double phi;
  double theta;
  double yi;
  double yr;

  printf ( "\n" );
  printf ( "SPHERICAL_HARMONIC_VALUES_TEST:\n" );
  printf ( "  SPHERICAL_HARMONIC_VALUES stores values of\n" );
  printf ( "  the spherical harmonic function.\n" );
  printf ( "\n" );
  printf ( "   L   M    THETA       PHI           Yr                    Yi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %2d  %2d  %8f  %8f  %24.16e  %24.16e\n",
    l, m, theta, phi, yr, yi );
  }
  return;
}
/******************************************************************************/

void sqrt_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SQRT_VALUES_TEST tests SQRT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SQRT_VALUES_TEST:\n" );
  printf ( "  SQRT_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X       Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sqrt_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void stirling1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STIRLING1_VALUES_TEST tests STIRLING1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int n_data;
  int s1;

  printf ( "\n" );
  printf ( "STIRLING1_VALUES_TEST:\n" );
  printf ( "  STIRLING1_VALUES returns values of\n" );
  printf ( "  the Stirling numbers of the first kind.\n" );
  printf ( "\n" );
  printf ( "     N     N        S1\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( &n_data, &n, &m, &s1 );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", n, m, s1 );
  }
  return;
}
/******************************************************************************/

void stirling2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STIRLING2_VALUES_TEST tests STIRLING2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int n_data;
  int s2;

  printf ( "\n" );
  printf ( "STIRLING2_VALUES_TEST:\n" );
  printf ( "  STIRLING2_VALUES returns values of\n" );
  printf ( "  the Stirling numbers of the second kind.\n" );
  printf ( "\n" );
  printf ( "     N     N        S2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( &n_data, &n, &m, &s2 );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", n, m, s2 );
  }
  return;
}
/******************************************************************************/

void stromgen_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STROMGEN_VALUES_TEST tests STROMGEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STROMGEN_VALUES_TEST:\n" );
  printf ( "  STROMGEN_VALUES stores values of \n" );
  printf ( "  the Stromgen function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stromgen_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void struve_h0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STRUVE_H0_VALUES_TEST tests STRUVE_H0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STRUVE_H0_VALUES_TEST:\n" );
  printf ( "  STRUVE_H0_VALUES stores values of\n" );
  printf ( "  the Struve H0 function.\n" );
  printf ( "\n" );
  printf ( "      X            H0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    struve_h0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void struve_h1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STRUVE_H1_VALUES_TEST tests STRUVE_H1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STRUVE_H1_VALUES_TEST:\n" );
  printf ( "  STRUVE_H1_VALUES stores values of\n" );
  printf ( "  the Struve H1 function.\n" );
  printf ( "\n" );
  printf ( "      X            H1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   struve_h1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void struve_l0_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STRUVE_L0_VALUES_TEST tests STRUVE_L0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STRUVE_L0_VALUES_TEST:\n" );
  printf ( "  STRUVE_L0_VALUES stores values of\n" );
  printf ( "  the Struve L0 function.\n" );
  printf ( "\n" );
  printf ( "      X            L0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    struve_l0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void struve_l1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STRUVE_L1_VALUES_TEST tests STRUVE_L1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STRUVE_L1_VALUES_TEST:\n" );
  printf ( "  STRUVE_L1_VALUES stores values of\n" );
  printf ( "  the Struve L1 function.\n" );
  printf ( "\n" );
  printf ( "      X            L1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   struve_l1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void student_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STUDENT_CDF_VALUES_TEST tests STUDENT_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2005

  Author:

    John Burkardt
*/
{
  double c;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STUDENT_CDF_VALUES_TEST:\n" );
  printf ( "  STUDENT_CDF_VALUES returns values of\n" );
  printf ( "  the Student T Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      C     X       CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   student_cdf_values ( &n_data, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  16f  %16f  %24.16e\n", c, x, fx );
  }
  return;
}
/******************************************************************************/

void student_noncentral_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    STUDENT_NONCENTRAL_CDF_VALUES_TEST tests STUDENT_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "STUDENT_NONCENTRAL_CDF_VALUES_TEST:\n" );
  printf ( "  STUDENT_NONCENTRAL_CDF_VALUES returns values of\n" );
  printf ( "  the noncentral Student T Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "    DF     LAMBDA        X        CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %8f  %24.16e\n", df, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void subfactorial_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SUBFACTORIAL_VALUES_TEST tests SUBFACTORIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 March 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "SUBFACTORIAL_VALUES_TEST:\n" );
  printf ( "  SUBFACTORIAL_VALUES returns values of\n" );
  printf ( "  the subfactorial function.\n" );
  printf ( "\n" );
  printf ( "      N       Subfactorial[N]\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    subfactorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void surten_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    SURTEN_VALUES_TEST tests SURTEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double sigma;
  double tc;

  printf ( "\n" );
  printf ( "SURTEN_VALUES_TEST:\n" );
  printf ( "  SURTEN_VALUES stores values of\n" );
  printf ( "  the surface tension of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            SIGMA(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   surten_values ( &n_data, &tc, &sigma );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", tc, sigma );
  }
  return;
}
/******************************************************************************/

void synch1_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SYNCH1_VALUES_TEST tests SYNCH1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SYNCH1_VALUES_TEST:\n" );
  printf ( "  SYNCH1_VALUES stores values of \n" );
  printf ( "  the Synchrotron function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    synch1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void synch2_values_test ( )

/******************************************************************************/
/*
  Purpose:

    SYNCH2_VALUES_TEST tests SYNCH2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "SYNCH2_VALUES_TEST:\n" );
  printf ( "  SYNCH2_VALUES stores values of \n" );
  printf ( "  the Synchrotron function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    synch2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tan_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TAN_VALUES_TEST tests TAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TAN_VALUES_TEST:\n" );
  printf ( "   TAN_VALUES stores values of the tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tanh_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TANH_VALUES_TEST tests TANH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TANH_VALUES_TEST:\n" );
  printf ( "   TANH_VALUES stores values of the hyperbolic tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tanh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tau_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TAU_VALUES_TEST tests TAU_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TAU_VALUES_TEST:\n" );
  printf ( "  TAU_VALUES returns values of\n" );
  printf ( "  the TAU function.\n" );
  printf ( "\n" );
  printf ( "     N         TAU(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   tau_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void thercon_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    THERCON_VALUES_TEST tests THERCON_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double lambda;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "THERCON_VALUES_TEST:\n" );
  printf ( "  THERCON_VALUES stores values of\n" );
  printf ( "  the thermal conductivity of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            LAMBDA(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   thercon_values ( &n_data, &tc, &p, &lambda );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, lambda );
  }
  return;
}
/******************************************************************************/

void three_j_values_test ( )

/******************************************************************************/
/*
  Purpose:

    THREE_J_VALUES_TEST tests THREE_J_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  printf ( "\n" );
  printf ( "THREE_J_VALUES_TEST:\n" );
  printf ( "  THREE_J_VALUES returns values of\n" );
  printf ( "  the Wigner 3J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      M1      M2      M3        THREE_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    three_j_values ( &n_data, &j1, &j2, &j3, &m1, &m2, &m3, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, m1, m2, m3, fx );
  }
  return;
}
/******************************************************************************/

void tran02_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN02_VALUES_TEST tests TRAN02_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN02_VALUES_TEST:\n" );
  printf ( "  TRAN02_VALUES stores values of \n" );
  printf ( "  the Transport function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran02_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran03_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN03_VALUES_TEST tests TRAN03_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN03_VALUES_TEST:\n" );
  printf ( "  TRAN03_VALUES stores values of \n" );
  printf ( "  the Transport function of order 3.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran03_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran04_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN04_VALUES_TEST tests TRAN04_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN04_VALUES_TEST:\n" );
  printf ( "  TRAN04_VALUES stores values of \n" );
  printf ( "  the Transport function of order 4.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran04_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran05_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN05_VALUES_TEST tests TRAN05_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN05_VALUES_TEST:\n" );
  printf ( "  TRAN05_VALUES stores values of \n" );
  printf ( "  the Transport function of order 5.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran05_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran06_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN06_VALUES_TEST tests TRAN06_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN06_VALUES_TEST:\n" );
  printf ( "  TRAN06_VALUES stores values of \n" );
  printf ( "  the Transport function of order 6.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran06_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran07_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN07_VALUES_TEST tests TRAN07_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN07_VALUES_TEST:\n" );
  printf ( "  TRAN07_VALUES stores values of \n" );
  printf ( "  the Transport function of order 7.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran07_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran08_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN08_VALUES_TEST tests TRAN08_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN08_VALUES_TEST:\n" );
  printf ( "  TRAN08_VALUES stores values of \n" );
  printf ( "  the Transport function of order 8.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran08_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void tran09_values_test ( )

/******************************************************************************/
/*
  Purpose:

    TRAN09_VALUES_TEST tests TRAN09_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRAN09_VALUES_TEST:\n" );
  printf ( "  TRAN09_VALUES stores values of \n" );
  printf ( "  the Transport function of order 9.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran09_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void trigamma_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRIGAMMA_VALUES_TEST tests TRIGAMMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TRIGAMMA_VALUES_TEST\n" );
  printf ( "  TRIGAMMA_VALUES stores values of\n" );
  printf ( "  the TriGamma function.\n" );
  printf ( "\n" );
  printf ( "      X            FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    trigamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_ab_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_AB_CDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_CDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_AB_CDF_VALUES stores values of\n" );
  printf ( "  the Truncated Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       A         B         X        CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_cdf_values ( &n_data, &mu, &sigma, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_ab_pdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_AB_PDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_PDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_AB_PDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_AB_PDF_VALUES stores values of\n" );
  printf ( "  the Truncated Normal Probability Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       A         B         X        PDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_pdf_values ( &n_data, &mu, &sigma, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_a_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_A_CDF_VALUES_TEST tests TRUNCATED_NORMAL_A_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_CDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_A_CDF_VALUES stores values of\n" );
  printf ( "  the Lower Truncated Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       A         X        CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_cdf_values ( &n_data, &mu, &sigma, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, a, x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_a_pdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_A_PDF_VALUES_TEST tests TRUNCATED_NORMAL_A_PDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_A_PDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_A_PDF_VALUES stores values of\n" );
  printf ( "  the Lower Truncated Normal Probability Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       A         X        PDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_pdf_values ( &n_data, &mu, &sigma, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, a, x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_b_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_B_CDF_VALUES_TEST tests TRUNCATED_NORMAL_B_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_CDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_B_CDF_VALUES stores values of\n" );
  printf ( "  the Upper Truncated Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       B         X        CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_cdf_values ( &n_data, &mu, &sigma, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, b, x, fx );
  }
  return;
}
/******************************************************************************/

void truncated_normal_b_pdf_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TRUNCATED_NORMAL_B_PDF_VALUES_TEST tests TRUNCATED_NORMAL_B_PDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2013

  Author:

    John Burkardt
*/
{
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TRUNCATED_NORMAL_B_PDF_VALUES_TEST:\n" );
  printf ( "  TRUNCATED_NORMAL_B_PDF_VALUES stores values of\n" );
  printf ( "  the Upper Truncated Normal Probability Density Function.\n" );
  printf ( "\n" );
  printf ( "        MU     SIGMA       B         X        PDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_pdf_values ( &n_data, &mu, &sigma, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8.1f  %8.1f  %8.1f  %8.1f  %24.16e\n", mu, sigma, b, x, fx );
  }
  return;
}
/******************************************************************************/

void tsat_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    TSAT_VALUES_TEST tests TSAT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TSAT_VALUES_TEST:\n" );
  printf ( "  TSAT_VALUES stores values of\n" );
  printf ( "  the saturation temperature\n" );
  printf ( "  as a function of pressure.\n" );
  printf ( "\n" );
  printf ( "      P           Tsat(P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   tsat_values ( &n_data, &p, &tc );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", p, tc );
  }
  return;
}
/******************************************************************************/

void van_der_corput_values_test ( )

/******************************************************************************/
/*
  Purpose:

    VAN_DER_CORPUT_VALUES_TEST tests VAN_DER_CORPUT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int base;
  int n_data;
  int seed;
  double value;

  printf ( "\n" );
  printf ( "VAN_DER_CORPUT_VALUES_TEST:\n" );
  printf ( "  VAN_DER_CORPUT_VALUES stores values of\n" );
  printf ( "  the van der Corput sequence in a given base.\n" );
  printf ( "\n" );
  printf ( "      BASE      SEED    VDC(BASE,SEED)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    van_der_corput_values ( &n_data, &base, &seed, &value );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %14f\n", base, seed, value );
  }

  return;
}
/******************************************************************************/

void viscosity_values_test ( )

/******************************************************************************/
/*
  Purpose: 

    VISCOSITY_VALUES_TEST tests VISCOSITY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double eta;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "VISCOSITY_VALUES_TEST:\n" );
  printf ( "  VISCOSITY_VALUES stores values of\n" );
  printf ( "  the viscosity of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            ETA(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   viscosity_values ( &n_data, &tc, &p, &eta );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, eta );
  }
  return;
}
/******************************************************************************/

void von_mises_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    VON_MISES_CDF_VALUES_TEST tests VON_MISES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "VON_MISES_CDF_VALUES_TEST:\n" );
  printf ( "  VON_MISES_CDF_VALUES stores values of\n" );
  printf ( "  the von Mises CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    von_mises_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void weekday_values_test ( )

/******************************************************************************/
/*
  Purpose:

    WEEKDAY_VALUES_TEST tests WEEKDAY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 January 2015

  Author:

    John Burkardt
*/
{
  int d;
  int m;
  int n_data;
  int w;
  int y;

  printf ( "\n" );
  printf ( "WEEKDAY_VALUES_TEST:\n" );
  printf ( "  WEEKDAY_VALUES returns values of \n" );
  printf ( "  the weekday for a given Y/M/D date\n" );
  printf ( "\n" );
  printf ( "     Y     M     D     W\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    weekday_values ( &n_data, &y, &m, &d, &w );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %4d  %4d  %4d  %4d\n", y, m, d, w );
  }
  return;
}
/******************************************************************************/

void weibull_cdf_values_test ( )

/******************************************************************************/
/*
  Purpose:

    WEIBULL_CDF_VALUES_TEST tests WEIBULL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "WEIBULL_CDF_VALUES_TEST:\n" );
  printf ( "  WEIBULL_CDF_VALUES returns values of \n" );
  printf ( "  the Weibull Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Alpha   Beta        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    weibull_cdf_values ( &n_data, &alpha, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", alpha, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void zeta_values_test ( )

/******************************************************************************/
/*
  Purpose:

    ZETA_VALUES_TEST tests ZETA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int n;
  int n_data;
  double zeta;

  printf ( "\n" );
  printf ( "ZETA_VALUES_TEST:\n" );
  printf ( "  ZETA_VALUES returns values of \n" );
  printf ( "  the Riemann Zeta function.\n" );
  printf ( "\n" );
  printf ( "     N        ZETA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &zeta );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24e\n", n, zeta );
  }
  return;
}
