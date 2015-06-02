# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <math.h>
# include <complex.h>

# include "polpak.h"

int main ( );

void agud_test ( );
void align_enum_test ( );
void bell_test ( );
void benford_test ( );
void bernoulli_number_test ( );
void bernoulli_number2_test ( );
void bernoulli_number3_test ( );
void bernoulli_poly_test ( );
void bernoulli_poly2_test ( );
void bernstein_poly_test ( );
void bpab_test ( );
void cardan_poly_test ( );
void cardan_poly_coef_test ( );
void cardinal_cos_test ( );
void cardinal_sin_test ( );
void catalan_test ( );
void catalan_row_next_test ( );
void charlier_test ( );
void cheby_t_poly_test ( );
void cheby_t_poly_coef_test ( );
void cheby_t_poly_zero_test ( );
void cheby_u_poly_test ( );
void cheby_u_poly_coef_test ( );
void cheby_u_poly_zero_test ( );
void chebyshev_discrete_test ( );
void collatz_count_test ( );
void collatz_count_max_test ( );
void comb_row_next_test ( );
void commul_test ( );
void complete_symmetric_poly_test ( );
void cos_power_int_test ( );
void euler_number_test ( );
void euler_number2_test ( );
void euler_poly_test ( );
void eulerian_test ( );
void f_hofstadter_test ( );
void fibonacci_direct_test ( );
void fibonacci_floor_test ( );
void fibonacci_recursive_test ( );
void g_hofstadter_test ( );
void gegenbauer_poly_test ( );
void gen_hermite_poly_test ( );
void gen_laguerre_poly_test ( );
void gud_test ( );
void h_hofstadter_test ( );
void hail_test ( );
void hermite_poly_phys_test ( );
void hermite_poly_phys_coef_test ( );
void i4_choose_test ( );
void i4_factor_test ( );
void i4_factorial_test ( );
void i4_factorial2_test ( );
void i4_is_triangular_test ( );
void i4_partition_distinct_count_test ( );
void i4_to_triangle_test ( );
void jacobi_poly_test ( );
void jacobi_symbol_test ( );
void krawtchouk_test ( );
void laguerre_associated_test ( );
void laguerre_poly_test ( );
void laguerre_poly_coef_test ( );
void legendre_associated_test ( );
void legendre_associated_normalized_test ( );
void legendre_function_q_test ( );
void legendre_poly_test ( );
void legendre_poly_coef_test ( );
void legendre_symbol_test ( );
void lerch_test ( );
void lgamma_test ( );
void lock_test ( );
void meixner_test ( );
void mertens_test ( );
void moebius_test ( );
void motzkin_test ( );
void normal_01_cdf_inverse_test ( );
void omega_test ( );
void pentagon_num_test ( );
void phi_test ( );
void plane_partition_num_test ( );
void poly_bernoulli_test ( );
void poly_coef_count_test ( );
void prime_test ( );
void pyramid_num_test ( );
void pyramid_square_num_test ( );
void r8_agm_test ( );
void r8_beta_test ( );
void r8_choose_test ( );
void r8_erf_test ( );
void r8_erf_inverse_test ( );
void r8_euler_constant_test ( );
void r8_factorial_test ( );
void r8_factorial_log_test ( );
void r8_hyper_2f1_test ( );
void r8_psi_test ( );
void r8poly_degree_test ( );
void r8poly_print_test ( );
void r8poly_value_horner_test ( );
void sigma_test ( );
void simplex_num_test ( );
void sin_power_int_test ( );
void slice_test ( );
void spherical_harmonic_test ( );
void stirling1_test ( );
void stirling2_test ( );
void tau_test ( );
void tetrahedron_num_test ( );
void triangle_num_test ( );
void triangle_to_i4_test ( );
void trinomial_test ( );
void v_hofstadter_test ( );
void vibonacci_test ( );
void zeckendorf_test ( );
void zernike_poly_test ( );
void zernike_poly_coef_test ( );
void zeta_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POLPAK_PRB.

  Discussion:

    POLPAK_PRB tests the POLPAK library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLPAK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the POLPAK library.\n" );

  agud_test ( );
  align_enum_test ( );
  bell_test ( );
  benford_test ( );
  bernoulli_number_test ( );
  bernoulli_number2_test ( );
  bernoulli_number3_test ( );
  bernoulli_poly_test ( );
  bernoulli_poly2_test ( );
  bernstein_poly_test ( );
  bpab_test ( );
  cardan_poly_test ( );
  cardan_poly_coef_test ( );
  cardinal_cos_test ( );
  cardinal_sin_test ( );
  catalan_test ( );
  catalan_row_next_test ( );
  charlier_test ( );
  cheby_t_poly_test ( );
  cheby_t_poly_coef_test ( );
  cheby_t_poly_zero_test ( );
  cheby_u_poly_test ( );
  cheby_u_poly_coef_test ( );
  cheby_u_poly_zero_test ( );
  chebyshev_discrete_test ( );
  collatz_count_test ( );
  collatz_count_max_test ( );
  comb_row_next_test ( );
  commul_test ( );
  complete_symmetric_poly_test ( );
  cos_power_int_test ( );
  euler_number_test ( );
  euler_number2_test ( );
  euler_poly_test ( );
  eulerian_test ( );
  f_hofstadter_test ( );
  fibonacci_direct_test ( );
  fibonacci_floor_test ( );
  fibonacci_recursive_test ( );
  g_hofstadter_test ( );
  gegenbauer_poly_test ( );
  gen_hermite_poly_test ( );
  gen_laguerre_poly_test ( );
  gud_test ( );
  h_hofstadter_test ( );
  hail_test ( );
  hermite_poly_phys_test ( );
  hermite_poly_phys_coef_test ( );
  i4_choose_test ( );
  i4_factor_test ( );
  i4_factorial_test ( );
  i4_factorial2_test ( );
  i4_is_triangular_test ( );
  i4_partition_distinct_count_test ( );
  i4_to_triangle_test ( );
  jacobi_poly_test ( );
  jacobi_symbol_test ( );
  krawtchouk_test ( );
  laguerre_associated_test ( );
  laguerre_poly_test ( );
  laguerre_poly_coef_test ( );
  legendre_associated_test ( );
  legendre_associated_normalized_test ( );
  legendre_function_q_test ( );
  legendre_poly_test ( );
  legendre_poly_coef_test ( );
  legendre_symbol_test ( );
  lerch_test ( );
  lgamma_test ( );
  lock_test ( );
  meixner_test ( );
  mertens_test ( );
  moebius_test ( );
  motzkin_test ( );
  normal_01_cdf_inverse_test ( );
  omega_test ( );
  pentagon_num_test ( );
  phi_test ( );
  plane_partition_num_test ( );
  poly_bernoulli_test ( );
  poly_coef_count_test ( );
  prime_test ( );
  pyramid_num_test ( );
  pyramid_square_num_test ( );
  r8_agm_test ( );
  r8_beta_test ( );
  r8_choose_test ( );
  r8_erf_test ( );
  r8_erf_inverse_test ( );
  r8_euler_constant_test ( );
  r8_factorial_test ( );
  r8_factorial_log_test ( );
  r8_hyper_2f1_test ( );
  r8_psi_test ( );
  r8poly_degree_test ( );
  r8poly_print_test ( );
  r8poly_value_horner_test ( );
  sigma_test ( );
  simplex_num_test ( );
  sin_power_int_test ( );
  slice_test ( );
  spherical_harmonic_test ( );
  stirling1_test ( );
  stirling2_test ( );
  tau_test ( );
  tetrahedron_num_test ( );
  triangle_num_test ( );
  triangle_to_i4_test ( );
  trinomial_test ( );
  v_hofstadter_test ( );
  vibonacci_test ( );
  zeckendorf_test ( );
  zernike_poly_test ( );
  zernike_poly_coef_test ( );
  zeta_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLPAK_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void agud_test ( )

/******************************************************************************/
/*
  Purpose:

    AGUD_TEST tests AGUD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2010

  Author:

    John Burkardt
*/
{
  double g;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "AGUD_TEST\n" );
  printf ( "  AGUD computes the inverse Gudermannian;\n" );
  printf ( "\n" );
  printf ( "         X     GUD(X)     AGUD(GUD(X))\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    g = gud ( x );
    x2 = agud ( g );

    printf ( "  %10f  %10f  %10f\n", x, g, x2 );
  }

  return;
}
/******************************************************************************/

void align_enum_test ( )

/******************************************************************************/
/*
  Purpose:

    ALIGN_ENUM_TEST tests ALIGN_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
# define M_MAX 10
# define N_MAX 10

  int i;
  int j;

  printf ( "\n" );
  printf ( "ALIGN_ENUM_TEST\n" );
  printf ( "  ALIGN_ENUM counts the number of possible\n" );
  printf ( "  alignments of two biological sequences.\n" );

  printf ( "\n" );
  printf ( "  Alignment enumeration table:\n" );
  printf ( "\n" );

  printf ( "      " );
  for ( j = 0; j <= 5; j++ )
  {
    printf ( "%8d  ", j );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( i = 0; i <= M_MAX; i++ )
  {
    printf ( "   %2d  ", i );
    for ( j = 0; j <= 5; j++ )
    {
      printf ( "%8d  ", align_enum ( i, j ) );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "      " );
  for ( j = 6; j <= N_MAX; j++ )
  {
    printf ( "%8d  ", j );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( i = 0; i <= M_MAX; i++ )
  {
    printf ( "  %2d  ", i );
    for ( j = 6; j <= N_MAX; j++ )
    {
      printf ( "%8d  ", align_enum ( i, j ) );
    }
    printf ( "\n" );
  }
  return;
# undef M_MAX
# undef N_MAX
}
/******************************************************************************/

void bell_test ( )

/******************************************************************************/
/*
  Purpose:

    BELL_TEST tests BELL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BELL_TEST\n" );
  printf ( "  BELL computes Bell numbers.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    bell ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void benford_test ( )

/******************************************************************************/
/*
  Purpose:

    BENFORD_TEST tests BENFORD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "BENFORD_TEST\n" );
  printf ( "  BENFORD(I) is the Benford probability of the\n" );
  printf ( "  initial digit sequence I.\n" );
  printf ( "\n" );
  printf ( "     I  BENFORD(I)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 9; i++ )
  {
    printf ( "  %4d  %10.4f\n", i, benford ( i ) );
  }

  return;
}
/******************************************************************************/

void bernoulli_number_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_NUMBER_TEST tests BERNOULLI_NUMBER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BERNOULLI_NUMBER_TEST\n" );
  printf ( "  BERNOULLI_NUMBER computes Bernoulli numbers;\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number ( n, c1 );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1[n] );
  }

  return;
}
/******************************************************************************/

void bernoulli_number2_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_NUMBER2_TEST tests BERNOULLI_NUMBER2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BERNOULLI_NUMBER2_TEST\n" );
  printf ( "  BERNOULLI_NUMBER2 computes Bernoulli numbers;\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number2 ( n, c1 );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1[n] );
  }

  return;
}
/******************************************************************************/

void bernoulli_number3_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_NUMBER3_TEST tests BERNOULLI_NUMBER3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "BERNOULLI_NUMBER3_TEST\n" );
  printf ( "  BERNOULLI_NUMBER3 computes Bernoulli numbers;\n" );
  printf ( "\n" );
  printf ( "   I      Exact     BERNOULLI_NUMBER3\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    c1 = bernoulli_number3 ( n );

    printf ( "  %4d  %10g  %10g\n", n, c0, c1 );
  }

  return;
}
/******************************************************************************/

void bernoulli_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_POLY_TEST tests BERNOULLI_POLY;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  printf ( "\n" );
  printf ( "BERNOULLI_POLY_TEST\n" );
  printf ( "  BERNOULLI_POLY evaluates Bernoulli polynomials;\n" );
  printf ( "\n" );
  printf ( "  X = %g\n", x );
  printf ( "\n" );
  printf ( "  I          BX\n" );
  printf ( "\n" );

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly ( i, x );

    printf ( "  %6d  %10g\n", i, bx );
  }

  return;
}
/******************************************************************************/

void bernoulli_poly2_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNOULLI_POLY2_TEST tests BERNOULLI_POLY2;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  printf ( "\n" );
  printf ( "BERNOULLI_POLY2_TEST\n" );
  printf ( "  BERNOULLI_POLY2 evaluates Bernoulli polynomials;\n" );
  printf ( "\n" );
  printf ( "  X = %g\n", x );
  printf ( "\n" );
  printf ( "  I          BX\n" );
  printf ( "\n" );

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly2 ( i, x );

    printf ( "  %6d  %10g\n", i, bx );
  }

  return;
}
/******************************************************************************/

void bernstein_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    BERNSTEIN_POLY_TEST tests BERNSTEIN_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double b;
  double bvec[11];
  int k;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "BERNSTEIN_POLY_TEST:\n" );
  printf ( "  BERNSTEIN_POLY evaluates the Bernstein polynomials.\n" );
  printf ( "\n" );
  printf ( "   N   K   X   Exact   B(N,K)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }

    bernstein_poly ( n, x, bvec );

    printf ( "  %4d  %4d  %7g  %14g  %14g\n", n, k, x, b, bvec[k] );
  }

  return;
}
/******************************************************************************/

void bpab_test ( )

/******************************************************************************/
/*
  Purpose:

    BPAB_TEST tests BPAB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double a;
  double b;
  double bern[N+1];
  int i;
  double x;

  printf ( "\n" );
  printf ( "BPAB_TEST\n" );
  printf ( "  BPAB evaluates Bernstein polynomials.\n" );
  printf ( "\n" );

  x = 0.3;
  a = 0.0;
  b = 1.0;

  bpab ( N, x, a, b, bern );

  printf ( "  The Bernstein polynomials of degree %d\n", N );
  printf ( "  based on the interval from %g\n", a );
  printf ( "  to %g\n", b );
  printf ( "  evaluated at X = %g\n", x );
  printf ( "\n" );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %14g\n", i, bern[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void cardan_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    CARDAN_POLY_TEST tests CARDAN_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  double *cx2;
  int i;
  int n;
  double s;
  double x;

  n = N_MAX;
  x = 0.25;
  s = 0.5;

  printf ( "\n" );
  printf ( "CARDAN_POLY_TEST\n" );
  printf ( "  CARDAN_POLY evaluates the Cardan polynomial.\n" );
  printf ( "\n" );
  printf ( "  Compare CARDAN_POLY_COEF + R8POLY_VALUE_HORNER\n" );
  printf ( "  versus CARDAN_POLY alone.\n" );
  printf ( "\n" );
  printf ( "  Evaluate polynomials at X = %g\n", x );
  printf ( "  We use the parameter S = %g\n", s );
  printf ( "\n" );
  printf ( "  Order       Horner          Direct\n" );
  printf ( "\n" );

  cx2 = cardan_poly ( n, x, s );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );

    cx1 = r8poly_value_horner ( n, c, x );

    printf ( "  %2d  %14g  %14g\n", n, cx1, cx2[n] );
  }
  free ( cx2 );

  return;
# undef N_MAX
}
/******************************************************************************/

void cardan_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    CARDAN_POLY_COEF_TEST tests CARDAN_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  int i;
  int n;
  double s;
  double x;

  s = 1.0;

  printf ( "\n" );
  printf ( "CARDAN_POLY_COEF_TEST\n" );
  printf ( "  CARDAN_POLY_COEF returns the coefficients of a\n" );
  printf ( "  Cardan polynomial.\n" );
  printf ( "\n" );
  printf ( "  We use the parameter S = %g\n", s );
  printf ( "\n" );
  printf ( "  Table of polynomial coefficients:\n" );
  printf ( "\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );
    printf ( "  %2d  ", n );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%5d  ", c[i] );
    }
    printf ( "\n" );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void cardinal_cos_test ( )

/******************************************************************************/
/*
  Purpose:

    CARDINAL_COS_TEST tests CARDINAL_COS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int m = 11;
  const double r8_pi = 3.141592653589793;
  double *t;

  printf ( "\n" );
  printf ( "CARDINAL_COS_TEST\n" );
  printf ( "  CARDINAL_COS evaluates cardinal cosine functions.\n" );
  printf ( "  Ci(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).\n" );
  printf ( "  A simple check of all pairs should form the identity matrix.\n" );

  printf ( "\n" );
  printf ( "  The CARDINAL_COS test matrix:\n" );
  printf ( "\n" );

  t = r8vec_linspace_new ( m + 2, 0.0, r8_pi );

  for ( j = 0; j <= m + 1; j++ )
  {
    c = cardinal_cos ( j, m, m + 2, t );
    for ( i = 0; i <= m + 1; i++ )
    {
      printf ( "  %4.1f", c[i] );
    }
    printf ( "\n" );
    free ( c );
  }

  free ( t );

  return;
}
/******************************************************************************/

void cardinal_sin_test ( )

/******************************************************************************/
/*
  Purpose:

    CARDINAL_SIN_TEST tests CARDINAL_SIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int m = 11;
  const double r8_pi = 3.141592653589793;
  double *s;
  double *t;

  printf ( "\n" );
  printf ( "CARDINAL_SIN_TEST\n" );
  printf ( "  CARDINAL_SIN evaluates cardinal sine functions.\n" );
  printf ( "  Si(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).\n" );
  printf ( "  A simple check of all pairs should form the identity matrix.\n" );

  printf ( "\n" );
  printf ( "  The CARDINAL_SIN test matrix:\n" );
  printf ( "\n" );

  t = r8vec_linspace_new ( m + 2, 0.0, r8_pi );

  for ( j = 0; j <= m + 1; j++ )
  {
    s = cardinal_sin ( j, m, m + 2, t );
    for ( i = 0; i <= m + 1; i++ )
    {
      printf ( "  %4.1f", s[i] );
    }
    printf ( "\n" );
    free ( s );
  }

  free ( t );

  return;
}
/******************************************************************************/

void catalan_test ( )

/******************************************************************************/
/*
  Purpose:

    CATALAN_TEST tests CATALAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "CATALAN_TEST\n" );
  printf ( "  CATALAN computes Catalan numbers.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    catalan ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void catalan_row_next_test ( )

/******************************************************************************/
/*
  Purpose:

    CATALAN_ROW_NEXT_TEST tests CATALAN_ROW_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  bool next;

  printf ( "\n" );
  printf ( "CATALAN_ROW_NEXT_TEST\n" );
  printf ( "  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.\n" );
  printf ( "\n" );
  printf ( "  First, compute row 7:\n" );
  printf ( "\n" );

  next = false;
  n = 7;
  catalan_row_next ( next, n, c );

  printf ( "%4d  ", n );
  for ( i = 0; i <= n; i++ )
  {
    printf ( "%8d  ", c[i] );
  }
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  Now compute rows consecutively, one at a time:\n" );
  printf ( "\n" );

  next = false;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = true;

    printf ( "%4d  ", i );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%8d  ", c[i] );
    }
    printf ( "\n" );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void charlier_test ( )

/******************************************************************************/
/*
  Purpose:

    CHARLIER_TEST tests CHARLIER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5
# define N 5

  double a;
  double a_test[TEST_NUM] = { 0.25, 0.5, 1.0, 2.0, 10.0 };
  int i;
  int j;
  int n;
  int test;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "CHARLIER_TEST:\n" );
  printf ( "  CHARLIER evaluates Charlier polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      A         X        P(N,A,X)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    a = a_test[test];

    printf ( "\n" );

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      charlier ( n, a, x, value );

      printf ( "\n" );
      for ( i = 0; i <= 5; i++ )
      {

        printf ( "  %6d  %8g  %8g  %14g\n", i, a, x, value[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void cheby_t_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_T_POLY_TEST tests CHEBY_T_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "CHEBY_T_POLY_TEST:\n" );
  printf ( "  CHEBY_T_POLY evaluates the Chebyshev T polynomial.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       T(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_t_poly ( 1, n, x_vec );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void cheby_t_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_T_POLY_COEF_TEST tests CHEBY_T_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 April 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 5;

  printf ( "\n" );
  printf ( "CHEBY_T_POLY_COEF_TEST\n" );
  printf ( "  CHEBY_T_POLY_COEF determines the  polynomial coefficients\n" );
  printf ( "  of the Chebyshev polynomial T(n,x).\n" );

  c = cheby_t_poly_coef ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  T(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(n+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(n+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(n+1)], j );
        }
      }
    }
  }
 
  free ( c );

  return;
}
/******************************************************************************/

void cheby_t_poly_zero_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_T_POLY_ZERO_TEST tests CHEBY_T_POLY_ZERO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  printf ( "\n" );
  printf ( "CHEBY_T_POLY_ZERO_TEST:\n" );
  printf ( "  CHEBY_T_POLY_ZERO returns zeroes of T(N,X).\n" );
  printf ( "\n" );
  printf ( "       N      X        T(N,X)\n" );
  printf ( "\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_t_poly_zero ( n );
    fx = cheby_t_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[i+n*n] );
    }
    printf ( "\n" );
    free ( fx );
    free ( z );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void cheby_u_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_U_POLY_TEST tests CHEBY_U_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 January 2015

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "CHEBY_U_POLY_TEST:\n" );
  printf ( "  CHEBY_U_POLY evaluates the Chebyshev U polynomial.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       U(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_u_poly ( 1, n, x_vec );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void cheby_u_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_U_POLY_COEF_TEST tests CHEBY_U_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "CHEBY_U_POLY_COEF_TEST\n" );
  printf ( "  CHEBY_U_POLY_COEF determines the polynomial coefficients\n" );
  printf ( "  of the Chebyshev polynomial U(n,x).\n" );

  cheby_u_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  U(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(N+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(N+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
        }
      }
    }
  }
 
  return;
# undef N
}
/******************************************************************************/

void cheby_u_poly_zero_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBY_U_POLY_ZERO_TEST tests CHEBY_U_POLY_ZERO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 January 2015

  Author:

    John Burkardt
*/
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  printf ( "\n" );
  printf ( "CHEBY_U_POLY_ZERO_TEST:\n" );
  printf ( "  CHEBY_U_POLY_ZERO returns zeroes of U(N,X).\n" );
  printf ( "\n" );
  printf ( "       N      X        U(N,X)\n" );
  printf ( "\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_u_poly_zero ( n );
    fx = cheby_u_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[i+n*n] );
    }
    printf ( "\n" );
    free ( fx );
    free ( z );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void chebyshev_discrete_test ( )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV_DISCRETE_TEST tests CHEBYSHEV_DISCRETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5
# define N 5

  int i;
  int j;
  int m;
  int n;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "CHEBYSHEV_DISCRETE_TEST:\n" );
  printf ( "  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      M         X        T(N,M,X)\n" );

  m = 5;
  n = N;

  for ( j = 0; j <= 5; j++ )
  {
    x = ( double ) ( j ) / 2.0;

    chebyshev_discrete ( n, m, x, value );

    printf ( "\n" );
    for ( i = 0; i <= 5; i++ )
    {
      printf ( "  %6d  %6d  %8g  %14g\n", i, m, x, value[i] );
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void collatz_count_test ( )

/******************************************************************************/
/*
  Purpose:

    COLLATZ_COUNT_TEST tests COLLATZ_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  int count;
  int count2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "COLLATZ_COUNT_TEST:\n" );
  printf ( "  COLLATZ_COUNT(N) counts the length of the\n" );
  printf ( "  Collatz sequence beginning with N.\n" );
  printf ( "\n" );
  printf ( "       N       COUNT(N)     COUNT(N)\n" );
  printf ( "              (computed)    (table)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }

    count2 = collatz_count ( n );

    printf ( "  %8d  %8d  %8d\n", n, count, count2 );
  }

  return;
}
/******************************************************************************/

void collatz_count_max_test ( )

/******************************************************************************/
/*
  Purpose:

    COLLATZ_COUNT_MAX_TEST tests COLLATZ_COUNT_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 May 2012

  Author:

    John Burkardt
*/
{
  int i_max;
  int j_max;
  int n;

  printf ( "\n" );
  printf ( "COLLATZ_COUNT_MAX_TEST:\n" );
  printf ( "  COLLATZ_COUNT_MAX(N) returns the length of the\n" );
  printf ( "  longest Collatz sequence from 1 to N.\n" );
  printf ( "\n" );
  printf ( "         N     I_MAX     J_MAX\n" );
  printf ( "\n" );

  n = 10;

  while ( n <= 100000 )
  {
    collatz_count_max ( n, &i_max, &j_max );

    printf ( "  %8d  %8d  %8d\n", n, i_max, j_max );

    n = n * 10;
  }

  return;
}
/******************************************************************************/

void comb_row_next_test ( )

/******************************************************************************/
/*
  Purpose:

    COMB_ROW_NEXT_TEST tests COMB_ROW_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 December 2014

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;

  printf ( "\n" );
  printf ( "COMB_ROW_NEXT_TEST\n" );
  printf ( "  COMB_ROW_NEXT computes the next row of Pascal's triangle.\n" );
  printf ( "\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    comb_row_next ( n, c );
    printf ( "  %2d  ", n );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%5d", c[i] );
    }
    printf ( "\n" );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void commul_test ( )

/******************************************************************************/
/*
  Purpose:

    COMMUL_TEST tests COMMUL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 November 2013

  Author:

    John Burkardt
*/
{
  int n;
  int factor[4];
  int i;
  int ncomb;
  int nfactor;

  printf ( "\n" );
  printf ( "COMMUL_TEST\n" );
  printf ( "  COMMUL computes a multinomial coefficient.\n" );
  printf ( "\n" );

  n = 8;
  nfactor = 2;
  factor[0] = 6;
  factor[1] = 2;
  ncomb = commul ( n, nfactor, factor );
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of factors = %d\n", nfactor );
  for ( i = 0; i < nfactor; i++ )
  {
    printf ( "  %2d  %8d\n", i, factor[i] );
  }
  printf ( "  Value of coefficient = %d\n", ncomb );

  n = 8;
  nfactor = 3;
  factor[0] = 2;
  factor[1] = 2;
  factor[2] = 4;
  ncomb = commul ( n, nfactor, factor );
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of factors = %d\n", nfactor );
  for ( i = 0; i < nfactor; i++ )
  {
    printf ( "  %2d  %8d\n", i, factor[i] );
  }
  printf ( "  Value of coefficient = %d\n", ncomb );

  n = 13;
  nfactor = 4;
  factor[0] = 5;
  factor[1] = 3;
  factor[2] = 3;
  factor[3] = 2;
  ncomb = commul ( n, nfactor, factor );
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of factors = %d\n", nfactor );
  for ( i = 0; i < nfactor; i++ )
  {
    printf ( "  %2d  %8d\n", i, factor[i] );
  }
  printf ( "  Value of coefficient = %d\n", ncomb );

  return;
}
/******************************************************************************/

void complete_symmetric_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    TEST02407 tests COMPLETE_SYMMETRIC_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 November 2013

  Author:

    John Burkardt
*/
{
  int n = 5;
  int nn;
  int r;
  int rr;
  double value;
  double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  printf ( "\n" );
  printf ( "COMPLETE_SYMMETRIC_POLY_TEST\n" );
  printf ( "  COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric.\n" );
  printf ( "  polynomial in a given set of variables X.\n" );
 
  r8vec_print ( n, x, "  Variable vector X:" );

  printf ( "\n" );
  printf ( "   N\\R     0       1       2       3       4       5\n" );
  printf ( "\n" );

  for ( nn = 0; nn <= n; nn++ )
  {
    printf ( "  %2d", nn );
    for ( rr = 0; rr <= 5; rr++ )
    {
      value = complete_symmetric_poly ( nn, rr, x );
      printf ( "  %6.0f", value );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void cos_power_int_test ( )

/******************************************************************************/
/*
  Purpose:

    COS_POWER_INT_TEST tests COS_POWER_INT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "COS_POWER_INT_TEST:\n" );
  printf ( "  COS_POWER_INT computes the integral of the N-th power\n" );
  printf ( "  of the cosine function.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        Exact    Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = cos_power_int ( a, b, n );

    printf ( "  %8g  %8g  %6d  %12g  %12g\n", a, b, n, fx, fx2 );
  }
  return;
}
/******************************************************************************/

void euler_number_test ( )

/******************************************************************************/
/*
  Purpose:

    EULER_NUMBER_TEST tests EULER_NUMBER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int c1;
  int c2[13];
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "EULER_NUMBER_TEST\n" );
  printf ( "  EULER_NUMBER computes Euler numbers.\n" );
  printf ( "\n" );
  printf ( "  N  exact   EULER_NUMBER\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    euler_number ( n, c2 );

    printf ( "  %4d  %12d  %12d\n", n, c1, c2[n] );

  }
 
  return;
}
/******************************************************************************/

void euler_number2_test ( )

/******************************************************************************/
/*
  Purpose:

    EULER_NUMBER2_TEST tests EULER_NUMBER2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int c1;
  int c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "EULER_NUMBER2_TEST\n" );
  printf ( "  EULER_NUMBER2 computes Euler numbers.\n" );
  printf ( "\n" );
  printf ( "  N  exact   EULER_NUMBER2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = euler_number2 ( n );

    printf ( "  %4d  %12d  %12d\n", n, c1, c2 );

  }
 
  return;
}
/******************************************************************************/

void euler_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    EULER_POLY_TEST tests EULER_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double f;
  int i;
  int n = 15;
  double x;

  x = 0.5;
 
  printf ( "\n" );
  printf ( "EULER_POLY_TEST\n" );
  printf ( "  EULER_POLY evaluates Euler polynomials.\n" );
  printf ( "\n" );
  printf ( "  N         X              F(X)\n" );
  printf ( "\n" );
   
  for ( i = 0; i <= n; i++ )
  {
    f = euler_poly ( i, x );

    printf ( "  %2d  %14g  %14g\n", i, x, f );
  }
 
  return;
}
/******************************************************************************/

void eulerian_test ( )

/******************************************************************************/
/*
  Purpose:

    EULERIAN_TEST tests EULERIAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
# define N 7

  int e[N*N];
  int i;
  int j;

  printf ( "\n" );
  printf ( "EULERIAN_TEST\n" );
  printf ( "  EULERIAN evaluates Eulerian numbers.\n" );
  printf ( "\n" );
 
  eulerian ( N, e );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      printf ( "%6d  ", e[i+j*N] );
    }
    printf ( "\n" );
  }
 
  return;
# undef N
}
/******************************************************************************/

void f_hofstadter_test ( )

/******************************************************************************/
/*
  Purpose:

    F_HOFSTADTER_TEST tests F_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;

  printf ( "\n" );
  printf ( "F_HOFSTADTER_TEST\n" );
  printf ( "  F_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  F function.\n" );
  printf ( "\n" );
  printf ( "     N   F(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    f = f_hofstadter ( i );

    printf ( "  %6d  %6d\n", i, f );
  }

  return;
}
/******************************************************************************/

void fibonacci_direct_test ( )

/******************************************************************************/
/*
  Purpose:

    FIBONACCI_DIRECT_TEST tests FIBONACCI_DIRECT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;
  int n = 20;

  printf ( "\n" );
  printf ( "FIBONACCI_DIRECT_TEST\n" );
  printf ( "  FIBONACCI_DIRECT evalutes a Fibonacci number directly.\n" );
  printf ( "\n" );
  
  for ( i = 1; i <= n; i++ )
  {
    f = fibonacci_direct ( i );

    printf ( "  %6d  %10d\n", i, f );
  }
 
  return;
}
/******************************************************************************/

void fibonacci_floor_test ( )

/******************************************************************************/
/*
  Purpose:

    FIBONACCI_FLOOR_TEST tests FIBONACCI_FLOOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int f;
  int i;
  int n;

  printf ( "\n" );
  printf ( "FIBONACCI_FLOOR_TEST\n" );
  printf ( "  FIBONACCI_FLOOR computes the largest Fibonacci number\n" );
  printf ( "  less than or equal to a given positive integer.\n" );
  printf ( "\n" );
  printf ( "     N  Fibonacci  Index\n" );
  printf ( "\n" );

  for ( n = 1; n <= 20; n++ )
  {
    fibonacci_floor ( n, &f, &i );

    printf ( "  %6d  %6d  %6d\n", n, f, i );
  }
 
  return;
}
/******************************************************************************/

void fibonacci_recursive_test ( )

/******************************************************************************/
/*
  Purpose:

    FIBONACCI_RECURSIVE_TEST tests FIBONACCI_RECURSIVE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
# define N 20

  int f[N];
  int i;

  printf ( "\n" );
  printf ( "FIBONACCI_RECURSIVE_TEST\n" );
  printf ( "  FIBONACCI_RECURSIVE computes the Fibonacci sequence.\n" );
  printf ( "\n" );
 
  fibonacci_recursive ( N, f );
 
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %10d\n", i, f[i] );
  }
 
  return;
# undef N
}
/******************************************************************************/

void g_hofstadter_test ( )

/******************************************************************************/
/*
  Purpose:

    G_HOFSTADTER_TEST tests G_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "G_HOFSTADTER_TEST\n" );
  printf ( "  G_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  G function.\n" );
  printf ( "\n" );
  printf ( "     N   G(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, g_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void gegenbauer_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    GEGENBAUER_POLY_TEST tests GEGENBAUER_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GEGENBAUER_POLY_TEST\n" );
  printf ( "  GEGENBAUER_POLY evaluates the Gegenbauer polynomials.\n" );
  printf ( "\n" );
  printf ( "        N       A       X       GPV      GEGENBAUER\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {

    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

    gegenbauer_poly ( n, a, x, c );
    fx2 = c[n];

    printf ( "  %6d  %10g  %10g  %14g  %14g\n", n, a, x, fx, fx2 );

    free ( c );
  }
 
  return;
}
/******************************************************************************/

void gen_hermite_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    GEN_HERMITE_POLY_TEST tests GEN_HERMITE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2015

  Author:

    John Burkardt
*/
{
# define N 10
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  double mu;
  double mu_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  printf ( "\n" );
  printf ( "GEN_HERMITE_POLY_TEST\n" );
  printf ( "  GEN_HERMITE_POLY evaluates the generalized Hermite\n" );
  printf ( "  polynomials.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    mu = mu_test[i];

    printf ( "\n" );
    printf ( "  Table of H(N,MU)(X) for\n" );
    printf ( "\n" );
    printf ( "    N(max) = %d\n", N );
    printf ( "    MU =     %g\n", mu );
    printf ( "    X =      %g\n", x );
    printf ( "\n" );
  
    gen_hermite_poly ( N, x, mu, c );
 
    for ( j = 0; j <= N; j++ )
    {
      printf ( "  %6d  %14g\n", j, c[j] );
    }
  }
 
  return;
# undef N
# undef N_TEST
}
/******************************************************************************/

void gen_laguerre_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    GEN_LAGUERRE_POLY_TEST tests GEN_LAGUERRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
# define N 10
# define N_TEST 6

  double alpha;
  double alpha_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double c[N+1];
  int i;
  int j;
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  printf ( "\n" );
  printf ( "GEN_LAGUERRE_POLY_TEST\n" );
  printf ( "  GEN_LAGUERRE_POLY evaluates the generalized Laguerre\n" );
  printf ( "  functions.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    alpha = alpha_test[i];

    printf ( "\n" );
    printf ( "  Table of L(N,ALPHA,X) for\n" );
    printf ( "\n" );
    printf ( "    N(max) = %d\n", N );
    printf ( "    ALPHA =  %g\n", alpha );
    printf ( "    X =      %g\n", x );
    printf ( "\n" );
  
    gen_laguerre_poly ( N, alpha, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      printf ( "  %6d  %14g\n", j, c[j] );
    }
  }
 
  return;
# undef N
# undef N_TEST
}
/******************************************************************************/

void gud_test ( )

/******************************************************************************/
/*
  Purpose:

    GUD_TEST tests GUD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "GUD_TEST:\n" );
  printf ( "  GUD evaluates the Gudermannian function.\n" );
  printf ( "\n" );
  printf ( "     X      Exact F       GUD(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gud ( x );

    printf ( "  %10.6g  %10.6g  %10.6g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void h_hofstadter_test ( )

/******************************************************************************/
/*
  Purpose:

    H_HOFSTADTER_TEST tests H_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "H_HOFSTADTER_TEST\n" );
  printf ( "  H_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  H function.\n" );

  printf ( "\n" );
  printf ( "     N   H(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, h_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void hail_test ( )

/******************************************************************************/
/*
  Purpose:

    HAIL_TEST tests HAIL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "HAIL_TEST\n" );
  printf ( "  HAIL(I) computes the length of the hail sequence\n" );
  printf ( "  for I, also known as the 3*N+1 sequence.\n" );
  printf ( "\n" );
  printf ( "  I,  HAIL(I)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    printf ( "  %4d  %6d\n", i,  hail ( i ) );
  }
 
  return;
}
/******************************************************************************/

void hermite_poly_phys_test ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLY_PHYS_TEST tests HERMITE_POLY_PHYS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "HERMITE_POLY_PHYS_TEST:\n" );
  printf ( "  HERMITE_POLY_PHYS evaluates the physicist's Hermite polynomial.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       H(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_phys_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    hermite_poly_phys ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void hermite_poly_phys_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLY_PHYS_COEF_TEST tests HERMITE_POLY_PHYS_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "HERMITE_POLY_PHYS_COEF_TEST\n" );
  printf ( "  HERMITE_POLY_PHYS_COEF determines physicist's Hermite polynomial coefficients.\n" );

  hermite_poly_phys_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  H(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }

  return;
# undef N
}
/******************************************************************************/

void i4_choose_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_CHOOSE_TEST tests I4_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  int cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "I4_CHOOSE_TEST\n" );
  printf ( "  I4_CHOOSE evaluates C(N,K).\n" );
  printf ( "\n" );
  printf ( "   N     K    CNK\n" );
  printf ( "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      printf ( "  %6d  %6d  %6d\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void i4_factor_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FACTOR_TEST tests I4_FACTOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 February 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int maxfactor = 10;
  int n;
  int n_test[3] = { 60, 664048, 8466763 };
  int nfactor;
  int nleft;
  int factor[10];
  int power[10];

  printf ( "\n" );
  printf ( "I4_FACTOR_TEST:\n" );
  printf ( "  I4_FACTOR tries to factor an I4\n" );

  for ( i = 0; i < 3; i++ )
  {
    n = n_test[i];
    i4_factor ( n, maxfactor, &nfactor, factor, power, &nleft );
    printf ( "\n" );
    printf ( "  Factors of N = %d\n", n );
    for ( j = 0; j < nfactor; j++ )
    {
      printf ( "    %d^%d\n", factor[j], power[j] );
    }
    if ( nleft != 1 )
    {
      printf ( "  Unresolved factor NLEFT = %d\n", nleft );
    }
  }

  return;
}
/******************************************************************************/

void i4_factorial_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FACTORIAL_TEST tests I4_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  int fn;
  int fn2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "I4_FACTORIAL_TEST:\n" );
  printf ( "  I4_FACTORIAL evaluates the factorial function.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       I4_FACTORIAL(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial ( n );

    printf ( "  %4d  %12d  %12d\n", n, fn, fn2 );

  }

  return;
}
/******************************************************************************/

void i4_factorial2_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_FACTORIAL2_TEST tests I4_FACTORIAL2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int fn;
  int fn2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "I4_FACTORIAL2_TEST:\n" );
  printf ( "  I4_FACTORIAL2 evaluates the double factorial function.\n" );
  printf ( "\n" );
  printf ( "   N   Exact  I4_FACTORIAL2(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial2 ( n );

    printf ( "  %4d  %8d  %8d\n", n, fn, fn2 );
  }

  return;
}
/******************************************************************************/

void i4_is_triangular_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_IS_TRIANGULAR_TEST tests I4_IS_TRIANGULAR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 December 2014

  Author:

    John Burkardt
*/
{
  int i;
  bool l;

  printf ( "\n" );
  printf ( "I4_IS_TRIANGULAR_TEST\n" );
  printf ( "  I4_IS_TRIANGULAR returns 0 or 1 depending on\n" );
  printf ( "  whether I is triangular.\n" );
  printf ( "\n" );
  printf ( "   I  =>   0/1\n" );
  printf ( "\n" );

  for ( i = 0; i <= 20; i++ )
  {
    l = i4_is_triangular ( i );

    printf ( "  %4d  %1d\n", i, l );
  }
 
  return;
}
/******************************************************************************/

void i4_partition_distinct_count_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_PARTITION_DISTINCT_COUNT_TEST tests I4_PARTITION_DISTINCT_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int c2;
  int n;
  int n_data;
  int n_max = 20;

  printf ( "\n" );
  printf ( "I4_PARTITION_DISTINCT_COUNT_TEST:\n" );
  printf ( "  For the number of partitions of an integer\n" );
  printf ( "  into distinct parts,\n" );
  printf ( "  I4_PARTITION_DISTINCT_COUNT computes any value.\n" );
  printf ( "\n" );
  printf ( "     N       Exact F    Q(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n_max < n )
    {
      continue;
    }

    c2 = i4_partition_distinct_count ( n );

    printf ( "  %10d  %10d  %10d\n", n, c, c2 );
  }

  return;
}
/******************************************************************************/

void i4_to_triangle_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_TO_TRIANGLE_TEST tests I4_TO_TRIANGLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "I4_TO_TRIANGLE_TEST\n" );
  printf ( "  I4_TO_TRIANGLE converts a linear index to a\n" );
  printf ( "  triangular one.\n" );
  printf ( "\n" );
  printf ( "     K  =>   I     J\n" );
  printf ( "\n" );

  for ( k = 0; k <= 20; k++ )
  {
    i4_to_triangle ( k, &i, &j );

    printf ( "  %4d    %4d  %4d\n", k, i, j );
  }
 
  return;
}
/******************************************************************************/

void jacobi_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_POLY_TEST tests JACOBI_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "JACOBI_POLY_TEST:\n" );
  printf ( "  JACOBI_POLY computes values of the Jacobi polynomial..\n" );
  printf ( "\n" );
  printf ( "       N       A       B      X       JPV      JACOBI\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_poly_values ( &n_data, &n, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = jacobi_poly ( n, a, b, x );
    fx2 = c[n];

    printf ( "  %8d  %8f  %8f  %10.4f  %14.6f  %14.6f\n", n, a, b, x, fx, fx2 );

    free ( c );
  }

  return;
}
/******************************************************************************/

void jacobi_symbol_test ( )

/******************************************************************************/
/*
  Purpose:

    JACOBI_SYMBOL_TEST tests JACOBI_SYMBOL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N_TEST 4

  int i;
  int p;
  int ptest[N_TEST] = { 3, 9, 10, 12 };
  int q;

  printf ( "\n" );
  printf ( "JACOBI_SYMBOL_TEST\n" );
  printf ( "  JACOBI_SYMBOL computes the Jacobi symbol\n" );
  printf ( "  (Q/P), which records if Q is a quadratic\n" );
  printf ( "  residue modulo the number P.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    printf ( "\n" );
    printf ( "Jacobi Symbols for P = %d\n", p );
    printf ( "\n" );
    for ( q = 0; q <= p; q++ )
    {
      printf ( "  %8d  %8d  %8d\n", p, q, jacobi_symbol ( q, p ) );
    }
  }

  return;
# undef N_TEST
}
/******************************************************************************/

void krawtchouk_test ( )

/******************************************************************************/
/*
  Purpose:

    KRAWTCHOUK_TEST tests KRAWTCHOUK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2
# define N 5

  int i;
  int j;
  int m;
  int n;
  double p;
  double p_test[TEST_NUM] = { 0.25, 0.5 };
  int test;
  double x;
  double value[N+1];

  printf ( "\n" );
  printf ( "KRAWTCHOUK_TEST:\n" );
  printf ( "  KRAWTCHOUK evaluates Krawtchouk polynomials.\n" );
  printf ( "\n" );
  printf ( "        N         P         X          M      K(N,P,X,M)\n" );
  printf ( "\n" );

  m = 5;
  n = N;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];

    printf ( "\n" );

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      krawtchouk ( n, p, x, m, value );

      printf ( "\n" );
      for ( i = 0; i <= 5; i++ )
      {

        printf ( "  %8d  %8g  %8g  %8d  %14g\n", i, p, x, m, value[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void laguerre_associated_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_ASSOCIATED_TEST tests LAGUERRE_ASSOCIATED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 6
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  int m;
  int m_test[N_TEST] = { 0, 0, 1, 2, 3, 1 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  printf ( "\n" );
  printf ( "LAGUERRE_ASSOCIATED_TEST\n" );
  printf ( "  LAGUERRE_ASSOCIATED evaluates the associated Laguerre\n" );
  printf ( "  polynomials.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    m = m_test[i];
    x = x_test[i];

    printf ( "\n" );
    printf ( "  Table of L(N,M,X) for\n" );
    printf ( "\n" );
    printf ( "  N(max) = %d\n", N );
    printf ( "  M      = %d\n", m );
    printf ( "  X =      %g\n", x );
    printf ( "\n" );
 
    laguerre_associated ( N, m, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      printf ( "  %6d  %14g\n", j, c[j] );
    }
  }

  return;
# undef N
# undef N_TEST
}
/******************************************************************************/

void laguerre_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLY_TEST tests LAGUERRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LAGUERRE_POLY_TEST:\n" );
  printf ( "  LAGUERRE_POLY evaluates the Laguerre polynomial.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       L(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    laguerre_poly ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void laguerre_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLY_COEF_TEST tests LAGUERRE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  double fact;
  int i;
  int j;

  printf ( "\n" );
  printf ( "LAGUERRE_POLY_COEF_TEST\n" );
  printf ( "  LAGUERRE_POLY_COEF determines Laguerre \n" );
  printf ( "  polynomial coefficients.\n" );

  laguerre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  L(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
 
  for ( i = 0; i <= N; i++ )
  {
    fact = r8_factorial ( i );
    printf ( "\n" );
    printf ( "  Factorially scaled L(%d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
  return;
# undef N
}
/******************************************************************************/

void legendre_associated_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ASSOCIATED_TEST tests LEGENDRE_ASSOCIATED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_ASSOCIATED_TEST:\n" );
  printf ( "  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.\n" );
  printf ( "\n" );
  printf ( "      N       M    X     Exact F     PNM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated ( n, m, x, fx2 );

    printf ( "  %8d  %8d  %8f  %14f  %14f\n", n, m, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void legendre_associated_normalized_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ASSOCIATED_NORMALIZED_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 February 2015

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_ASSOCIATED_NORMALIZED_TEST:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED evaluates \n" );
  printf ( "  normalized associated Legendre functions.\n" );
  printf ( "\n" );
  printf ( "      N       M    X     Exact F     PNM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_sphere_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated_normalized ( n, m, x, fx2 );

    printf ( "  %8d  %8d  %8f  %14f  %14f\n", n, m, x, fx, fx2[n] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void legendre_function_q_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_FUNCTION_Q_TEST tests LEGENDRE_FUNCTION_Q.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_FUNCTION_Q_TEST:\n" );
  printf ( "  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       Q(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_function_q ( n, x, fx2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void legendre_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_POLY_TEST tests LEGENDRE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 12

  double fx;
  double fp2[N_MAX+1];
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LEGENDRE_POLY_TEST:\n" );
  printf ( "  LEGENDRE_POLY evaluates the Legendre PN function.\n" );
  printf ( "\n" );
  printf ( "     N      X        Exact F       P(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_poly ( n, x, fx2, fp2 );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

  }

  return;
# undef N_MAX
}
/******************************************************************************/

void legendre_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_POLY_COEF_TEST tests LEGENDRE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  printf ( "\n" );
  printf ( "LEGENDRE_POLY_COEF_TEST\n" );
  printf ( "  LEGENDRE_POLY_COEF determines the Legendre P \n" );
  printf ( "  polynomial coefficients.\n" );

  legendre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  P(%8d)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
      }
    }
  }
 
  return;
# undef N
}
/******************************************************************************/

void legendre_symbol_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_SYMBOL_TEST tests LEGENDRE_SYMBOL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N_TEST 4

  int i;
  int l;
  int p;
  int ptest[N_TEST] = { 7, 11, 13, 17 };
  int q;

  printf ( "\n" );
  printf ( "LEGENDRE_SYMBOL_TEST\n" );
  printf ( "  LEGENDRE_SYMBOL computes the Legendre\n" );
  printf ( "  symbol (Q/P) which records whether Q is \n" );
  printf ( "  a quadratic residue modulo the prime P.\n" );

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    printf ( "\n" );
    printf ( "  Legendre Symbols for P = %d\n", p );
    printf ( "\n" );
    for ( q = 0; q <= p; q++ )
    {
      printf ( "  %8d  %8d  %8d\n", p, q, legendre_symbol ( q, p ) );
    }
  }

  return;
# undef N_TEST
}
/******************************************************************************/

void lerch_test ( )

/******************************************************************************/
/*
  Purpose:

    LERCH_TEST tests LERCH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int s;
  double z;

  printf ( "\n" );
  printf ( "LERCH_TEST:\n" );
  printf ( "  LERCH evaluates the Lerch function.\n" );
  printf ( "\n" );
  printf ( "       Z       S       A         Lerch           Lerch\n" );
  printf ( "                             Tabulated        Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lerch ( z, s, a );

    printf ( "  %8g  %4d  %8g  %14g  %14g\n", z, s, a, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void lgamma_test ( )

/******************************************************************************/
/*
  Purpose:

    LGAMMA_TEST tests LGAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "LGAMMA_TEST:\n" );
  printf ( "  LGAMMA is a C math library function which evaluates\n" );
  printf ( "  the logarithm of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       LGAMMA(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lgamma ( x );

    printf ( "  %8g  %10g  %10g\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void lock_test ( )

/******************************************************************************/
/*
  Purpose:

    LOCK_TEST tests LOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int a[N+1];
  int i;

  printf ( "\n" );
  printf ( "LOCK_TEST\n" );
  printf ( "  LOCK counts the combinations on a button lock.\n" );
  printf ( "\n" );
  printf ( "     I      LOCK(I)\n" );
  printf ( "\n" );

  lock ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %10d\n", i, a[i] );
  }
 
  return;
# undef N
}
/******************************************************************************/

void meixner_test ( )

/******************************************************************************/
/*
  Purpose:

    MEIXNER_TEST tests MEIXNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
# define N 5
# define TEST_NUM 3

  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.0 };
  double c;
  double c_test[TEST_NUM] = { 0.125, 0.25, 0.5 };
  int i;
  int j;
  int n;
  int test;
  double v[N+1];
  double x;

  printf ( "\n" );
  printf ( "MEIXNER_TEST:\n" );
  printf ( "  MEIXNER evaluates Meixner polynomials.\n" );
  printf ( "\n" );
  printf ( "       N      BETA         C         X        M(N,BETA,C,X)\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    beta = beta_test[test];
    c = c_test[test];

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      meixner ( n, beta, c, x, v );

      printf ( "\n" );

      for ( i = 0; i <= n; i++ )
      {
        printf ( "  %8d  %8g  %8g  %8g  %14g\n", i, beta, c, x, v[i] );
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void mertens_test ( )

/******************************************************************************/
/*
  Purpose:

    MERTENS_TEST tests MERTENS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "MERTENS_TEST\n" );
  printf ( "  MERTENS computes the Mertens function.\n" );
  printf ( "\n" );
  printf ( "      N   Exact   MERTENS(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
     mertens_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %10d  %10d\n", n, c, mertens ( n ) );
  }
 
  return;
}
/******************************************************************************/

void moebius_test ( )

/******************************************************************************/
/*
  Purpose:

    MOEBIUS_TEST tests MOEBIUS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "MOEBIUS_TEST\n" );
  printf ( "  MOEBIUS computes the Moebius function.\n" );
  printf ( "\n" );
  printf ( "      N   Exact   MOEBIUS(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
     moebius_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %8d  %10d  %10d\n", n, c, moebius ( n ) );
  }
 
  return;
}
/******************************************************************************/

void motzkin_test ( )

/******************************************************************************/
/*
  Purpose:

    MOTZKIN_TEST tests MOTZKIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int a[N+1];
  int i;

  printf ( "\n" );
  printf ( "MOTKIN_TEST\n" );
  printf ( "  MOTZKIN computes the Motzkin numbers A(0:N).\n" );
  printf ( "  A(N) counts the paths from (0,0) to (N,0).\n" );
  printf ( "\n" );
  printf ( "  I,  A(I)\n" );
  printf ( "\n" );

  motzkin ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "  %4d  %10d\n", i, a[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void normal_01_cdf_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_CDF_INVERSE_TEST tests NORMAL_01_CDF_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 February 2015

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "NORMAL_01_CDF_INVERSE_TEST:\n" );
  printf ( "  NORMAL_01_CDF_INVERSE inverts the Normal 01 CDF.\n" );
  printf ( "\n" );
  printf ( "    FX           X1           X2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x1, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = normal_01_cdf_inverse ( fx );

    printf ( "  %8f  %14f  %14f\n", fx, x1, x2 );
  }

  return;
}
/******************************************************************************/

void omega_test ( )

/******************************************************************************/
/*
  Purpose:

    OMEGA_TEST tests OMEGA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "OMEGA_TEST\n" );
  printf ( "  OMEGA computes the OMEGA function.\n" );
  printf ( "\n" );
  printf ( "          N   Exact   OMEGA(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %10d  %10d\n", n, c, omega ( n ) );
  }
 
  return;
}
/******************************************************************************/

void pentagon_num_test ( )

/******************************************************************************/
/*
  Purpose:

    PENTAGON_NUM_TEST tests PENTAGON_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "PENTAGON_NUM_TEST\n" );
  printf ( "  PENTAGON_NUM computes the pentagonal numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, pentagon_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void phi_test ( )

/******************************************************************************/
/*
  Purpose:

    PHI_TEST tests PHI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "PHI_TEST\n" );
  printf ( "  PHI computes the PHI function.\n" );
  printf ( "\n" );
  printf ( "  N   Exact   PHI(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, phi ( n ) );

  }
 
  return;
}
/******************************************************************************/

void plane_partition_num_test ( )

/******************************************************************************/
/*
  Purpose:

    PLANE_PARTITION_NUM_TEST tests PLANE_PARTITION_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 February 2015

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "PLANE_PARTITION_NUM_TEST\n" );
  printf ( "  PLANE_PARTITION_NUM counts the number of plane\n" );
  printf ( "  partitions of an integer.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, plane_partition_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void poly_bernoulli_test ( )

/******************************************************************************/
/*
  Purpose:

    POLY_BERNOULLI_TEST tests POLY_BERNOULLI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int b;
  int k;
  int n;

  printf ( "\n" );
  printf ( "POLY_BERNOULLI_TEST\n" );
  printf ( "  POLY_BERNOULLI computes the poly-Bernoulli numbers\n" );
  printf ( "  of negative index, B_n^(-k)\n" );
  printf ( "\n" );
  printf ( "   N   K    B_N^(-K)\n" );
  printf ( "\n" );

  for ( k = 0; k <= 6; k++ )
  {
    printf ( "\n" );
    for ( n = 0; n <= 6; n++ )
    {
      b = poly_bernoulli ( n, k );

      printf ( "  %2d  %2d  %12d\n", n, k, b );
    }
  }

  return;
}
/******************************************************************************/

void poly_coef_count_test ( )

/******************************************************************************/
/*
  Purpose:

    POLY_COEF_COUNT_TEST tests POLY_COEF_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int degree;
  int dim;
  int n;

  printf ( "\n" );
  printf ( "POLY_COEF_COUNT_TEST\n" );
  printf ( "  POLY_COEF_COUNT counts the number of coefficients\n" );
  printf ( "  in a polynomial of degree DEGREE and dimension DIM.\n" );
  printf ( "\n" );
  printf ( " Dimension    Degree     Count\n" );

  for ( dim = 1; dim <= 10; dim = dim + 3 )
  {
    printf ( "\n" );
    for ( degree = 0; degree <= 5; degree++ )
    {
      printf ( "  %8d  %8d  %8d\n", dim, degree, poly_coef_count ( dim, degree ) );
    }
  }

  return;
}
/******************************************************************************/

void prime_test ( )

/******************************************************************************/
/*
  Purpose:

    PRIME_TEST tests PRIME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2014

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int prime_max;

  printf ( "\n" );
  printf ( "PRIME_TEST\n" );
  printf ( "  PRIME returns primes from a table.\n" );

  n = -1;
  prime_max = prime ( n );
  printf ( "\n" );
  printf ( "  Number of primes stored is %d\n", prime_max );
  printf ( "\n" );
  printf ( "     I    Prime(I)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %4d  %6d\n", i, prime ( i ) );
  }
  printf ( "\n" );
  for ( i = prime_max - 10; i <= prime_max; i++ )
  {
    printf ( "  %4d  %6d\n", i, prime ( i ) );
  }
  
  return;
}
/******************************************************************************/

void pyramid_num_test ( )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_NUM_TEST tests PYRAMID_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "PYRAMID_NUM_TEST\n" );
  printf ( "  PYRAMID_NUM computes the pyramidal numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, pyramid_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void pyramid_square_num_test ( )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_SQUARE_NUM_TEST tests PYRAMID_SQUARE_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 December 2014

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "PYRAMID_SQUARE_NUM_TEST\n" );
  printf ( "  PYRAMID_SQARE_NUM computes the pyramidal square numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, pyramid_square_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void r8_agm_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_AGM_TEST tests R8_AGM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 December 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;

  printf ( "\n" );
  printf ( "R8_AGM_TEST\n" );
  printf ( "  R8_AGM computes the arithmetic geometric mean.\n" );
  printf ( "\n" );
  printf ( "      A           B          " );
  printf ( "   AGM                       AGM                   Diff" );
  printf ( "                             " );
  printf ( "  (Tabulated)             R8_AGM(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( &n_data, &a, &b, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_agm ( a, b );

    printf ( "  %10.6f  %10.6f  %24.16f  %24.16f  %10.6e\n",
      a, b, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void r8_beta_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_BETA_TEST tests R8_BETA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 January 2015

  Author:

    John Burkardt
*/
{
  double fxy;
  double fxy2;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "R8_BETA_TEST:\n" );
  printf ( "  R8_BETA evaluates the Beta function.\n" );
  printf ( "\n" );
  printf ( "     X      Y        Exact F       R8_BETA(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = r8_beta ( x, y );

    printf ( "  %10f  %10f  %10g  %10g\n", x, y, fxy, fxy2 );
  }

  return;
}
/******************************************************************************/

void r8_choose_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE_TEST tests R8_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
  double cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "R8_CHOOSE_TEST\n" );
  printf ( "  R8_CHOOSE evaluates C(N,K) using real arithmetic.\n" );
  printf ( "\n" );
  printf ( "       N       K          CNK\n" );
  printf ( "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );

      printf ( "  %6d  %6d  %10g\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void r8_erf_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ERF_TEST tests R8_ERF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 August 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_ERF_TEST:\n" );
  printf ( "  R8_ERF evaluates the error function.\n" );
  printf ( "\n" );
  printf ( "     X      Exact F       R8_ERF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_erf ( x );

    printf ( "  %8f  %14f  %14f\n", x, fx, fx2 );
  }

  return;
}
/******************************************************************************/

void r8_erf_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ERF_INVERSE_TEST tests R8_ERF_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2010

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "R8_ERF_INVERSE_TEST:\n" );
  printf ( "  R8_ERF_INVERSE inverts the error function.\n" );
  printf ( "\n" );
  printf ( "    FX           X1           X2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x1, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = r8_erf_inverse ( fx );

    printf ( "  %8f  %14f  %14f\n", fx, x1, x2 );
  }

  return;
}
/******************************************************************************/

void r8_euler_constant_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_EULER_CONSTANT_TEST tests R8_EULER_CONSTANT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt
*/
{
  double g;
  double g_approx;
  int i;
  int n;
  double n_r8;
  int test;

  g = r8_euler_constant ( );

  printf ( "\n" );
  printf ( "R8_EULER_CONSTANT_TEST:\n" );
  printf ( "  R8_EULER_CONSTANT returns the Euler-Mascheroni constant\n" );
  printf ( "  sometimes denoted by 'gamma'.\n" );
  printf ( "\n" );
  printf ( "  gamma = limit ( N -> oo ) ( sum ( 1 <= I <= N ) 1 / I ) - log ( N )\n" );
  printf ( "\n" );
  printf ( "  Numerically, g = %g\n", g );
  printf ( "\n" );
  printf ( "         N      Partial Sum    |gamma - partial sum|\n" );
  printf ( "\n" );

  n = 1;
  for ( test = 0; test <= 20; test++ )
  {
    n_r8 = ( double ) ( n );
    g_approx = - log ( n_r8 );
    for ( i = 1; i <= n; i++ )    
    {
      g_approx = g_approx + 1.0 / ( double ) ( i );
    }
    printf ( "  %8d  %14.6g  %14.6g\n", n, g_approx, fabs ( g_approx - g ) );
    n = n * 2;
  }

  return;
}
/******************************************************************************/

void r8_factorial_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL_TEST tests R8_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  double fn;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "R8_FACTORIAL_TEST:\n" );
  printf ( "  R8_FACTORIAL evaluates the factorial function.\n" );
  printf ( "\n" );
  printf ( "     N       Exact F       R8_FACTORIAL(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %14g  %14g\n", n, fn, r8_factorial ( n ) );
  }

  return;
}
/******************************************************************************/

void r8_factorial_log_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL_LOG_TEST tests R8_FACTORIAL_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
  double fn;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "R8_FACTORIAL_LOG_TEST:\n" );
  printf ( "  R8_FACTORIAL_LOG evaluates the logarithm of the\n" );
  printf ( "  factorial function.\n" );
  printf ( "\n" );
  printf ( "     N	   Exact F	 R8_FACTORIAL_LOG(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %5d  %14g  %14g\n", n, fn, r8_factorial_log ( n ) );
  }

  return;
}
/******************************************************************************/

void r8_hyper_2f1_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_HYPER_2F1_TEST tests R8_HYPER_2F1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_HYPER_2F1_TEST:\n" );
  printf ( "  R8_HYPER_2F1 evaluates the hypergeometric function 2F1.\n" );
  printf ( "\n" );
  printf ( "      A       B       C       X      " );
  printf ( " 2F1                       2F1                     DIFF\n" );
  printf ( "                                     " );
  printf ( "(tabulated)               (computed)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_hyper_2f1 ( a, b, c, x );

    printf ( "  %6f  %6f  %6f  %6f  %24.16g  %24.16g  %10.4g\n",
      a, b, c, x, fx, fx2, fabs ( fx - fx2 )  );
  }
  return;
}
/******************************************************************************/

void r8_psi_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_PSI_TEST tests R8_PSI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_PSI_TEST:\n" );
  printf ( "  R8_PSI evaluates the Psi function.\n" );
  printf ( "\n" );
  printf ( "         X                  Psi(X)           " );
  printf ( "         Psi(X)          DIFF\n" );
  printf ( "                         (Tabulated)         " );
  printf ( "       (R8_PSI)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_psi ( x );

    printf ( "  %8.2g  %24.16g  %24.16g  %10.4g\n", x, fx, fx2, fabs ( fx - fx2 ) );

  }

  return;
}
/******************************************************************************/

void r8poly_degree_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_DEGREE_TEST tests R8POLY_DEGREE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2015

  Author:

    John Burkardt
*/
{
  double c1[4] = { 1.0, 2.0, 3.0, 4.0 }; 
  double c2[4] = { 1.0, 2.0, 3.0, 0.0 };
  double c3[4] = { 1.0, 2.0, 0.0, 4.0 };
  double c4[4] = { 1.0, 0.0, 0.0, 0.0 };
  double c5[4] = { 0.0, 0.0, 0.0, 0.0 };
  int d;
  int m;
 
  printf ( "\n" );
  printf ( "R8POLY_DEGREE_TEST\n" );
  printf ( "  R8POLY_DEGREE determines the degree of an R8POLY.\n" );

  m = 3;

  r8poly_print ( m, c1, "  The R8POLY:" );
  d = r8poly_degree ( m, c1 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c2, "  The R8POLY:" );
  d = r8poly_degree ( m, c2 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c3, "  The R8POLY:" );
  d = r8poly_degree ( m, c3 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c4, "  The R8POLY:" );
  d = r8poly_degree ( m, c4 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c5, "  The R8POLY:" );
  d = r8poly_degree ( m, c5 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  return;
}
/******************************************************************************/

void r8poly_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_PRINT_TEST tests R8POLY_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 January 2015

  Author:

    John Burkardt
*/
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  printf ( "\n" );
  printf ( "R8POLY_PRINT_TEST\n" );
  printf ( "  R8POLY_PRINT prints an R8POLY.\n" );

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
/******************************************************************************/

void r8poly_value_horner_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2015

  Author:

    John Burkardt
*/
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  printf ( "\n" );
  printf ( "R8POLY_VALUE_HORNER_TEST\n" );
  printf ( "  R8POLY_VALUE_HORNER evaluates a polynomial at\n" );
  printf ( "  one point, using Horner's method.\n" );

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  printf ( "\n" );
  printf ( "   I    X    P(X)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    printf ( "  %2d  %8.4f  %14.6g\n", i, x[i], p );
  }

  free ( x );

  return;
}
/******************************************************************************/

void sigma_test ( )

/******************************************************************************/
/*
  Purpose:

    SIGMA_TEST tests SIGMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "SIGMA_TEST\n" );
  printf ( "  SIGMA computes the SIGMA function.\n" );
  printf ( "\n" );
  printf ( "  N   Exact   SIGMA(N)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    sigma_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, sigma ( n ) );
  }
 
  return;
}
/******************************************************************************/

void simplex_num_test ( )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX_NUM_TEST tests SIMPLEX_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 February 2015

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int value;

  printf ( "\n" );
  printf ( "SIMPLEX_NUM_TEST\n" );
  printf ( "  SIMPLEX_NUM computes the N-th simplex number\n" );
  printf ( "  in M dimensions.\n" );
  printf ( "\n" );
  printf ( "      M: 0     1     2     3     4     5\n" );
  printf ( "   N\n" );
 
  for ( n = 0; n <= 10; n++ )
  {
    printf ( "  %2d", n );
    for ( m = 0; m <= 5; m++ )
    {
      value = simplex_num ( m, n );
      printf ( "  %4d", value );
    }
    printf ( "\n" );
  } 
  return;
}
/******************************************************************************/

void sin_power_int_test ( )

/******************************************************************************/
/*
  Purpose:

    SIN_POWER_INT_TEST tests SIN_POWER_INT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "SIN_POWER_INT_TEST:\n" );
  printf ( "  SIN_POWER_INT computes the integral of the N-th power\n" );
  printf ( "  of the sine function.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        Exact    Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = sin_power_int ( a, b, n );

    printf ( "  %8g  %8g  %6d  %12g  %12g\n", a, b, n, fx, fx2 );
  }
  return;
}
/******************************************************************************/

void slice_test ( )

/******************************************************************************/
/*
  Purpose:

    SLICE_TEST tests SLICE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 August 2011

  Author:

    John Burkardt
*/
{
# define DIM_MAX 5
# define SLICE_MAX 8

  int dim_max = DIM_MAX;
  int dim_num;
  int p[DIM_MAX*SLICE_MAX];
  int piece_num;
  int slice_max = SLICE_MAX;
  int slice_num;

  printf ( "\n" );
  printf ( "SLICE_TEST:\n" );
  printf ( "  SLICE determines the maximum number of pieces created\n" );
  printf ( "  by SLICE_NUM slices in a DIM_NUM space.\n" );

  for ( dim_num = 1; dim_num <= dim_max; dim_num++ )
  {
    for ( slice_num = 1; slice_num <= slice_max; slice_num++ )
    {
      piece_num = slice ( dim_num, slice_num );
      p[dim_num-1+(slice_num-1)*dim_max] = piece_num;
    }
  }

  i4mat_print ( dim_max, slice_max, p, "  Slice Array:" );

  return;
# undef DIM_MAX
# undef SLICE_MAX
}
/******************************************************************************/

void spherical_harmonic_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERICAL_HARMONIC_TEST tests SPHERICAL_HARMONIC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double c[N_MAX+1];
  int l;
  int m;
  int n_data;
  double phi;
  double s[N_MAX+1];
  double theta;
  double yi;
  double yi2;
  double yr;
  double yr2;

  printf ( "\n" );
  printf ( "SPHERICAL_HARMONIC_TEST:\n" );
  printf ( "  SPHERICAL_HARMONIC evaluates spherical harmonic functions.\n" );
  printf ( "\n" );
  printf ( "         N         M    THETA      PHI            YR            YI\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }

    spherical_harmonic ( l, m, theta, phi, c, s );

    yr2 = c[l];
    yi2 = s[l];

    printf ( "  %8d  %8d  %8g  %8g  %14g  %14g\n", l, m, theta, phi, yr, yi );

    printf ( "                                          %14g  %14g\n", yr2, yi2 );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void stirling1_test ( )

/******************************************************************************/
/*
  Purpose:

    STIRLING1_TEST tests STIRLING1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m = 8;
  int n = 8;
  int *s1;

  printf ( "\n" );
  printf ( "STIRLING1_TEST\n" );
  printf ( "  STIRLING1: Stirling numbers of first kind.\n" );
  printf ( "  Get rows 1 through %d\n", m );
  printf ( "\n" );
 
  s1 = stirling1 ( m, n );
 
  for ( i = 0; i < m; i++ )
  {
    printf ( "%6d  ", i + 1 );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%6d  ", s1[i+j*m] );
    }
    printf ( "\n" );
  }

  free ( s1 );
 
  return;
}
/******************************************************************************/

void stirling2_test ( )

/******************************************************************************/
/*
  Purpose:

    STIRLING2_TEST tests STIRLING2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m = 8;
  int n = 8;
  int *s2;

  printf ( "\n" );
  printf ( "STIRLING2_TEST\n" );
  printf ( "  STIRLING2: Stirling numbers of second kind.\n" );
  printf ( "  Get rows 1 through %d\n", m );
  printf ( "\n" );
 
  s2 = stirling2 ( m, n );
 
  for ( i = 0; i < m; i++ )
  {
    printf ( "%6d  ", i + 1 );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%6d  ", s2[i+j*m] );
    }
    printf ( "\n" );
  }
 
  free ( s2 );

  return;
}
/******************************************************************************/

void tau_test ( )

/******************************************************************************/
/*
  Purpose:

    TAU_TEST tests TAU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TAU_TEST\n" );
  printf ( "  TAU computes the Tau function.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );
 
  n_data = 0;

  for ( ; ; )
  {
    tau_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %4d  %10d  %10d\n", n, c, tau ( n ) );
  }
 
  return;
}
/******************************************************************************/

void tetrahedron_num_test ( )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NUM_TEST tests TETRAHEDRON_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TETRAHEDRON_NUM_TEST\n" );
  printf ( "  TETRAHEDRON_NUM computes the tetrahedron numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, tetrahedron_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void triangle_num_test ( )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_NUM_TEST tests TRIANGLE_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TRIANGLE_NUM_TEST\n" );
  printf ( "  TRIANGLE_NUM computes the triangular numbers.\n" );
  printf ( "\n" );
 
  for ( n = 1; n <= 10; n++ )
  {
    printf ( "  %4d  %6d\n", n, triangle_num ( n ) );
  }
 
  return;
}
/******************************************************************************/

void triangle_to_i4_test ( )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_TO_I4_TEST tests TRIANGLE_TO_I4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "TRIANGLE_TO_I4_TEST\n" );
  printf ( "  TRIANGLE_TO_I4 converts a triangular index to a\n" );
  printf ( "  linear one.\n" );
  printf ( "\n" );
  printf ( "     I     J   ==> K\n" );
  printf ( "\n" );

  for ( i = 0; i <= 4; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      k = triangle_to_i4 ( i, j );
      printf ( "  %4d  %4d    %4d\n", i, j, k );
    }
  }
 
  return;
}
/******************************************************************************/

void trinomial_test ( )

/******************************************************************************/
/*
  Purpose:

    TRINOMIAL_TEST tests TRINOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int t;

  printf ( "\n" );
  printf ( "TRINOMIAL_TEST\n" );
  printf ( "  TRINOMIAL evaluates the trinomial coefficient:\n" );
  printf ( "\n" );
  printf ( "  T(I,J,K) = (I+J+K)! / I! / J! / K!\n" );
  printf ( "\n" );
  printf ( "     I     J     K    T(I,J,K)\n" );
  printf ( "\n" );
 
  for ( k = 0; k <= 4; k++ )
  {
    for ( j = 0; j <= 4; j++ )
    {
      for ( i = 0; i <= 4; i++ )
      {
        t = trinomial ( i, j, k );
        printf ( "  %4d  %4d  %4d  %8d\n", i, j, k, t );
      }
    }
  }
 
  return;
}
/******************************************************************************/

void v_hofstadter_test ( )

/******************************************************************************/
/*
  Purpose:

    V_HOFSTADTER_TEST tests V_HOFSTADTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int i;
  int v;

  printf ( "\n" );
  printf ( "V_HOFSTADTER_TEST\n" );
  printf ( "  V_HOFSTADTER evaluates Hofstadter's recursive\n" );
  printf ( "  V function.\n" );
  printf ( "\n" );
  printf ( "     N   V(N)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 30; i++ )
  {
    printf ( "  %6d  %6d\n", i, v_hofstadter ( i ) );
  }

  return;
}
/******************************************************************************/

void vibonacci_test ( )

/******************************************************************************/
/*
  Purpose:

    VIBONACCI_TEST tests VIBONACCI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
# define N 20
  int i;
  int j;
  int seed;
  int v1[N];
  int v2[N];
  int v3[N];

  printf ( "\n" );
  printf ( "VIBONACCI_TEST\n" );
  printf ( "  VIBONACCI computes a Vibonacci sequence.\n" );
  printf ( "\n" );
  printf ( "  We compute the series 3 times.\n" );
  printf ( "\n" );
  printf ( "     I      V1      V2      V3\n" );
  printf ( "\n" );

  seed = 123456789;

  vibonacci ( N, &seed, v1 );
  vibonacci ( N, &seed, v2 );
  vibonacci ( N, &seed, v3 );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6d  %6d  %6d\n", i, v1[i], v2[i], v3[i] );
  } 

  return;
# undef N
}
/******************************************************************************/

void zeckendorf_test ( )

/******************************************************************************/
/*
  Purpose:

    ZECKENDORF_TEST tests ZECKENDORF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 May 2012

  Author:

    John Burkardt
*/
{
# define M_MAX 20

  int i;
  int i_list[M_MAX];
  int j;
  int f_list[M_MAX];
  int f_sum;
  int m;
  int n;

  printf ( "\n" );
  printf ( "ZECKENDORF_TEST\n" );
  printf ( "  ZECKENDORF computes the Zeckendorf decomposition of\n" );
  printf ( "  an integer N into nonconsecutive Fibonacci numbers.\n" );
  printf ( "\n" );
  printf ( "   N Sum M Parts\n" );
  printf ( "\n" );

  for ( n = 1; n <= 100; n++ )
  {
    zeckendorf ( n, M_MAX, &m, i_list, f_list );

    printf ( "%4d  ", n );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%4d  ", f_list[j] );
    }
    printf ( "\n" );

  }

  return;
# undef M_MAX
}
/******************************************************************************/

void zernike_poly_test ( )

/******************************************************************************/
/*
  Purpose:

    ZERNIKE_POLY_TEST tests ZERNIKE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int m;
  int n;
  double rho;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "ZERNIKE_POLY_TEST\n" );
  printf ( "  ZERNIKE_POLY evaluates a Zernike polynomial directly.\n" );
  printf ( "\n" );
  printf ( "  Table of polynomial coefficients:\n" );
  printf ( "\n" );
  printf ( "   N   M\n" );
  printf ( "\n" );

  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      printf ( "  %2d  %2d", n, m );
      for ( i = 0; i <= n; i++ )
      {
        printf ( "  %7g", c[i] );
      }
      printf ( "\n" );
      free ( c );
    }
  }

  rho = 0.987654321;

  printf ( "\n" );
  printf ( "  Z1: Compute polynomial coefficients,\n" );
  printf ( "  then evaluate by Horner's method;\n" );
  printf ( "  Z2: Evaluate directly by recursion.\n" );
  printf ( "\n" );
  printf ( "   N   M       Z1              Z2\n" );
  printf ( "\n" );

  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      z1 = r8poly_value_horner ( n, c, rho );

      z2 = zernike_poly ( m, n, rho );
      printf ( "  %2d  %2d  %16g  %16g\n", n, m, z1, z2 );

      free ( c );
    }
  }

  return;
}
/******************************************************************************/

void zernike_poly_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    ZERNIKE_POLY_COEF_TEST tests ZERNIKE_POLY_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  double *c;
  int m;
  int n;

  printf ( "\n" );
  printf ( "ZERNIKE_POLY_COEF_TEST\n" );
  printf ( "  ZERNIKE_POLY_COEF determines the Zernike\n" );
  printf ( "  polynomial coefficients.\n" );

  n = 5;

  for ( m = 0; m <= n; m++ )
  {
    c = zernike_poly_coef ( m, n );
    r8poly_print ( n, c, "  Zernike polynomial" );
    free ( c );
  }

  return;
}
/******************************************************************************/

void zeta_test ( )

/******************************************************************************/
/*
  Purpose:

    ZETA_TEST tests ZETA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;
  double n_real;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "ZETA_TEST\n" );
  printf ( "  ZETA computes the Zeta function.\n" );
  printf ( "\n" );
  printf ( "       N            exact Zeta         computed Zeta\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &z1 );

    if ( n_data == 0 )
    {
      break;
    }

    n_real = ( double ) n;

    z2 = zeta ( n_real );

    printf ( "  %6d  %20.14g  %20.14g\n", n, z1, z2 );
  }

  return;
}
