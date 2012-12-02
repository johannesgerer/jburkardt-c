void bashforth_set ( int order, double xtab[], double weight[] );
void bdf_set ( int order, double alpha[], double *beta, double *gamma );
void bdfc_set ( int order, double xtab[], double weight[] );
void bdfp_set ( int order, double xtab[], double weight[] );
double bdf_sum ( double func ( double x ), int order, double xtab[], 
  double weight[] );
char ch_cap ( char ch );
void cheb_set ( int order, double xtab[], double weight[] );
void chebyshev1_compute ( int order, double xtab[], double weight[] );
double chebyshev1_integral ( int expon );
void chebyshev2_compute ( int order, double xtab[], double weight[] );
double chebyshev2_integral ( int expon );
void chebyshev3_compute ( int order, double xtab[], double weight[] );
void clenshaw_curtis_compute ( int n, double x[], double w[] );
void clenshaw_curtis_set ( int order, double xtab[], double weight[] );
void fejer1_compute ( int n, double x[], double w[] );
void fejer1_set ( int order, double xtab[], double weight[] );
void fejer2_compute ( int n, double x[], double w[] );
void fejer2_set ( int order, double xtab[], double weight[] );
void gegenbauer_compute ( int order, double alpha, double xtab[], 
  double weight[] );
double gegenbauer_integral ( int expon, double alpha );
void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double c[] );
void gegenbauer_root ( double *x, int order, double alpha,  double *dp2, 
  double *p1, double c[] );
void gen_hermite_compute ( int order, double alpha, double x[], double w[] );
double gen_hermite_integral ( int expon, double alpha );
void gen_laguerre_compute ( int order, double alpha, double xtab[], 
  double weight[] );
double gen_laguerre_integral ( int expon, double alpha );
void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] );
void gen_laguerre_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] );
void hermite_ek_compute ( int n, double x[], double w[] );
void hermite_genz_keister_set ( int n, double x[], double w[] );
double hermite_integral ( int n );
double hermite_integral2 ( double a );
void hermite_set ( int order, double xtab[], double weight[] );
void hermite_ss_compute ( int order, double xtab[], double weight[] );
void hermite_ss_recur ( double *p2, double *dp2, double *p1, double x, int order );
void hermite_ss_root ( double *x, int order, double *dp2, double *p1 );
int i4_factorial2 ( int n );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
void imtqlx ( int n, double d[], double e[], double z[] );
void jacobi_compute ( int order, double alpha, double beta, double xtab[], 
  double weight[] );
double jacobi_integral ( int expon, double alpha, double beta );
void jacobi_root ( double *x, int order, double alpha, double beta, 
  double *dp2, double *p1, double b[], double c[] );
void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
  double alpha, double beta, double b[], double c[] );
void kronrod_set ( int order, double xtab[], double weight[] );
void laguerre_compute_dr ( int order, double xtab[], double weight[] );
double laguerre_integral ( int expon );
void laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double b[], double c[] );
void laguerre_root ( double *x, int order, double *dp2, 
  double *p1, double b[], double c[] );
void laguerre_set ( int order, double xtab[], double weight[] );
double laguerre_sum ( double func ( double x ), double a, int order, 
  double xtab[], double weight[] );
void legendre_compute ( int order, double xtab[], double weight[] );
double legendre_integral ( int expon );
void legendre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order );
void legendre_set ( int order, double xtab[], double weight[] );
void legendre_set_cos ( int order, double xtab[], double weight[] );
void legendre_set_cos2 ( int order, double xtab[], double weight[] );
void legendre_set_log ( int order, double xtab[], double weight[] );
void legendre_set_sqrtx_01 ( int order, double xtab[], double weight[] );
void legendre_set_sqrtx2_01 ( int order, double xtab[], double weight[] );
void legendre_set_x0_01 ( int order, double xtab[], double weight[] );
void legendre_set_x1 ( int order, double xtab[], double weight[] );
void legendre_set_x1_01 ( int order, double xtab[], double weight[] );
void legendre_set_x2 ( int order, double xtab[], double weight[] );
void legendre_set_x2_01 ( int order, double xtab[], double weight[] );
void lobatto_compute ( int n, double x[], double w[] );
void lobatto_set ( int order, double xtab[], double weight[] );
double log_gamma ( double x );
void moulton_set ( int order, double xtab[], double weight[] );
void nc_compute ( int order, double a, double b, double xtab[], double weight[] );
void ncc_compute ( int order, double xtab[], double weight[] );
void ncc_compute_points ( int n, double x[] );
void ncc_compute_weights ( int n, double w[] );
void ncc_set ( int order, double xtab[], double weight[] );
void nco_compute ( int n, double x[], double w[] );
void nco_compute_points ( int n, double x[] );
void nco_compute_weights ( int n, double w[] );
void nco_set ( int order, double xtab[], double weight[] );
void ncoh_compute ( int order, double xtab[], double weight[] );
void ncoh_set ( int order, double xtab[], double weight[] );
void patterson_set ( int order, double xtab[], double weight[] );
double r8_abs ( double x );
double r8_epsilon ( void );
double r8_factorial ( int n );
double r8_factorial2 ( int n );
double r8_gamma ( double x );
double r8_huge ( void );
double r8_hyper_2f1 ( double a, double b, double c, double x );
double r8_max ( double x, double y );
double r8_psi ( double xx );
double r8_sign ( double x );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_reverse ( int n, double x[] );
void radau_compute ( int n, double x[], double w[] );
void radau_set ( int order, double xtab[], double weight[] );
void rule_adjust ( double a, double b, double c, double d, int order, 
  double x[], double w[] );
int s_eqi ( char *s1, char *s2 );
double sum_sub ( double func ( double x ), double a, double b, int nsub, 
  int order, double xlo, double xhi, double xtab[], double weight[] );
double summer ( double func ( double x ), int order, double xtab[], 
  double weight[] );
void summer_gk ( double func ( double x ), int orderg, double weightg[],
  double *resultg, int orderk, double xtabk[], double weightk[], 
  double *resultk );
void sum_sub_gk ( double func ( double x ), double a, double b, int nsub, 
  int orderg, double weightg[], double *resultg, int orderk, double xtabk[], 
  double weightk[], double *resultk, double *error );
void timestamp ( void );
