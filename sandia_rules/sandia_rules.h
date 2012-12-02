void binary_vector_next ( int n, int bvec[] );

void chebyshev1_compute ( int order, double x[], double w[] );
void chebyshev1_compute_np ( int order, int np, double p[], double x[], double w[] );
void chebyshev1_compute_points ( int order, double x[] );
void chebyshev1_compute_points_np ( int order, int np, double p[], double x[] );
void chebyshev1_compute_weights ( int order, double w[] );
void chebyshev1_compute_weights_np ( int order, int np, double p[], double w[] );
double chebyshev1_integral ( int expon );

void chebyshev2_compute ( int order, double x[], double w[] );
void chebyshev2_compute_np ( int order, int np, double p[], double x[], double w[] );
void chebyshev2_compute_points ( int order, double x[] );
void chebyshev2_compute_points_np ( int order, int np, double p[], double x[] );
void chebyshev2_compute_weights ( int order, double w[] );
void chebyshev2_compute_weights_np ( int order, int np, double p[], double w[] );
double chebyshev2_integral ( int expon );

void clenshaw_curtis_compute ( int order, double x[], double w[] );
void clenshaw_curtis_compute_np ( int order, int np, double p[], double x[], double w[] );
void clenshaw_curtis_compute_points ( int order, double x[] );
void clenshaw_curtis_compute_points_np ( int order, int np, double p[], double x[] );
void clenshaw_curtis_compute_weights ( int order, double w[] );
void clenshaw_curtis_compute_weights_np ( int order, int np, double p[], double w[] );

void comp_next ( int n, int k, int a[], int *more, int *h, int *t );

double cpu_time ( void );

void fejer2_compute ( int order, double x[], double w[] );
void fejer2_compute_np ( int order, int np, double p[], double x[], double w[] );
void fejer2_compute_points ( int order, double x[] );
void fejer2_compute_points_np ( int order, int np, double p[], double x[] );
void fejer2_compute_weights ( int order, double w[] );
void fejer2_compute_weights_np ( int order, int np, double p[], double w[] );

void gegenbauer_compute ( int order, double alpha, double x[], double w[] );
void gegenbauer_compute_np ( int order, int np, double p[], double x[], double w[] );
void gegenbauer_compute_points ( int order, double alpha, double x[] );
void gegenbauer_compute_points_np ( int order, int np, double p[], double x[] );
void gegenbauer_compute_weights ( int order, double alpha, double w[] );
void gegenbauer_compute_weights_np ( int order, int np, double p[], double w[] );
double gegenbauer_integral ( int expon, double alpha );
void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, 
int order, double alpha, double c[] );
void gegenbauer_root ( double *x, int order, double alpha,double *dp2, 
double *p1, double c[] );

void gen_hermite_compute ( int order, double alpha, double x[], double w[] );
void gen_hermite_compute_np ( int order, int np, double p[], double x[], double w[] );
void gen_hermite_compute_points ( int order, double alpha, double x[] );
void gen_hermite_compute_points_np ( int order, int np, double p[], double x[] );
void gen_hermite_compute_weights ( int order, double alpha, double w[] );
void gen_hermite_compute_weights_np ( int order, int np, double p[], double w[] );
double gen_hermite_integral ( int expon, double alpha );

void gen_laguerre_compute ( int order, double alpha, double x[], double w[] );
void gen_laguerre_compute_np ( int order, int np, double p[], double x[], double w[] );
void gen_laguerre_compute_points ( int order, double alpha, double x[] );
void gen_laguerre_compute_points_np ( int order, int np, double p[], double x[] );
void gen_laguerre_compute_weights ( int order, double alpha, double w[] );
void gen_laguerre_compute_weights_np ( int order, int np, double p[], double w[] );
double gen_laguerre_integral ( int expon, double alpha );
void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
int order, double alpha, double b[], double c[] );
void gen_laguerre_root ( double *x, int order, double alpha, double *dp2, 
double *p1, double b[], double c[] );

void hermite_compute ( int order, double x[], double w[] );
void hermite_compute_np ( int order, int np, double p[], double x[], double w[] );
void hermite_compute_points ( int order, double x[] );
void hermite_compute_points_np ( int order, int np, double p[], double x[] );
void hermite_compute_weights ( int order, double w[] );
void hermite_compute_weights_np ( int order, int np, double p[], double w[] );

void hermite_genz_keister_lookup ( int n, double x[], double w[] );
void hermite_genz_keister_lookup_points ( int n, double x[] );
void hermite_genz_keister_lookup_points_np ( int n, int np, double p[], 
  double x[] );
void hermite_genz_keister_lookup_weights ( int n, double w[] );
void hermite_genz_keister_lookup_weights_np ( int n, int np, double p[], 
  double w[] );

double hermite_integral ( int n );
void hermite_recur ( double *p2, double *dp2, double *p1, double x, 
int order );
void hermite_root ( double *x, int order, double *dp2, double *p1 );

int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );

void i4mat_write ( char *output_filename, int m, int n, int table[] );

void i4vec_print ( int n, int a[], char *title );
int i4vec_product ( int n, int a[] );
int i4vec_sum ( int n, int a[] );

void jacobi_compute ( int order, double alpha, double beta, double x[], 
double w[] );
void jacobi_compute_np ( int order, int np, double p[], double x[], double w[] );
void jacobi_compute_points ( int order, double alpha, double beta, 
double x[] );
void jacobi_compute_points_np ( int order, int np, double p[], double x[] );
void jacobi_compute_weights ( int order, double alpha, double beta, 
double w[] );
void jacobi_compute_weights_np ( int order, int np, double p[], double w[] );
double jacobi_integral ( int expon, double alpha, double beta );
void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
double alpha, double beta, double b[], double c[] );
void jacobi_root ( double *x, int order, double alpha, double beta, 
double *dp2, double *p1, double b[], double c[] );

void laguerre_compute ( int order, double x[], double w[] );
void laguerre_compute_np ( int order, int np, double p[], double x[], double w[] );
void laguerre_compute_points ( int order, double x[] );
void laguerre_compute_points_np ( int order, int np, double p[], double x[] );
void laguerre_compute_weights ( int order, double w[] );
void laguerre_compute_weights_np ( int order, int np, double p[], double w[] );
double laguerre_integral ( int expon );
void laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
int order, double b[], double c[] );
void laguerre_root ( double *x, int order, double *dp2, double *p1, 
double b[], double c[] );

void legendre_compute ( int order, double x[], double w[] );
void legendre_compute_np ( int order, int np, double p[], double x[], double w[] );
void legendre_compute_points ( int order, double x[] );
void legendre_compute_points_np ( int order, int np, double p[], double x[] );
void legendre_compute_weights ( int order, double w[] );
void legendre_compute_weights_np ( int order, int np, double p[], double w[] );
double legendre_integral ( int expon );

void level_growth_to_order ( int dim_num, int level[], int rule[], int growth[],
  int order[] );
void level_to_order_default ( int dim_num, int level[], int rule[], 
  int order[] );
void level_to_order_exponential ( int dim_num, int level[], int rule[], 
  int order[] );
void level_to_order_exponential_slow ( int dim_num, int level[], int rule[], 
  int order[] );
void level_to_order_linear ( int dim_num, int level[], int rule[], 
  int order[] );

void nc_compute ( int n, double x_min, double x_max, double x[], double w[] );

void ncc_compute_points ( int n, double x[] );
void ncc_compute_weights ( int n, double w[] );

void nco_compute_points ( int n, double x[] );
void nco_compute_weights ( int n, double w[] );

void patterson_lookup ( int n, double x[], double w[] );
void patterson_lookup_points ( int order, double x[] );
void patterson_lookup_points_np ( int order, int np, double p[], double x[] );
void patterson_lookup_weights ( int order, double w[] );
void patterson_lookup_weights_np ( int order, int np, double p[], double w[] );

int point_radial_tol_unique_count ( int m, int n, double a[], double tol, 
  int *seed );
int point_radial_tol_unique_index ( int m, int n, double a[], double tol, 
  int *seed, int undx[], int xdnu[] );

void point_unique_index ( int m, int n, double a[], int unique_num, int undx[], 
  int xdnu[] );

void product_mixed_weight ( int dim_num, int order_1d[], int order_nd, 
  int rule[], double alpha[], double beta[], double weight_nd[] );

double r8_abs ( double x );
double r8_ceiling ( double x );
double r8_choose ( int n, int k );
double r8_epsilon ( void );
double r8_factorial ( int n );
double r8_factorial2 ( int n );
double r8_floor ( double x );
double r8_gamma ( double x );
double r8_huge ( void );
double r8_hyper_2f1 ( double a, double b, double c, double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_mop ( int i );
double r8_psi ( double xx );

int r8col_compare ( int m, int n, double a[], int i, int j );
void r8col_sort_heap_a ( int m, int n, double a[] );
int *r8col_sort_heap_index_a ( int m, int n, int base, double a[] );
int r8col_sorted_unique_count ( int m, int n, double a[], double tol );
void r8col_swap ( int m, int n, double a[], int j1, int j2 );
void r8col_tol_undex ( int x_dim, int x_num, double x_val[], int x_unique_num, 
  double tol, int undx[], int xdnu[] );
int r8col_tol_unique_count ( int m, int n, double a[], double tol );
void r8col_undex ( int x_dim, int x_num, double x_val[], int x_unique_num, 
  double tol, int undx[], int xdnu[] );
int *r8col_unique_index ( int m, int n, double a[], double tol );

void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void r8mat_write ( char *output_filename, int m, int n, double table[] );

double *r8vec_chebyshev_new ( int n, double a_first, double a_last );
int r8vec_compare ( int n, double a[], double b[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void r8vec_index_sorted_range ( int n, double r[], int indx[], double r_lo, 
  double r_hi, int *i_lo, int *i_hi );
double r8vec_min ( int n, double r8vec[] );
double r8vec_min_pos ( int n, double a[] );
void r8vec_print ( int n, double a[], char *title );
int *r8vec_sort_heap_index_a ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );

void timestamp ( void );

void vec_colex_next3 ( int dim_num, int base[], int a[], int *more );
