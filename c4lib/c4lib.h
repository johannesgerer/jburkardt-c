float          c4_abs ( float complex c );
float complex  c4_acos ( float complex c1 );
float complex  c4_acosh ( float complex c1 );
float complex  c4_add ( float complex c1, float complex c2 );
float          c4_arg ( float complex c );
float complex  c4_asin ( float complex c1 );
float complex  c4_asinh ( float complex c1 );
float complex  c4_atan ( float complex c1 );
float complex  c4_atanh ( float complex c1 );
float complex  c4_conj ( float complex c1 );
void           c4_copy ( float complex c1, float complex c2 );
float complex  c4_cos ( float complex c1 );
float complex  c4_cosh ( float complex c1 );
float complex  c4_cube_root ( float complex c1 );
float complex  c4_div ( float complex c1, float complex c2 );
float complex  c4_div_r4 ( float complex c1, float r );
float complex  c4_exp ( float complex c1 );
float complex  c4_i ( );
float          c4_imag ( float complex c );
float complex  c4_inv ( float complex c1 );
int            c4_le_l1 ( float complex x, float complex y );
int            c4_le_l2 ( float complex x, float complex y );
int            c4_le_li ( float complex x, float complex y );
float complex  c4_log ( float complex c1 );
float          c4_mag ( float complex c );
float complex  c4_mul ( float complex c1, float complex c2 );
float complex  c4_neg ( float complex c1 );
float complex  c4_nint ( float complex c1 );
float          c8_norm_l1 ( float complex x );
float          c8_norm_l2 ( float complex x );
float          c8_norm_li ( float complex x );
float complex  c4_normal_01 ( int *seed );
float complex  c4_one ( );
void           c4_print ( float complex a, char *title );
float          c4_real ( float complex c );
float complex  c4_sin ( float complex c1 );
float complex  c4_sinh ( float complex c1 );
float complex  c4_sqrt ( float complex c1 );
float complex  c4_sub ( float complex c1, float complex c2 );
void           c4_swap ( float complex c1, float complex c2 );
float complex  c4_tan ( float complex c1 );
float complex  c4_tanh ( float complex c1 );
void           c4_to_cartesian ( float complex c, float *x, float *y );
void           c4_to_polar ( float complex c, float *r, float *theta );
float complex  c4_uniform_01 ( int *seed );
float complex  c4_zero ( );
void           c4mat_add ( int m, int n, float complex alpha, float complex a[],
               float complex beta, float complex b[], float complex c[] );
void           c4mat_add_r4 ( int m, int n, float alpha, float complex a[],
               float beta, float complex b[], float complex c[] );
void           c4mat_copy ( int m, int n, float complex a[], float complex b[] );
float complex *c4mat_copy_new ( int m, int n, float complex a1[] );
void           c4mat_fss ( int n, float complex a[], int nb, float complex x[] );
float complex *c4mat_fss_new ( int n, float complex a[], int nb, 
               float complex b[] );
float complex *c4mat_identity_new ( int n );
float complex *c4mat_indicator_new ( int m, int n );
void           c4mat_minvm ( int n1, int n2, float complex a[], 
               float complex b[], float complex c[] );
float complex *c4mat_minvm_new ( int n1, int n2, float complex a[], 
               float complex b[] );
void           c4mat_mm ( int n1, int n2, int n3, float complex a[], float complex b[], 
               float complex c[] );
float complex *c4mat_mm_new ( int n1, int n2, int n3, float complex a[], 
               float complex b[] );
void           c4mat_nint ( int m, int n, float complex a[] );
float          c4mat_norm_fro ( int m, int n, float complex a[] );
float          c4mat_norm_l1 ( int m, int n, float complex a[] );
float          c4mat_norm_li ( int m, int n, float complex a[] );
void           c4mat_print ( int m, int n, float complex a[], char *title );
void           c4mat_print_some ( int m, int n, float complex a[], int ilo, int jlo, 
               int ihi, int jhi, char *title );
void           c4mat_scale ( int m, int n, float complex alpha, float complex a[] );
void           c4mat_scale_r4 ( int m, int n, float alpha, float complex a[] );
float complex *c4mat_uniform_01 ( int m, int n, int *seed );
float complex *c4mat_zero_new ( int m, int n );
void           c4vec_copy ( int n, float complex a[], float complex b[] );
float complex *c4vec_copy_new ( int n, float complex a1[] );
void           c4vec_nint ( int n, float complex a[] );
float          c4vec_norm_l2 ( int n, float complex a[] );
void           c4vec_print ( int n, float complex a[], char *title );
void           c4vec_print_part ( int n, float complex a[], int max_print, char *title );
void           c4vec_print_some ( int n, float complex a[], int i_lo, int i_hi, 
               char *title );
void           c4vec_sort_a_l2 ( int n, float complex x[] );
float complex *c4vec_spiral ( int n, int m, float complex c1, 
               float complex c2 );
float complex *c4vec_uniform_01_new ( int n, int *seed );
float complex *c4vec_unity ( int n );
float complex  cartesian_to_c4 ( float x, float y );
float complex  polar_to_c4 ( float r, float theta );
float complex  r4_csqrt ( float x );
void           r4poly2_root ( float a, float b, float c, float complex *r1,
               float complex *r2 );
void           r4poly3_root ( float a, float b, float c, float d,
               float complex *r1, float complex *r2, float complex *r3 );
void           r4poly4_root ( float a, float b, float c, float d, float e,
               float complex *r1, float complex *r2, float complex *r3,
               float complex *r4 );
void           sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
