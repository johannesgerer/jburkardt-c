/*
  Utility to get a seed for the random number generator.
*/
int get_seed ( );
/*
  Utility for one of the tests.
*/
double *hilbert_inverse ( int n );
/*
  Utilities for integers.
*/
int i4_huge ( void );
int i4_log_10 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_uniform ( int a, int b, int *seed );
/*
  Utilities for integer vectors.
*/
void i4vec_print ( int n, int a[], char *title );
int i4vec_search_binary_a ( int n, int a[], int b );
/*
  R4 Utilities.
*/
float r4_abs ( float x );
int r4_nint ( float x );
float r4_uniform ( float b, float c, int *seed );
float r4_uniform_01 ( int *seed );
/*
  R8 Utilities.
*/
double r8_abs ( double x );
int r8_is_int ( double r );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
double r8_sign2 ( double x, double y );
void r8_swap ( double *x, double *y );
double r8_uniform ( double rlo, double rhi, int *seed );
double r8_uniform_01 ( int *seed );
/*
  Real double precision Triadiagonal.
*/
double *r83_cr_fa ( int n, double a[] );
double *r83_cr_sl ( int n, double a_cr[], double b[] );
double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] );
void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, int job );
double *r83_indicator ( int n );
void r83_jac_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job );
double *r83_mxv ( int n, double a[], double x[] );
double r83_np_det ( int n, double a_lu[] );
int r83_np_fa ( int n, double a[] );
double *r83_np_fs ( int n, double a[], double b[] );
double *r83_np_ml ( int n, double a_lu[], double x[], int job );
double *r83_np_sl ( int n, double a_lu[], double b[], int job );
void r83_print ( int n, double a[], char *title );
void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  char *title );
double *r83_random ( int n, int *seed );
double *r83_to_r8ge ( int n, double a[] );
double *r83_vxm ( int n, double a[], double x[] );
double *r83_zero ( int n );
/*
  Real double precision Tridiagonal No Pivoting.
*/
double *r83np_fs ( int n, double a[], double b[] );
/*
  Real double precision Tridiagonal Periodic.
*/
double r83p_det ( int n, double a[], double work4 );
int r83p_fa ( int n, double a[], double work2[], double work3[], double *work4 );
double *r83p_indicator ( int n );
double *r83p_ml ( int n, double a[], double x[], int job );
double *r83p_mxv ( int n, double a[], double x[] );
void r83p_print ( int n, double a[], char *title );
void r83p_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
double *r83p_random ( int n, int *seed );
double *r83p_sl ( int n, double a[], double b[], int job, double work2[], 
  double work3[], double work4 );
double *r83p_to_r8ge ( int n, double a[] );
double *r83p_vxm ( int n, double a[], double x[] );
double *r83p_zero ( int n );
/*
  Real double precision Pentagonal.
*/
double *r85_indicator ( int n );
double *r85_mxv ( int n, double a[], double x[] );
double *r85_np_fs ( int n, double a[], double b[] );
void r85_print ( int n, double a[], char *title );
void r85_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
double *r85_random ( int n, int *seed );
double *r85_to_r8ge ( int n, double a[] );
double *r85_vxm ( int n, double a[], double x[] );
double *r85_zero ( int n );
/*
  Real double precision general.
*/
double r8ge_co ( int n, double a[], int pivot[] );
double r8ge_det ( int n, double a[], int pivot[] );
double *r8ge_dilu ( int m, int n, double a[] );
int r8ge_fa ( int n, double a[], int pivot[] );
void r8ge_fs ( int n, double a[], double x[] );
double *r8ge_fs_new ( int n, double a[], double b[] );
void r8ge_fss ( int n, double a[], int nb, double x[] );
double *r8ge_fss_new ( int n, double a[], int nb, double b[] );
double *r8ge_identity ( int n );
void r8ge_ilu ( int m, int n, double a[], double l[], double u[] );
double *r8ge_indicator ( int m, int n );
double *r8ge_inverse ( int n, double a[], int pivot[] );
double *r8ge_ml ( int n, double a[], int pivot[], double x[], int job );
double *r8ge_mu ( int m, int n, double a[], char trans, int pivot[], double x[] );
double *r8ge_mxm ( int n, double a[], double b[] );
double *r8ge_mxv ( int m, int n, double a[], double x[] );
double r8ge_np_det ( int n, double a[] );
int r8ge_np_fa ( int n, double a[] );
double *r8ge_np_inverse ( int n, double a[] );
double *r8ge_np_ml ( int n, double a[], double x[], int job );
double *r8ge_np_sl ( int n, double a[], double b[], int job );
double *r8ge_np_trm ( int m, int n, double a[], double x[], int job );
int r8ge_np_trf ( int m, int n, double a[] );
double *r8ge_np_trs ( int n, int nrhs, char trans, double a[], double b[] );
void r8ge_plu ( int m, int n, double a[], double p[], double l[], double u[] );
double *r8ge_poly ( int n, double a[] );
void r8ge_print ( int m, int n, double a[], char *title );
void r8ge_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
double *r8ge_random ( int m, int n, int *seed );
double *r8ge_res ( int m, int n, double a[], double x[], double b[] );
double *r8ge_sl ( int n, double a[], int pivot[], double b[], int job );
double *r8ge_sl_it ( int n, double a[], double alu[], int pivot[], double b[], 
  int job, double x[] );
double *r8ge_to_r8gb ( int m, int n, int ml, int mu, double a[] );
double *r8ge_to_r8vec ( int m, int n, double *a );
int r8ge_trf ( int m, int n, double a[], int pivot[] );
double *r8ge_trs ( int n, int nrhs, char trans, double a[], int pivot[], double b[] );
double *r8ge_vxm ( int m, int n, double a[], double x[] );
double *r8ge_zero ( int m, int n );
/*
  Real double precision Lower Triangular, Full Storage.
*/
double r8lt_det ( int n, double a[] );
double *r8lt_indicator ( int m, int n );
double *r8lt_inverse ( int n, double a[] );
double *r8lt_mxm ( int n, double a[], double b[] );
double *r8lt_mxv ( int m, int n, double a[], double x[] );
void r8lt_print ( int m, int n, double a[], char *title );
void r8lt_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
double *r8lt_random ( int m, int n, int *seed );
double *r8lt_sl ( int n, double a[], double b[], int job );
double *r8lt_vxm ( int m, int n, double a[], double x[] );
double *r8lt_zero ( int m, int n );
/*
  R8MAT Utilities.
*/
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
/*
  R8VEC Utilities.
*/
double *r8vec_indicator_new ( int n );
void r8vec_print ( int n, double a[], char *title );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title );
double *r8vec_uniform ( int n, double b, double c, int *seed );

int s_len_trim ( char *s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( void );
