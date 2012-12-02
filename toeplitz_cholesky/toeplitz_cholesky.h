int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *r8mat_copy_new ( int m, int n, double a1[] );
double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void timestamp ( void );
double *toep_cholesky_lower ( int n, double g[] );
double *toep_cholesky_upper ( int n, double g[] );
double *toeplitz_cholesky_lower ( int n, double a[] );
double *toeplitz_cholesky_upper ( int n, double a[] );

