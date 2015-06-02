double cholesky_upper_error ( int n, double a[], double c[] );
double eigen_error ( int n, int k, double a[], double x[], double lambda[] );
double inverse_error ( int n, double a[], double b[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );

double *l1dd_apply ( int n, double h, double u[] );
double *l1dd_cholesky ( int n, double h );
void l1dd_eigen ( int n, double h, double v[], double lambda[] );
double *l1dd ( int n, double h );
double *l1dd_inverse ( int n, double h );
void l1dd_lu ( int n, double h, double l[], double u[] );

double *l1dn_apply ( int n, double h, double u[] );
double *l1dn_cholesky ( int n, double h );
void l1dn_eigen ( int n, double h, double v[], double lambda[] );
double *l1dn ( int n, double h );
double *l1dn_inverse ( int n, double h );
void l1dn_lu ( int n, double h, double l[], double u[] );

double *l1nd_apply ( int n, double h, double u[] );
double *l1nd_cholesky ( int n, double h );
void l1nd_eigen ( int n, double h, double v[], double lambda[] );
double *l1nd ( int n, double h );
double *l1nd_inverse ( int n, double h );
void l1nd_lu ( int n, double h, double l[], double u[] );

double *l1nn_apply ( int n, double h, double u[] );
double *l1nn_cholesky ( int n, double h );
void l1nn_eigen ( int n, double h, double v[], double lambda[] );
double *l1nn ( int n, double h );
void l1nn_lu ( int n, double h, double l[], double u[] );

double *l1pp_apply ( int n, double h, double u[] );
double *l1pp_cholesky ( int n, double h );
void l1pp_eigen ( int n, double h, double v[], double lambda[] );
double *l1pp ( int n, double h );
void l1pp_lu ( int n, double h, double l[], double u[] );

double lu_error ( int n, double a[], double l[], double u[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] );
double r8mat_norm_fro ( int m, int n, double a[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
double *r8mat_sub_new ( int m, int n, double a[], double b[] );
void r8vec_print ( int n, double a[], char *title );
void timestamp ( void );

