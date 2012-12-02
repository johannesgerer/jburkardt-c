double cpu_time ( void );
double *fibonacci2 ( int n );
void power_method ( int n, double a[], double y[], int it_max, double tol,
  double *lambda, int *it_num );
double r8_abs ( double x );
double r8_epsilon ( void );
void r8mat_mv ( int m, int n, double a[], double x[], double ax[] );
double r8vec_dot ( int n, double a1[], double a2[] );
double r8vec_norm_l2 ( int n, double a[] );
double *r8vec_uniform_01 ( int n, int *seed );
void timestamp ( void );
