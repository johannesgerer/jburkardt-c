double root_rc ( double x, double fx, double *ferr, double *xerr, double q[9] );
void roots_rc ( int n, double x[], double fx[], double *ferr, double xnew[], 
  double q[] );
double r8_abs ( double x );
double r8_epsilon ( void );
double r8_huge ( void );
double r8_sign ( double x );
void r8mat_fs ( int lda, int n, double a[], double x[] );
void timestamp ( void );
