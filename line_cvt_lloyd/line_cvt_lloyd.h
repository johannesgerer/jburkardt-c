void energy_plot ( int it_num, double e_plot[], char *header );
void evolution_plot ( int n, int it_num, double x_plot[], char *header );
void line_ccvt_lloyd ( int n, double a, double b, int it_num, char *header, 
  double x[] );
void line_ccvt_lloyd_step ( int n, double a, double b, double x[] );
double line_cvt_energy ( int n, double a, double b, double x[] );
void line_cvt_lloyd ( int n, double a, double b, int it_num, char *header, 
  double x[] );
void line_cvt_lloyd_step ( int n, double a, double b, double x[] );
void motion_plot ( int it_num, double xm_plot[], char *header );
void r8vec_copy ( int n, double a1[], double a2[] );
double *r8vec_copy_new ( int n, double a1[] );
void r8vec_print ( int n, double a[], char *title );
void r8vec_sort_insert_a ( int n, double a[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int *seed );
void timestamp ( );

