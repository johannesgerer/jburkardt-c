int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void r83vec_print_part ( int n, double a[], int max_print, char *title );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
double *tetrahedron_grid ( int n, double t[], int ng );
int tetrahedron_grid_count ( int n );
void timestamp ( void );
