int polygon_grid_count ( int n, int nv );
void polygon_grid_display ( int n, int nv, double v[], int ng, double xg[], 
  char *prefix );
double *polygon_grid_points ( int n, int nv, double v[], int ng );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void timestamp ( );

