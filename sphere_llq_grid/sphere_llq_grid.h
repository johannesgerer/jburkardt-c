void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8vec_print ( int n, double a[], char *title );
void sphere_llq_grid_display ( int ng, double xg[], int line_num, 
  int line_data[], char *prefix );
int sphere_llq_grid_line_count ( int lat_num, int long_num );
int *sphere_llq_grid_lines ( int nlat, int nlong, int line_num );
int sphere_llq_grid_point_count ( int lat_num, int long_num );
double *sphere_llq_grid_points ( double r, double pc[3], int lat_num, int lon_num,
  int point_num );
void timestamp ( );
