double *cc_compute_points ( int n );
int i4vec_product ( int n, int a[] );
double *lagrange_basis_1d ( int nd, double xd[], int ni, double xi[] );
double *lagrange_interp_nd_grid ( int m, int n_1d[], double a[], double b[], int nd );
double *lagrange_interp_nd_grid2 ( int m, int ind[], double a[], double b[], int nd );
int lagrange_interp_nd_size ( int m, int ind[] );
int lagrange_interp_nd_size2 ( int m, int ind[] );
double *lagrange_interp_nd_value ( int m, int n_1d[], double a[], double b[], int nd, 
  double zd[], int ni, double xi[] );
double *lagrange_interp_nd_value2 ( int m, int ind[], double a[], double b[], int nd, 
  double zd[], int ni, double xi[] );
int order_from_level_135 ( int l );
