void alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *alpha_min, double *alpha_ave,
  double *alpha_area );
double angle_rad_2d ( double p1[2], double p2[2], double p3[2] );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_sign ( int i );
int i4_wrap ( int ival, int ilo, int ihi );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void i4vec_heap_d ( int n, int a[] );
int *i4vec_indicator_new ( int n );
int i4vec_min ( int n, int a[] );
void i4vec_sort_heap_a ( int n, int a[] );
int i4vec_sorted_unique ( int n, int a[] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2, double yv2,
  double dv );
int perm_check2 ( int n, int p[], int base );
void perm_inverse ( int n, int p[] );
int *points_delaunay_naive_2d ( int n, double p[], int *ntri );
void points_hull_2d ( int node_num, double node_xy[], int *hull_num,
  int hull[] );
void quad_convex_random ( int *seed, double xy[] );
double r8_abs ( double x );
double r8_acos ( double c );
double r8_epsilon ( void );
double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
void r82vec_part_quick_a ( int n, double a[], int *l, int *r );
void r82vec_permute ( int n, int p[], int base, double a[] );
int *r82vec_sort_heap_index_a ( int n, int base, double a[] );
void r82vec_sort_quick_a ( int n, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
int r8tris2 ( int point_num, double point_xy[], int *tri_num,
  int tri_vert[], int tri_nabe[] );
int r8vec_eq ( int n, double a1[], double a2[] );
int r8vec_gt ( int n, double a1[], double a2[] );
int r8vec_lt ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], char *title );
void r8vec_swap ( int n, double a1[], double a2[] );
int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
  int stack[] );
void timestamp ( void );
double *triangle_circumcenter_2d ( double t[] );
void triangulation_order3_plot ( char *file_out_name, int g_num, double g_xy[],
  int tri_num, int nod_tri[], int node_show, int triangle_show );
void triangulation_order3_print ( int point_num, int tri_num, double xc[], 
  int tri_vert[], int tri_nabe[] );
void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );

