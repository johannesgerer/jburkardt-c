int between ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
int collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
int diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
int diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
int in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
int intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
int intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
int l4_xor ( int l1, int l2 );
double *monomial_value ( int m, int n, int e[], double x[] );
double polygon_area ( int nv, double v[] );
double polygon_monomial_integral ( int nv, double v[], int e[] );
double *polygon_sample ( int nv, double v[], int n, int *seed );
int *polygon_triangulate ( int n, double x[], double y[] );
double r8_choose ( int n, int k );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int *seed );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
void timestamp ( );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc );

