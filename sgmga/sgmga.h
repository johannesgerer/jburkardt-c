void sgmga_aniso_normalize ( int option, int dim_num, double level_weight[] );

void sgmga_importance_to_aniso ( int dim_num, double importance[], 
  double level_weight[] );

void sgmga_index ( int dim_num, double level_weight[], int level_max, 
  int rule[], int point_num, int point_total_num, int sparse_unique_index[],
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  int sparse_order[], int sparse_index[] );

void sgmga_point ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  int point_num, int sparse_order[], int sparse_index[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  double sparse_point[] );

void sgmga_product_weight ( int dim_num, int order_1d[], int order_nd, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  double weight_nd[] );

int sgmga_size ( int dim_num, double level_weight[], int level_max, int rule[], 
  int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol,
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ) );

int sgmga_size_total ( int dim_num, double level_weight[], int level_max, 
  int rule[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ) );

void sgmga_unique_index ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol, int point_num, int point_total_num, 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ), 
  int sparse_unique_index[] );

void sgmga_vcn ( int n, double level_weight[], int x_max[], int x[], 
  double q_min, double q_max, int *more );

double sgmga_vcn_coef ( int n, double level_weight[], int x_max[], int x[], 
  double q_min, double q_max );

void sgmga_vcn_ordered ( int dim_num, double level_weight[], int x_max[], int x[], 
  double q_min, double q_max, int *more );

void sgmga_weight ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  int point_num, int point_total_num, int sparse_unique_index[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  double sparse_weight[] );

void sgmga_write ( int dim_num, double level_weight[], int rule[], int np[],
  double p[], int point_num, double sparse_weight[], double sparse_point[],
  char *file_name );
