int between ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
int diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
int diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
void i4mat_print ( int m, int n, int a[], char *title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void i4mat_write ( char *output_filename, int m, int n, int table[] );
int in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
int intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
int intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
int l4_xor ( int l1, int l2 );
void l4vec_print ( int n, int a[], char *title );
int *polygon_triangulate ( int n, double x[], double y[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
int s_len_trim ( char *s );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc );

