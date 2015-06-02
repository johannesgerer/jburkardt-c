char ch_cap ( char c );
int ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
void fem_data_read ( char *node_coord_file_name, char *element_file_name, 
  char *node_data_file_name, int dim_num, int node_num, int element_num,
  int element_order, int node_data_num, double **node_coord, int **element_node,
  double **node_data );
void fem_header_print ( int dim_num, int node_num, int element_num, 
  int element_order, int node_data_num );
void fem_header_read ( char *node_coord_file_name, char *element_file_name, 
  char *node_data_file_name, int *dim_num, int *node_num, int *element_num,
  int *element_order, int *node_data_num );
void fem_write ( char *node_coord_file_name, char *element_file_name, 
  char *node_data_file_name, int dim_num, int node_num, int element_num, 
  int element_order, int node_data_num, double node_coord[], 
  int element_node[], double node_data[] );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( char *input_filename, int m, int n );
void i4mat_header_read ( char *input_filename, int *m, int *n );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void i4mat_write ( char *output_filename, int m, int n, int table[] );
double r8_epsilon ( );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
int s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( );
