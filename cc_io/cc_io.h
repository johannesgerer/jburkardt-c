void cc_data_read ( char *prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
void cc_header_read ( char *prefix, int *ncc, int *n );
double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  char *title );
void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], char *title );
void cc_write ( char *prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
int file_row_count ( char *input_filename );
void i4vec_data_read ( char *input_filename, int n, int a[] );
void i4vec_dec ( int n, int a[] );
void i4vec_inc ( int n, int a[] );
void i4vec_write ( char *output_filename, int n, int table[] );
void r8vec_data_read ( char *input_filename, int n, double x[] );
void r8vec_write ( char *output_filename, int n, double x[] );
int s_len_trim ( char *s );
void timestamp ( );
