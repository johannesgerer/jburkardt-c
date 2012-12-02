char ch_cap ( char c );
int ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int file_column_count ( char *input_filename );
char *file_name_ext_swap ( char *filename, char *ext );
int file_row_count ( char *input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
void r8mat_latinize ( int m, int n, double table[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
int *r8vec_sort_heap_index_a_new ( int n, double a[] );
int s_index_last_c ( char *s, char c );
int s_len_trim ( char *s );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( void );

