char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void ffmsh_2d_data_example ( int v_num, int e_num, int t_num, double v_xy[], 
  int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_data_print ( char *title, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_data_read ( char *ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_size_example ( int *v_num, int *e_num, int *t_num );
void ffmsh_2d_size_print ( char *title, int v_num, int e_num, int t_num );
void ffmsh_2d_size_read ( char *ffmsh_filename, int *v_num, int *e_num, 
  int *t_num );
void ffmsh_2d_write ( char *ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void i4mat_copy ( int m, int n, int a1[], int a2[] );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void i4vec_copy ( int n, int a1[], int a2[] );
void i4vec_print ( int n, int a[], char *title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
double s_to_r8 ( char *s, int *lchar, int *error );
void timestamp ( );
