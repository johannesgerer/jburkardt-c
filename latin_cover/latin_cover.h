int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_uniform ( int a, int b, int *seed );
int i4_wrap ( int ival, int ilo, int ihi );
void i4block_print ( int l, int m, int n, int a[], char *title );
void i4mat_print ( int m, int n, int a[], char *title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void i4vec_print ( int n, int a[], char *title );
int *latin_cover ( int n, int p[] );
int *latin_cover_2d ( int n, int p1[], int p2[] );
int *latin_cover_3d ( int n, int p1[], int p2[], int p3[] );
int perm_check ( int n, int p[] );
void perm_print ( int n, int p[], char *title );
int *perm_uniform_new ( int n, int *seed );
int s_len_trim ( char *s );
void timestamp ( void );
