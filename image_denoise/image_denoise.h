int *gray_median_news ( int m, int n, int gray[] );
int i4vec_frac ( int n, int a[], int k );
int i4vec_median ( int n, int a[] );
void pgma_read_data ( FILE *file_in, int xsize, int ysize, int *g );
void pgma_read_header ( FILE *file_in, int *xsize, int *ysize, int *maxg );
void pgma_write ( char *file_out_name, int xsize, int ysize, int *g );
void pgma_write_data ( FILE *file_out, int xsize, int ysize, int *g );
void pgma_write_header ( FILE *file_out, char *file_out_name, int xsize, 
  int ysize, int maxg );
void timestamp ( void );
