int i4_huge ( void );
int *i4mat_histogram ( int m, int n, int a[], int histo_num );
int i4mat_max ( int m, int n, int a[] );
int *news ( int m, int n, int a[] );
void pbma_write ( char *file_out_name, int xsize, int ysize, int *b );
void pbma_write_data ( FILE *file_out, int xsize, int ysize, int *b );
void pbma_write_header ( FILE *file_out, char *file_out_name, int xsize, 
  int ysize );
void pgma_read_data ( FILE *file_in, int xsize, int ysize, int *g );
void pgma_read_header ( FILE *file_in, int *xsize, int *ysize, int *maxg );
void timestamp ( void );
