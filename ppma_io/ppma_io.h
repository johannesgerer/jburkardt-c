char ch_cap ( char ch );
int i4_max ( int i1, int i2 );

int ppma_check_data ( int xsize, int ysize, int rgb_max, int *r, int *g, 
  int *b );

void ppma_example ( int xsize, int ysize, int *r, int *g, int *b );

void ppma_read ( char *input_name, int *xsize, int *ysize, int *rgb_max,
  int **r, int **g, int **b );
void ppma_read_data ( FILE *input, int xsize, int ysize, int *r,
  int *g, int *b );
void ppma_read_header ( FILE *input, int *xsize, int *ysize, int *rgb_max );
void ppma_read_test ( char *input_name );

int ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
  int *g, int *b );
int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r,
  int *g, int *b );
int ppma_write_header ( FILE *file_out, char *file_out_name, int xsize,
  int ysize, int rgb_max );
int ppma_write_test ( char *file_out_name );

void timestamp ( void );
