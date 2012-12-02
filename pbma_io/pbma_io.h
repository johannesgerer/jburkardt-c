void pbma_check_data ( int xsize, int ysize, int *barray );
void pbma_example ( int xsize, int ysize, int *barray );

void pbma_read ( char *file_in_name, int *xsize, int *ysize, int **barrary );
void pbma_read_data ( FILE *file_in, int xsize, int ysize, int *barray );
void pbma_read_header ( FILE *file_in, int *xsize, int *ysize );
void pbma_read_test ( char *file_in_name );

void pbma_write ( char *file_out_name, int xsize, int ysize, int *barray );
void pbma_write_data ( FILE *file_out, int xsize, int ysize, int *barray );
void pbma_write_header ( FILE *file_out, char *file_out_name, int xsize, 
       int ysize );
void pbma_write_test ( char *file_out_name );
