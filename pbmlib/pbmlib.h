void bitchr75 ( char c, int *pattern );

int  pbm_check_data ( int xsize, int ysize, int *barray );
int  pbm_example ( int xsize, int ysize, int *barray );

int  pbma_read ( char *filein_name, int *xsize, int *ysize, int **barray );
int  pbma_read_data ( FILE *filein, int xsize, int ysize, int *barray );
int  pbma_read_header ( FILE *filein, int *xsize, int *ysize );
int  pbma_read_test ( char *filein_name );

int  pbma_write ( char *fileout_name, int xsize, int ysize, int *barray );
int  pbma_write_data ( FILE *fileout, int xsize, int ysize, int *barray );
int  pbma_write_header ( FILE *fileout, char *fileout_name, int xsize, 
       int ysize );
int  pbma_write_test ( char *fileout_name );

int  pbmb_read ( char *filein_name, int *xsize, int *ysize, int **barray );
int  pbmb_read_data ( FILE *filein, int xsize, int ysize, int *barray );
int  pbmb_read_header ( FILE *filein, int *xsize, int *ysize );
int  pbmb_read_test ( char *filein_name );

int  pbmb_write ( char *fileout_name, int xsize, int ysize, int *barray );
int  pbmb_write_data ( FILE *fileout, int xsize, int ysize, int *barray );
int  pbmb_write_header ( FILE *fileout, int xsize, int ysize );
int  pbmb_write_test ( char *fileout_name );

int  pgm_check_data ( int xsize, int ysize, int maxgray, int *garray );
int  pgm_example ( int xsize, int ysize, int *garray );

int  pgma_read ( char *filein_name, int *xsize, int *ysize, int *maxgray,
       int **garray );
int  pgma_read_data ( FILE *filein, int xsize, int ysize, int *garray );
int  pgma_read_header ( FILE *filein, int *xsize, int *ysize, int *maxgray );
int  pgma_read_test ( char *filein_name );

int  pgma_write ( char *fileout_name, int xsize, int ysize, int *barray );
int  pgma_write_data ( FILE *fileout, int xsize, int ysize, int *barray );
int  pgma_write_header ( FILE *fileout, char *fileout_name, int xsize, 
       int ysize, int maxgray );
int  pgma_write_test ( char *fileout_name );

int  pgmb_read ( char *filein_name, int *xsize, int *ysize, int *maxgray,
      int **garray );
int  pgmb_read_data ( FILE *filein, int xsize, int ysize, int *garray );
int  pgmb_read_header ( FILE *filein, int *xsize, int *ysize, int *maxgray );
int  pgmb_read_test ( char *filein_name );

int  pgmb_write ( char *fileout_name, int xsize, int ysize, int *barray );
int  pgmb_write_data ( FILE *fileout, int xsize, int ysize, int *barray );
int  pgmb_write_header ( FILE *fileout, int xsize, int ysize, int maxgray );
int  pgmb_write_test ( char *fileout_name );

int  ppm_check_data ( int xsize, int ysize, int maxrgb, int *rarray,
       int *garray, int *barray );
int  ppm_example ( int xsize, int ysize, int *rarray, int *garray, int *barray );

int  ppma_read ( char *filein_name, int *xsize, int *ysize, int *maxrgb,
       int **rarrary, int **garray, int **barray );
int  ppma_read_data ( FILE *filein, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
int  ppma_read_header ( FILE *filein, int *xsize, int *ysize, int *maxrgb );
int  ppma_read_test ( char *filein_name );

int  ppma_write ( char *fileout_name, int xsize, int ysize, int *rarray, 
      int *garray, int *barray );
int  ppma_write_data ( FILE *fileout, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
int  ppma_write_header ( FILE *fileout, char *fileout_name, int xsize, 
       int ysize, int maxrgb );
int  ppma_write_test ( char *fileout_name );

int  ppmb_read ( char *filein_name, int *xsize, int *ysize, int *maxrgb,
       int **rarray, int **garray, int **barray );
int  ppmb_read_data ( FILE *filein, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
int  ppmb_read_header ( FILE *filein, int *xsize, int *ysize, int *maxrgb );
int  ppmb_read_test ( char *filein_name );

int  ppmb_write ( char *fileout_name, int xsize, int ysize, int *rarray,
      int *garray, int *barray );
int  ppmb_write_data ( FILE *fileout, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
int  ppmb_write_header ( FILE *fileout, int xsize, int ysize, int maxrgb );
int  ppmb_write_test ( char *fileout_name );
void timestamp ( void );
