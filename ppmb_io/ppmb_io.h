bool ppmb_check_data ( int xsize, int ysize, int maxrgb, unsigned char *rarray,
  unsigned char *garray, unsigned char *barray );
bool ppmb_example ( int xsize, int ysize, unsigned char *rarray, 
  unsigned char *garray, unsigned char *barray );
bool ppmb_read ( char *file_name, int *xsize, int *ysize, int *maxrgb,
  unsigned char **rarray, unsigned char **garray, unsigned char **barray );
bool ppmb_read_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *rarray, unsigned char *garray, unsigned char *barray );
bool ppmb_read_header ( FILE *file_pointer, int *xsize, int *ysize, int *maxrgb );
bool ppmb_read_test ( char *file_name );
bool ppmb_write ( char *file_name, int xsize, int ysize, unsigned char *rarray, 
  unsigned char *garray, unsigned char *barray );
bool ppmb_write_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *rarray, unsigned char *garray, unsigned char *barray );
bool ppmb_write_header ( FILE *file_pointer, int xsize, int ysize, int maxrgb );
bool ppmb_write_test ( char *file_name );
