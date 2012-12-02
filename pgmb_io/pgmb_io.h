char ch_cap ( char ch );
bool pgmb_check_data ( int xsize, int ysize, unsigned char maxgray, 
  unsigned char *garray );
bool pgmb_example ( int xsize, int ysize, unsigned char *garray );
bool pgmb_read ( char *file_name, int *xsize, int *ysize, unsigned char *maxgray,
  unsigned char **garray );
bool pgmb_read_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *garray );
bool pgmb_read_header ( FILE *file_pointer, int *xsize, int *ysize, 
  unsigned char *maxgray );
bool pgmb_read_test ( char *file_name );
bool pgmb_write ( char *file_name, int xsize, int ysize, unsigned char *garray );
bool pgmb_write_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *garray );
bool pgmb_write_header ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char maxgray );
bool pgmb_write_test ( char *file_name );
void s_adjustl ( char *s );
int s_eqi ( char *s1, char *s2 );
int s_len_trim ( char *s );
char *s_word_extract_first ( char *s );
