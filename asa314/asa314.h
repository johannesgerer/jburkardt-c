int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4mat_print ( int m, int n, int a[], char *title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void invmod ( int mat[], int imat[], int rmod[], int cmod[], int nrow, 
  int *ifault );
void msort ( int mat[], int imat[], int rmod[], int cmod[], int rsort[], 
  int csort[], int nrow );
void musort ( int mat[], int imat[], int rmod[], int cmod[], int rsort[], 
  int csort[], int nrow );
void timestamp ( );

