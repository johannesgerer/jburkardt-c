float cabs1 ( _Complex float z );
float cabs2 ( _Complex float z );
void caxpy ( int n, _Complex float ca, _Complex float cx[], 
  int incx, _Complex float cy[], int incy );
void ccopy ( int n, _Complex float cx[], int incx, _Complex float cy[], 
  int incy );
_Complex float cdotc ( int n, _Complex float cx[], int incx, 
  _Complex float cy[], int incy );
_Complex float cdotu ( int n, _Complex float cx[], int incx, 
  _Complex float cy[], int incy );
float cmach ( int job );
void crotg ( _Complex float *ca, _Complex float cb, float *c, 
  _Complex float *s );
void cscal ( int n, _Complex float ca, _Complex float cx[], int incx );
_Complex float csign1 ( _Complex float z1, _Complex float z2 );
_Complex float csign2 ( _Complex float z1, _Complex float z2 );
void csrot ( int n, _Complex float cx[], int incx, _Complex float cy[], 
  int incy, float c, float s );
void csscal ( int n, float sa, _Complex float cx[], int incx );
void cswap ( int n, _Complex float cx[], int incx, _Complex float cy[], 
  int incy );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int icamax ( int n, _Complex float x[], int incx );
int lsame ( char ca, char cb );
float r4_abs ( float x );
float r4_sign ( float x );
float scasum ( int n, _Complex float x[], int incx );
float scnrm2 ( int n, _Complex float x[], int incx );
void xerbla ( char *srname, int info );
