double dznrm2 ( int n, _Complex double x[], int incx );
double dzasum ( int n, _Complex double x[], int incx );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int izamax ( int n, _Complex double x[], int incx );
int lsame ( char ca, char cb );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_sign ( double x );
void xerbla ( char *srname, int info );
double zabs1 ( _Complex double z );
double zabs2 ( _Complex double z );
void zaxpy ( int n, _Complex double ca, _Complex double cx[], 
  int incx, _Complex double cy[], int incy );
void zcopy ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy );
_Complex double zdotc ( int n, _Complex double cx[], int incx, 
  _Complex double cy[], int incy );
_Complex double zdotu ( int n, _Complex double cx[], int incx, 
  _Complex double cy[], int incy );
void zdrot ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy, double c, double s );
void zdscal ( int n, double sa, _Complex double cx[], int incx );
double zmach ( int job );
void zrotg ( _Complex double *ca, _Complex double cb, double *c, 
  _Complex double *s );
void zscal ( int n, _Complex double ca, _Complex double cx[], int incx );
_Complex double zsign1 ( _Complex double z1, _Complex double z2 );
_Complex double zsign2 ( _Complex double z1, _Complex double z2 );
void zswap ( int n, _Complex double cx[], int incx, _Complex double cy[], 
  int incy );
