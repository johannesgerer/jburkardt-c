double dznrm2 ( int n, double complex x[], int incx );
double dzasum ( int n, double complex x[], int incx );
int izamax ( int n, double complex x[], int incx );
void zaxpy ( int n, double complex ca, double complex cx[], 
  int incx, double complex cy[], int incy );
void zcopy ( int n, double complex cx[], int incx, double complex cy[], 
  int incy );
double complex zdotc ( int n, double complex cx[], int incx, 
  double complex cy[], int incy );
double complex zdotu ( int n, double complex cx[], int incx, 
  double complex cy[], int incy );
void zdrot ( int n, double complex cx[], int incx, double complex cy[], 
  int incy, double c, double s );
void zdscal ( int n, double sa, double complex cx[], int incx );
void zrotg ( double complex *ca, double complex cb, double *c, 
  double complex *s );
void zscal ( int n, double complex ca, double complex cx[], int incx );
void zswap ( int n, double complex cx[], int incx, double complex cy[], 
  int incy );
