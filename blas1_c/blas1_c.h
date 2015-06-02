void caxpy ( int n, float complex ca, float complex cx[], 
  int incx, float complex cy[], int incy );
void ccopy ( int n, float complex cx[], int incx, float complex cy[], 
  int incy );
float complex cdotc ( int n, float complex cx[], int incx, 
  float complex cy[], int incy );
float complex cdotu ( int n, float complex cx[], int incx, 
  float complex cy[], int incy );
void crotg ( float complex *ca, float complex cb, float *c, 
  float complex *s );
void cscal ( int n, float complex ca, float complex cx[], int incx );
void csrot ( int n, float complex cx[], int incx, float complex cy[], 
  int incy, float c, float s );
void csscal ( int n, float sa, float complex cx[], int incx );
void cswap ( int n, float complex cx[], int incx, float complex cy[], 
  int incy );
int icamax ( int n, float complex x[], int incx );
float scasum ( int n, float complex x[], int incx );
float scnrm2 ( int n, float complex x[], int incx );

