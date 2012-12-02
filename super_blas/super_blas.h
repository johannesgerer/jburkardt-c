/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 1, 1997
 *
 */
#ifndef __SUPERLU_CNAMES /* allow multiple inclusions */
#define __SUPERLU_CNAMES

/*
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 */

#define ADD_       0
#define NOCHANGE   1
#define UPCASE     2
#define C_CALL     3

#ifdef UpCase
#define F77_CALL_C UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C ADD_
#endif

#ifndef F77_CALL_C
#define F77_CALL_C ADD_
#endif

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm_(...)
 *
 * This is the default.
 */

#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void DGEMM(...)
 */
#define isamax_   ISAMAX
#define sasum_    SASUM
#define saxpy_    SAXPY
#define scopy_    SCOPY
#define sdot_     SDOT
#define sgemm_    SGEMM
#define sgemv_    SGEMV
#define sger_     SGER
#define snrm2_    SNRM2
#define srot_     SROT
#define sscal_    SSCAL
#define ssymv_    SSYMV
#define ssyr2_    SSYR2
#define strsm_    STRSM
#define strsv_    STRSV

#define dasum_    SASUM
#define daxpy_    SAXPY
#define dcopy_    SCOPY
#define ddot_     SDOT
#define dgemm_    SGEMM
#define dgemv_    SGEMV
#define dger_     SGER
#define dnrm2_    SNRM2
#define drot_     SROT
#define dscal_    SSCAL
#define dsymv_    SSYMV
#define dsyr2_    SSYR2
#define dtrsm_    STRSM
#define dtrsv_    STRSV
#define idamax_   ISAMAX

#define caxpy_    CAXPY
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2
#define icamax_   ICAMAX
#define scasum_   SCASUM
#define scnrm2_   SCNRM2

#define dzasum_   SCASUM
#define izamax_   ICAMAX
#define zcopy_    CCOPY
#define zscal_    CSCAL
#define dznrm2_   SCNRM2
#define zaxpy_    CAXPY
#define zgemv_    CGEMV
#define ztrsv_    CTRSV
#define zgemm_    CGEMM
#define ztrsm_    CTRSM
#define zgerc_    CGERC
#define zhemv_    CHEMV
#define zher2_    CHER2

#define c_bridge_dgssv_ C_BRIDGE_DGSSV
#define c_fortran_dgssv_ C_FORTRAN_DGSSV
#endif

#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm(...)
 */
#define isamax_   isamax
#define sasum_    sasum
#define saxpy_    saxpy
#define scopy_    scopy
#define sdot_     sdot
#define sgemm_    sgemm
#define sgemv_    sgemv
#define sger_     sger
#define snrm2_    snrm2
#define srot_     srot
#define sscal_    sscal
#define ssymv_    ssymv
#define ssyr2_    ssyr2
#define strsm_    strsm
#define strsv_    strsv

#define dasum_    dasum
#define daxpy_    daxpy
#define dcopy_    dcopy
#define ddot_     ddot
#define dgemm_    dgemm
#define dgemv_    dgemv
#define dger_     dger
#define dnrm2_    dnrm2
#define drot_     drot
#define dscal_    dscal
#define dsymv_    dsymv
#define dsyr2_    dsyr2
#define dtrsm_    dtrsm
#define dtrsv_    dtrsv
#define idamax_   idamax

#define scasum_   scasum
#define icamax_   icamax
#define ccopy_    ccopy
#define cscal_    cscal
#define scnrm2_   scnrm2
#define caxpy_    caxpy
#define cgemv_    cgemv
#define ctrsv_    ctrsv
#define cgemm_    cgemm
#define ctrsm_    ctrsm
#define cgerc_    cgerc
#define chemv_    chemv
#define cher2_    cher2

#define dzasum_   dzasum
#define izamax_   izamax
#define zcopy_    zcopy
#define zscal_    zscal
#define dznrm2_   dznrm2
#define zaxpy_    zaxpy
#define zgemv_    zgemv
#define ztrsv_    ztrsv
#define zgemm_    zgemm
#define ztrsm_    ztrsm
#define zgerc_    zgerc
#define zhemv_    zhemv
#define zher2_    zher2

#define c_bridge_dgssv_ c_bridge_dgssv
#define c_fortran_dgssv_ c_fortran_dgssv
#endif

#endif /* __SUPERLU_CNAMES */


#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#if 0
typedef long int integer; /* 64 on 64-bit machine */
typedef long int logical;
#endif

typedef int integer;
typedef int logical;

typedef char *address;
typedef struct { float r, i; } complex;
typedef short int shortint;
typedef float real;
typedef struct { double r, i; } doublecomplex;
typedef double doublereal;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

#define VOID void

#endif

/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/* Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )

/* Complex Addition c = a + b */
#define c_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define c_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define cs_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define cc_mult(c, a, b) { \
	float cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

/* Complex equality testing */
#define c_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

double c_abs ( complex *z );
double c_abs1 ( complex *z );
void c_div ( complex *c, complex *a, complex *b );
void c_exp ( complex *r, complex *z );

int caxpy_ ( integer *n, complex *ca, complex *cx, integer *incx, 
  complex *cy, integer *incy );
int ccopy_ ( integer *n, complex *cx, integer *incx, complex *cy, 
  integer *incy );
VOID cdotc_ ( complex * ret_val, integer *n, complex *cx, 
  integer *incx, complex *cy, integer *incy );
int cgemv_ ( char *trans, integer *m, integer *n, complex *alpha, complex *a, 
  integer *lda, complex *x, integer *incx, complex *beta, complex *y, 
  integer *incy );
int cgerc_ ( integer *m, integer *n, complex *alpha, complex *x, 
  integer *incx, complex *y, integer *incy, complex *a, integer *lda );
int chemv_ ( char *uplo, integer *n, complex *alpha, complex *a, 
  integer *lda, complex *x, integer *incx, complex *beta, complex *y,
  integer *incy );
int cher2_ ( char *uplo, integer *n, complex *alpha, complex *x, 
  integer *incx, complex *y, integer *incy, complex *a, integer *lda );
void clsolve ( int ldm, int ncol, complex *M, complex *rhs );
void cmatvec ( int ldm, int nrow, int ncol, complex *M, complex *vec,
  complex *Mxvec );
int cscal_ ( integer *n, complex *ca, complex *cx, integer *incx );
int ctrsv_ ( char *uplo, char *trans, char *diag, integer *n, 
  complex *a, integer *lda, complex *x, integer *incx );
void cusolve ( int ldm, int ncol, complex *M, complex *rhs );

void d_cnjg(doublecomplex *r, doublecomplex *z);
double d_imag(doublecomplex *z);

doublereal dasum_ ( integer *n, doublereal *dx, integer *incx );
int daxpy_ ( integer *n, doublereal *da, doublereal *dx, 
  integer *incx, doublereal *dy, integer *incy );
doublereal dcabs1_ ( doublecomplex *z );
int dcopy_ ( integer *n, doublereal *dx, integer *incx, 
  doublereal *dy, integer *incy );
doublereal ddot_ ( integer *n, doublereal *dx, integer *incx, doublereal *dy, 
  integer *incy );
int dgemv_ ( char *trans, integer *m, integer *n, doublereal *alpha, 
  doublereal *a, integer *lda, doublereal *x, integer *incx, 
  doublereal *beta, doublereal *y, integer *incy );
int dger_ ( integer *m, integer *n, doublereal *alpha, 
  doublereal *x, integer *incx, doublereal *y, integer *incy, 
  doublereal *a, integer *lda );
void dlsolve ( int ldm, int ncol, double *M, double *rhs );
void dmatvec ( int ldm, int nrow, int ncol, double *M, double *vec, 
  double *Mxvec );
doublereal dnrm2_ ( integer *n, doublereal *x, integer *incx );
int drot_ ( integer *n, doublereal *dx, integer *incx, 
  doublereal *dy, integer *incy, doublereal *c, doublereal *s );
int dscal_ ( integer *n, doublereal *da, doublereal *dx, 
  integer *incx );
int dsymv_ ( char *uplo, integer *n, doublereal *alpha, 
  doublereal *a, integer *lda, doublereal *x, integer *incx, 
  doublereal *beta, doublereal *y, integer *incy );
int dsyr2_ ( char *uplo, integer *n, doublereal *alpha, 
  doublereal *x, integer *incx, doublereal *y, integer *incy, 
  doublereal *a, integer *lda );
int dtrsv_ ( char *uplo, char *trans, char *diag, integer *n, 
  doublereal *a, integer *lda, doublereal *x, integer *incx );
void dusolve ( int ldm, int ncol, double *M, double *rhs );
doublereal dzasum_ ( integer *n, doublecomplex *zx, integer *incx );
doublereal dznrm2_ ( integer *n, doublecomplex *x, integer *incx );

integer icamax_ ( integer *n, complex *cx, integer *incx );
integer idamax_ ( integer *n, doublereal *dx, integer *incx );
integer isamax_ ( integer *n, real *sx, integer *incx );
integer izamax_ ( integer *n, doublecomplex *zx, integer *incx );

int lsame_ ( char *ca, char *cb );

void r_cnjg ( complex *r, complex *z );
double r_imag ( complex *z );

real sasum_ ( integer *n, real *sx, integer *incx );
int saxpy_ ( integer *n, real *sa, real *sx, integer *incx, 
  real *sy, integer *incy );
real scasum_ ( integer *n, complex *cx, integer *incx );
real scnrm2_ ( integer *n, complex *x, integer *incx );
int scopy_ ( integer *n, real *sx, integer *incx, real *sy, 
  integer *incy );
real sdot_ ( integer *n, real *sx, integer *incx, real *sy, integer *incy );
int sgemv_ ( char *trans, integer *m, integer *n, real *alpha, 
  real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
  integer *incy );
int sger_ ( integer *m, integer *n, real *alpha, real *x, 
  integer *incx, real *y, integer *incy, real *a, integer *lda );
void slsolve ( int ldm, int ncol, float *M, float *rhs );
void smatvec ( int ldm, int nrow, int ncol, float *M, float *vec, 
  float *Mxvec );
real snrm2_ ( integer *n, real *x, integer *incx );
int sp_ienv ( int ispec );
int srot_ ( integer *n, real *sx, integer *incx, real *sy, 
  integer *incy, real *c, real *s );
int sscal_ ( integer *n, real *sa, real *sx, integer *incx );
int ssymv_( char *uplo, integer *n, real *alpha, real *a, 
  integer *lda, real *x, integer *incx, real *beta, real *y, 
  integer *incy );
int ssyr2_ ( char *uplo, integer *n, real *alpha, real *x, 
  integer *incx, real *y, integer *incy, real *a, integer *lda );
int strsv_ ( char *uplo, char *trans, char *diag, integer *n, 
  real *a, integer *lda, real *x, integer *incx );
void susolve ( int ldm, int ncol, float *M, float *rhs );

int xerbla_ ( char *srname, int *info );

double z_abs(doublecomplex *z);
double z_abs1(doublecomplex *z);
void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b);
void z_exp(doublecomplex *r, doublecomplex *z);

int zaxpy_ ( integer *n, doublecomplex *za, doublecomplex *zx, 
  integer *incx, doublecomplex *zy, integer *incy );
int zcopy_ ( integer *n, doublecomplex *zx, integer *incx, 
  doublecomplex *zy, integer *incy );
VOID zdotc_ ( doublecomplex * ret_val, integer *n, 
  doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy );
int zgemv_ ( char *trans, integer *m, integer *n, 
  doublecomplex *alpha, doublecomplex *a, integer *lda, 
  doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, 
  integer *incy );
int zgerc_ ( integer *m, integer *n, doublecomplex *alpha, 
  doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
  doublecomplex *a, integer *lda );
int zhemv_ ( char *uplo, integer *n, doublecomplex *alpha, 
  doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
  doublecomplex *beta, doublecomplex *y, integer *incy );
int zher2_ ( char *uplo, integer *n, doublecomplex *alpha, 
  doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
  doublecomplex *a, integer *lda );
void zlsolve ( int ldm, int ncol, doublecomplex *M, doublecomplex *rhs );
void zmatvec ( int ldm, int nrow, int ncol, doublecomplex *M, doublecomplex *vec, 
  doublecomplex *Mxvec );
int zscal_ ( integer *n, doublecomplex *za, doublecomplex *zx, 
  integer *incx );
int ztrsv_ ( char *uplo, char *trans, char *diag, integer *n, 
  doublecomplex *a, integer *lda, doublecomplex *x, integer *incx );
void zusolve ( int ldm, int ncol, doublecomplex *M, doublecomplex *rhs );
