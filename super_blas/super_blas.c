# include <stdlib.h>
# include <stdio.h>

# include "super_blas.h"

/*******************************************************************************/

double c_abs ( complex *z )

/*******************************************************************************/
/*
  Purpose:

    C_ABS returns the absolute value of a complex number.

  Parameters:

    Input, complex *Z, the number whose absolute value is desired.

    Output, double C_ABS, the absolute value of *Z.
*/
{
  float temp;
  float real = z->r;
  float imag = z->i;

  if ( real < 0 ) 
  {
    real = -real;
  }

  if ( imag < 0 ) 
  {
    imag = -imag;
  }

  if ( real < imag )
  {
    temp = real;
    real = imag;
    imag = temp;
  }

  if ( ( real + imag ) == real ) 
  { 
    return ( real );
  }
  
  temp = imag / real;
  temp = real * sqrt ( 1.0 + temp * temp );  /*overflow!!*/
  return ( temp );
}
/*******************************************************************************/

double c_abs1 ( complex *z )

/*******************************************************************************/
/*
  Purpose:

    C_ABS1 returns the L1 norm of a complex number.

  Parameters:

    Input, complex *Z, the complex number.

    Output, double C_ABS1, the L1 norm of *Z.
*/
{
  float real = z->r;
  float imag = z->i;
  
  if (real < 0) real = -real;
  if (imag < 0) imag = -imag;

  return (real + imag);
}
/*******************************************************************************/

void c_div ( complex *c, complex *a, complex *b )

/*******************************************************************************/
/*
  Purpose:

    C_DIV carries out complex division.

  Parameters:

    Output, complex *C, the result of the division.

    Input, complex *A, *B, the numerator and denominator.
*/
{
  float abi;
  float abr;
  float ci;
  float cr;
  float den;
  float ratio;
  
  if( (abr = b->r) < 0.0 )
	abr = -abr;

  if( (abi = b->i) < 0.0 )
	abi = -abi;

  if( abr <= abi ) 
  {
    if (abi == 0.0 ) 
    {
      fprintf(stderr, "z_div.c: division by zero");
      exit (-1);
    }	  
    ratio = b->r / b->i ;
    den = b->i * (1 + ratio*ratio);
    cr = (a->r*ratio + a->i) / den;
    ci = (a->i*ratio - a->r) / den;
  } 
  else
  {
    ratio = b->i / b->r ;
    den = b->r * (1 + ratio*ratio);
    cr = (a->r + a->i*ratio) / den;
    ci = (a->i - a->r*ratio) / den;
  }

  c->r = cr;
  c->i = ci;
}
/*******************************************************************************/

void c_exp ( complex *r, complex *z )

/*******************************************************************************/
/*
  Purpose:

    C_EXP returns the result of complex exponentiation.

  Parameters:

    Output, complex *R, the value of EXP ( *Z ).

    Input, complex *Z, the number to be exponentiated.
*/
{
  float expx;

  expx = exp ( z->r );
  r->r = expx * cos ( z->i );
  r->i = expx * sin ( z->i );
}
/*******************************************************************************/

int caxpy_ ( integer *n, complex *ca, complex *cx, integer *incx, 
  complex *cy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    CAXPY adds a multiple of one complex vector to another.

  Modified:

    15 May 2004

  Author:

    Jack Dongarra

  Parameters:
*/
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    static integer i, ix, iy;


/*
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if ((r__1 = ca->r, dabs(r__1)) + (r__2 = r_imag(ca), dabs(r__2)) == 0.f) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = iy;
	i__3 = iy;
	i__4 = ix;
	q__2.r = ca->r * CX(ix).r - ca->i * CX(ix).i, q__2.i = ca->r * CX(
		ix).i + ca->i * CX(ix).r;
	q__1.r = CY(iy).r + q__2.r, q__1.i = CY(iy).i + q__2.i;
	CY(iy).r = q__1.r, CY(iy).i = q__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	i__4 = i;
	q__2.r = ca->r * CX(i).r - ca->i * CX(i).i, q__2.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	q__1.r = CY(i).r + q__2.r, q__1.i = CY(i).i + q__2.i;
	CY(i).r = q__1.r, CY(i).i = q__1.i;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

int ccopy_ ( integer *n, complex *cx, integer *incx, complex *cy, 
  integer *incy )

/*******************************************************************************/
/*
  Purpose:

    CCOPY copies one complex vector to another.

  Modified:

    15 May 2004

  Author:

    Jack Dongarra

  Parameters:
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i, ix, iy;


/*    
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = iy;
	i__3 = ix;
	CY(iy).r = CX(ix).r, CY(iy).i = CX(ix).i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	CY(i).r = CX(i).r, CY(i).i = CX(i).i;
/* L30: */
    }
    return 0;
} /* ccopy_ */
/*******************************************************************************/

VOID cdotc_ ( complex * ret_val, integer *n, complex *cx, 
  integer *incx, complex *cy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    CDOTC computes the conjugated dot product of two complex vectors.

  Modified:

    15 May 2004

  Author:

    Jack Dongarra

  Parameters:
*/
{
    /* System generated locals */
    integer i__1, i__2;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i;
    static complex ctemp;
    static integer ix, iy;


/*
   Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    ctemp.r = 0.f, ctemp.i = 0.f;
     ret_val->r = 0.f,  ret_val->i = 0.f;
    if (*n <= 0) {
	return ;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	r_cnjg(&q__3, &cx[ix]);
	i__2 = iy;
	q__2.r = q__3.r * cy[iy].r - q__3.i * cy[iy].i, q__2.i = q__3.r * 
		cy[iy].i + q__3.i * cy[iy].r;
	q__1.r = ctemp.r + q__2.r, q__1.i = ctemp.i + q__2.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	r_cnjg(&q__3, &cx[i]);
	i__2 = i;
	q__2.r = q__3.r * cy[i].r - q__3.i * cy[i].i, q__2.i = q__3.r * 
		cy[i].i + q__3.i * cy[i].r;
	q__1.r = ctemp.r + q__2.r, q__1.i = ctemp.i + q__2.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
/* L30: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;
} /* cdotc_ */
/*******************************************************************************/

int cgemv_ ( char *trans, integer *m, integer *n, complex *
	alpha, complex *a, integer *lda, complex *x, integer *incx, complex *
	beta, complex *y, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    CGEMV computes the product of a general matrix with a vector.

  Discussion:

    CGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or   

       y := alpha*conjg( A' )*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

  Modified:

    15 May 2004

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs. 

  Parameters:

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - COMPLEX          array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX         .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX          array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer lenx, leny, i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj;


/*
    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("CGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r 
	    == 1.f && beta->i == 0.f)) {
	return 0;
    }

    noconj = lsame_(trans, "T");

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (beta->r != 1.f || beta->i != 0.f) {
	if (*incy == 1) {
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    Y(i).r = 0.f, Y(i).i = 0.f;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    i__3 = i;
		    q__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    q__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0.f, Y(iy).i = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    q__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    q__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0.f || X(jx).i != 0.f) {
		    i__2 = jx;
		    q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i;
			i__4 = i;
			i__5 = i + j * a_dim1;
			q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				q__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + 
				q__2.i;
			Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0.f || X(jx).i != 0.f) {
		    i__2 = jx;
		    q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = iy;
			i__4 = iy;
			i__5 = i + j * a_dim1;
			q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				q__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + 
				q__2.i;
			Y(iy).r = q__1.r, Y(iy).i = q__1.i;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
 */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0.f, temp.i = 0.f;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i;
			q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(i)
				.i, q__2.i = A(i,j).r * X(i).i + A(i,j)
				.i * X(i).r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L90: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			r_cnjg(&q__3, &A(i,j));
			i__3 = i;
			q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				q__2.i = q__3.r * X(i).i + q__3.i * X(i)
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jy += *incy;
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0.f, temp.i = 0.f;
		ix = kx;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = ix;
			q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(ix)
				.i, q__2.i = A(i,j).r * X(ix).i + A(i,j)
				.i * X(ix).r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L120: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			r_cnjg(&q__3, &A(i,j));
			i__3 = ix;
			q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				q__2.i = q__3.r * X(ix).i + q__3.i * X(ix)
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L130: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jy += *incy;
/* L140: */
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

int cgerc_ ( integer *m, integer *n, complex *alpha, 
  complex *x, integer *incx, complex *y, integer *incy, 
  complex *a, integer *lda )

/*******************************************************************************/
/*
  Purpose:   

    CGERC performs the rank 1 operation  A := alpha*x*conjg( y' ) + A,   

  Discussion:

    alpha is a scalar, x is an m element vector, y is an n element 
    vector and A is an m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("CGERC ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0.f || Y(jy).i != 0.f) {
		r_cnjg(&q__2, &Y(jy));
		q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			alpha->r * q__2.i + alpha->i * q__2.r;
		temp.r = q__1.r, temp.i = q__1.i;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = i;
		    q__2.r = X(i).r * temp.r - X(i).i * temp.i, q__2.i =
			     X(i).r * temp.i + X(i).i * temp.r;
		    q__1.r = A(i,j).r + q__2.r, q__1.i = A(i,j).i + q__2.i;
		    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0.f || Y(jy).i != 0.f) {
		r_cnjg(&q__2, &Y(jy));
		q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			alpha->r * q__2.i + alpha->i * q__2.r;
		temp.r = q__1.r, temp.i = q__1.i;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = ix;
		    q__2.r = X(ix).r * temp.r - X(ix).i * temp.i, q__2.i =
			     X(ix).r * temp.i + X(ix).i * temp.r;
		    q__1.r = A(i,j).r + q__2.r, q__1.i = A(i,j).i + q__2.i;
		    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;
}
/*******************************************************************************/

int chemv_(char *uplo, integer *n, complex *alpha, complex *a, 
  integer *lda, complex *x, integer *incx, complex *beta, complex *y,
  integer *incy )

/*******************************************************************************/
/*
  Purpose: 

    CHEMV performs the matrix-vector operation y := alpha*A*x + beta*y,   

  Discussion:

    alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n hermitian matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the hermitian matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the hermitian matrix and the strictly   
             upper triangular part of A is not referenced.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set and are assumed to be zero.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX         .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("CHEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r == 1.f && 
	    beta->i == 0.f)) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (beta->r != 1.f || beta->i != 0.f) {
	if (*incy == 1) {
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    Y(i).r = 0.f, Y(i).i = 0.f;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    i__3 = i;
		    q__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    q__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0.f && beta->i == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0.f, Y(iy).i = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    q__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    q__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }
    if (lsame_(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		q__1.r = alpha->r * X(j).r - alpha->i * X(j).i, q__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + q__2.i;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
		    r_cnjg(&q__3, &A(i,j));
		    i__3 = i;
		    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, q__2.i =
			     q__3.r * X(i).i + q__3.i * X(i).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L50: */
		}
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		q__3.r = d__1 * temp1.r, q__3.i = d__1 * temp1.i;
		q__2.r = Y(j).r + q__3.r, q__2.i = Y(j).i + q__3.i;
		q__4.r = alpha->r * temp2.r - alpha->i * temp2.i, q__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, q__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + q__2.i;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    r_cnjg(&q__3, &A(i,j));
		    i__3 = ix;
		    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, q__2.i =
			     q__3.r * X(ix).i + q__3.i * X(ix).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		q__3.r = d__1 * temp1.r, q__3.i = d__1 * temp1.i;
		q__2.r = Y(jy).r + q__3.r, q__2.i = Y(jy).i + q__3.i;
		q__4.r = alpha->r * temp2.r - alpha->i * temp2.i, q__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		q__1.r = alpha->r * X(j).r - alpha->i * X(j).i, q__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		q__2.r = d__1 * temp1.r, q__2.i = d__1 * temp1.i;
		q__1.r = Y(j).r + q__2.r, q__1.i = Y(j).i + q__2.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(i).r + q__2.r, q__1.i = Y(i).i + q__2.i;
		    Y(i).r = q__1.r, Y(i).i = q__1.i;
		    r_cnjg(&q__3, &A(i,j));
		    i__3 = i;
		    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, q__2.i =
			     q__3.r * X(i).i + q__3.i * X(i).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L90: */
		}
		i__2 = j;
		i__3 = j;
		q__2.r = alpha->r * temp2.r - alpha->i * temp2.i, q__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = Y(j).r + q__2.r, q__1.i = Y(j).i + q__2.i;
		Y(j).r = q__1.r, Y(j).i = q__1.i;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		q__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, q__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = q__1.r, temp1.i = q__1.i;
		temp2.r = 0.f, temp2.i = 0.f;
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		q__2.r = d__1 * temp1.r, q__2.i = d__1 * temp1.i;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    q__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    q__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    q__1.r = Y(iy).r + q__2.r, q__1.i = Y(iy).i + q__2.i;
		    Y(iy).r = q__1.r, Y(iy).i = q__1.i;
		    r_cnjg(&q__3, &A(i,j));
		    i__3 = ix;
		    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, q__2.i =
			     q__3.r * X(ix).i + q__3.i * X(ix).r;
		    q__1.r = temp2.r + q__2.r, q__1.i = temp2.i + q__2.i;
		    temp2.r = q__1.r, temp2.i = q__1.i;
/* L110: */
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp2.r - alpha->i * temp2.i, q__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		q__1.r = Y(jy).r + q__2.r, q__1.i = Y(jy).i + q__2.i;
		Y(jy).r = q__1.r, Y(jy).i = q__1.i;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

int cher2_(char *uplo, integer *n, complex *alpha, complex *x, 
  integer *incx, complex *y, integer *incy, complex *a, integer *lda )

/*******************************************************************************/
/*
  Purpose:

    CHER2 performs A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,   

  Discussion:

    This is a Hermitian rank 2 operation.  

    alpha is a scalar, x and y are n element vectors and A is an n 
    by n hermitian matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX         .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the hermitian matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the hermitian matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set, they are assumed to be zero, and on exit they   
             are set to zero.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("CHER2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0.f && alpha->i == 0.f) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0.f || X(j).i != 0.f || (Y(j).r != 0.f 
			|| Y(j).i != 0.f)) {
		    r_cnjg(&q__2, &Y(j));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = j;
		    q__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    q__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			q__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				q__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = i;
			q__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				q__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L10: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    q__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    q__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    q__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    q__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0.f || X(jx).i != 0.f || (Y(jy).r != 0.f 
			|| Y(jy).i != 0.f)) {
		    r_cnjg(&q__2, &Y(jy));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = jx;
		    q__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    ix = kx;
		    iy = ky;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			q__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				q__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = iy;
			q__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				q__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    q__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    q__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    q__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    q__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0.f || X(j).i != 0.f || (Y(j).r != 0.f 
			|| Y(j).i != 0.f)) {
		    r_cnjg(&q__2, &Y(j));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = j;
		    q__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    q__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    q__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    q__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    q__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    q__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			q__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				q__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = i;
			q__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				q__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L50: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0.f || X(jx).i != 0.f || (Y(jy).r != 0.f 
			|| Y(jy).i != 0.f)) {
		    r_cnjg(&q__2, &Y(jy));
		    q__1.r = alpha->r * q__2.r - alpha->i * q__2.i, q__1.i = 
			    alpha->r * q__2.i + alpha->i * q__2.r;
		    temp1.r = q__1.r, temp1.i = q__1.i;
		    i__2 = jx;
		    q__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    q__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    r_cnjg(&q__1, &q__2);
		    temp2.r = q__1.r, temp2.i = q__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    q__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    q__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    q__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    q__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		    d__1 = A(j,j).r + q__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			ix += *incx;
			iy += *incy;
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			q__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				q__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			q__2.r = A(i,j).r + q__3.r, q__2.i = A(i,j).i + 
				q__3.i;
			i__6 = iy;
			q__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				q__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
			A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L70: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.f;
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

void clsolve ( int ldm, int ncol, complex *M, complex *rhs )

/*******************************************************************************/
/*
  Purpose:

    CLSOLVE solves a dense UNIT lower triangular system.

  Discussion:

    The unit lower triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
    The solution will be returned in the rhs vector.

 */
{
    int k;
    complex x0, x1, x2, x3, temp;
    complex *M0;
    complex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;

    M0 = &M[0];


    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      	Mki0 = M0 + 1;
      	Mki1 = Mki0 + ldm + 1;
      	Mki2 = Mki1 + ldm + 1;
      	Mki3 = Mki2 + ldm + 1;

      	x0 = rhs[firstcol];
      	cc_mult(&temp, &x0, Mki0); Mki0++;
      	c_sub(&x1, &rhs[firstcol+1], &temp);
      	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x2, &rhs[firstcol+2], &temp);
	cc_mult(&temp, &x1, Mki1); Mki1++;
	c_sub(&x2, &x2, &temp);
      	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x3, &rhs[firstcol+3], &temp);
	cc_mult(&temp, &x1, Mki1); Mki1++;
	c_sub(&x3, &x3, &temp);
	cc_mult(&temp, &x2, Mki2); Mki2++;
	c_sub(&x3, &x3, &temp);

 	rhs[++firstcol] = x1;
      	rhs[++firstcol] = x2;
      	rhs[++firstcol] = x3;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    cc_mult(&temp, &x0, Mki0); Mki0++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x1, Mki1); Mki1++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x2, Mki2); Mki2++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x3, Mki3); Mki3++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	}

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol];
	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x1, &rhs[firstcol+1], &temp);

      	rhs[++firstcol] = x1;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    cc_mult(&temp, &x0, Mki0); Mki0++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x1, Mki1); Mki1++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	} 
    }
    
}
/*******************************************************************************/

void cmatvec ( int ldm, int nrow, int ncol, complex *M, complex *vec, 
  complex *Mxvec )

/*******************************************************************************/
/*
  Purpose:

    CMATVEC performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.

  Discussion:

    The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
 */
{
    complex vi0, vi1, vi2, vi3;
    complex *M0, temp;
    complex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;
    int k;

    M0 = &M[0];

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */
	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) {
	    cc_mult(&temp, &vi0, Mki0); Mki0++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi1, Mki1); Mki1++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi2, Mki2); Mki2++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi3, Mki3); Mki3++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	}

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */
 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++) {
	    cc_mult(&temp, &vi0, Mki0); Mki0++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	}
	M0 += ldm;
    }
	
}
/*******************************************************************************/

int cscal_ ( integer *n, complex *ca, complex *cx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    CSCAL scales a vector by a constant.

  Author:

    Jack Dongarra

  Parameters:
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1;

    /* Local variables */
    static integer i, nincx;


/*
   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	i__3 = i;
	i__4 = i;
	q__1.r = ca->r * CX(i).r - ca->i * CX(i).i, q__1.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	CX(i).r = q__1.r, CX(i).i = q__1.i;
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */

L20:
    i__2 = *n;
    for (i = 1; i <= *n; ++i) {
	i__1 = i;
	i__3 = i;
	q__1.r = ca->r * CX(i).r - ca->i * CX(i).i, q__1.i = ca->r * CX(
		i).i + ca->i * CX(i).r;
	CX(i).r = q__1.r, CX(i).i = q__1.i;
/* L30: */
    }
    return 0;
} /* cscal_ */

/*******************************************************************************/

int ctrsv_ ( char *uplo, char *trans, char *diag, integer *n, 
  complex *a, integer *lda, complex *x, integer *incx )

/*******************************************************************************/
/*  
  Purpose: 

    CTRSV solves a unit/nonunit upper or lower triangular linear system.
  
  Discussion:

    CTRSV solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   conjg( A' )*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj, nounit;


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("CTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    noconj = lsame_(trans, "T");
    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    if (X(j).r != 0.f || X(j).i != 0.f) {
			if (nounit) {
			    i__1 = j;
			    c_div(&q__1, &X(j), &A(j,j));
			    X(j).r = q__1.r, X(j).i = q__1.i;
			}
			i__1 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			for (i = j - 1; i >= 1; --i) {
			    i__1 = i;
			    i__2 = i;
			    i__3 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(i).r - q__2.r, q__1.i = X(i).i - 
				    q__2.i;
			    X(i).r = q__1.r, X(i).i = q__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    i__1 = jx;
		    if (X(jx).r != 0.f || X(jx).i != 0.f) {
			if (nounit) {
			    i__1 = jx;
			    c_div(&q__1, &X(jx), &A(j,j));
			    X(jx).r = q__1.r, X(jx).i = q__1.i;
			}
			i__1 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    i__1 = ix;
			    i__2 = ix;
			    i__3 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(ix).r - q__2.r, q__1.i = X(ix).i - 
				    q__2.i;
			    X(ix).r = q__1.r, X(ix).i = q__1.i;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    if (X(j).r != 0.f || X(j).i != 0.f) {
			if (nounit) {
			    i__2 = j;
			    c_div(&q__1, &X(j), &A(j,j));
			    X(j).r = q__1.r, X(j).i = q__1.i;
			}
			i__2 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    i__3 = i;
			    i__4 = i;
			    i__5 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(i).r - q__2.r, q__1.i = X(i).i - 
				    q__2.i;
			    X(i).r = q__1.r, X(i).i = q__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = jx;
		    if (X(jx).r != 0.f || X(jx).i != 0.f) {
			if (nounit) {
			    i__2 = jx;
			    c_div(&q__1, &X(jx), &A(j,j));
			    X(jx).r = q__1.r, X(jx).i = q__1.i;
			}
			i__2 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    i__3 = ix;
			    i__4 = ix;
			    i__5 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(ix).r - q__2.r, q__1.i = X(ix).i - 
				    q__2.i;
			    X(ix).r = q__1.r, X(ix).i = q__1.i;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = i;
			    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, q__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L90: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__3 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__2 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L110: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    ix = kx;
		    i__2 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = ix;
			    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, q__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix += *incx;
/* L120: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__3 = ix;
			    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				    q__2.i = q__3.r * X(ix).i + q__3.i * X(
				    ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix += *incx;
/* L130: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__2 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx += *incx;
/* L140: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = i;
			    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, q__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L150: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__2 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L160: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__1 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L170: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    ix = kx;
		    i__1 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = ix;
			    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, q__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix -= *incx;
/* L180: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__2 = ix;
			    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				    q__2.i = q__3.r * X(ix).i + q__3.i * X(
				    ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix -= *incx;
/* L190: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__1 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx -= *incx;
/* L200: */
		}
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

void cusolve ( int ldm, int ncol, complex *M, complex *rhs )

/*******************************************************************************/
/*
  Purpose:

    CUSOLVE solves a dense upper triangular system. 

    The upper triangular matrix is stored in a 2-dim array M(1:ldm,1:ncol). 
    The solution will be returned in the rhs vector.
 */
{
    complex xj, temp;
    int jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	c_div(&xj, &rhs[jcol], &M[jcol + jcol*ldm]); /* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++) {
	    cc_mult(&temp, &xj, &M[irow+jcol*ldm]); /* M(irow, jcol) */
	    c_sub(&rhs[irow], &rhs[irow], &temp);
	}

	jcol--;

    }
}
/*******************************************************************************/

void d_cnjg ( doublecomplex *r, doublecomplex *z )

/*******************************************************************************/
/* 
  Purpose:

    D_CNJG returns the complex conjugate of a complex number.
*/
{
    r->r = z->r;
    r->i = -z->i;
}
/*******************************************************************************/

double d_imag ( doublecomplex *z )

/*******************************************************************************/
/* 
  Purpose:

    D_IMAG returns the imaginary part of a complex number.
*/
{
    return (z->i);
}
/*******************************************************************************/

doublereal dasum_ ( integer *n, doublereal *dx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DASUM returns the sum of the absolute values of the entries of a vector.

  Author:

    Jack Dongarra

  Parameters:

*/
{


    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer nincx, mp1;


/*

   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	dtemp += (d__1 = DX(i), abs(d__1));
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	dtemp += (d__1 = DX(i), abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 6) {
	dtemp = dtemp + (d__1 = DX(i), abs(d__1)) + (d__2 = DX(i + 1), abs(
		d__2)) + (d__3 = DX(i + 2), abs(d__3)) + (d__4 = DX(i + 3), 
		abs(d__4)) + (d__5 = DX(i + 4), abs(d__5)) + (d__6 = DX(i + 5)
		, abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
}
/*******************************************************************************/

int daxpy_( integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    DAXPY adds a multiple of a vector to another vector.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*
    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DY(iy) += *da * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
	DY(i) += *da * DX(i);
	DY(i + 1) += *da * DX(i + 1);
	DY(i + 2) += *da * DX(i + 2);
	DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    return 0;
} /* daxpy_ */

/*******************************************************************************/

doublereal dcabs1_ ( doublecomplex *z )

/*******************************************************************************/
/*
  Purpose:

    DCABS1 returns the L1 norm of a complex number.
*/
{
/*  

       System generated locals */
    doublereal ret_val;
    static doublecomplex equiv_0[1];

    /* Local variables */
#define t ((doublereal *)equiv_0)
#define zz (equiv_0)

    zz->r = z->r, zz->i = z->i;
    ret_val = abs(t[0]) + abs(t[1]);
    return ret_val;
} /* dcabs1_ */

#undef zz
#undef t
/*******************************************************************************/

int dcopy_ ( integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    DCOPY copies one double precision vector to another.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*

    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DY(iy) = DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) = DX(i);
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
	DY(i) = DX(i);
	DY(i + 1) = DX(i + 1);
	DY(i + 2) = DX(i + 2);
	DY(i + 3) = DX(i + 3);
	DY(i + 4) = DX(i + 4);
	DY(i + 5) = DX(i + 5);
	DY(i + 6) = DX(i + 6);
/* L50: */
    }
    return 0;
}
/*******************************************************************************/

doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)

/*******************************************************************************/
/*
  Purpose:

    DDOT computes the dot product of two double precision vectors.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*
 
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp += DX(ix) * DY(iy);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	dtemp += DX(i) * DY(i);
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
	dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
		DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
}
/*******************************************************************************/

int dgemv_(char *trans, integer *m, integer *n, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy)

/*******************************************************************************/
/*
  Purpose:

    DGEMV adds a general matrix-vector product to a vector.

  Discussion:

    DGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   



*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			Y(i) += temp * A(i,j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			Y(iy) += temp * A(i,j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(i);
/* L90: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(ix);
		    ix += *incx;
/* L110: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

int dger_ ( integer *m, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda )

/*******************************************************************************/
/*  
  Purpose:

    DGER performs the rank 1 operation A := alpha*x*y' + A,   

  Discussion:

    alpha is a scalar, x is an m element vector, y is an n element 
    vector and A is an m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*

       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DGER  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.) {
		temp = *alpha * Y(jy);
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(i) * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.) {
		temp = *alpha * Y(jy);
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(ix) * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;
}
/*******************************************************************************/

void dlsolve ( int ldm, int ncol, double *M, double *rhs )

/*******************************************************************************/
/*
  Purpose:

    DLSOLVE solves a dense UNIT lower triangular system. 

  Discussion:

    The unit lower triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
    The solution will be returned in the rhs vector.

*/
{
    int k;
    double x0, x1, x2, x3, x4, x5, x6, x7;
    double *M0;
    register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;

    M0 = &M[0];

    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;
      Mki4 = Mki3 + ldm + 1;
      Mki5 = Mki4 + ldm + 1;
      Mki6 = Mki5 + ldm + 1;
      Mki7 = Mki6 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;
      x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
      x4 = rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++;
      x5 = rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++;
      x6 = rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++;
      x7 = rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
			   - x6 * *Mki6++;

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      rhs[++firstcol] = x4;
      rhs[++firstcol] = x5;
      rhs[++firstcol] = x6;
      rhs[++firstcol] = x7;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	                - x2 * *Mki2++ - x3 * *Mki3++
                        - x4 * *Mki4++ - x5 * *Mki5++
			- x6 * *Mki6++ - x7 * *Mki7++;
 
      M0 += 8 * ldm + 8;
    }

    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;
      x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	                - x2 * *Mki2++ - x3 * *Mki3++;
 
      M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;

      rhs[++firstcol] = x1;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
 
    }
    
}
/*******************************************************************************/

void dmatvec ( int ldm, int nrow, int ncol, double *M, double *vec, double *Mxvec )

/*******************************************************************************/
/*
  Purpose:

    DMATVEC performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.

  Discussion:

    The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
*/
{
    double vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
    double *M0;
    register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;
    int k;

    M0 = &M[0];
    while ( firstcol < ncol - 7 ) {	/* Do 8 columns */

	Mki0 = M0;
	Mki1 = Mki0 + ldm;
        Mki2 = Mki1 + ldm;
        Mki3 = Mki2 + ldm;
	Mki4 = Mki3 + ldm;
	Mki5 = Mki4 + ldm;
	Mki6 = Mki5 + ldm;
	Mki7 = Mki6 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	vi4 = vec[firstcol++];
	vi5 = vec[firstcol++];
	vi6 = vec[firstcol++];
	vi7 = vec[firstcol++];	

	for (k = 0; k < nrow; k++) 
	    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
		      + vi2 * *Mki2++ + vi3 * *Mki3++ 
		      + vi4 * *Mki4++ + vi5 * *Mki5++
		      + vi6 * *Mki6++ + vi7 * *Mki7++;

	M0 += 8 * ldm;
    }

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */

	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) 
	    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
		      + vi2 * *Mki2++ + vi3 * *Mki3++ ;

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */

 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++)
	    Mxvec[k] += vi0 * *Mki0++;

	M0 += ldm;
    }
	
}
/*******************************************************************************/

doublereal dnrm2_ ( integer *n, doublereal *x, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DNRM2 returns the euclidean norm of a double precision vector.
*/
{


    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal norm, scale, absxi;
    static integer ix;
    static doublereal ssq;


/*  

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = abs(X(1));
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.) {
		absxi = (d__1 = X(ix), abs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;
}
/*******************************************************************************/

int drot_ ( integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c, doublereal *s )

/*******************************************************************************/
/*
  Purpose:

    DROT applies a plane rotation.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    static doublereal dtemp;
    static integer ix, iy;


/* 
  
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(ix) + *s * DY(iy);
	DY(iy) = *c * DY(iy) - *s * DX(ix);
	DX(ix) = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(i) + *s * DY(i);
	DY(i) = *c * DY(i) - *s * DX(i);
	DX(i) = dtemp;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

int dscal_ ( integer *n, doublereal *da, doublereal *dx, 
	integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DSCAL multiplies a double precision vector by a constant.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
	DX(i) = *da * DX(i);
	DX(i + 1) = *da * DX(i + 1);
	DX(i + 2) = *da * DX(i + 2);
	DX(i + 3) = *da * DX(i + 3);
	DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    return 0;
}
/*******************************************************************************/

int dsymv_ ( char *uplo, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal 
	*beta, doublereal *y, integer *incy )

/*******************************************************************************/
/*  
  Purpose:

    DSYMV performs the symmetric matrix-vector operation y := alpha*A*x + beta*y,   

  Discussion:

    alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters: 

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("DSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L50: */
		}
		Y(j) = Y(j) + temp1 * A(j,j) + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		Y(jy) = Y(jy) + temp1 * A(j,j) + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		Y(j) += temp1 * A(j,j);
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L90: */
		}
		Y(j) += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		Y(jy) += temp1 * A(j,j);
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
/* L110: */
		}
		Y(jy) += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

int dsyr2_ ( char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda )

/*******************************************************************************/
/*
  Purpose:

    DSYR2 performs the symmetric rank 2 operation A := alpha*x*y' + alpha*y*x' + A,   

  Discussion:

    alpha is a scalar, x and y are n element vectors and A is an n 
    by n symmetric matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DSYR2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

int dtrsv_ ( char *uplo, char *trans, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *x, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DTRSV solves a unit or nonunit upper or lower triangular linear system.

  Discussion:

    The routine solves

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			for (i = j - 1; i >= 1; --i) {
			    X(i) -= temp * A(i,j);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    X(ix) -= temp * A(i,j);
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    X(i) -= temp * A(i,j);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    X(ix) -= temp * A(i,j);
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(i);
/* L90: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = kx;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(ix);
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(i);
/* L130: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = kx;
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(ix);
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

void dusolve ( int ldm, int ncol, double *M, double *rhs )

/*******************************************************************************/
/*
  Purpose:

    DUSOLVE solves a dense upper triangular system. 

  Discussion:

    The upper triangular matrix is stored in a 2-dim array M(1:ldm,1:ncol). 
    The solution will be returned in the rhs vector.
*/
{
    double xj;
    int jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	xj = rhs[jcol] / M[jcol + jcol*ldm]; 		/* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++)
	    rhs[irow] -= xj * M[irow + jcol*ldm];	/* M(irow, jcol) */

	jcol--;

    }
}
/*******************************************************************************/

doublereal dzasum_ ( integer *n, doublecomplex *zx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DZASUM sums the absolute values of the entries of a double precision complex vector.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i;
    static doublereal stemp;
    extern doublereal dcabs1_(doublecomplex *);
    static integer ix;


/*

   Parameter adjustments   
       Function Body */
#define ZX(I) zx[(I)-1]


    ret_val = 0.;
    stemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp += dcabs1_(&ZX(ix));
	ix += *incx;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp += dcabs1_(&ZX(i));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
}
/*******************************************************************************/

doublereal dznrm2_ ( integer *n, doublecomplex *x, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    DZNRM2 returns the euclidean norm of a double precision complex vector.
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static doublereal temp, norm, scale;
    static integer ix;
    static doublereal ssq;


/*

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL ZLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    i__3 = ix;
	    if (X(ix).r != 0.) {
		i__3 = ix;
		temp = (d__1 = X(ix).r, abs(d__1));
		if (scale < temp) {
/* Computing 2nd power */
		    d__1 = scale / temp;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    d__1 = temp / scale;
		    ssq += d__1 * d__1;
		}
	    }
	    if (d_imag(&X(ix)) != 0.) {
		temp = (d__1 = d_imag(&X(ix)), abs(d__1));
		if (scale < temp) {
/* Computing 2nd power */
		    d__1 = scale / temp;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    d__1 = temp / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;
}
/*******************************************************************************/

integer icamax_ ( integer *n, complex *cx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    ICAMAX returns the index of the entry of maximum absolute value in a complex vector.

  Author:

    Jack Dongarra
*/
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static real smax;
    static integer i, ix;
/* 

   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }
/*        code for increment not equal to 1 */
    ix = 1;
    smax = (r__1 = CX(1).r, dabs(r__1)) + (r__2 = r_imag(&CX(1)), dabs(r__2));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = ix;
	if ((r__1 = CX(ix).r, dabs(r__1)) + (r__2 = r_imag(&CX(ix)), dabs(
		r__2)) <= smax) {
	    goto L5;
	}
	ret_val = i;
	i__2 = ix;
	smax = (r__1 = CX(ix).r, dabs(r__1)) + (r__2 = r_imag(&CX(ix)), 
		dabs(r__2));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;
/*        code for increment equal to 1 */
L20:
    smax = (r__1 = CX(1).r, dabs(r__1)) + (r__2 = r_imag(&CX(1)), dabs(r__2));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = i;
	if ((r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(i)), dabs(
		r__2)) <= smax) {
	    goto L30;
	}
	ret_val = i;
	i__2 = i;
	smax = (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(i)), dabs(
		r__2));
L30:
	;
    }
    return ret_val;
}
/*******************************************************************************/

integer idamax_ ( integer *n, doublereal *dx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    IDAMAX returns the index of the entry of maximum absolute value in a double precision vector.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal dmax__;
    static integer i, ix;


/*  
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax__ = abs(DX(1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(ix), abs(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(ix), abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax__ = abs(DX(1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(i), abs(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(i), abs(d__1));
L30:
	;
    }
    return ret_val;
}
/*******************************************************************************/

integer isamax_ ( integer *n, real *sx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    ISAMAX returns the index of the entry of maximum absolute value in a real vector.

  Author:

    Jack Dongarra
*/
{
/* 
  System generated locals 
*/
  integer ret_val;
  integer i__1;
  real r__1;

/* 
  Local variables 
*/
  static real smax;
  static integer i, ix;

#define SX(I) sx[(I)-1]

  ret_val = 0;

  if ( *n < 1 || *incx <= 0 ) 
  {
    return ret_val;
  }

  ret_val = 1;

  if ( *n == 1 ) 
  {
    return ret_val;
  }
  if ( *incx == 1 )
  {
    goto L20;
  }
/*
  Code for increment not equal to 1.
*/
  ix = 1;
  smax = dabs( SX(1) );
  ix += *incx;
  i__1 = *n;
  for (i = 2; i <= *n; ++i) 
  {
    if ( ( r__1 = SX(ix), dabs ( r__1 ) ) <= smax) 
    {
      goto L5;
    }
    ret_val = i;
    smax = (r__1 = SX(ix), dabs ( r__1 ) );
L5:
    ix += *incx;
/* L10: */
  }
  return ret_val;
/*
  Code for increment equal to 1.
*/
L20:
  smax = dabs ( SX(1) );
  i__1 = *n;
  for ( i = 2; i <= *n; ++i ) 
  {
    if ( (r__1 = SX(i), dabs ( r__1 ) ) <= smax ) 
    {
      goto L30;
    }
    ret_val = i;
    smax = (r__1 = SX(i), dabs ( r__1 ) );
L30:
  ;
  }
  return ret_val;
}
/*******************************************************************************/

integer izamax_ ( integer *n, doublecomplex *zx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    IZAMAX returns the index of the entry of maximum absolute value in a double precision complex vector.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static doublereal smax;
    static integer i;
    extern doublereal dcabs1_(doublecomplex *);
    static integer ix;


/*
    
   Parameter adjustments   
       Function Body */
#define ZX(I) zx[(I)-1]


    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smax = dcabs1_(&ZX(1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if (dcabs1_(&ZX(ix)) <= smax) {
	    goto L5;
	}
	ret_val = i;
	smax = dcabs1_(&ZX(ix));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smax = dcabs1_(&ZX(1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if (dcabs1_(&ZX(i)) <= smax) {
	    goto L30;
	}
	ret_val = i;
	smax = dcabs1_(&ZX(i));
L30:
	;
    }
    return ret_val;
}
/*******************************************************************************/

int lsame_ ( char *ca, char *cb )

/*******************************************************************************/
/*
  Purpose:

    LSAME returns TRUE if CA is the same letter as CB regardless of case.

  Parameters:

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   
*/
{

    /* System generated locals */
    int ret_val;
    
    /* Local variables */
    int inta, intb, zcode;

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

    /* Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

    /* Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {
	/* ASCII is assumed - ZCODE is the ASCII code of either lower or   
          upper case 'Z'. */
	if (inta >= 97 && inta <= 122) inta += -32;
	if (intb >= 97 && intb <= 122) intb += -32;

    } else if (zcode == 233 || zcode == 169) {
	/* EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or   
          upper case 'Z'. */
	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169)
	    inta += 64;
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169)
	    intb += 64;
    } else if (zcode == 218 || zcode == 250) {
	/* ASCII is assumed, on Prime machines - ZCODE is the ASCII code   
          plus 128 of either lower or upper case 'Z'. */
	if (inta >= 225 && inta <= 250) inta += -32;
	if (intb >= 225 && intb <= 250) intb += -32;
    }
    ret_val = inta == intb;
    return ret_val;
}
/*******************************************************************************/

void r_cnjg ( complex *r, complex *z )

/*******************************************************************************/
/*
  Purpose:

    R_CNJG returns the complex conjugate of a complex number

  Modified:

    15 May 2004

  Parameters:

    Output, complex *R, the complex conjugate of Z.

    Input, complex *Z, the number whose complex conjugate is to be computed.
*/
{
  r->r = z->r;
  r->i = -z->i;
}
/*******************************************************************************/

double r_imag ( complex *z )

/*******************************************************************************/
/*
  Purpose:

    R_IMAG returns the imaginary part of a complex number.

  Modified:

    15 May 2004

  Parameters:

    Input, complex *Z, the number whose imaginary part is desired.

    Output, double R_IMAG, the imaginary part of Z.
*/
{
  return ( z->i );
}
/*******************************************************************************/

real sasum_ ( integer *n, real *sx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    SASUM sums the absolute values of the entries of a real vector.

  Author:

    Jack Dongarra

*/
{


    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static integer i, m, nincx;
    static real stemp;
    static integer mp1;


/*

   Parameter adjustments   
       Function Body */
#define SX(I) sx[(I)-1]


    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	stemp += (r__1 = SX(i), dabs(r__1));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	stemp += (r__1 = SX(i), dabs(r__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 6) {
	stemp = stemp + (r__1 = SX(i), dabs(r__1)) + (r__2 = SX(i + 1), dabs(
		r__2)) + (r__3 = SX(i + 2), dabs(r__3)) + (r__4 = SX(i + 3), 
		dabs(r__4)) + (r__5 = SX(i + 4), dabs(r__5)) + (r__6 = SX(i + 
		5), dabs(r__6));
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
}
/*******************************************************************************/

int saxpy_ ( integer *n, real *sa, real *sx, integer *incx, 
	real *sy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    SAXPY adds a multiple of one vector to another.

  Author:

    Jack Dongarra

*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*
    
   Parameter adjustments   
       Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*sa == 0.f) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	SY(iy) += *sa * SX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	SY(i) += *sa * SX(i);
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
	SY(i) += *sa * SX(i);
	SY(i + 1) += *sa * SX(i + 1);
	SY(i + 2) += *sa * SX(i + 2);
	SY(i + 3) += *sa * SX(i + 3);
/* L50: */
    }
    return 0;
}
/*******************************************************************************/

real scasum_ ( integer *n, complex *cx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    SCASUM sums the absolute values of the entries of a complex vector.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    static integer i, nincx;
    static real stemp;


/* 
    
   Parameter adjustments   
       Function Body */
#define CX(I) cx[(I)-1]


    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	i__3 = i;
	stemp = stemp + (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(
		i)), dabs(r__2));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */

L20:
    i__2 = *n;
    for (i = 1; i <= *n; ++i) {
	i__1 = i;
	stemp = stemp + (r__1 = CX(i).r, dabs(r__1)) + (r__2 = r_imag(&CX(
		i)), dabs(r__2));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
}
/*******************************************************************************/

real scnrm2_ ( integer *n, complex *x, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    SCNRM2 returns the euclidean norm of a complex vector.
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    real ret_val, r__1;

    /* Builtin functions */
    double r_imag(complex *), sqrt(doublereal);

    /* Local variables */
    static real temp, norm, scale;
    static integer ix;
    static real ssq;


/*      
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.f;
    } else {
	scale = 0.f;
	ssq = 1.f;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL CLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    i__3 = ix;
	    if (X(ix).r != 0.f) {
		i__3 = ix;
		temp = (r__1 = X(ix).r, dabs(r__1));
		if (scale < temp) {
/* Computing 2nd power */
		    r__1 = scale / temp;
		    ssq = ssq * (r__1 * r__1) + 1.f;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    r__1 = temp / scale;
		    ssq += r__1 * r__1;
		}
	    }
	    if (r_imag(&X(ix)) != 0.f) {
		temp = (r__1 = r_imag(&X(ix)), dabs(r__1));
		if (scale < temp) {
/* Computing 2nd power */
		    r__1 = scale / temp;
		    ssq = ssq * (r__1 * r__1) + 1.f;
		    scale = temp;
		} else {
/* Computing 2nd power */
		    r__1 = temp / scale;
		    ssq += r__1 * r__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

}
/*******************************************************************************/

int scopy_ ( integer *n, real *sx, integer *incx, real *sy, 
	integer *incy )

/*******************************************************************************/
/*
  Purpose:

    SCOPY copies one real vector to another.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*         
   Parameter adjustments   
       Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	SY(iy) = SX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	SY(i) = SX(i);
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
	SY(i) = SX(i);
	SY(i + 1) = SX(i + 1);
	SY(i + 2) = SX(i + 2);
	SY(i + 3) = SX(i + 3);
	SY(i + 4) = SX(i + 4);
	SY(i + 5) = SX(i + 5);
	SY(i + 6) = SX(i + 6);
/* L50: */
    }
    return 0;
}
/*******************************************************************************/

real sdot_ ( integer *n, real *sx, integer *incx, real *sy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    SDOT computes the dot product of two real vectors.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i, m;
    static real stemp;
    static integer ix, iy, mp1;


/* 
   Parameter adjustments   
       Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]


    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp += SX(ix) * SY(iy);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	stemp += SX(i) * SY(i);
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
	stemp = stemp + SX(i) * SY(i) + SX(i + 1) * SY(i + 1) + SX(i + 2) * 
		SY(i + 2) + SX(i + 3) * SY(i + 3) + SX(i + 4) * SY(i + 4);
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
}
/*******************************************************************************/

int sgemv_ ( char *trans, integer *m, integer *n, real *alpha, 
	real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy )

/*******************************************************************************/
/*  
  Purpose:

    SGEMV adds a general matrix-vector product to a vector.

  Discussion:

    The routine performs one of the operations:

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - REAL             array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - REAL            .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - REAL             array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static real temp;
    static integer lenx, leny, i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*

       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("SGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = 0.f;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.f) {
	return 0;
    }
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.f) {
		    temp = *alpha * X(jx);
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			Y(i) += temp * A(i,j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.f) {
		    temp = *alpha * X(jx);
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			Y(iy) += temp * A(i,j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.f;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(i);
/* L90: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.f;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(ix);
		    ix += *incx;
/* L110: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

int sger_ ( integer *m, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda )

/*******************************************************************************/
/*  
  Purpose:

    SGER performs the rank 1 operation A := alpha*x*y' + A,   

  Discussion:

    alpha is a scalar, x is an m element vector, y is an n element 
    vector and A is an m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - REAL             array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static real temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*

       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SGER  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.f) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.f) {
		temp = *alpha * Y(jy);
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(i) * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.f) {
		temp = *alpha * Y(jy);
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(ix) * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

}
/*******************************************************************************/

void slsolve ( int ldm, int ncol, float *M, float *rhs )

/*******************************************************************************/
/*
  Purpose:

    SLSOLVE solves a dense UNIT lower triangular system. 

  Discussion:

    The unit lower triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
    The solution will be returned in the rhs vector.
 */
{
    int k;
    float x0, x1, x2, x3, x4, x5, x6, x7;
    float *M0;
    register float *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;

    M0 = &M[0];

    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;
      Mki4 = Mki3 + ldm + 1;
      Mki5 = Mki4 + ldm + 1;
      Mki6 = Mki5 + ldm + 1;
      Mki7 = Mki6 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;
      x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
      x4 = rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++;
      x5 = rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++;
      x6 = rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++;
      x7 = rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
	                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
			   - x6 * *Mki6++;

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      rhs[++firstcol] = x4;
      rhs[++firstcol] = x5;
      rhs[++firstcol] = x6;
      rhs[++firstcol] = x7;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	                - x2 * *Mki2++ - x3 * *Mki3++
                        - x4 * *Mki4++ - x5 * *Mki5++
			- x6 * *Mki6++ - x7 * *Mki7++;
 
      M0 += 8 * ldm + 8;
    }

    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;
      x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	                - x2 * *Mki2++ - x3 * *Mki3++;
 
      M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;

      x0 = rhs[firstcol];
      x1 = rhs[firstcol+1] - x0 * *Mki0++;

      rhs[++firstcol] = x1;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++)
	rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
 
    }
    
}
/*******************************************************************************/

void smatvec ( int ldm, int nrow, int ncol, float *M, float *vec, float *Mxvec )

/*******************************************************************************/
/* 
  Purpose:

    SMATVEC performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.

  Discussion:

    The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
*/
{
    float vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
    float *M0;
    register float *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;
    int k;

    M0 = &M[0];
    while ( firstcol < ncol - 7 ) {	/* Do 8 columns */

	Mki0 = M0;
	Mki1 = Mki0 + ldm;
        Mki2 = Mki1 + ldm;
        Mki3 = Mki2 + ldm;
	Mki4 = Mki3 + ldm;
	Mki5 = Mki4 + ldm;
	Mki6 = Mki5 + ldm;
	Mki7 = Mki6 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	vi4 = vec[firstcol++];
	vi5 = vec[firstcol++];
	vi6 = vec[firstcol++];
	vi7 = vec[firstcol++];	

	for (k = 0; k < nrow; k++) 
	    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
		      + vi2 * *Mki2++ + vi3 * *Mki3++ 
		      + vi4 * *Mki4++ + vi5 * *Mki5++
		      + vi6 * *Mki6++ + vi7 * *Mki7++;

	M0 += 8 * ldm;
    }

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */

	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) 
	    Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
		      + vi2 * *Mki2++ + vi3 * *Mki3++ ;

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */

 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++)
	    Mxvec[k] += vi0 * *Mki0++;

	M0 += ldm;
    }
	
}
/*******************************************************************************/

real snrm2_ ( integer *n, real *x, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    SNRM2 returns the euclidean norm of a real vector.
*/
{


    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real norm, scale, absxi;
    static integer ix;
    static real ssq;


/*      
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.f;
    } else if (*n == 1) {
	norm = dabs(X(1));
    } else {
	scale = 0.f;
	ssq = 1.f;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.f) {
		absxi = (r__1 = X(ix), dabs(r__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    r__1 = scale / absxi;
		    ssq = ssq * (r__1 * r__1) + 1.f;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    r__1 = absxi / scale;
		    ssq += r__1 * r__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

}
/*******************************************************************************/

int sp_ienv ( int ispec )

/*******************************************************************************/
/*
  Purpose:  

    SP_IENV chooses machine-dependent parameters for the local environment. 

  Discussion:

    See ISPEC for a description of the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

  Parameters:

    ISPEC   (input) int
            Specifies the parameter to be returned as the value of SP_IENV.   
            = 1: the panel size w; a panel consists of w consecutive
	         columns of matrix A in the process of Gaussian elimination.
		 The best value depends on machine's cache characters.
            = 2: the relaxation parameter relax; if the number of
	         nodes (columns) in a subtree of the elimination tree is less
		 than relax, this subtree is considered as one supernode,
		 regardless of their row structures.
            = 3: the maximum size for a supernode;
	    = 4: the minimum row dimension for 2-D blocking to be used;
	    = 5: the minimum column dimension for 2-D blocking to be used;
	    = 6: the estimated fills factor for L and U, compared with A;
	    
   (SP_IENV) (output) int
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if SP_IENV = -k, the k-th argument had an illegal value. 
*/
{
    int i;

    switch (ispec) {
	case 1: return (8);
	case 2: return (1);
	case 3: return (100);
	case 4: return (200);
	case 5: return (40);
        case 6: return (20);
    }

    /* Invalid value for ISPEC */
    i = 1;
    xerbla_("sp_ienv", &i);
    return 0;

}
/*******************************************************************************/

int srot_ ( integer *n, real *sx, integer *incx, real *sy, 
	integer *incy, real *c, real *s )

/*******************************************************************************/
/*
  Purpose:

    SROT applies a Givens plane rotation.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    static real stemp;
    static integer ix, iy;


/* 
    
   Parameter adjustments   
       Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp = *c * SX(ix) + *s * SY(iy);
	SY(iy) = *c * SY(iy) - *s * SX(ix);
	SX(ix) = stemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	stemp = *c * SX(i) + *s * SY(i);
	SY(i) = *c * SY(i) - *s * SX(i);
	SX(i) = stemp;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

int sscal_ ( integer *n, real *sa, real *sx, integer *incx )

/*******************************************************************************/
/*
  Purpose:

    SSCAL scales a vector by a constant.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*        
   Parameter adjustments   
       Function Body */
#define SX(I) sx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	SX(i) = *sa * SX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	SX(i) = *sa * SX(i);
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
	SX(i) = *sa * SX(i);
	SX(i + 1) = *sa * SX(i + 1);
	SX(i + 2) = *sa * SX(i + 2);
	SX(i + 3) = *sa * SX(i + 3);
	SX(i + 4) = *sa * SX(i + 4);
/* L50: */
    }
    return 0;
}
/*******************************************************************************/

int ssymv_ ( char *uplo, integer *n, real *alpha, real *a, 
	integer *lda, real *x, integer *incx, real *beta, real *y, integer *
	incy )

/*******************************************************************************/
/*  
  Purpose:

    SSYMV performs the symmetric matrix-vector operation y := alpha*A*x + beta*y.

  Discussion:

    alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - REAL            .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static real temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);



/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("SSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = 0.f;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = 0.f;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.f) {
	return 0;
    }
    if (lsame_(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.f;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L50: */
		}
		Y(j) = Y(j) + temp1 * A(j,j) + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.f;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		Y(jy) = Y(jy) + temp1 * A(j,j) + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.f;
		Y(j) += temp1 * A(j,j);
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    Y(i) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(i);
/* L90: */
		}
		Y(j) += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.f;
		Y(jy) += temp1 * A(j,j);
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    Y(iy) += temp1 * A(i,j);
		    temp2 += A(i,j) * X(ix);
/* L110: */
		}
		Y(jy) += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

int ssyr2_ ( char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *y, integer *incy, real *a, integer *lda )

/*******************************************************************************/
/*  
  Purpose:  

    SSYR2 performs the symmetric rank 2 operation A := alpha*x*y' + alpha*y*x' + A,   

  Discussion:

    alpha is a scalar, x and y are n element vectors and A is an n 
    by n symmetric matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - REAL            .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

*/

    
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static real temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("SSYR2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.f) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0.f || Y(j) != 0.f) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.f || Y(jy) != 0.f) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0.f || Y(j) != 0.f) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.f || Y(jy) != 0.f) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

int strsv_ ( char *uplo, char *trans, char *diag, integer *n, 
	real *a, integer *lda, real *x, integer *incx )

/*******************************************************************************/
/*  
  Purpose:  

    STRSV solves a unit/nonunit upper or lower triangular linear system.

  Discussion:

    The routine solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - REAL             array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - REAL             array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static real temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("STRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.f) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			for (i = j - 1; i >= 1; --i) {
			    X(i) -= temp * A(i,j);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.f) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    X(ix) -= temp * A(i,j);
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.f) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    X(i) -= temp * A(i,j);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.f) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    X(ix) -= temp * A(i,j);
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(i);
/* L90: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = kx;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(ix);
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(i);
/* L130: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = kx;
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(ix);
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

void susolve ( int ldm, int ncol, float *M, float *rhs )

/*******************************************************************************/
/*
  Purpose:

    SUSOLVE solves a dense upper triangular system. 

  Discussion:

    The upper triangular matrix is stored in a 2-dim array M(1:ldm,1:ncol). 
    The solution will be returned in the rhs vector.
*/
{
    float xj;
    int jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	xj = rhs[jcol] / M[jcol + jcol*ldm]; 		/* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++)
	    rhs[irow] -= xj * M[irow + jcol*ldm];	/* M(irow, jcol) */

	jcol--;

    }
}
/*******************************************************************************/

int xerbla_ ( char *srname, int *info )

/*******************************************************************************/
/*
  Purpose:

    XERBLA  is an error handler for the LAPACK routines.   

  Discussion:

    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

  Parameters:

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INT   
            The position of the invalid parameter in the parameter list   
            of the calling routine.   
*/
{


    printf("** On entry to %6s, parameter number %2d had an illegal value\n",
		srname, *info);

    return 0;
}
/*******************************************************************************/

double z_abs(doublecomplex *z)

/*******************************************************************************/
/* 
  Purpose:

    Z_ABS returns the absolute value of a double precision complex number.
*/
{
    double temp;
    double real = z->r;
    double imag = z->i;

    if (real < 0) real = -real;
    if (imag < 0) imag = -imag;
    if (imag > real) {
	temp = real;
	real = imag;
	imag = temp;
    }
    if ((real+imag) == real) return(real);
  
    temp = imag/real;
    temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
    return (temp);
}
/*******************************************************************************/

double z_abs1 ( doublecomplex *z )

/*******************************************************************************/
/*
  Purpose:

    Z_ABS1 returns the L1 norm of a double precision complex number.
*/
{
    double real = z->r;
    double imag = z->i;
  
    if (real < 0) real = -real;
    if (imag < 0) imag = -imag;

    return (real + imag);
}
/*******************************************************************************/

void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b)

/*******************************************************************************/
/*
  Purpose:

    Z_DIV carries out complex double precision division.
*/
{
    double ratio, den;
    double abr, abi, cr, ci;
  
    if( (abr = b->r) < 0.)
	abr = - abr;
    if( (abi = b->i) < 0.)
	abi = - abi;
    if( abr <= abi ) {
	if (abi == 0) {
	    fprintf(stderr, "z_div.c: division by zero");
	    exit (-1);
	}	  
	ratio = b->r / b->i ;
	den = b->i * (1 + ratio*ratio);
	cr = (a->r*ratio + a->i) / den;
	ci = (a->i*ratio - a->r) / den;
    } else {
	ratio = b->i / b->r ;
	den = b->r * (1 + ratio*ratio);
	cr = (a->r + a->i*ratio) / den;
	ci = (a->i - a->r*ratio) / den;
    }
    c->r = cr;
    c->i = ci;
}
/*******************************************************************************/

void z_exp(doublecomplex *r, doublecomplex *z)

/*******************************************************************************/
/*
  Purpose:

    Z_EXP carries out complex double precision exponentiation.
*/
{
    double expx;

    expx = exp(z->r);
    r->r = expx * cos(z->i);
    r->i = expx * sin(z->i);
}
/*******************************************************************************/

int zaxpy_ ( integer *n, doublecomplex *za, doublecomplex *zx, 
	integer *incx, doublecomplex *zy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    ZAXPY adds a multiple of a complex double precision vector to another.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i;
    extern doublereal dcabs1_(doublecomplex *);
    static integer ix, iy;


/*
   Parameter adjustments   
       Function Body */
#define ZY(I) zy[(I)-1]
#define ZX(I) zx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (dcabs1_(za) == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = iy;
	i__3 = iy;
	i__4 = ix;
	z__2.r = za->r * ZX(ix).r - za->i * ZX(ix).i, z__2.i = za->r * ZX(
		ix).i + za->i * ZX(ix).r;
	z__1.r = ZY(iy).r + z__2.r, z__1.i = ZY(iy).i + z__2.i;
	ZY(iy).r = z__1.r, ZY(iy).i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	i__4 = i;
	z__2.r = za->r * ZX(i).r - za->i * ZX(i).i, z__2.i = za->r * ZX(
		i).i + za->i * ZX(i).r;
	z__1.r = ZY(i).r + z__2.r, z__1.i = ZY(i).i + z__2.i;
	ZY(i).r = z__1.r, ZY(i).i = z__1.i;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

int zcopy_ ( integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    ZCOPY copies one double precision complex vector to another.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i, ix, iy;


/* 
    
   Parameter adjustments   
       Function Body */
#define ZY(I) zy[(I)-1]
#define ZX(I) zx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = iy;
	i__3 = ix;
	ZY(iy).r = ZX(ix).r, ZY(iy).i = ZX(ix).i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	ZY(i).r = ZX(i).r, ZY(i).i = ZX(i).i;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

VOID zdotc_ ( doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy )

/*******************************************************************************/
/*
  Purpose:

    ZDOTC computes the conjugated dot product of two complex double precision vectors.

  Author:

    Jack Dongarra
*/
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i;
    static doublecomplex ztemp;
    static integer ix, iy;


/*
    
   Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    ztemp.r = 0., ztemp.i = 0.;
     ret_val->r = 0.,  ret_val->i = 0.;
    if (*n <= 0) {
	return ;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	d_cnjg(&z__3, &zx[ix]);
	i__2 = iy;
	z__2.r = z__3.r * zy[iy].r - z__3.i * zy[iy].i, z__2.i = z__3.r * 
		zy[iy].i + z__3.i * zy[iy].r;
	z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
	ztemp.r = z__1.r, ztemp.i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
     ret_val->r = ztemp.r,  ret_val->i = ztemp.i;
    return ;

/*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	d_cnjg(&z__3, &zx[i]);
	i__2 = i;
	z__2.r = z__3.r * zy[i].r - z__3.i * zy[i].i, z__2.i = z__3.r * 
		zy[i].i + z__3.i * zy[i].r;
	z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
	ztemp.r = z__1.r, ztemp.i = z__1.i;
/* L30: */
    }
     ret_val->r = ztemp.r,  ret_val->i = ztemp.i;
    return ;
}
/*******************************************************************************/

int zgemv_(char *trans, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *
	incy)

/*******************************************************************************/
/*  
  Purpose:

    ZGEMV adds a general matrix-vector product to a vector.

  Discussion:

    The routine carries out of one the operations:

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or   

       y := alpha*conjg( A' )*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX*16      .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - COMPLEX*16       array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX*16      .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX*16       array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


*/

{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer lenx, leny, i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj;


/*

       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("ZGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
	return 0;
    }

    noconj = lsame_(trans, "T");

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (beta->r != 1. || beta->i != 0.) {
	if (*incy == 1) {
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    Y(i).r = 0., Y(i).i = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = i;
		    i__3 = i;
		    z__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    z__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = z__1.r, Y(i).i = z__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0., Y(iy).i = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    z__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    z__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = z__1.r, Y(iy).i = z__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0. || X(jx).i != 0.) {
		    i__2 = jx;
		    z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i;
			i__4 = i;
			i__5 = i + j * a_dim1;
			z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				z__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			z__1.r = Y(i).r + z__2.r, z__1.i = Y(i).i + 
				z__2.i;
			Y(i).r = z__1.r, Y(i).i = z__1.i;
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		if (X(jx).r != 0. || X(jx).i != 0.) {
		    i__2 = jx;
		    z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__1.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    iy = ky;
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = iy;
			i__4 = iy;
			i__5 = i + j * a_dim1;
			z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				z__2.i = temp.r * A(i,j).i + temp.i * A(i,j)
				.r;
			z__1.r = Y(iy).r + z__2.r, z__1.i = Y(iy).i + 
				z__2.i;
			Y(iy).r = z__1.r, Y(iy).i = z__1.i;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
 */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0., temp.i = 0.;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i;
			z__2.r = A(i,j).r * X(i).r - A(i,j).i * X(i)
				.i, z__2.i = A(i,j).r * X(i).i + A(i,j)
				.i * X(i).r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
/* L90: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			d_cnjg(&z__3, &A(i,j));
			i__3 = i;
			z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				z__2.i = z__3.r * X(i).i + z__3.i * X(i)
				.r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
/* L100: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		z__1.r = Y(jy).r + z__2.r, z__1.i = Y(jy).i + z__2.i;
		Y(jy).r = z__1.r, Y(jy).i = z__1.i;
		jy += *incy;
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp.r = 0., temp.i = 0.;
		ix = kx;
		if (noconj) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = ix;
			z__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(ix)
				.i, z__2.i = A(i,j).r * X(ix).i + A(i,j)
				.i * X(ix).r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
			ix += *incx;
/* L120: */
		    }
		} else {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			d_cnjg(&z__3, &A(i,j));
			i__3 = ix;
			z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				z__2.i = z__3.r * X(ix).i + z__3.i * X(ix)
				.r;
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
			temp.r = z__1.r, temp.i = z__1.i;
			ix += *incx;
/* L130: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		z__1.r = Y(jy).r + z__2.r, z__1.i = Y(jy).i + z__2.i;
		Y(jy).r = z__1.r, Y(jy).i = z__1.i;
		jy += *incy;
/* L140: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

int zgerc_(integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *a, integer *lda)

/*******************************************************************************/
/*  
  Purpose:   

    ZGERC performs the rank 1 operation A := alpha*x*conjg( y' ) + A,   

  Discussion:

    alpha is a scalar, x is an m element vector, y is an n element 
    vector and A is an m by n matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:  

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX*16      .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   


*/

{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*

       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("ZGERC ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0. || Y(jy).i != 0.) {
		d_cnjg(&z__2, &Y(jy));
		z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			alpha->r * z__2.i + alpha->i * z__2.r;
		temp.r = z__1.r, temp.i = z__1.i;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = i;
		    z__2.r = X(i).r * temp.r - X(i).i * temp.i, z__2.i =
			     X(i).r * temp.i + X(i).i * temp.r;
		    z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + z__2.i;
		    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = jy;
	    if (Y(jy).r != 0. || Y(jy).i != 0.) {
		d_cnjg(&z__2, &Y(jy));
		z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			alpha->r * z__2.i + alpha->i * z__2.r;
		temp.r = z__1.r, temp.i = z__1.i;
		ix = kx;
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * a_dim1;
		    i__4 = i + j * a_dim1;
		    i__5 = ix;
		    z__2.r = X(ix).r * temp.r - X(ix).i * temp.i, z__2.i =
			     X(ix).r * temp.i + X(ix).i * temp.r;
		    z__1.r = A(i,j).r + z__2.r, z__1.i = A(i,j).i + z__2.i;
		    A(i,j).r = z__1.r, A(i,j).i = z__1.i;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

}
/*******************************************************************************/

int zhemv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy)

/*******************************************************************************/
/*  
  Purpose:

    ZHEMV adds a Hermitian matrix-vector product to a vector.

  Discussion:

    The routine performs the matrix-vector operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n hermitian matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX*16      .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the hermitian matrix and the strictly   
             lower triangular part of A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the hermitian matrix and the strictly   
             upper triangular part of A is not referenced.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set and are assumed to be zero.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - COMPLEX*16      .   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("ZHEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A.   

       First form  y := beta*y. */

    if (beta->r != 1. || beta->i != 0.) {
	if (*incy == 1) {
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    Y(i).r = 0., Y(i).i = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = i;
		    i__3 = i;
		    z__1.r = beta->r * Y(i).r - beta->i * Y(i).i, 
			    z__1.i = beta->r * Y(i).i + beta->i * Y(i)
			    .r;
		    Y(i).r = z__1.r, Y(i).i = z__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    Y(iy).r = 0., Y(iy).i = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    i__2 = iy;
		    i__3 = iy;
		    z__1.r = beta->r * Y(iy).r - beta->i * Y(iy).i, 
			    z__1.i = beta->r * Y(iy).i + beta->i * Y(iy)
			    .r;
		    Y(iy).r = z__1.r, Y(iy).i = z__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }
    if (lsame_(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		z__1.r = alpha->r * X(j).r - alpha->i * X(j).i, z__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = z__1.r, temp1.i = z__1.i;
		temp2.r = 0., temp2.i = 0.;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    z__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    z__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    z__1.r = Y(i).r + z__2.r, z__1.i = Y(i).i + z__2.i;
		    Y(i).r = z__1.r, Y(i).i = z__1.i;
		    d_cnjg(&z__3, &A(i,j));
		    i__3 = i;
		    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, z__2.i =
			     z__3.r * X(i).i + z__3.i * X(i).r;
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
		    temp2.r = z__1.r, temp2.i = z__1.i;
/* L50: */
		}
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
		z__2.r = Y(j).r + z__3.r, z__2.i = Y(j).i + z__3.i;
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
		Y(j).r = z__1.r, Y(j).i = z__1.i;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, z__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = z__1.r, temp1.i = z__1.i;
		temp2.r = 0., temp2.i = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    z__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    z__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    z__1.r = Y(iy).r + z__2.r, z__1.i = Y(iy).i + z__2.i;
		    Y(iy).r = z__1.r, Y(iy).i = z__1.i;
		    d_cnjg(&z__3, &A(i,j));
		    i__3 = ix;
		    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, z__2.i =
			     z__3.r * X(ix).i + z__3.i * X(ix).r;
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
		    temp2.r = z__1.r, temp2.i = z__1.i;
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
		z__2.r = Y(jy).r + z__3.r, z__2.i = Y(jy).i + z__3.i;
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
		Y(jy).r = z__1.r, Y(jy).i = z__1.i;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		z__1.r = alpha->r * X(j).r - alpha->i * X(j).i, z__1.i =
			 alpha->r * X(j).i + alpha->i * X(j).r;
		temp1.r = z__1.r, temp1.i = z__1.i;
		temp2.r = 0., temp2.i = 0.;
		i__2 = j;
		i__3 = j;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
		z__1.r = Y(j).r + z__2.r, z__1.i = Y(j).i + z__2.i;
		Y(j).r = z__1.r, Y(j).i = z__1.i;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    i__3 = i;
		    i__4 = i;
		    i__5 = i + j * a_dim1;
		    z__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    z__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    z__1.r = Y(i).r + z__2.r, z__1.i = Y(i).i + z__2.i;
		    Y(i).r = z__1.r, Y(i).i = z__1.i;
		    d_cnjg(&z__3, &A(i,j));
		    i__3 = i;
		    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, z__2.i =
			     z__3.r * X(i).i + z__3.i * X(i).r;
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
		    temp2.r = z__1.r, temp2.i = z__1.i;
/* L90: */
		}
		i__2 = j;
		i__3 = j;
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		z__1.r = Y(j).r + z__2.r, z__1.i = Y(j).i + z__2.i;
		Y(j).r = z__1.r, Y(j).i = z__1.i;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		z__1.r = alpha->r * X(jx).r - alpha->i * X(jx).i, z__1.i =
			 alpha->r * X(jx).i + alpha->i * X(jx).r;
		temp1.r = z__1.r, temp1.i = z__1.i;
		temp2.r = 0., temp2.i = 0.;
		i__2 = jy;
		i__3 = jy;
		i__4 = j + j * a_dim1;
		d__1 = A(j,j).r;
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
		z__1.r = Y(jy).r + z__2.r, z__1.i = Y(jy).i + z__2.i;
		Y(jy).r = z__1.r, Y(jy).i = z__1.i;
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    ix += *incx;
		    iy += *incy;
		    i__3 = iy;
		    i__4 = iy;
		    i__5 = i + j * a_dim1;
		    z__2.r = temp1.r * A(i,j).r - temp1.i * A(i,j).i, 
			    z__2.i = temp1.r * A(i,j).i + temp1.i * A(i,j)
			    .r;
		    z__1.r = Y(iy).r + z__2.r, z__1.i = Y(iy).i + z__2.i;
		    Y(iy).r = z__1.r, Y(iy).i = z__1.i;
		    d_cnjg(&z__3, &A(i,j));
		    i__3 = ix;
		    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, z__2.i =
			     z__3.r * X(ix).i + z__3.i * X(ix).r;
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
		    temp2.r = z__1.r, temp2.i = z__1.i;
/* L110: */
		}
		i__2 = jy;
		i__3 = jy;
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
		z__1.r = Y(jy).r + z__2.r, z__1.i = Y(jy).i + z__2.i;
		Y(jy).r = z__1.r, Y(jy).i = z__1.i;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

}
/*******************************************************************************/

int zher2_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *a, integer *lda)

/*******************************************************************************/
/*  
  Purpose:

    ZHER2 performs a hermitian rank 2 update.

  Discussion:

    The routine carries out the operation   

       A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,   

    where alpha is a scalar, x and y are n element vectors and A is an n 
    by n hermitian matrix.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - COMPLEX*16      .   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the hermitian matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the hermitian matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   
             Note that the imaginary parts of the diagonal elements need 
  
             not be set, they are assumed to be zero, and on exit they   
             are set to zero.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp1, temp2;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);


/*
     Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("ZHER2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0. || X(j).i != 0. || (Y(j).r != 0. || 
			Y(j).i != 0.)) {
		    d_cnjg(&z__2, &Y(j));
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
		    temp1.r = z__1.r, temp1.i = z__1.i;
		    i__2 = j;
		    z__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    z__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    d_cnjg(&z__1, &z__2);
		    temp2.r = z__1.r, temp2.i = z__1.i;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			z__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				z__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			z__2.r = A(i,j).r + z__3.r, z__2.i = A(i,j).i + 
				z__3.i;
			i__6 = i;
			z__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				z__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L10: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    z__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    z__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    z__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    z__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    d__1 = A(j,j).r + z__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0. || X(jx).i != 0. || (Y(jy).r != 0. || 
			Y(jy).i != 0.)) {
		    d_cnjg(&z__2, &Y(jy));
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
		    temp1.r = z__1.r, temp1.i = z__1.i;
		    i__2 = jx;
		    z__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    d_cnjg(&z__1, &z__2);
		    temp2.r = z__1.r, temp2.i = z__1.i;
		    ix = kx;
		    iy = ky;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			z__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				z__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			z__2.r = A(i,j).r + z__3.r, z__2.i = A(i,j).i + 
				z__3.i;
			i__6 = iy;
			z__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				z__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    z__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    z__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    z__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    z__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    d__1 = A(j,j).r + z__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		i__3 = j;
		if (X(j).r != 0. || X(j).i != 0. || (Y(j).r != 0. || 
			Y(j).i != 0.)) {
		    d_cnjg(&z__2, &Y(j));
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
		    temp1.r = z__1.r, temp1.i = z__1.i;
		    i__2 = j;
		    z__2.r = alpha->r * X(j).r - alpha->i * X(j).i, 
			    z__2.i = alpha->r * X(j).i + alpha->i * X(j)
			    .r;
		    d_cnjg(&z__1, &z__2);
		    temp2.r = z__1.r, temp2.i = z__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = j;
		    z__2.r = X(j).r * temp1.r - X(j).i * temp1.i, 
			    z__2.i = X(j).r * temp1.i + X(j).i * 
			    temp1.r;
		    i__5 = j;
		    z__3.r = Y(j).r * temp2.r - Y(j).i * temp2.i, 
			    z__3.i = Y(j).r * temp2.i + Y(j).i * 
			    temp2.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    d__1 = A(j,j).r + z__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = i;
			z__3.r = X(i).r * temp1.r - X(i).i * temp1.i, 
				z__3.i = X(i).r * temp1.i + X(i).i * 
				temp1.r;
			z__2.r = A(i,j).r + z__3.r, z__2.i = A(i,j).i + 
				z__3.i;
			i__6 = i;
			z__4.r = Y(i).r * temp2.r - Y(i).i * temp2.i, 
				z__4.i = Y(i).r * temp2.i + Y(i).i * 
				temp2.r;
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L50: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = jx;
		i__3 = jy;
		if (X(jx).r != 0. || X(jx).i != 0. || (Y(jy).r != 0. || 
			Y(jy).i != 0.)) {
		    d_cnjg(&z__2, &Y(jy));
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
		    temp1.r = z__1.r, temp1.i = z__1.i;
		    i__2 = jx;
		    z__2.r = alpha->r * X(jx).r - alpha->i * X(jx).i, 
			    z__2.i = alpha->r * X(jx).i + alpha->i * X(jx)
			    .r;
		    d_cnjg(&z__1, &z__2);
		    temp2.r = z__1.r, temp2.i = z__1.i;
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    i__4 = jx;
		    z__2.r = X(jx).r * temp1.r - X(jx).i * temp1.i, 
			    z__2.i = X(jx).r * temp1.i + X(jx).i * 
			    temp1.r;
		    i__5 = jy;
		    z__3.r = Y(jy).r * temp2.r - Y(jy).i * temp2.i, 
			    z__3.i = Y(jy).r * temp2.i + Y(jy).i * 
			    temp2.r;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    d__1 = A(j,j).r + z__1.r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			ix += *incx;
			iy += *incy;
			i__3 = i + j * a_dim1;
			i__4 = i + j * a_dim1;
			i__5 = ix;
			z__3.r = X(ix).r * temp1.r - X(ix).i * temp1.i, 
				z__3.i = X(ix).r * temp1.i + X(ix).i * 
				temp1.r;
			z__2.r = A(i,j).r + z__3.r, z__2.i = A(i,j).i + 
				z__3.i;
			i__6 = iy;
			z__4.r = Y(iy).r * temp2.r - Y(iy).i * temp2.i, 
				z__4.i = Y(iy).r * temp2.i + Y(iy).i * 
				temp2.r;
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
			A(i,j).r = z__1.r, A(i,j).i = z__1.i;
/* L70: */
		    }
		} else {
		    i__2 = j + j * a_dim1;
		    i__3 = j + j * a_dim1;
		    d__1 = A(j,j).r;
		    A(j,j).r = d__1, A(j,j).i = 0.;
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

void zlsolve ( int ldm, int ncol, doublecomplex *M, doublecomplex *rhs )

/*******************************************************************************/
/*
  Purpose:

    ZLSOLVE solves a dense UNIT lower triangular system. 

  Discussion:

    The unit lower triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
    The solution will be returned in the rhs vector.
*/
{
    int k;
    doublecomplex x0, x1, x2, x3, temp;
    doublecomplex *M0;
    doublecomplex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;

    M0 = &M[0];


    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      	Mki0 = M0 + 1;
      	Mki1 = Mki0 + ldm + 1;
      	Mki2 = Mki1 + ldm + 1;
      	Mki3 = Mki2 + ldm + 1;

      	x0 = rhs[firstcol];
      	zz_mult(&temp, &x0, Mki0); Mki0++;
      	z_sub(&x1, &rhs[firstcol+1], &temp);
      	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x2, &rhs[firstcol+2], &temp);
	zz_mult(&temp, &x1, Mki1); Mki1++;
	z_sub(&x2, &x2, &temp);
      	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x3, &rhs[firstcol+3], &temp);
	zz_mult(&temp, &x1, Mki1); Mki1++;
	z_sub(&x3, &x3, &temp);
	zz_mult(&temp, &x2, Mki2); Mki2++;
	z_sub(&x3, &x3, &temp);

 	rhs[++firstcol] = x1;
      	rhs[++firstcol] = x2;
      	rhs[++firstcol] = x3;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    zz_mult(&temp, &x0, Mki0); Mki0++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x1, Mki1); Mki1++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x2, Mki2); Mki2++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x3, Mki3); Mki3++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	}

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol];
	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x1, &rhs[firstcol+1], &temp);

      	rhs[++firstcol] = x1;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    zz_mult(&temp, &x0, Mki0); Mki0++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x1, Mki1); Mki1++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	} 
    }
    
}
/*******************************************************************************/

void zmatvec ( int ldm, int nrow, int ncol, doublecomplex *M, doublecomplex *vec, 
  doublecomplex *Mxvec )

/*******************************************************************************/
/*
  Purpose:

    ZMATVEC performs a dense matrix-vector multiply: 

  Discussion:

    Mxvec = Mxvec + M * vec.

    The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
 */
{
    doublecomplex vi0, vi1, vi2, vi3;
    doublecomplex *M0, temp;
    doublecomplex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;
    int k;

    M0 = &M[0];

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */
	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) {
	    zz_mult(&temp, &vi0, Mki0); Mki0++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi1, Mki1); Mki1++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi2, Mki2); Mki2++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi3, Mki3); Mki3++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	}

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */
 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++) {
	    zz_mult(&temp, &vi0, Mki0); Mki0++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	}
	M0 += ldm;
    }
	
}
/*******************************************************************************/

int zscal_(integer *n, doublecomplex *za, doublecomplex *zx, 
	integer *incx)

/*******************************************************************************/
/*
  Purpose:

    ZSCAL scales a double precision complex vector by a constant.

  Author:

    Jack Dongarra
*/
{


    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i, ix;


/*
    
   Parameter adjustments   
       Function Body */
#define ZX(I) zx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	i__3 = ix;
	z__1.r = za->r * ZX(ix).r - za->i * ZX(ix).i, z__1.i = za->r * ZX(
		ix).i + za->i * ZX(ix).r;
	ZX(ix).r = z__1.r, ZX(ix).i = z__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	i__3 = i;
	z__1.r = za->r * ZX(i).r - za->i * ZX(i).i, z__1.i = za->r * ZX(
		i).i + za->i * ZX(i).r;
	ZX(i).r = z__1.r, ZX(i).i = z__1.i;
/* L30: */
    }
    return 0;
}
/*******************************************************************************/

int ztrsv_(char *uplo, char *trans, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx)

/*******************************************************************************/
/*
  Purpose: 

    ZTRSV solves a unit or nonunit, upper or lower triangular linear system.

  Discussion:

    The routine solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

  Author:

    Jack Dongarra, Argonne National Lab.   
    Jeremy Du Croz, Nag Central Office.   
    Sven Hammarling, Nag Central Office.   
    Richard Hanson, Sandia National Labs.   

  Parameters:

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   conjg( A' )*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
*/
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj, nounit;


/*
       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("ZTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    noconj = lsame_(trans, "T");
    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			if (nounit) {
			    i__1 = j;
			    z_div(&z__1, &X(j), &A(j,j));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__1 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			for (i = j - 1; i >= 1; --i) {
			    i__1 = i;
			    i__2 = i;
			    i__3 = i + j * a_dim1;
			    z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    z__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    i__1 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			if (nounit) {
			    i__1 = jx;
			    z_div(&z__1, &X(jx), &A(j,j));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__1 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    i__1 = ix;
			    i__2 = ix;
			    i__3 = i + j * a_dim1;
			    z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    z__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			if (nounit) {
			    i__2 = j;
			    z_div(&z__1, &X(j), &A(j,j));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__2 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    i__3 = i;
			    i__4 = i;
			    i__5 = i + j * a_dim1;
			    z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    z__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			if (nounit) {
			    i__2 = jx;
			    z_div(&z__1, &X(jx), &A(j,j));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__2 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    i__3 = ix;
			    i__4 = ix;
			    i__5 = i + j * a_dim1;
			    z__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    z__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = i;
			    z__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, z__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L90: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(j,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    d_cnjg(&z__3, &A(i,j));
			    i__3 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L100: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(j,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__2 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L110: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    ix = kx;
		    i__2 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = ix;
			    z__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, z__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L120: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(j,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    d_cnjg(&z__3, &A(i,j));
			    i__3 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L130: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(j,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__2 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx += *incx;
/* L140: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = i;
			    z__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, z__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L150: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(j,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    d_cnjg(&z__3, &A(i,j));
			    i__2 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L160: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(j,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__1 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L170: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    ix = kx;
		    i__1 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = ix;
			    z__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, z__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L180: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(j,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    d_cnjg(&z__3, &A(i,j));
			    i__2 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L190: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(j,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__1 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx -= *incx;
/* L200: */
		}
	    }
	}
    }

    return 0;
}
/*******************************************************************************/

void zusolve ( int ldm, int ncol, doublecomplex *M, doublecomplex *rhs )

/*******************************************************************************/
/*
  Purpose:

    ZUSOLVE solves a dense upper triangular system. 

  Discussion:

    The upper triangular matrix is stored in a 2-dim array M(1:ldm,1:ncol). 

    The solution will be returned in the rhs vector.
*/
{
    doublecomplex xj, temp;
    int jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	z_div(&xj, &rhs[jcol], &M[jcol + jcol*ldm]); /* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++) {
	    zz_mult(&temp, &xj, &M[irow+jcol*ldm]); /* M(irow, jcol) */
	    z_sub(&rhs[irow], &rhs[irow], &temp);
	}

	jcol--;

    }
}


