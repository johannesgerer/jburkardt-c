/*

Translated to C by Bonnie Toy 5/88
  (modified on 2/25/94  to fix a problem with daxpy  for
   unequal increments or equal increments not equal to 1.
     Jack Dongarra)

To compile single precision version for Sun-4:

	cc -DSP -O4 -fsingle -fsingle2 clinpack.c -lm

To compile double precision version for Sun-4:

	cc -DDP -O4 clinpack.c -lm

To obtain rolled source BLAS, add -DROLL to the command lines.
To obtain unrolled source BLAS, add -DUNROLL to the command lines.

You must specify one of -DSP or -DDP to compile correctly.

You must specify one of -DROLL or -DUNROLL to compile correctly.

*/

#ifdef SP
#define REAL float
#define ZERO 0.0
#define ONE 1.0
#define PREC "Single "
#endif

#ifdef DP
#define REAL double
#define ZERO 0.0e0
#define ONE 1.0e0
#define PREC "Double "
#endif

#define NTIMES 10

#ifdef ROLL
#define ROLLING "Rolled "
#endif
#ifdef UNROLL
#define ROLLING "Unrolled "
#endif

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

static REAL times[9][9];

int main ( void );
double cpu_time ( void );
void daxpy ( int n, REAL da, REAL dx[], int incx, REAL dy[], int incy );
REAL ddot ( int n, REAL dx[], int incx, REAL dy[], int incy );
void dgefa ( REAL a[], int lda, int n, int ipvt[], int *info );
void dgesl ( REAL a[], int lda, int n, int ipvt[], REAL b[], int job );
void dmxpy ( int n1, REAL y[], int n2, int ldm, REAL x[], REAL m[] );
void dscal ( int n, REAL da, REAL dx[], int incx );
REAL epslon ( REAL x );
int idamax ( int n, REAL dx[], int incx );
void matgen ( REAL a[], int lda, int n, REAL b[], REAL *norma );
void print_time ( int row );

/*******************************************************************************/

int main ( void )

/*******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINPACK_BENCH.

  Modified:

    07 March 2008
*/
{
	static REAL aa[1000*1000],a[1000*1001],b[1000],x[1000];
	REAL cray,ops,total,norma,normx;
	REAL resid,residn,eps,t1,tm,tm2;
	REAL kf;
	static int ipvt[1000],n,i,ntimes,info,lda,ldaa,mflops;

	lda = 1001;
	ldaa = 1000;
	cray = .056; 
	n = 1000;

  printf(ROLLING);
  printf(PREC);
  printf("Precision Linpack\n\n");

        ops = (2.0e0*(n*n*n))/3.0 + 2.0*(n*n);

        matgen(a,lda,n,b,&norma);

        t1 = cpu_time ( );
        dgefa(a,lda,n,ipvt,&info);
        times[0][0] = cpu_time ( ) - t1;

        t1 = cpu_time ( );
        dgesl(a,lda,n,ipvt,b,0);
        times[1][0] = cpu_time ( ) - t1;

        total = times[0][0] + times[1][0];

/*
  compute a residual to verify results.  
*/ 

        for (i = 0; i < n; i++) {
            	x[i] = b[i];
	}
        matgen(a,lda,n,b,&norma);
        for (i = 0; i < n; i++) {
            	b[i] = -b[i];
	}
        dmxpy(n,b,n,lda,x,a);
        resid = 0.0;
        normx = 0.0;
        for (i = 0; i < n; i++) {
            	resid = (resid > fabs((double)b[i])) 
			? resid : fabs((double)b[i]);
            	normx = (normx > fabs((double)x[i])) 
			? normx : fabs((double)x[i]);
	}
        eps = epslon((REAL)ONE);
        residn = resid/( n*norma*normx*eps );
	
   	printf("     norm. resid      resid           machep");
    printf("         x[0]-1        x[n-1]-1\n");
	printf("  %8.1f      %16.8e%16.8e%16.8e%16.8e\n",
	       (double)residn, (double)resid, (double)eps, 
               (double)x[0]-1, (double)x[n-1]-1);

    printf("\n" );
   	printf("    times are reported for matrices of order %5d\n",n);
	printf("      dgefa      dgesl      total       Mflops     unit");
	printf("      ratio\n");

        times[2][0] = total;
        times[3][0] = ops/(1.0e6*total);
        times[4][0] = 2.0e6/times[3][0];
        times[5][0] = total/cray;

    printf ( "\n" );
   	printf(" times for array with leading dimension of %5d\n",lda);
    printf ( "\n" );

	print_time(0);

        matgen(a,lda,n,b,&norma);
        t1 = cpu_time ( );
        dgefa(a,lda,n,ipvt,&info);
        times[0][1] = cpu_time() - t1;

        t1 = cpu_time ( );
        dgesl(a,lda,n,ipvt,b,0);
        times[1][1] = cpu_time() - t1;

        total = times[0][1] + times[1][1];
        times[2][1] = total;
        times[3][1] = ops/(1.0e6*total);
        times[4][1] = 2.0e6/times[3][1];
        times[5][1] = total/cray;

        matgen(a,lda,n,b,&norma);
        t1 = cpu_time ( );
        dgefa(a,lda,n,ipvt,&info);
        times[0][2] = cpu_time ( ) - t1;
        t1 = cpu_time ( );
        dgesl(a,lda,n,ipvt,b,0);
        times[1][2] = cpu_time ( ) - t1;
        total = times[0][2] + times[1][2];
        times[2][2] = total;
        times[3][2] = ops/(1.0e6*total);
        times[4][2] = 2.0e6/times[3][2];
        times[5][2] = total/cray;

        ntimes = NTIMES;
        tm2 = 0.0;
        t1 = cpu_time ( );

	for (i = 0; i < ntimes; i++) {
            	tm = cpu_time ( );
		matgen(a,lda,n,b,&norma);
		tm2 = tm2 + cpu_time ( ) - tm;
		dgefa(a,lda,n,ipvt,&info);
	}

        times[0][3] = (cpu_time ( ) - t1 - tm2)/ntimes;
        t1 = cpu_time ( );

	for (i = 0; i < ntimes; i++) {
            	dgesl(a,lda,n,ipvt,b,0);
	}

        times[1][3] = (cpu_time ( ) - t1)/ntimes;
        total = times[0][3] + times[1][3];
        times[2][3] = total;
        times[3][3] = ops/(1.0e6*total);
        times[4][3] = 2.0e6/times[3][3];
        times[5][3] = total/cray;

	print_time(1);
	print_time(2);
	print_time(3);

        matgen(aa,ldaa,n,b,&norma);
        t1 = cpu_time ( );
        dgefa(aa,ldaa,n,ipvt,&info);
        times[0][4] = cpu_time ( ) - t1;
        t1 = cpu_time ( );
        dgesl(aa,ldaa,n,ipvt,b,0);
        times[1][4] = cpu_time ( ) - t1;
        total = times[0][4] + times[1][4];
        times[2][4] = total;
        times[3][4] = ops/(1.0e6*total);
        times[4][4] = 2.0e6/times[3][4];
        times[5][4] = total/cray;

        matgen(aa,ldaa,n,b,&norma);
        t1 = cpu_time ( );
        dgefa(aa,ldaa,n,ipvt,&info);
        times[0][5] = cpu_time ( ) - t1;
        t1 = cpu_time ( );
        dgesl(aa,ldaa,n,ipvt,b,0);
        times[1][5] = cpu_time ( ) - t1;
        total = times[0][5] + times[1][5];
        times[2][5] = total;
        times[3][5] = ops/(1.0e6*total);
        times[4][5] = 2.0e6/times[3][5];
        times[5][5] = total/cray;

	matgen(aa,ldaa,n,b,&norma);
	t1 = cpu_time ( );
	dgefa(aa,ldaa,n,ipvt,&info);
	times[0][6] = cpu_time ( ) - t1;
	t1 = cpu_time ( );
	dgesl(aa,ldaa,n,ipvt,b,0);
	times[1][6] = cpu_time ( ) - t1;
	total = times[0][6] + times[1][6];
	times[2][6] = total;
	times[3][6] = ops/(1.0e6*total);
	times[4][6] = 2.0e6/times[3][6];
	times[5][6] = total/cray;

	ntimes = NTIMES;
	tm2 = 0;
	t1 = cpu_time ( );
	for (i = 0; i < ntimes; i++) {
		tm = cpu_time ( );
		matgen(aa,ldaa,n,b,&norma);
		tm2 = tm2 + cpu_time ( ) - tm;
		dgefa(aa,ldaa,n,ipvt,&info);
	}
	times[0][7] = (cpu_time ( ) - t1 - tm2)/ntimes;
	t1 = cpu_time ( );
	for (i = 0; i < ntimes; i++) {
		dgesl(aa,ldaa,n,ipvt,b,0);
	}
	times[1][7] = (cpu_time ( ) - t1)/ntimes;
	total = times[0][7] + times[1][7];
	times[2][7] = total;
	times[3][7] = ops/(1.0e6*total);
	times[4][7] = 2.0e6/times[3][7];
	times[5][7] = total/cray;

	/* the following code sequence implements the semantics of
	   the Fortran intrinsics "nint(min(times[3][3],times[3][7]))"	*/

	kf = (times[3][3] < times[3][7]) ? times[3][3] : times[3][7];
	kf = (kf > ZERO) ? (kf + .5) : (kf - .5);
	if (fabs((double)kf) < ONE) 
		mflops = 0;
	else {
		mflops = floor(fabs((double)kf));
		if (kf < ZERO) mflops = -mflops;
	}

    printf ( "\n" );
	printf(" times for array with leading dimension of%4d\n",ldaa);
    printf ( "\n " );
	print_time(4);
	print_time(5);
	print_time(6);
	print_time(7);
	printf(ROLLING);
    printf(PREC);
	printf(" Precision %5d Mflops ; %d Reps \n",mflops,NTIMES);

  return 0;
}
/*******************************************************************************/

double cpu_time ( void )

/*******************************************************************************/
/*
  Purpose:
 
    CPU_TIME reports the total CPU time for a program.

  Modified:

    27 September 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/*******************************************************************************/

void daxpy ( int n, REAL da, REAL dx[], int incx, REAL dy[], int incy )

/*******************************************************************************/
/*
  Purpose:

    DAXPY computes a constant times a vector plus a vector.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Jack Dongarra

  Reference:

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Algorithm 539: 
    Basic Linear Algebra Subprograms for Fortran Usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
*/
{
	int i,ix,iy,m,mp1;

	if(n <= 0) return;
	if (da == ZERO) return;

	if(incx != 1 || incy != 1) {

		/* code for unequal increments or equal increments
		   not equal to 1 					*/

		ix = 0;
		iy = 0;
		if(incx < 0) ix = (-n+1)*incx;
		if(incy < 0)iy = (-n+1)*incy;
		for (i = 0;i < n; i++) {
			dy[iy] = dy[iy] + da*dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
      		return;
	}

	/* code for both increments equal to 1 */

#ifdef ROLL
	for (i = 0;i < n; i++) {
		dy[i] = dy[i] + da*dx[i];
	}
#endif
#ifdef UNROLL

	m = n % 4;
	if ( m != 0) {
		for (i = 0; i < m; i++) 
			dy[i] = dy[i] + da*dx[i];
		if (n < 4) return;
	}
	for (i = m; i < n; i = i + 4) {
		dy[i] = dy[i] + da*dx[i];
		dy[i+1] = dy[i+1] + da*dx[i+1];
		dy[i+2] = dy[i+2] + da*dx[i+2];
		dy[i+3] = dy[i+3] + da*dx[i+3];
	}
#endif
  return;
}
/*******************************************************************************/

REAL ddot ( int n, REAL dx[], int incx, REAL dy[], int incy )

/*******************************************************************************/
/*
  Purpose:

    DDOT forms the dot product of two vectors.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Jack Dongarra

  Reference:

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Algorithm 539: 
    Basic Linear Algebra Subprograms for Fortran Usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
*/
{
	REAL dtemp;
	int i,ix,iy,m,mp1;

	dtemp = ZERO;

	if(n <= 0) return(ZERO);

	if(incx != 1 || incy != 1) {

		/* code for unequal increments or equal increments
		   not equal to 1					*/

		ix = 0;
		iy = 0;
		if (incx < 0) ix = (-n+1)*incx;
		if (incy < 0) iy = (-n+1)*incy;
		for (i = 0;i < n; i++) {
			dtemp = dtemp + dx[ix]*dy[iy];
			ix = ix + incx;
			iy = iy + incy;
		}
		return(dtemp);
	}

	/* code for both increments equal to 1 */

#ifdef ROLL
	for (i=0;i < n; i++)
		dtemp = dtemp + dx[i]*dy[i];
	return(dtemp);
#endif
#ifdef UNROLL

	m = n % 5;
	if (m != 0) {
		for (i = 0; i < m; i++)
			dtemp = dtemp + dx[i]*dy[i];
		if (n < 5) return(dtemp);
	}
	for (i = m; i < n; i = i + 5) {
		dtemp = dtemp + dx[i]*dy[i] +
		dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] +
		dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
	}
	return(dtemp);
#endif
}
/*******************************************************************************/

void dgefa ( REAL a[], int lda, int n, int ipvt[], int *info )

/*******************************************************************************/
/*
  Purpose:

    DGEFA factors a double precision matrix by gaussian elimination.

  Discussion:

    We would like to declare a[][lda], but c does not allow it.  In this
    function, references to a[i][j] are written a[lda*i+j].

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Cleve Moler.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

     on entry

        a       REAL precision[n][lda]
                the matrix to be factored.

        lda     integer
                the leading dimension of the array  a .

        n       integer
                the order of the matrix  a .

     on return

        a       an upper triangular matrix and the multipliers
                which were used to obtain it.
                the factorization can be written  a = l*u  where
                l  is a product of permutation and unit lower
                triangular matrices and  u  is upper triangular.

        ipvt    integer[n]
                an integer vector of pivot indices.

        info    integer
                = 0  normal value.
                = k  if  u[k][k] .eq. 0.0 .  this is not an error
                     condition for this subroutine, but it does
                     indicate that dgesl or dgedi will divide by zero
                     if called.  use  rcond  in dgeco for a reliable
                     indication of singularity.
*/
{
REAL t;
int j,k,kp1,l,nm1;


/*     gaussian elimination with partial pivoting	*/

	*info = 0;
	nm1 = n - 1;
	if (nm1 >=  0) {
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;

          		/* find l = pivot index	*/

			l = idamax(n-k,&a[lda*k+k],1) + k;
			ipvt[k] = l;

			/* zero pivot implies this column already 
			   triangularized */

			if (a[lda*k+l] != ZERO) {

				/* interchange if necessary */

				if (l != k) {
					t = a[lda*k+l];
					a[lda*k+l] = a[lda*k+k];
					a[lda*k+k] = t; 
				}

				/* compute multipliers */

				t = -ONE/a[lda*k+k];
				dscal(n-(k+1),t,&a[lda*k+k+1],1);

				/* row elimination with column indexing */

				for (j = kp1; j < n; j++) {
					t = a[lda*j+l];
					if (l != k) {
						a[lda*j+l] = a[lda*j+k];
						a[lda*j+k] = t;
					}
					daxpy(n-(k+1),t,&a[lda*k+k+1],1,
					      &a[lda*j+k+1],1);
  				} 
  			}
			else { 
            			*info = k;
			}
		} 
	}
	ipvt[n-1] = n-1;
	if (a[lda*(n-1)+(n-1)] == ZERO) *info = n-1;

  return;
}
/*******************************************************************************/

void dgesl ( REAL a[], int lda, int n, int ipvt[], REAL b[], int job )

/*******************************************************************************/
/* We would like to declare a[][lda], but c does not allow it.  In this
function, references to a[i][j] are written a[lda*i+j].  */

/*
  Purpose:

    DGESL solves A*x=b or A'*x=b, after A has been factored.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Cleve Moler

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

     on entry

        a       double precision[n][lda]
                the output from dgeco or dgefa.

        lda     integer
                the leading dimension of the array  a .

        n       integer
                the order of the matrix  a .

        ipvt    integer[n]
                the pivot vector from dgeco or dgefa.

        b       double precision[n]
                the right hand side vector.

        job     integer
                = 0         to solve  a*x = b ,
                = nonzero   to solve  trans(a)*x = b  where
                            trans(a)  is the transpose.

    on return

        b       the solution vector  x .

     error condition

        a division by zero will occur if the input factor contains a
        zero on the diagonal.  technically this indicates singularity
        but it is often caused by improper arguments or improper
        setting of lda .  it will not occur if the subroutines are
        called correctly and if dgeco has set rcond .gt. 0.0
        or dgefa has set info .eq. 0 .

     to compute  inverse(a) * c  where  c  is a matrix
     with  p  columns
           dgeco(a,lda,n,ipvt,rcond,z)
           if (!rcond is too small){
           	for (j=0,j<p,j++)
              		dgesl(a,lda,n,ipvt,c[j][0],0);
	   }

     linpack. this version dated 08/14/78 .
     cleve moler, university of new mexico, argonne national lab.

     functions

     blas daxpy,ddot
*/
{
/*     internal variables	*/

	REAL t;
	int k,kb,l,nm1;

	nm1 = n - 1;
	if (job == 0) {

		/* job = 0 , solve  a * x = b
		   first solve  l*y = b    	*/

		if (nm1 >= 1) {
			for (k = 0; k < nm1; k++) {
				l = ipvt[k];
				t = b[l];
				if (l != k){ 
					b[l] = b[k];
					b[k] = t;
				}	
				daxpy(n-(k+1),t,&a[lda*k+k+1],1,&b[k+1],1);
			}
		} 

		/* now solve  u*x = y */

		for (kb = 0; kb < n; kb++) {
		    k = n - (kb + 1);
		    b[k] = b[k]/a[lda*k+k];
		    t = -b[k];
		    daxpy(k,t,&a[lda*k+0],1,&b[0],1);
		}
	}
	else { 

		/* job = nonzero, solve  trans(a) * x = b
		   first solve  trans(u)*y = b 			*/

		for (k = 0; k < n; k++) {
			t = ddot(k,&a[lda*k+0],1,&b[0],1);
			b[k] = (b[k] - t)/a[lda*k+k];
		}

		/* now solve trans(l)*x = y	*/

		if (nm1 >= 1) {
			for (kb = 1; kb < nm1; kb++) {
				k = n - (kb+1);
				b[k] = b[k] + ddot(n-(k+1),&a[lda*k+k+1],1,&b[k+1],1);
				l = ipvt[k];
				if (l != k) {
					t = b[l];
					b[l] = b[k];
					b[k] = t;
				}
			}
		}
	}
  return;
}
/*******************************************************************************/

void dmxpy ( int n1, REAL y[], int n2, int ldm, REAL x[], REAL m[] )

/*******************************************************************************/
/*
  Purpose:

    DMXPY computes y = y + M * x.

  Discussion:

    We would like to declare m[][ldm], but c does not allow it.  In this
    function, references to m[i][j] are written m[ldm*i+j].

  Modified:

    07 March 2008

  Parameters:

     n1 integer, number of elements in vector y, and number of rows in
         matrix m

     y double [n1], vector of length n1 to which is added 
         the product m*x

     n2 integer, number of elements in vector x, and number of columns
         in matrix m

     ldm integer, leading dimension of array m

     x double [n2], vector of length n2

     m double [ldm][n2], matrix of n1 rows and n2 columns
*/
{
	int j,i,jmin;
	/* cleanup odd vector */

	j = n2 % 2;
	if (j >= 1) {
		j = j - 1;
		for (i = 0; i < n1; i++) 
            		y[i] = (y[i]) + x[j]*m[ldm*j+i];
	} 

	/* cleanup odd group of two vectors */

	j = n2 % 4;
	if (j >= 2) {
		j = j - 1;
		for (i = 0; i < n1; i++)
            		y[i] = ( (y[i])
                  	       + x[j-1]*m[ldm*(j-1)+i]) + x[j]*m[ldm*j+i];
	} 

	/* cleanup odd group of four vectors */

	j = n2 % 8;
	if (j >= 4) {
		j = j - 1;
		for (i = 0; i < n1; i++)
			y[i] = ((( (y[i])
			       + x[j-3]*m[ldm*(j-3)+i]) 
			       + x[j-2]*m[ldm*(j-2)+i])
			       + x[j-1]*m[ldm*(j-1)+i]) + x[j]*m[ldm*j+i];
	} 

	/* cleanup odd group of eight vectors */

	j = n2 % 16;
	if (j >= 8) {
		j = j - 1;
		for (i = 0; i < n1; i++)
			y[i] = ((((((( (y[i])
			       + x[j-7]*m[ldm*(j-7)+i]) + x[j-6]*m[ldm*(j-6)+i])
		  	       + x[j-5]*m[ldm*(j-5)+i]) + x[j-4]*m[ldm*(j-4)+i])
			       + x[j-3]*m[ldm*(j-3)+i]) + x[j-2]*m[ldm*(j-2)+i])
			       + x[j-1]*m[ldm*(j-1)+i]) + x[j]  *m[ldm*j+i];
	} 
	
	/* main loop - groups of sixteen vectors */

	jmin = (n2%16)+16;
	for (j = jmin-1; j < n2; j = j + 16) {
		for (i = 0; i < n1; i++) 
			y[i] = ((((((((((((((( (y[i])
			       	+ x[j-15]*m[ldm*(j-15)+i]) 
				+ x[j-14]*m[ldm*(j-14)+i])
			        + x[j-13]*m[ldm*(j-13)+i]) 
				+ x[j-12]*m[ldm*(j-12)+i])
			        + x[j-11]*m[ldm*(j-11)+i]) 
				+ x[j-10]*m[ldm*(j-10)+i])
			        + x[j- 9]*m[ldm*(j- 9)+i]) 
				+ x[j- 8]*m[ldm*(j- 8)+i])
			        + x[j- 7]*m[ldm*(j- 7)+i]) 
				+ x[j- 6]*m[ldm*(j- 6)+i])
			        + x[j- 5]*m[ldm*(j- 5)+i]) 
				+ x[j- 4]*m[ldm*(j- 4)+i])
			        + x[j- 3]*m[ldm*(j- 3)+i]) 
				+ x[j- 2]*m[ldm*(j- 2)+i])
			        + x[j- 1]*m[ldm*(j- 1)+i]) 
				+ x[j]   *m[ldm*j+i];
	}
  return;
} 
/*******************************************************************************/

void dscal ( int n, REAL da, REAL dx[], int incx )

/*******************************************************************************/
/*
  Purpose:

    DSCAL scales a vector by a constant.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Jack Dongarra

  Reference:

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Algorithm 539: 
    Basic Linear Algebra Subprograms for Fortran Usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
*/
{
	int i,m,mp1,nincx;

	if(n <= 0)return;
	if(incx != 1) {

		/* code for increment not equal to 1 */

		nincx = n*incx;
		for (i = 0; i < nincx; i = i + incx)
			dx[i] = da*dx[i];
		return;
	}

	/* code for increment equal to 1 */

#ifdef ROLL
	for (i = 0; i < n; i++)
		dx[i] = da*dx[i];
#endif
#ifdef UNROLL

	m = n % 5;
	if (m != 0) {
		for (i = 0; i < m; i++)
			dx[i] = da*dx[i];
		if (n < 5) return;
	}
	for (i = m; i < n; i = i + 5){
		dx[i] = da*dx[i];
		dx[i+1] = da*dx[i+1];
		dx[i+2] = da*dx[i+2];
		dx[i+3] = da*dx[i+3];
		dx[i+4] = da*dx[i+4];
	}
#endif

  return;
}
/*******************************************************************************/

REAL epslon ( REAL x )

/*******************************************************************************/
/*
  Purpose:

    EPSLON estimates the unit roundoff in quantities of size X.

  Discussion:

     this program should function properly on all systems
     satisfying the following two assumptions,
        1.  the base used in representing dfloating point
            numbers is not a power of three.
        2.  the quantity  a  in statement 10 is represented to 
            the accuracy used in dfloating point variables
            that are stored in memory.
     the statement number 10 and the go to 10 are intended to
     force optimizing compilers to generate code satisfying 
     assumption 2.
     under these assumptions, it should be true that,
            a  is not exactly equal to four-thirds,
            b  has a zero for its last bit or digit,
            c  is not exactly equal to one,
            eps  measures the separation of 1.0 from
                 the next larger dfloating point number.
     the developers of eispack would appreciate being informed
     about any systems where these assumptions do not hold.

  Modified:

    07 March 2008
*/
{
	REAL a,b,c,eps;


	a = 4.0e0/3.0e0;
	eps = ZERO;
	while (eps == ZERO) {
		b = a - ONE;
		c = b + b + b;
		eps = fabs((double)(c-ONE));
	}
	return(eps*fabs((double)x));
}
/*******************************************************************************/

int idamax ( int n, REAL dx[], int incx )

/*******************************************************************************/
/*
  Purpose:

    IDAMAX indexes the vector element of maximum absolute value.

  Modified:

    07 March 2008

  Author:

    FORTRAN77 original version by Jack Dongarra

  Reference:

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Algorithm 539: 
    Basic Linear Algebra Subprograms for Fortran Usage,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
*/
{
	REAL dmax;
	int i, ix, itemp;

	if( n < 1 ) return(-1);
	if(n ==1 ) return(0);
	if(incx != 1) {

		/* code for increment not equal to 1 */

		ix = 1;
		dmax = fabs((double)dx[0]);
		ix = ix + incx;
		for (i = 1; i < n; i++) {
			if(fabs((double)dx[ix]) > dmax)  {
				itemp = i;
				dmax = fabs((double)dx[ix]);
			}
			ix = ix + incx;
		}
	}
	else {

		/* code for increment equal to 1 */

		itemp = 0;
		dmax = fabs((double)dx[0]);
		for (i = 1; i < n; i++) {
			if(fabs((double)dx[i]) > dmax) {
				itemp = i;
				dmax = fabs((double)dx[i]);
			}
		}
	}
	return (itemp);
}
/*******************************************************************************/

void matgen ( REAL a[], int lda, int n, REAL b[], REAL *norma )

/*******************************************************************************/
/* 
  Purpose:

    MATGEN generates a "random" matrix for testing.

  Discussion:

    We would like to declare a[][lda], but c does not allow it.  In this
    function, references to a[i][j] are written a[lda*i+j].  

  Modified:

    07 March 2008

*/
{
  int i;
  int init;
  int j;

  init = 1325;
  *norma = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      init = 3125 * init % 65536;
      a[lda*j+i] = ( init - 32768.0 ) / 16384.0;
      if ( *norma < a[lda*j+i] )
      {
        *norma = a[lda*j+i];
      }
    }
  }
  for ( i = 0; i < n; i++ ) 
  {
    b[i] = 0.0;
  }
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = b[i] + a[lda*j+i];
    }
  }
  return;
}
/*******************************************************************************/

void print_time ( int row )

/*******************************************************************************/
/*
  Purpose:

    PRINT_TIME prints a row of the time table.

  Modified:

    07 March 2008
*/
{
  printf("%11.2f%11.2f%11.2f%11.0f%11.2f%11.2f\n",   (double)times[0][row],
        (double)times[1][row], (double)times[2][row], (double)times[3][row], 
        (double)times[4][row], (double)times[5][row]);

  return;
}
