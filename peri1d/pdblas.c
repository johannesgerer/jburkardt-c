/*
* This code isdistrubuted under GPLv3
* 
* Code Author: Miroslav Stoyanov, Nov 2011
* 
* Copyright (C) 2011  Miroslav Stoyanov
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* Since the GNU General Public License is longer than this entire code, 
* a copy of it can be obtained separately at <http://www.gnu.org/licenses/>
*/

#include "pdblas.h"

#ifdef PDBLAS_CBLAS
	// CBLAS doesn't require defining any functions, everythign is done with macros in facml.h

#else
// BLAS Level 1
double ddot( int N, double *a, int inca, double *b, int incb ){
	return ddot_( &N, a, &inca, b, &incb );
};
double dnrm2( int N, double *a, int inca ){
	return dnrm2_( &N, a, &inca );
};
void dcopy(  int N, double *a, int inca, double *b, int incb ){
	dcopy_( &N, a, &inca, b, &incb );
};
void daxpy( int N, double alpha, double *a, int inca, double *b, int incb ){
	daxpy_( &N, &alpha, a, &inca, b, &incb );
};
void dscal( int N, double alpha, double *a, int inca ){
	dscal_( &N, &alpha, a, &inca );
};

// BLAS Level 2
void dtrsv( char uplo, char trans, char diag, int N, double *A, int lda, double *x, int incx ){
    dtrsv_( &uplo, &trans, &diag, &N, A, &lda, x, &incx );
};

// BLAS Level 3
void dgemv( char trans, int m, int n, double alpha, double *A, int lda, double *x, int incx, double beta, double *y, int incy ){
	dgemv_( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
};


#endif

double ddoti( int N, double *a, int *indx, double *y ){
	int i;
	double count = 0;
	for( i=0; i<N; i++ ){ count += a[i] * y[indx[i]-1]; };
	return count;
};
