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

/*
 * This library creates an interface from either a Fortran implementarion of Blas or C implementation of CBlas
 * to function calls similar to ACML (AMD Math Core LIbrary). Fortran uses "pass-by-reference" convention,
 * for example, if we want to call dnrm2 from within C, we need to add a dummy variable for incx:
 * 
 * int incx = 1;
 * double norm = dnrm2_( N, x, &incx );
 * 
 * CBlas adds the proper C convention and allows for calls of the form:
 * double norm = cblas_dnrm2( N, x, 1 );
 * however, the naming conventino requires cblas_ to be added before all names.
 * 
 * ACML creates a proper C interface, where we can call funcitons of the form:
 * double norm = dnrm2( N, x, 1 );
 * This is the most natural way to do things in C. Since ACML isn't always available, this interface allows
 * to use the ACML function calls and link to wither Fortran Blas or CBlas
 * 
 * By defualt, this assumes that we are using Fortran Blas. If you want to use CBlas, then add
 * -DFACML_CBLAS to the compiler variables. Note that I have seen buggy implementations of CBlas and hence
 * I prefer to link to either ACML or Fortran Blas.
 * 
 * Also note that SuiteSparse requires Fortran Blas, thus if you plan on using with together with SMPACK and
 * SuiteSparse, then you should probably use Fortran Blas here as well (but you don't have to).
 * 
 * Most of my programs use only a small number of Blas functions, hence this creates only a limited interface.
 * 
 */

#ifdef PDBLAS_CBLAS
	#define ddot cblas_ddot
	#define dnrm2 cblas_dnrm2
	#define dcopy cblas_dcopy
	#define daxpy cblas_daxpy
	#define dscal cblas_dscal
	
	#define dtrsv cblas_dtrsv
	
	#define dgemv cblas_dgemv
	
	double ddoti( int N, double *a, int *indx, double *y ); // the sparse function is not included in regular BLAS
#else

	// BLAS Level 1
	double ddot( int N, double *a, int inca, double *b, int incb );
	double ddoti( int N, double *a, int *indx, double *y ); // the sparse function is not included in regular BLAS
	double dnrm2( int N, double *a, int inca );
	void dcopy(  int N, double *a, int inca, double *b, int incb );
	void daxpy( int N, double alpha, double *a, int inca, double *b, int incb );
	void dscal( int N, double alpha, double *a, int inca );

	// BLAS Level 2
	void dtrsv( char uplo, char trans, char diag, int N, double *A, int lda, double *x, int incx );

	// BLAS Level 3
	void dgemv( char trans, int m, int n, double alpha, double *A, int lda, double *x, int incx, double beta, double *y, int incy );

		// Fortran adds underscore _ after the function definitions
		// Fortran BLAS Level 1
		double ddot_( const int *N, const double *a, const int *inca, const double *b, const int *incb );
		double dnrm2_( const int *N, const double *a, const int *inca );
		void dcopy_(  const int *N, const double *a, const int *inca, const double *b, const int *incb );
		void daxpy_( const int *N, const double *alpha, const double *a, const int *inca, const double *b, const int *incb );
		void dscal_( const int *N, const double *alpha, const double *a, const int *inca );

		// Fortran BLAS Level 2
		void dtrsv_( const char *uplo, const char *trans, const char *diag, const int *N, const double *A, const int *lda, const double *x, const int *incx );

		// Fortran BLAS Level 3
		void dgemv_( const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *lda, const double *x, const int *incx, const double *beta, const double *y, const int *incy );


#endif

	
