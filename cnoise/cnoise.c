#ifndef __CNOISE_C
#define __CNOISE_C

/*
* This code isdistrubuted under GPLv3
* 
* This code generates corelated or colored noise with 1/f^alpha power spectrum distribution. 
* It uses the algorithm by:
* 
* Jeremy Kasdin
* Discrete Simulation of Colored Noise and Stochastic Processes and $ 1/f^\alpha $ Power Law Noise Generation
* Proceedings of the IEEE
* Volume 83, Number 5, 1995, pages 802-827.
* 
* This code uses GSL fast Fourier transform gsl_fft_complex_forward(...) and the GCC rand() functions
* 
* Code Author: Miroslav Stoyanov, Jan 2011
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
* 
* updated June 2011: fixed a small bug causing "nan" to be returned sometimes
* 
*/

#include "cnoise.h"

#define _CNOISE_PI 3.14159265358979323846 // pi

inline double cnoise_uniform01(){
	return ( ( (double) rand() + 1.0 ) / ( (double) RAND_MAX + 1.0 ) );
};

inline void cnoise_two_gaussian_truncated( double *a, double *b, double std, double range ){
	double gauss_u = cnoise_uniform01();
	double gauss_v = cnoise_uniform01();
	*a = std * sqrt( - 2 * log( gauss_u ) ) * cos( 2 * _CNOISE_PI * gauss_v );
	*b = std * sqrt( - 2 * log( gauss_u ) ) * sin( 2 * _CNOISE_PI * gauss_v );
	while ( (*a < -range) || (*a > range) ){
		gauss_u = cnoise_uniform01(); gauss_v = cnoise_uniform01();
		*a = std * sqrt( - 2 * log( gauss_u ) ) * cos( 2 * _CNOISE_PI * gauss_v );
	};
	while ( (*b < -range) || (*b > range) ){
		gauss_u = cnoise_uniform01(); gauss_v = cnoise_uniform01();
		*b = std * sqrt( - 2 * log( gauss_u ) ) * cos( 2 * _CNOISE_PI * gauss_v );
	};
};

void cnoise_generate_colored_noise( double *x, int N, double alpha, double std ){
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	
	int i, j;
	double tmp;
	double gauss_second, gauss_u, gauss_v;
	double * fh = malloc( 4*N*sizeof(double) );
	double * fw = malloc( 4*N*sizeof(double) );
	
	fh[0] = 1.0; fh[1] = 0.0;
	for( i=1; i<N; i++ ){
		fh[2*i] = fh[2*(i-1)] * ( 0.5 * alpha + (double)(i-1) ) / ((double) i );
		fh[2*i+1] = 0.0;
	};
	for( i=0; i<N; i+=2 ){
		gauss_u = cnoise_uniform01(); gauss_v = cnoise_uniform01();
		fw[2*i] = std * sqrt( - 2 * log( gauss_u ) ) * cos( 2 * _CNOISE_PI * gauss_v ); fw[2*i+1] = 0.0;
		fw[2*i+2] = std * sqrt( - 2 * log( gauss_u ) ) * sin( 2 * _CNOISE_PI * gauss_v ); fw[2*i+3] = 0.0;
	};
	for( i=2*N; i<4*N; i++ ){ fh[i] = 0.0; fw[i] = 0.0; };
	
	wavetable = gsl_fft_complex_wavetable_alloc(2*N);
	workspace = gsl_fft_complex_workspace_alloc(2*N);

	gsl_fft_complex_forward (fh, 1, 2*N, wavetable, workspace);
	gsl_fft_complex_forward (fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N+1; i++ ){
		tmp = fw[2*i];
		fw[2*i]   = tmp*fh[2*i] - fw[2*i+1]*fh[2*i+1];
		fw[2*i+1] = tmp*fh[2*i+1] + fw[2*i+1]*fh[2*i];
	};

	fw[0] /= 2.0; fw[1] /= 2.0;
	fw[2*N] /= 2.0; fw[2*N+1] /= 2.0;
	
	for( i=N+1; i<2*N; i++ ){
		fw[2*i] = 0.0; fw[2*i+1] = 0.0;
	};
	
	gsl_fft_complex_inverse( fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N; i++ ){
		x[i] = 2.0*fw[2*i];
	};

	free( fh );
	free( fw );

	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
};

void cnoise_generate_colored_noise_uniform( double *x, int N, double alpha, double range ){
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	
	int i, j;
	double tmp;
	double gauss_second, gauss_u, gauss_v;
	double * fh = malloc( 4*N*sizeof(double) );
	double * fw = malloc( 4*N*sizeof(double) );
	
	fh[0] = 1.0; fh[1] = 0.0;
	fw[0] = 2.0*range*cnoise_uniform01()-range; fw[1] = 0.0;
	for( i=1; i<N; i++ ){
		fh[2*i] = fh[2*(i-1)] * ( 0.5 * alpha + (double)(i-1) ) / ((double) i );
		fh[2*i+1] = 0.0;
		fw[2*i] = 2.0*range*cnoise_uniform01()-range; fw[2*i+1] = 0.0;
	};
	for( i=2*N; i<4*N; i++ ){ fh[i] = 0.0; fw[i] = 0.0; };
	
	wavetable = gsl_fft_complex_wavetable_alloc(2*N);
	workspace = gsl_fft_complex_workspace_alloc(2*N);

	gsl_fft_complex_forward (fh, 1, 2*N, wavetable, workspace);
	gsl_fft_complex_forward (fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N+1; i++ ){
		tmp = fw[2*i];
		fw[2*i]   = tmp*fh[2*i] - fw[2*i+1]*fh[2*i+1];
		fw[2*i+1] = tmp*fh[2*i+1] + fw[2*i+1]*fh[2*i];
	};

	fw[0] /= 2.0; fw[1] /= 2.0;
	fw[2*N] /= 2.0; fw[2*N+1] /= 2.0;
	
	for( i=N+1; i<2*N; i++ ){
		fw[2*i] = 0.0; fw[2*i+1] = 0.0;
	};
	
	gsl_fft_complex_inverse( fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N; i++ ){
		x[i] = 2.0*fw[2*i];
	};

	free( fh );
	free( fw );

	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
};

void cnoise_generate_colored_noise_truncated( double *x, int N, double alpha, double std, double range ){
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	
	int i, j;
	double tmp;
	double gauss_second, gauss_u, gauss_v;
	double * fh = malloc( 4*N*sizeof(double) );
	double * fw = malloc( 4*N*sizeof(double) );
	
	fh[0] = 1.0; fh[1] = 0.0;
	for( i=1; i<N; i++ ){
		fh[2*i] = fh[2*(i-1)] * ( 0.5 * alpha + (double)(i-1) ) / ((double) i );
		fh[2*i+1] = 0.0;
	};
	for( i=0; i<N; i+=2 ){
		cnoise_two_gaussian_truncated( &fw[2*i], &fw[2*i+2], std, range );
	};
	for( i=2*N; i<4*N; i++ ){ fh[i] = 0.0; fw[i] = 0.0; };
	
	wavetable = gsl_fft_complex_wavetable_alloc(2*N);
	workspace = gsl_fft_complex_workspace_alloc(2*N);

	gsl_fft_complex_forward (fh, 1, 2*N, wavetable, workspace);
	gsl_fft_complex_forward (fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N+1; i++ ){
		tmp = fw[2*i];
		fw[2*i]   = tmp*fh[2*i] - fw[2*i+1]*fh[2*i+1];
		fw[2*i+1] = tmp*fh[2*i+1] + fw[2*i+1]*fh[2*i];
	};

	fw[0] /= 2.0; fw[1] /= 2.0;
	fw[2*N] /= 2.0; fw[2*N+1] /= 2.0;
	
	for( i=N+1; i<2*N; i++ ){
		fw[2*i] = 0.0; fw[2*i+1] = 0.0;
	};
	
	gsl_fft_complex_inverse( fw, 1, 2*N, wavetable, workspace);
	
	for( i=0; i<N; i++ ){
		x[i] = 2.0*fw[2*i];
	};

	free( fh );
	free( fw );

	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
};

#endif
