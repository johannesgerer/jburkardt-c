#ifndef __CNOISE_H
#define __CNOISE_H

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
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


void cnoise_generate_colored_noise( double *x, int N, double alpha, double std );
// generates a vector x of size N with a 1/f^alpha frequency distribution
// std is the standard deviation of the underlying gaussian distribution
// Variables: double *x - Input: allocated vector of size N
//                        Output: a realization of 1/f^alpha noise vector of size N, using an undelying Gaussian (0,std)
//             int N    - Input: the size of the vector
//            double alpha - Input: the decay rate for the power spectrum (i.e. 1/f^alpha
//            double std - Input: the standard deviation of the underlying Gaussian distribution
// NOTE: you should call srand( seed ) before you call dcnoise_generate_colored_noise

void cnoise_generate_colored_noise_uniform( double *x, int N, double alpha, double range );
// same as cnoise_generate_colored_noise_uniform, except the white noise vector comes from
// uniform distribution on (-range,+range)

void cnoise_generate_colored_noise_truncated( double *x, int N, double alpha, double std, double range );
// same as cnoise_generate_colored_noise_uniform, except the white noise vector comes from
// truncated Gaussian distribution with mean 0, standard deviation std and truncated to (-range,range)

#endif
