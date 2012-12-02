/*
* This code is distributed under GPLv3
*
* This code generates correlated or colored noise with 1/f^alpha power
* spectrum distribution.
*
* It uses the algorithm by:
*
* Jeremy Kasdin
* Discrete Simulation of Colored Noise and Stochastic Processes and
* $ 1/f^\alpha $ Power Law Noise Generation
* Proceedings of the IEEE
* Volume 83, Number 5, 1995, pages 802-827.
*
* This code uses GSL fast Fourier transform gsl_fft_complex_forward(...)
* and the GCC rand() functions
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

Installation:

Edit Config.mk and setup the variables for your compiler, flags
and the path to your GSL installation.

Afterwards, use

make

You should also go to
cd examples/
make
./test

to make sure everything is working properly. See test.c for details.
