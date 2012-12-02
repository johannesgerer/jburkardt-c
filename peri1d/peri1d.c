/*
* This code is distrubuted under GPLv3
* 
* Code Author: Miroslav Stoyanov, Jan 2012
* 
* Copyright (C) 2012 Miroslav Stoyanov
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

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include "pdblas.h"


const int N = 10;
const double dt = 1.0/ ( 100*1024.0 );
const double T = 500;


// for output



// ----------------------------

//#define __KERNEL_SMOOTH // the kernel is smoothed by exp( 1 - delta / ( delta - dist) )
#define __KERNEL_JUMP   // the kernel jump vanishes at dist > delta

//#define __USE_OPENMP

#define _SMPACK_OMP_READ_SPLIT( me, nthreads, N, start, end, length ) \
	(me) = omp_get_thread_num(); \
	(nthreads) = omp_get_num_threads(); \
	(start) = (me) * (N) / (nthreads); \
	(end) = (me+1) * (N) / (nthreads) - 1; \
	(length) = (end) - (start) + 1; 

struct elem1d
{
  double xl, xr;
  int il, ir; // global index
  int iul, iur; // index of the unknown
  double ul, ur; // values for the solution
};

struct mesh1d
{
  int ne, nu; // number of elements
  struct elem1d *el; // elements
};

int main ( int argc, char** argv );
double bnd_u ( double x );
void external_force( struct mesh1d *mesh, double *B, double (*b)(double) );
void force_action_quad3( struct mesh1d *mesh, double *F, double delta, 
  double (*kernel)(double,double,double,double,double) );
double force_b ( double x );
void generateMeshDL( struct mesh1d *mesh, int N, double delta, int bL, int bR );
double int_over_interval( double delta, 
  double (*kernel)(double,double,double,double,double), double ul, double ur, 
  double xl, double xr, double u, double x );
void invert_massDL( struct mesh1d *mesh, double *u );
void load_uvalues( struct mesh1d *mesh, double *u,  double (*bnd_u)(double) );
void read_solution( struct mesh1d *mesh, const char * filename );
double smooth_kernel ( double delta, double up, double u, double xp, double x );
void timestamp ( void );
void write_solution( struct mesh1d *mesh, const char * filename );
void wwrite_vector( double *x, int n, const char * filename );

/******************************************************************************/

int main ( int argc, char** argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for the 1D peridynamics solver.
*/
{
  char buffer[50];
  const double delta = 0.005;
  int i;
  int itr = 0;
  struct mesh1d mesh;
  int nu;
  double out_time;
  double t = 0.0;
  const double Tint = 10;

  timestamp ( );
  printf ( "\n" );
  printf ( "PERI1D:\n" );
  printf ( "  Solve a 1D peridynamics problem\n" );
  printf ( "\n" );
  printf ( "  Initial time t = %g\n", t );
  printf ( "  Output interval Tint = %g\n", Tint );
  printf ( "  Final time T = %g\n", T );
  printf ( "  Number of elements N = %d\n", N );
  printf ( "  Horizon, delta = %g\n", delta );

  generateMeshDL( &mesh, N, delta, 1, 1 );
	
  nu = mesh.nu;

  printf ( "  MESH.NU = %d\n", mesh.nu );
  printf ( "  MESH.NE = %d\n", mesh.ne );

  double *u = malloc( mesh.nu * sizeof(double) );
  double *u_old = malloc( mesh.nu * sizeof(double) );
  double *u_new = malloc( mesh.nu * sizeof(double) );
  double *b = malloc( mesh.nu * sizeof(double) );
  double *F = malloc( mesh.nu * sizeof(double) );
  double *tmp;
	
  for ( i = 0; i < nu; i++ )
  {
    u[i] = 0.0;
    u_old[i] = 0.0;
  }
//
// Removed the integration part, there is no need for that now
//
  out_time = t + Tint;

  while ( t < T )
  {
    load_uvalues( &mesh, u, &bnd_u );

    force_action_quad3( &mesh, F, delta, &smooth_kernel );	
    external_force( &mesh, b, &force_b );
//
//  F = F + b.
//
    daxpy ( nu, 1.0, b, 1, F, 1 );
//
//  Solve M * F, storing result in F;
//
    invert_massDL( &mesh, F );
//
//  Save current solution U as U_NEW.
//
    dcopy ( nu, u, 1, u_new, 1 );
    dscal ( nu, 2.0, u, 1 );
    daxpy ( nu, dt*dt, F, 1, u, 1 );
    daxpy ( nu, -1.0, u_old, 1, u, 1 );
//
//  Interchange U_OLD and U_NEW.
//	
    tmp = u_old;
    u_old = u_new;
    u_new = tmp;
    tmp = NULL;
				
    t = t + dt;
//
//  Occasionally, write a solution to a file.
//
    if ( out_time <= t )
    {
      printf ( "  Saving solution at time %g\n", t );
      out_time = t + Tint;
      load_uvalues( &mesh, u, &bnd_u );
      itr++;
      sprintf ( buffer, "outfile%d", itr );
      write_solution( &mesh, buffer );
    };

  };
	
  load_uvalues( &mesh, u, &bnd_u );
  write_solution( &mesh, "outfile" );
/*
  Free memory.
*/
  free ( u );
  free ( u_old );
  free ( u_new );
  free ( b );
  free ( F );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PERI1D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
};
/******************************************************************************/

double bnd_u ( double x )

/******************************************************************************/
{
  return 0;
};
/******************************************************************************/

void external_force( struct mesh1d *mesh, double *B, double (*b)(double) )

/******************************************************************************/
/*
  Purpose:

    EXTERNAL_FORCE computes the body force acting on the material.
*/
{
	int i;
	double q3x = sqrt(15.0)/5.0, w3x = 5.0/9.0, w3z = 8.0/9.0; // quadrature
	double integl, integr, val;
	double xl, xr;
	
#ifdef __USE_OPENMP
	int start, end, length, me, nthreads;
	#pragma omp parallel private(start, end, length, me, nthreads,i,integl,integr,val,xl,xr)
	{
		_SMPACK_OMP_READ_SPLIT( me, nthreads, mesh ->ne, start, end, length );
		for ( i = start; i <= end; i++ ){
			xl = mesh -> el[i].xl; xr = mesh -> el[i].xr;
			val = b( -q3x*(xr-xl)/2 + (xr+xl)/2 );
			integl  = w3x * ( 0.5 + q3x/2 ) * val;
			integr  = w3x * ( 0.5 - q3x/2 ) * val;
			
			val = w3z * 0.5 * b( (xr+xl)/2 );
			integl += val;
			integr += val;
			
			val = b( q3x*(xr-xl)/2 + (xr+xl)/2 );
			integl += w3x * ( 0.5 - q3x/2 ) * val;
			integr += w3x * ( 0.5 + q3x/2 ) * val;
			
			if ( mesh ->el[i].iul != -1 ){ // if we are using this test function
				B[mesh ->el[i].iul] = 0.5 * ( xr - xl ) * integl;
			};
			if ( mesh ->el[i].iur != -1 ){ // if we are using this test function
				B[mesh ->el[i].iur] = 0.5 * ( xr - xl ) * integr;
			};
		};
	};
#else		
	for ( i = 0; i < mesh ->ne; i++ ){
		xl = mesh -> el[i].xl; xr = mesh -> el[i].xr;
		val = b( -q3x*(xr-xl)/2 + (xr+xl)/2 );
		integl  = w3x * ( 0.5 + q3x/2 ) * val;
		integr  = w3x * ( 0.5 - q3x/2 ) * val;
		
		val = w3z * 0.5 * b( (xr+xl)/2 );
		integl += val;
		integr += val;
		
		val = b( q3x*(xr-xl)/2 + (xr+xl)/2 );
		integl += w3x * ( 0.5 - q3x/2 ) * val;
		integr += w3x * ( 0.5 + q3x/2 ) * val;
		
		if ( mesh ->el[i].iul != -1 ){ // if we are using this test function
			B[mesh ->el[i].iul] = 0.5 * ( xr - xl ) * integl;
		};
		if ( mesh ->el[i].iur != -1 ){ // if we are using this test function
			B[mesh ->el[i].iur] = 0.5 * ( xr - xl ) * integr;
		};
	};
#endif
};
/******************************************************************************/

void force_action_quad3( struct mesh1d *mesh, double *F, double delta, 
  double (*kernel)(double,double,double,double,double) )

/******************************************************************************/
/*
  Purpose:

    FORCE_ACTION_QUAD3 computes the forces and stores them in F.

  Discussion:

    We need the kernel density "kernel" and the distance of interaction "delta" 
    so that we don't integrate over unrelated elements.

    Note that kernel takes input of the form kernel( delta, u', u, x', x );
*/
{
	int i, j, e;
	//double q3x = sqrt(15.0)/5.0, w3x = 5.0/9.0, w3z = 8.0/9.0;
	double xl, xr, ul, ur;
	double xlj, xrj, ulj, urj;
	//double qu, qx; // quad points
	double integl, integr, val;

	for( e=0; e < mesh ->ne; e++ ){
		xl = mesh ->el[e].xl; xr = mesh ->el[e].xr;
		ul = mesh ->el[e].ul; ur = mesh ->el[e].ur;

		// int for the right test function
		if ( e + 1 < mesh ->ne ){
			xlj = mesh ->el[e+1].xl; xrj = mesh ->el[e+1].xr;
			ulj = mesh ->el[e+1].ul; urj = mesh ->el[e+1].ur;
			val = int_over_interval( delta, kernel, ulj, ulj + delta*(urj-ulj)/(xrj-xlj), xlj, xlj+delta, ur, xr );
		}else{
			val = 0.0;
		};
		val += int_over_interval( delta, kernel, ur - delta*(ur-ul)/(xr-xl), ur, xr-delta, xr, ur, xr );
		
		integr = val;
		
		// int for the left test function
		if ( e > 0 ){
			xlj = mesh ->el[e-1].xl; xrj = mesh ->el[e-1].xr;
			ulj = mesh ->el[e-1].ul; urj = mesh ->el[e-1].ur;
			val = int_over_interval( delta, kernel, urj - delta*(urj-ulj)/(xrj-xlj), urj, xrj-delta, xrj, ul, xl );
		}else{
			val = 0.0;
		};
		val += int_over_interval( delta, kernel, ul, ul + delta*(ur-ul)/(xr-xl), xl, xl+delta, ul, xl );
		
		integl = val;
		
		// int for the mid-point
		val  = int_over_interval( delta, kernel, 0.5*(ul+ur), 0.5*(ul+ur) + delta*(ur-ul)/(xr-xl), 0.5*(xl+xr), 0.5*(xl+xr)+delta, 0.5*(ul+ur), 0.5*(xl+xr) );
		val += int_over_interval( delta, kernel, 0.5*(ul+ur) - delta*(ur-ul)/(xr-xl), 0.5*(ul+ur), 0.5*(xl+xr)-delta, 0.5*(xl+xr), 0.5*(ul+ur), 0.5*(xl+xr) );
		
		integr += 2.0 * val;
		integl += 2.0 * val;
		
		integr *= (xr-xl)/6.0;
		integl *= (xr-xl)/6.0;
		
		if ( mesh ->el[e].iul != -1 ){ // if we are using this test function
			F[mesh ->el[e].iul] = integl;
		};
		if ( mesh ->el[e].iur != -1 ){ // if we are using this test function
			F[mesh ->el[e].iur] = integr;
		};
	};
};
/******************************************************************************/

double force_b ( double x )

/******************************************************************************/
/*
  Purpose:

    FORCE_B ...
*/
{
  double dx = 1.0 / N;
  if ( x < dx+1.E-7 )
  {
    return -1.0/dx;
  }
  else
  {
    return 0.0;
  };
};
/******************************************************************************/

void generateMeshDL( struct mesh1d *mesh, int N, double delta, int bL, int bR )

/******************************************************************************/
/*
  Purpose:

    GENERATEMESHDL generates the mesh and related quantities.

  Discussion:

    This function computes the mesh, based on the approximate number of 
    elements, the value of horizon, and the left/right boundary indicators 
    (1 boundary, 0 no boundary).
*/
{
	int i, count, indx, indxu;
	double dx = 1.0 / N;
	double x, x_max;
	// do a one pass to get the number of elements
	if ( bL == 1 ){ // left boundary
		count = 1; // expand to account for the boundary region
		x = -delta;
	}else{
		count = 0;
		x = 0;
	};
	if ( bR == 1 ){ // right boundary
		count++; // add one element on the right
		x_max = 1.0 + delta;
	}else{
		x_max = 1.0;
	};
	while( x < x_max ){
		count ++;
		x = x + dx;
		if ( x > x_max ){
			x = x_max;
		}else if ( x + dx/2 > x_max ){
			x = x_max;
		};
	};
	//printf(" Number of elements: %d \n",count);
	mesh ->ne = count;
	mesh ->el = malloc( count * sizeof( struct elem1d ) );
	// set the number of elements
	indx = 0; indxu = 0;
	if ( bL == 1 ){ // left boundary
		count = 1;
		x = 0.0;
		mesh ->el[0].xl = -delta;
		mesh ->el[0].xr = 0.0;
		mesh ->el[0].il = indx;
		indx++;
		mesh ->el[0].ir = indx;
		indx++;
		mesh ->el[0].iul = -1;
		mesh ->el[0].iur = -1;
	}else{
		count = 0;
		x = 0;
	};
	if ( bR == 1 ){ // right boundary
		// one element on the right
		x_max = 1.0 + delta;
	}else{
		x_max = 1.0;
	};
	while( x < x_max ){
		mesh ->el[count].xl = x;
		x = x + dx;
		if ( x > x_max ){
			x = x_max;
		};
		mesh ->el[count].xr = x;
		mesh ->el[count].il = indx;
		indx++;
		mesh ->el[count].ir = indx;
		indx++;
		mesh ->el[count].iul = indxu;
		indxu++;
		mesh ->el[count].iur = indxu;
		indxu++;
		count++;
	};
	if ( bR == 1 ){ // right boundary
		// one element on the right
		mesh ->el[count-1].xl = 1.0;
		mesh ->el[count-1].xr = 1.0 + delta;
		mesh ->el[count-1].il = indx;
		indx++;
		mesh ->el[count-1].ir = indx;
		indx++;
		mesh ->el[count-1].iul = -1;
		mesh ->el[count-1].iur = -1;
	}
	mesh ->nu = indxu;
	/*for( i = 0; i<mesh->ne; i++ ){
		printf(" Element %d  (xl,xr) %f %f   (il,ir) %d %d   (iul,iur) %d %d\n",i,mesh->el[i].xl,mesh->el[i].xr,mesh->el[i].il,mesh->el[i].ir,mesh->el[i].iul,mesh->el[i].iur);
	};*/
};
/******************************************************************************/

double int_over_interval( double delta, 
  double (*kernel)(double,double,double,double,double), double ul, double ur, 
  double xl, double xr, double u, double x )

/******************************************************************************/
/*
  Purpose:

    INT_OVER_INTERVAL integrates the kernel over (xl, xr) assuming we are over one element.

  Discussion:

    The 4 point quadrature rule is used.
*/
{
  double intg;	
  double q4x1 = sqrt( (3.0 - 2.0 * sqrt( 6.0/5.0 )) / 7.0 ), q4x2 = sqrt( (3.0 + 2.0 * sqrt( 6.0/5.0 )) / 7.0 );
  double w4x1 = ( 18 + sqrt(30) )/36, w4x2 = ( 18 - sqrt(30) )/36;

  intg  = w4x1*kernel( delta, -q4x1*(ur-ul)/2 + (ur+ul)/2, u, -q4x1*(xr-xl)/2 + (xr+xl)/2, x );
  intg += w4x2*kernel( delta, -q4x2*(ur-ul)/2 + (ur+ul)/2, u, -q4x2*(xr-xl)/2 + (xr+xl)/2, x );
  intg += w4x2*kernel( delta,  q4x2*(ur-ul)/2 + (ur+ul)/2, u,  q4x2*(xr-xl)/2 + (xr+xl)/2, x );
  intg += w4x1*kernel( delta,  q4x1*(ur-ul)/2 + (ur+ul)/2, u,  q4x1*(xr-xl)/2 + (xr+xl)/2, x );

  return 0.5 * (xr - xl) * intg;
};
/******************************************************************************/

void invert_massDL( struct mesh1d *mesh, double *u )

/******************************************************************************/
/*
  Purpose:

    INVERT_MASSDL inverts the mass matrix on U.
*/
{
	int i;
	double dx, utmp;
#ifdef __USE_OPENMP
	int start, end, length, me, nthreads;
	#pragma omp parallel private(start, end, length, me, nthreads, i, dx, utmp)
	{
		_SMPACK_OMP_READ_SPLIT( me, nthreads, mesh ->ne, start, end, length );
		for( i=start; i <= end; i++ ){
			dx = mesh ->el[i].xr - mesh ->el[i].xl;
			if (mesh -> el[i].iul != -1){
				if ( mesh -> el[i].iur != -1 ){
					utmp = u[ mesh ->el[i].iul ];
					u[ mesh ->el[i].iul ] = (1/dx) * (4 * utmp -2* u[ mesh ->el[i].iur ]);
					u[ mesh ->el[i].iur ] = (1/dx) * (-2* utmp +4* u[ mesh ->el[i].iur ]);
				}else{
					u[ mesh ->el[i].iul ] = (3/dx) * u[ mesh ->el[i].iul ];
				};
			}else{
				if ( mesh -> el[i].iur != -1 ){
					u[ mesh ->el[i].iur ] = (3/dx) * u[ mesh ->el[i].iur ];
				};
			};
		};
	};
#else
	for( i=0; i < mesh->ne; i++ ){
		dx = mesh ->el[i].xr - mesh ->el[i].xl;
		if (mesh -> el[i].iul != -1){
			if ( mesh -> el[i].iur != -1 ){
				utmp = u[ mesh ->el[i].iul ];
				u[ mesh ->el[i].iul ] = (1/dx) * (4 * utmp -2* u[ mesh ->el[i].iur ]);
				u[ mesh ->el[i].iur ] = (1/dx) * (-2* utmp +4* u[ mesh ->el[i].iur ]);
			}else{
				u[ mesh ->el[i].iul ] = (3/dx) * u[ mesh ->el[i].iul ];
			};
		}else{
			if ( mesh -> el[i].iur != -1 ){
				u[ mesh ->el[i].iur ] = (3/dx) * u[ mesh ->el[i].iur ];
			};
		};
	};
#endif
};
/******************************************************************************/

void load_uvalues ( struct mesh1d *mesh, double *u,  double (*bnd_u)(double) )

/******************************************************************************/
/*
  Purpose:

    LOAD_UVALUES sets the left and right values of U in each element.
*/
{
  int i;

#ifdef __USE_OPENMP
  int start, end, length, me, nthreads;
#pragma omp parallel private(start, end, length, me, nthreads,i)
  {
    _SMPACK_OMP_READ_SPLIT( me, nthreads, mesh ->ne, start, end, length );
    for( i=start; i<=end; i++ )
    {
      mesh->el[i].ul = ( mesh ->el[i].iul == -1 ) ? bnd_u( mesh->el[i].xl ) : u[ mesh ->el[i].iul ];
      mesh->el[i].ur = ( mesh ->el[i].iur == -1 ) ? bnd_u( mesh->el[i].xr ) : u[ mesh ->el[i].iur ];
    };
  };
#else
  for ( i = 0; i < mesh->ne; i++ )
  {
    mesh->el[i].ul = ( mesh ->el[i].iul == -1 ) ? bnd_u( mesh->el[i].xl ) : u[ mesh ->el[i].iul ];
    mesh->el[i].ur = ( mesh ->el[i].iur == -1 ) ? bnd_u( mesh->el[i].xr ) : u[ mesh ->el[i].iur ];
  };
#endif
};
/******************************************************************************/

void read_solution( struct mesh1d *mesh, const char * filename )

/******************************************************************************/
/*
  Purpose:

    READ_SOLUTION reads the solution from a file.
*/
{
	int i, dummy;
	double xl, xr, ul, ur;
	FILE * fout = fopen( filename, "r");
	dummy = fscanf(fout,"%d",&dummy);
	for( i=0; i<mesh ->ne; i++ ){
		dummy = fscanf(fout,"%lf %lf",&xl,&ul );
		dummy = fscanf(fout,"%lf %lf",&xr,&ur );
		mesh ->el[i].ul = ul;
		mesh ->el[i].ur = ur;
	};
	fclose( fout );
};
/******************************************************************************/

double smooth_kernel ( double delta, double up, double u, double xp, double x )

/******************************************************************************/
/*
  Purpose:

    SMOOTH_KERNEL...
*/
{
  double dist = fabs( xp - x );
#ifdef __KERNEL_SMOOTH
  return ( dist <= delta ) ? exp( 1 - delta / ( delta - dist) ) * 2 * (up - u)/( delta*delta*dist) : 0;
#endif
#ifdef __KERNEL_JUMP
  return ( dist <= delta ) ? 2 * (up - u)/( delta*delta*dist) : 0;
#endif
};
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void write_solution( struct mesh1d *mesh, const char * filename )

/******************************************************************************/
/*
  Purpose:

    WRITE_SOLUTION writes the solution to a file.
*/
{
	int i;
	FILE * fout = fopen( filename, "w");
	fprintf(fout,"%d\n",2*mesh->ne);
	for( i=0; i<mesh ->ne; i++ ){
		fprintf(fout,"%f %f\n",mesh->el[i].xl,mesh->el[i].ul );
		fprintf(fout,"%f %f\n",mesh->el[i].xr,mesh->el[i].ur );
	};
	fclose( fout );
};
/******************************************************************************/

void wwrite_vector( double *x, int n, const char * filename )

/******************************************************************************/
/*
  Purpose:

    WWRITE_VECTOR writes a vector to a file.
*/
{
	int dummy;
	FILE * fout = fopen( filename, "wb");
	dummy = fwrite( &n, sizeof(int),1,fout );
	dummy = fwrite( x, sizeof(double), n, fout );
	printf(" Written %d, out of %d \n",dummy,n );
	fclose( fout );
};
