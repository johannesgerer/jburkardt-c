# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>

# include "dream.h"
# include "dream_user.h"
# include "pdflib.h"
# include "problem1_covariance.h"
/*
  This is the definition, or "home" of the variables.

  Because it is outside of any function in this file, all functions
  in this file can "see" (read or write) the variables.

  The EXTERN statement in "problem1_covariance.h" allows this
  data to be accessed by functions external to this file.
*/
double *c = NULL;
double c_det = 0.0;
double *c_factor = NULL;
int c_factored = 0;
double *zp_mean = NULL;

void covariance_initialize ( int par_num );
double *covariance_set ( int par_num );

/******************************************************************************/

void covariance_initialize ( int par_num )

/******************************************************************************/
/*
  Purpose:

    COVARIANCE_INITIALIZE sets covariance matrix, Cholesky factor, determinant.

  Discussion:

    Note that C_FACTOR is the upper triangular Cholesky factor of the
    covariance matrix C, so that 

      C = C_FACTOR' * C_FACTOR

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int PAR_NUM, the number of parameters.
*/
{
  int i;

  c = covariance_set ( par_num );

  c_factor = r8mat_pofac ( par_num, c );

  c_det = r8mat_podet ( par_num, c_factor );

  zp_mean = ( double * ) malloc ( par_num * sizeof ( double ) );
  for ( i = 0; i < par_num; i++ )
  {
    zp_mean[i] = 0.0;
  }

  c_factored = 1;

  return;
}
/******************************************************************************/

double *covariance_set ( int par_num )

/******************************************************************************/
/*
  Purpose:

    COVARIANCE_SET sets the covariance matrix.

  Discussion:

    This is a multivariate normal distribution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Reference:

    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
    Dave Higdon,
    Accelerating Markov Chain Monte Carlo Simulation by Differential 
    Evolution with Self-Adaptive Randomized Subspace Sampling,
    International Journal of Nonlinear Sciences and Numerical Simulation,
    Volume 10, Number 3, March 2009, pages 271-288.

  Parameters:

    Input, int PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.

    Output, double C[PAR_NUM*PAR_NUM], the covariance matrix.
*/
{
  int i;
  int j;

  c = ( double * ) malloc ( par_num * par_num * sizeof ( double ) );

  for ( j = 0; j < par_num; j++ )
  {
    for ( i = 0; i < par_num; i++ )
    {
      c[i+j*par_num] = 0.5;
    }
  }

  for ( i = 0; i < par_num; i++ )
  {
    c[i+i*par_num] = ( double ) ( i + 1 );
  }

  return c;
}
/******************************************************************************/

void problem_size ( int *chain_num, int *cr_num, int *gen_num, int *pair_num, 
  int *par_num )

/******************************************************************************/
/*
  Purpose:

    PROBLEM_SIZE sets information having to do with dimensions.

  Discussion:

    In the Vrugt paper, PAR_NUM is 100.  For testing, it is reasonable
    to try a tiny value like PAR_NUM = 5 instead.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Reference:

    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
    Dave Higdon,
    Accelerating Markov Chain Monte Carlo Simulation by Differential 
    Evolution with Self-Adaptive Randomized Subspace Sampling,
    International Journal of Nonlinear Sciences and Numerical Simulation,
    Volume 10, Number 3, March 2009, pages 271-288.

  Parameters:

    Output, int *CHAIN_NUM, the total number of chains.
    3 <= CHAIN_NUM.

    Output, int *CR_NUM, the total number of CR values.
    1 <= CR_NUM.

    Output, int *GEN_NUM, the total number of generations.
    2 <= GEN_NUM.

    Output, int *PAIR_NUM, the number of pairs of 
    crossover chains.
    0 <= PAIR_NUM.

    Output, int *PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.
*/
{
  *chain_num = 10;
  *cr_num = 3;
  *gen_num = 10;
  *pair_num = 3;
/*
  *par_num = 100;
*/
  *par_num = 5;
/*
  Initialize the covariance information.
*/
  covariance_initialize ( *par_num );

  return;
}
/******************************************************************************/

void problem_value ( char **chain_filename, char **gr_filename, 
  double *gr_threshold, int *jumpstep, double limits[], int par_num,
  int *printstep, char **restart_read_filename, char **restart_write_filename )

/******************************************************************************/
/*
  Purpose:

    PROBLEM_VALUE sets information, including numeric data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Output, char *CHAIN_FILENAME, the "base" filename
    to be used for the chain files.  If this is NULL
    then the chain files will not be written.  This name should 
    include a string of 0's which will be replaced by the chain 
    indices.  For example, "chain000.txt" would work as long as the
    number of chains was 1000 or less.

    Output, char **GR_FILENAME, the name of the file
    in which values of the Gelman-Rubin statistic will be recorded,
    or NULL if this file is not to be written.

    Output, double *GR_THRESHOLD, the convergence tolerance for
    the Gelman-Rubin statistic.

    Output, int *JUMPSTEP, forces a "long jump" every
    JUMPSTEP generations.

    Output, double LIMITS[2*PAR_NUM], lower and upper bounds
    for each parameter.

    Input, int PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.

    Output, int *PRINTSTEP, the interval between generations on 
    which the Gelman-Rubin statistic will be computed and written to a file.

    Output, char *RESTART_READ_FILENAME, the name of the file
    containing restart information.  If this calculation is not a restart,
    then this should be NULL.

    Output, char *RESTART_WRITE_FILENAME, the name of the file
    to be written, containing restart information.  If a restart file is not
    to be written, this should be NULL.
*/
{
  int j;

  *chain_filename = ( char * ) malloc ( 22 * sizeof ( char ) );
  strcpy ( *chain_filename, "problem1_chain00.txt" );
  *gr_filename = ( char * ) malloc ( 16 * sizeof ( char ) );
  strcpy ( *gr_filename, "problem1_gr.txt" );
  *gr_threshold = 1.2;
  *jumpstep = 5;
  for ( j = 0; j < par_num; j++ )
  {
    limits[0+j*2] =   9.9;
    limits[1+j*2] = +10.0;
  }
  *printstep = 10;
  *restart_read_filename = NULL;
  *restart_write_filename = ( char * ) malloc ( 21 * sizeof ( char ) );
  strcpy ( *restart_write_filename, "problem1_restart.txt" );

  return;
}
/******************************************************************************/

double prior_density ( int par_num, double zp[] )

/******************************************************************************/
/*
  Purpose:

    PRIOR_DENSITY evaluates the prior density function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.

    Input, double ZP[PAR_NUM], the argument of the density
    function.

    Output, real PRIOR_DENSITY, the value of the prior density function.
*/
{
  double value;

  value = r8vec_multinormal_pdf ( par_num, zp_mean, c_factor, c_det, zp );

  return value;
}
/******************************************************************************/

double *prior_sample ( int par_num )

/******************************************************************************/
/*
  Purpose:

    PRIOR_SAMPLE samples from the prior distribution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.

    Output, double PRIOR_SAMPLE[PAR_NUM], the sample from the distribution.
*/
{
  int i;
  double *x;
  double *zp;

  x = ( double * ) malloc ( par_num * sizeof ( double ) );

  for ( i = 0; i < par_num; i++ )
  {
    x[i] = r8_normal_01_sample ( );
  }

  zp = r8mat_mtv_new ( par_num, par_num, c_factor, x );

  for ( i = 0; i < par_num; i++ )
  {
    zp[i] = zp[i] + zp_mean[i];
  }

  free ( x );

  return zp;
}
/******************************************************************************/

double sample_likelihood ( int par_num, double zp[] )

/******************************************************************************/
/*
  Purpose:

    SAMPLE_LIKELIHOOD computes the log likelihood function.

  Discussion:

    This is a one mode Gaussian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int PAR_NUM, the total number of parameters.
    1 <= PAR_NUM.

    Input, double ZP[PAR_NUM], a sample.

    Output, double SAMPLE_LIKELIHOOD, the log likelihood function 
    for the sample.
*/
{
  int i;
  const double pi = 3.141592653589793;
  double value;
  double *x;
  double xcx;
  double *y;

  x = ( double * ) malloc ( par_num * sizeof ( double ) );
  for ( i = 0; i < par_num; i++ )
  {
    x[i] = zp[i] - zp_mean[i];
  }

  y = r8mat_utsol ( par_num, c_factor, x );
/*
  Compute:
    (x-mu)' * inv(C)          * (x-mu)
  = (x-mu)' * inv(R'*R)       * (x-mu)
  = (x-mu)' * inv(R) * inv(R) * (x-mu)
  = y' * y.
*/
  xcx = r8vec_dot_product ( par_num, y, y );

  value = - 0.5 * ( double ) ( par_num ) * log ( 2.0 * pi ) 
    - 0.5 * log ( c_det ) 
    - 0.5 * xcx;

  free ( x );
  free ( y );
    
  return value;
}

