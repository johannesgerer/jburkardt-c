# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( );
int absorb ( int *seed );
double cross ( double e );
double dist2c ( double e, int *seed );
double energy ( int *seed );
void output ( int na, double ea, double sa, int nr, double er, double sr, 
  int nt, double et, double st, int ntot );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_uniform_01 ( int *seed );
void scatter ( int *seed, double *e, double *mu, double *azm );
void source ( int *seed, double *e, double *mu, double *azm, double *x, 
  double *y, double *z );
void timestamp ( void );
void update ( double mu, double azm, double d, double *x, double *y, double *z );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for the reactor shielding simulation.

  Discussion:

    This is a Monte Carlo simulation, using
    uniform random numbers, which investigates the
    effectiveness of a shield intended to absorb the
    neutrons emitted from a nuclear reactor.
   
    The reactor is modeled as a point source,
    located at (0,0,0).
   
    A particle emitted from the reactor has a random
    initial direction, and an energy selected from
    [Emin,Emax] with a 1/Sqrt(E) distribution.
   
    The shield is modeled as a wall of thickness THICK,
    extending from 0 to THICK in the X direction, and
    extending forever in the Y and Z directions.
   
    Based on the particle energy, a distance D is computed
    which measures how far the particle could travel through
    the shield before colliding.
   
    Based on the particle direction, the position is updated
    by D units.
   
    If the particle is now to the left of the shield, it is
    counted as being REFLECTED.
   
    If the particle is to the right of the shield, it is 
    counted as being ABSORBED.
   
    If the particle is inside the shield, it has COLLIDED.
    A particle that collides is either absorbed (end of story)
    or SCATTERED with a new random direction and a new (lower)
    energy.
   
    Every particle is followed from origin to its final fate,
    which is reflection, transmission, or absorption.
    At the end, a summary is printed, giving the number of
    particles with each fate, and the average energy of each
    group of particles.
   
    Increasing NTOT, the number of particles used, will improve the
    expected reliability of the results.
   
    Increasing THICK, the thickness of the shield, should 
    result in more absorptions and reflections.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Local Parameters:

    Local, double AZM, the azimuthal angle of the particle's
    direction.

    Local, double D, the distance that the particle can
    travel through the slab, given its current energy.

    Local, double E, the energy of the particle.

    Local, double EA, energy absorbed by the slab.

    Local, double ER, energy reflected by the slab.

    Local, double ET, energy transmitted through the slab.

    Local, double MU, the cosine of the angle between the
    particle's direction and the X axis.

    Local, int NA, number of particles absorbed by the slab.

    Local, int NPART, the index of the current particle.

    Local, int NR, number of particles reflected by the slab.

    Local, int NT, number of particles transmitted by the slab.

    Local, int NTOT, the total number of particles to be
    emitted from the neutron source.

    Local, double SA, standard deviation of absorbed energy.

    Local, double SR, standard deviation of reflected energy.

    Local, double ST, standard deviation of transmitted energy.

    Local, double THICK, the thickness of the slab that is
    intended to absorb most of the particles.

    Local, double X, Y, Z, the current position of the particle.
*/
{
  double azm;
  double d;
  double e;
  double ea;
  double er;
  double et;
  int i;
  double mu;
  int na;
  int npart;
  int nr;
  int nt;
  int ntot = 100000;
  int part;
  double sa;
  int seed;
  double sr;
  double st;
  int test;
  int test_num = 5;
  double thick = 2.0;
  double x;
  double y;
  double z;

  timestamp ( );

  printf ( "\n" );
  printf ( "REACTOR_SIMULATION\n" );
  printf ( "  C version\n" );
  printf ( "  The reactor shielding simulation.\n" );
  printf ( "\n" );
  printf ( "  Shield thickness is THICK = %g\n", thick );
  printf ( "  Number of simulated particles is NTOT = %d\n", ntot );
  printf ( "  Number of tests TEST_NUM = %d\n", test_num );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test # %d\n", test );
    printf ( "  SEED = %d\n", seed );
/*
  Initialize.
*/
    ea = 0.0;
    er = 0.0;
    et = 0.0;
    na = 0;
    nr = 0;
    nt = 0;
    sa = 0.0;
    sr = 0.0;
    st = 0.0;
/*
  Loop over the particles.
*/
    for ( part = 1; part <= ntot; part++ )
    {
/*
  Generate a new particle.
*/
      source ( &seed, &e, &mu, &azm, &x, &y, &z );

      while ( 1 )
      {
/*
  Compute the distance that the particle can travel through the slab,
  based on its current energy.
*/
        d = dist2c ( e, &seed );
/*
  Update the particle's position by D units.
*/
        update ( mu, azm, d, &x, &y, &z );
/*
  The particle was reflected by the shield, and this path is complete.
*/
        if ( x < 0.0 )
        {
          nr = nr + 1;
          er = er + e;
          sr = sr + e * e;
          break;
        }
/*
  The particle was transmitted through the shield, and this path is complete.
*/
        else if ( thick < x )
        {
          nt = nt + 1;
          et = et + e;
          st = st + e * e;
          break;
        }
/*
  The particle collided with the shield, and was absorbed.  This path is done.
*/
        else if ( absorb ( &seed ) )
        {
          na = na + 1;
          ea = ea + e;
          sa = sa + e * e;
          break;
        }
/*
  The particle collided with the shield and was scattered.
  Find the scattering angle and energy, and continue along the new path.
*/
        else
        {
          scatter ( &seed, &e, &mu, &azm );
        }
      }
    }
/*
  Print the results of the simulation.
*/
    output ( na, ea, sa, nr, er, sr, nt, et, st, ntot );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "REACTOR_SIMULATION:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int absorb ( int *seed )

/******************************************************************************/
/*
  Purpose:

    ABSORB determines if a colliding particle is absorbed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, logical ABSORB, is TRUE if the particle is absorbed.

  Local parameters:

    Local, double PA, the probability of absorption.
*/
{
  double pa = 0.1;
  double u;
  int value;

  u = r8_uniform_01 ( seed );

  if ( u <= pa )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

double cross ( double e )

/******************************************************************************/
/*
  Purpose:

    CROSS returns the "cross section" of a particle based on its energy.

  Discussion:

    The particle's cross section is a measure of its likelihood to collide
    with the material of the slab.  This quantity typically depends on both
    the particle's energy and the kind of medium through which it is traveling.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input, double E, the energy of the particle.

    Output, double CROSS, the cross section.
*/
{
  double s;
  double value;
  double y;

  s = r8_abs ( sin ( 100.0 * ( exp ( e ) - 1.0 ) ) 
    + sin ( 18.81 * ( exp ( e ) - 1.0 ) ) );

  y = r8_max ( 0.02, s );

  value = 10.0 * exp ( -0.1 / y );

  return value;
}
/******************************************************************************/

double dist2c ( double e, int *seed )

/******************************************************************************/
/*
  Purpose:

    DIST2C returns the distance to collision.

  Discussion:

    Assuming the particle has a given energy, and assuming it is currently
    somewhere inside the shield, it is possible to determine a typical distance
    which the particle can travel before it collides with the material of
    the shield.

    The computation of the collision distance is made by estimating a
    "cross section" (as though having more energy made the particle "bigger"
    and hence more likely to collide) and then randomly selecting a distance
    that is logarithmically distributed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input, double E, the energy of the particle.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double DIST2C, the distance the particle can travel
    through the slab before colliding.
*/
{
  double u;
  double value;

  u = r8_uniform_01 ( seed );

  value = - log ( u ) / cross ( e );

  return value;
}
/******************************************************************************/

double energy ( int *seed )

/******************************************************************************/
/*
  Purpose:

    ENERGY assigns an energy to an emitted particle.

  Discussion:

    The energy E is in the range [EMIN,EMAX], with distribution
    const/sqrt(energy).

    An inverse function approach is used to compute this.

    The energies are measured in MeV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double ENERGY, a randomly chosen energy that is
    distributed as described above.

  Local parameters:

    Local, double EMIN, EMAX, the minimum and maximum
    energies.
*/
{
  double c;
  double emax = 2.5;
  double emin = 1.0E-03;
  double u;
  double value;

  u = r8_uniform_01 ( seed );

  c = 1.0 / ( 2.0 * ( sqrt ( emax ) - sqrt ( emin ) ) );

  value = ( u / ( 2.0 * c ) + sqrt ( emin ) );
  value = value * value;

  return value;
}
/******************************************************************************/

void output ( int na, double ea, double sa, int nr, double er, double sr, 
  int nt, double et, double st, int ntot )

/******************************************************************************/
/*
  Purpose:

    OUTPUT prints the results of the reactor shielding simulation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input, int NA, number of particles absorbed by the slab.

    Input, double EA, energy absorbed by the slab.

    Input, double SA, the sum of the squares of the 
    absorbed energies.

    Input, int NR, number of particles reflected by the slab.

    Input, double ER, energy reflected by the slab.

    Input, double SR, the sum of the squares of the 
    reflected energies.

    Input, int NT, number of particles transmitted by the slab.

    Input, double ET, energy transmitted through the slab.

    Input, double ST, the sum of the squares of the 
    transmitted energies.

    Input, int NTOT, the total number of particles.
*/
{
  double ea_ave;
  double er_ave;
  double et_ave;
  double etot;
  double pa;
  double pr;
  double pt;
  double ptot;

  printf ( "\n" );
  printf ( "  The Reactor Shielding Problem:\n" );
  printf ( "\n" );
  printf ( "                           Total                   Average\n" );
  printf ( "                    #      Energy      " );
  printf ( "Percent     Energy         StDev\n" );
  printf ( "\n" );

  etot = ea + er + et;

  if ( 0 < na )
  {
    ea_ave = ea / ( double ) ( na );
    sa = sqrt ( sa / ( double ) ( na ) - ea_ave * ea_ave );
  }
  else
  {
    ea_ave = 0.0;
  }

  pa = ( double ) ( na * 100 ) / ( double ) ( ntot );

  printf ( "Absorbed     %8d  %14g  %6.2f  %14g  %14g\n",
    na, ea, pa, ea_ave, sa );

  if ( 0 < nr )
  {
    er_ave = er / ( double ) ( nr );
    sr = sqrt ( sr / ( double ) ( nr ) - er_ave * er_ave );
  }
  else
  {
    er_ave = 0.0;
  }

  pr = ( double ) ( nr * 100 ) / ( double ) ( ntot );

  printf ( "Reflected    %8d  %14g  %6.2f  %14g  %14g\n",
    nr, er, pr, er_ave, sr );

  if ( 0 < nt )
  {
    et_ave = et / ( double ) ( nt );
    st = sqrt ( st / ( double ) ( nt ) - et_ave * et_ave );
  }
  else
  {
    et_ave = 0.0;
  }

  pt = ( double ) ( nt * 100 ) / ( double ) ( ntot );

  printf ( "Transmitted  %8d  %14g  %6.2f  %14g  %14g\n",
    nt, et, pt, et_ave, st );

  ptot = 100.0;

  printf ( "\n" );
  printf ( "Total        %8d  %14g  %6.2f\n",
    ntot, etot, ptot );

  return;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void scatter ( int *seed, double *e, double *mu, double *azm )

/******************************************************************************/
/*
  Purpose:

    SCATTER returns the new direction and energy of a particle that is scattered.

  Discussion:

    The scattering direction is chosen uniformly on the sphere.

    The energy of the scattered particle is chosen uniformly in
    [ 0.3*E, E ].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input/output, int *SEED, a seed for the random
    number generator.

    Input/output, double *E.  On input, the particle energy
    before collision.  On output, the particle energy after collision
    and scattering.

    Output, double *MU, the cosine of the angle between the
    particle's direction and the X axis.

    Output, double *AZM, the azimuthal angle of the particle's
    direction.
*/
{
  double pi = 3.141592653589793;
  double u;

  u = r8_uniform_01 ( seed );
  *mu = - 1.0 + 2.0 * u;

  u = r8_uniform_01 ( seed );
  *azm = u * 2.0 * pi;

  u = r8_uniform_01 ( seed );
  *e = ( u * 0.7 + 0.3 ) * ( *e );

  return;
}
/******************************************************************************/

void source ( int *seed, double *e, double *mu, double *azm, double *x, 
  double *y, double *z )

/******************************************************************************/
/*
  Purpose:

    SOURCE generates a new particle from the neutron source.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double *E, the initial energy of the particle.

    Output, double *MU, the cosine of the angle between the
    particle's direction and the X axis.

    Output, double *AZM, the azimuthal angle of the particle's
    direction.

    Output, double *X, *Y, *Z, the initial coordinates of the particle.
*/
{
  double pi = 3.141592653589793;
  double u;

  u = r8_uniform_01 ( seed );
  *mu = u;

  u = r8_uniform_01 ( seed );
  *azm = u * 2.0 * pi;

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;

  *e = energy ( seed );

  return;
}
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void update ( double mu, double azm, double d, double *x, double *y, double *z )

/******************************************************************************/
/*
  Purpose:

    UPDATE determines the position of the particle after it has traveled D units.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2012

  Author:

    Original FORTRAN77 version by Kahaner, Moler, Nash.
    C version by John Burkardt.

  Reference:

    David Kahaner, Cleve Moler, Steven Nash,
    Numerical Methods and Software,
    Prentice Hall, 1989,
    ISBN: 0-13-627258-4,
    LC: TA345.K34.

  Parameters:

    Input, double MU, the cosine of the angle between the
    particle's direction and the X axis.

    Input, double AZM, the azimuthal angle of the particle's
    direction.

    Input, double D, the distance the particle traveled.

    Input/output, double *X, *Y, *Z.  On input, the previous
    coordinates of the particle.  On output, the updated coordinates of the
    particle.
*/
{
  double s;

  s = sqrt ( 1.0 - mu * mu );

  *x = *x + d * mu;
  *y = *y + d * s * cos ( azm );
  *z = *z + d * s * sin ( azm );

  return;
}
