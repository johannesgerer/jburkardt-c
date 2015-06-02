# include <math.h>
# include <stdio.h>

# include "toms722.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS722_PRB.

  Discussion:

    TOMS722_PRB_PRB tests the TOMS722 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TOMS722_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS722 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOMS722_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void test01 ( void ) 

/******************************************************************************/
{
  float eps;
  float epsneg;
  long int i;
  long int ibeta;
  long int iexp;
  static long int nan_f  = { 0x7f800004L } ;
  long int irnd;
  long int it;
  long int machep;
  long int maxexp;
  long int minexp;
  float mone;
  long int n;
  long int negep;
  long int ngrd;
  float one;
  float tt;
  float *x;
  float xmax;
  float xmin;
  float xx;
  float yy;
  float zz;

  union wjc {
    long int jj;
    float xbig;
  } uval;

  machar_s ( &ibeta, &it, &irnd, &ngrd, &machep, &negep,
    &iexp, &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax );

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  MACHAR_S determines constants for floating point arithmetic.\n");
  printf ( "\n" );

  printf ( "ibeta =  %d\n", ibeta );
  printf ( "it =     %d\n", it );
  printf ( "irnd =   %d\n", irnd );
  printf ( "ngrd =   %d\n", ngrd );
  printf ( "machep = %d\n", machep );
  printf ( "negep =  %d\n", negep );
  printf ( "iexp =   %d\n", iexp );
  printf ( "minexp = %d\n", minexp );
  printf ( "maxexp = %d\n", maxexp );

#define DISPLAY_F(s,x) { \
        uval.xbig = x ; \
        printf(s); \
        printf(" %24.16e ",(double)(float) x) ; \
        printf(" %8lX ",uval.jj) ; \
        printf("\n"); \
        }

  DISPLAY_F("eps   ",eps);
  DISPLAY_F("epsneg",epsneg);
  DISPLAY_F("xmin  ",xmin);
  DISPLAY_F("xmax  ",xmax);

  printf ( "\n" );
  printf ( "Tests with moderate numbers:\n");
  printf ( "\n");

  yy = eps;
  DISPLAY_F("yy                   =", yy );
  zz = epsneg;
  DISPLAY_F("zz                   =", zz );
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =", xx );
  xx = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =", xx );

  one = 1.0;
  mone = -1.0;
  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  tt = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",tt);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  tt = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",tt);

  printf(" \n");
  yy = 2.0e0 - eps;
  DISPLAY_F("yy                   =",yy);
  zz = 4.0e0;
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  tt = nextafter_f(-yy,-zz);
  DISPLAY_F("nextafter_f(-yy,-zz) =",tt);

  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  tt = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",tt);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  tt = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",tt);

  printf(" \n");
  yy = xx;
  DISPLAY_F("yy                   =",yy);
  zz = 0.0e0;
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  xx = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",xx);

  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  tt = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",tt);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  tt = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",tt);

  printf ( "\n");
  printf ( "Tests near smallest positive number \n");
  printf ( "\n");

  yy = 0.0e0;
  DISPLAY_F("yy                   =",yy);
  zz = 1.0e0;
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  tt = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",tt);
  tt = nextafter_f(-yy,-zz);
  DISPLAY_F("nextafter_f(-yy,-zz) =",tt);

  printf(" \n");
  yy = xx;
  DISPLAY_F("yy                   =",yy);
  zz = -yy;
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  xx = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",xx);
  tt = nextafter_f(-yy,-zz);
  DISPLAY_F("nextafter_f(-yy,-zz) =",tt);
  tt = copysign_f(tt,one);
  DISPLAY_F("copysign_f(-0.0,1.0) =",tt);
  tt = copysign_f(tt,mone);
  DISPLAY_F("copysign_f(0.0,-1.0) =",tt);

  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  tt = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",tt);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  tt = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",tt);

  if ( VER == 1 ) {
    printf(" \n");
    tt = 11.0e0 * yy;
    DISPLAY_F("tt                   =",tt);
    n = -1;
    xx = scalb_f(tt,n);
    DISPLAY_F("scalb_f(tt,-1)       =",xx);
    xx = tt * 0.5;
    DISPLAY_F("tt * 0.5             =",xx);
    n = -2;
    xx = scalb_f(tt,n);
    DISPLAY_F("scalb_f(tt,-2)       =",xx);
    xx = tt * 0.25;
    DISPLAY_F("tt * 0.25            =",xx);
    n = -3;
    xx = scalb_f(tt,n);
    DISPLAY_F("scalb_f(tt,-3)       =",xx);
    xx = tt * 0.125;
    DISPLAY_F("tt * 0.125           =",xx);
    n = 3;
    xx = scalb_f(tt,n);
    DISPLAY_F("scalb_f(tt,3)        =",xx);
  }

  printf(" \n");
  DISPLAY_F("yy                   =",yy);
  zz = 1.0e0;
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  tt = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",tt);

  printf(" \n");
  yy = xx;
  DISPLAY_F("yy                   =",yy);
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  xx = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",xx);

  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  tt = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",tt);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  tt = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",tt);

  printf ( "\n");
  printf ( "Test near largest positive number \n");
  printf ( "\n");

  yy = xmax;
  DISPLAY_F("yy                   =",yy);
  DISPLAY_F("zz                   =",zz);
  xx = nextafter_f(yy,zz);
  DISPLAY_F("nextafter_f(yy,zz)   =",xx);
  xx = nextafter_f(-yy,zz);
  DISPLAY_F("nextafter_f(-yy,zz)  =",xx);
  tt = copysign_f(yy,mone);
  DISPLAY_F("copysign_f(yy,-1.0)  =",tt);
  tt = copysign_f(tt,one);
  DISPLAY_F("copysign_f(-yy,1.0)  =",tt);

  printf(" \n");
  n = (long int) logb_f(yy);
  printf("n = (long int) logb_f(yy)  = %d \n",n);
  xx = scalb_f(yy,-n);
  DISPLAY_F("scalb_f(yy,-n)       =",xx);
  n = (long int) logb_f(-yy);
  printf("n = (long int) logb_f(-yy) = %d \n",n);
  xx = scalb_f(-yy,-n);
  DISPLAY_F("scalb_f(-yy,-n)      =",xx);

  if (irnd != 5) 
    printf(" \n No tests with infinity and NaN \n");

  if (irnd == 5) {

    printf ( "\n" );
    printf ( "Tests with infinity \n");
    printf ( "\n");

    yy = xmax;
    DISPLAY_F("yy                   =",yy);
    zz = scalb_f(yy,-machep);
    DISPLAY_F("zz                   =",zz);
    tt = copysign_f(zz,mone);
    DISPLAY_F("copysign_f(Inf,-1.0) =",tt);
    tt = copysign_f(tt,zz);
    DISPLAY_F("copysign_f(-Inf,Inf) =",tt);
    xx = nextafter_f(yy,zz);
    DISPLAY_F("nextafter_f(yy,zz)   =",xx);
    xx = nextafter_f(zz,yy);
    DISPLAY_F("nextafter_f(zz,yy)   =",xx);
    xx = nextafter_f(-zz,yy);
    DISPLAY_F("nextafter_f(-zz,yy)  =",xx);
    tt = logb_f(zz);
    DISPLAY_F("logb_f(zz)           =",tt);

    printf ( "\n" );
    printf ( "Tests with NaN \n" );
    printf ( "\n" );

    x = (float *) &nan_f;
    yy = *x;
    DISPLAY_F("yy                   =",yy);
    zz = xmax;
    DISPLAY_F("zz                   =",zz);
    tt = test_f(yy);
    DISPLAY_F("Using a NaN returns  ",tt);
    tt = copysign_f(yy,mone);
    DISPLAY_F("copysign_f(NaN,-1.0) =",tt);
    xx = copysign_f(tt,one);
    DISPLAY_F("copysign_f(-NaN,1.0) =",xx);
    xx = nextafter_f(yy,zz);
    DISPLAY_F("nextafter_f(yy,zz)   =",xx);
    xx = nextafter_f(-yy,zz);
    DISPLAY_F("nextafter_f(-yy,zz)  =",xx);
    xx = nextafter_f(zz,yy);
    DISPLAY_F("nextafter_f(zz,yy)   =",xx);
    xx = nextafter_f(-zz,yy);
    DISPLAY_F("nextafter_f(-zz,yy)  =",xx);

    xx = logb_f(yy);
    DISPLAY_F("logb_f(yy)           =",xx);
    i = 10;
    xx = scalb_f(yy,i);
    DISPLAY_F("scalb_f(yy,10)       =",xx);
  }

  printf ( "\n" );
  printf ( "Special test with 0.0:\n" );
  printf ( "\n" );

  yy = 0.0e0;
  DISPLAY_F("yy                   =",yy);
  i = 10;
  xx = scalb_f(yy,i);
  DISPLAY_F("scalb_f(yy,10)       =",xx);
  xx = logb_f(yy);
  DISPLAY_F("logb_f(yy)           =",xx);

  return;
}
/******************************************************************************/

void test02 ( void ) 

/******************************************************************************/
 /*
     Test driver for logb, scalb_d and nextafter functions.

     Stored double values may be thought of as an array of two
     long int values, lx[0] and lx[1].  When the exponent field
     is stored in lx[0], thusly

        ........ xx ........
        7ff00000    00000000
        ..lx[0].    ..lx[1].
 
     we designate the storage scheme as MN.  When the exponent
     field is stored in lx[1]

        ........ xx ........
        00000000    7ff00000
        ..lx[0].    ..lx[1].
 
     we designate the storage scheme as NM.  

     Systems that support the IEEE graceful underflow are
     designated with G; those that support only flush-to-zero
     underflow are designated with F.  We assume that G systems
     also support the IEEE infinity and NaN arithmetic, but
     that F systems support neither.

 */
{
  double eps;
  double epsneg;
  long int i;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  long int machep;
  long int maxexp;
  long int minexp;
  double mone;
  long int n;
  long int negep;
  long int ngrd;
  double one;
  double tt;
  double *x;
  double xmax;
  double xmin;
  double xx;
  double yy;
  double zz;

  union wjc{
    long int jj[2];
    double xbig;
  } uval;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  MACHAR_D determines constants for floating point arithmetic.\n");
  printf ( "\n" );

  machar_d(&ibeta,&it,&irnd,&ngrd,&machep,&negep,
    &iexp,&minexp,&maxexp,&eps,&epsneg,&xmin,&xmax);

  printf("ibeta = %d\n",ibeta);
  printf("it    = %d\n",it);
  printf("irnd  = %d\n",irnd);
  printf("ngrd  = %d\n",ngrd);
  printf("machep = %d\n",machep);
  printf("negep = %d\n",negep);
  printf("iexp = %d\n",iexp);
  printf("minexp = %d\n",minexp);
  printf("maxexp = %d\n",maxexp);

#define DISPLAY_D(s,x) { \
        uval.xbig = x ; \
        printf(s); \
        printf(" %24.16e ",(double) x) ; \
        for(i=0;i<2;i++) printf(" %8lX ",uval.jj[i]) ; \
        printf("\n"); \
        }

      DISPLAY_D("eps   ",eps);
      DISPLAY_D("epsneg",epsneg);
      DISPLAY_D("xmin  ",xmin);
      DISPLAY_D("xmax  ",xmax);

      printf(" \n Tests with moderate numbers \n");
      printf(" \n");
      yy = eps;
      DISPLAY_D("yy                   =",yy);
      zz = epsneg;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      xx = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",xx);

      one = 1.0e0;
      mone = -1.0e0;
      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

      printf(" \n");
      yy = 2.0e0 - eps;
      DISPLAY_D("yy                   =",yy);
      zz = 4.0e0;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,-zz);
      DISPLAY_D("nextafter_d(-yy,-zz) =",tt);

      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

      printf(" \n");
      yy = xx;
      DISPLAY_D("yy                   =",yy);
      zz = 0.0e0;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);

      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

      printf(" \n Tests near smallest positive number \n");
      printf(" \n");
      yy = 0.0e0;
      DISPLAY_D("yy                   =",yy);
      zz = 1.0e0;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);
      tt = nextafter_d(-yy,-zz);
      DISPLAY_D("nextafter_d(-yy,-zz) =",tt);

      printf(" \n");
      yy = xx;
      DISPLAY_D("yy                   =",yy);
      zz = -yy;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);
      tt = nextafter_d(-yy,-zz);
      DISPLAY_D("nextafter_d(-yy,-zz) =",tt);
      tt = copysign_d(tt,one);
      DISPLAY_D("copysign_d(-0.0,1.0) =",tt);
      tt = copysign_d(tt,mone);
      DISPLAY_D("copysign_d(0.0,-1.0) =",tt);

      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

   if (VER == 1) {
      printf(" \n");
      tt = 11.0e0 * yy;
      DISPLAY_D("tt                   =",tt);
      n = -1;
      xx = scalb_d(tt,n);
      DISPLAY_D("scalb_d(tt,-1)       =",xx);
      xx = tt * 0.5;
      DISPLAY_D("tt * 0.5             =",xx);
      n = -2;
      xx = scalb_d(tt,n);
      DISPLAY_D("scalb_d(tt,-2)       =",xx);
      xx = tt * 0.25;
      DISPLAY_D("tt * 0.25            =",xx);
      n = -3;
      xx = scalb_d(tt,n);
      DISPLAY_D("scalb_d(tt,-3)       =",xx);
      xx = tt * 0.125;
      DISPLAY_D("tt * 0.125           =",xx);
      n = 3;
      xx = scalb_d(tt,n);
      DISPLAY_D("scalb_d(tt,3)        =",xx);
     }

      printf(" \n");
      DISPLAY_D("yy                   =",yy);
      zz = 1.0e0;
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);

      printf(" \n");
      yy = xx;
      DISPLAY_D("yy                   =",yy);
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);

      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

      printf(" \n Test near largest positive number \n");
      printf(" \n");
      yy = xmax;
      DISPLAY_D("yy                   =",yy);
      DISPLAY_D("zz                   =",zz);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);
      tt = copysign_d(yy,mone);
      DISPLAY_D("copysign_d(yy,-1.0)  =",tt);
      tt = copysign_d(-tt,one);
      DISPLAY_D("copysign_d(-yy,1.0)  =",tt);

      printf(" \n");
      n = (long int) logb_d(yy);
      printf("n = (long int) logb_d(yy)  = %d \n",n);
      tt = scalb_d(yy,-n);
      DISPLAY_D("scalb_d(yy,-n)       =",tt);
      n = (long int) logb_d(-yy);
      printf("n = (long int) logb_d(-yy) = %d \n",n);
      tt = scalb_d(-yy,-n);
      DISPLAY_D("scalb_d(-yy,-n)      =",tt);

    if (irnd != 5) 
      printf(" \n No tests with infinity and NaN \n");

    if (irnd == 5) {
      printf(" \n Tests with infinity \n");
      printf(" \n");
      yy = xmax;
      DISPLAY_D("yy                   =",yy);
      zz = scalb_d(yy,-machep);
      DISPLAY_D("zz                   =",zz);
      tt = copysign_d(zz,mone);
      DISPLAY_D("copysign_d(Inf,-1.0) =",tt);
      tt = copysign_d(tt,zz);
      DISPLAY_D("copysign_d(-Inf,Inf) =",tt);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      xx = nextafter_d(zz,yy);
      DISPLAY_D("nextafter_d(zz,yy)   =",xx);
      xx = nextafter_d(-zz,yy);
      DISPLAY_D("nextafter_d(-zz,yy)  =",xx);
      xx = logb_d(zz);
      DISPLAY_D("logb_d(zz)           =",xx);

      printf(" \n Tests with NaN \n");
      printf(" \n");
      x = (double *) DNAN;
      yy = *x;
      DISPLAY_D("yy                   =",yy);
      zz = xmax;
      DISPLAY_D("zz                   =",zz);
      tt = test_d(yy);
      DISPLAY_D("Using a NaN returns ",tt);
      tt = copysign_d(yy,mone);
      DISPLAY_D("copysign_d(NaN,-1.0) =",tt);
      tt = copysign_d(tt,yy);
      DISPLAY_D("copysign_d(-NaN,NaN) =",tt);
      xx = nextafter_d(yy,zz);
      DISPLAY_D("nextafter_d(yy,zz)   =",xx);
      tt = nextafter_d(-yy,zz);
      DISPLAY_D("nextafter_d(-yy,zz)  =",tt);
      xx = nextafter_d(zz,yy);
      DISPLAY_D("nextafter_d(zz,yy)   =",xx);
      tt = nextafter_d(-zz,yy);
      DISPLAY_D("nextafter_d(-zz,yy)  =",tt);

      tt = logb_d(yy);
      DISPLAY_D("logb_d(yy)           =",tt);
      i = 10;
      tt = scalb_d(yy,i);
      DISPLAY_D("scalb_d(yy,10)       =",tt);
     }

      printf(" \n Special test with 0.0 \n");
      printf(" \n");
      yy = 0.0e0;
      DISPLAY_D("yy                   =",yy);
      i = 10;
      tt = scalb_d(yy,i);
      DISPLAY_D("scalb_d(yy,10)       =",tt);
      tt = logb_d(yy);
      DISPLAY_D("logb_d(yy)           =",tt);

  return;

}



