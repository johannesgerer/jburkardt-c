 /*

     The compiler directives BEG, BEF, LEG, and LEF designate
     four possible floating-point systems as follows.

     Stored double values may be thought of as an array of two
     long values, lx[0] and lx[1].  When the exponent field
     is stored in lx[0], thusly

        ........ xx ........
        7ff00000    00000000
        ..lx[0].    ..lx[1].

     we designate the storage scheme as BE for Big Endian.  When
     the exponent field is stored in lx[1]

        ........ xx ........
        00000000    7ff00000
        ..lx[0].    ..lx[1].

     we designate the storage scheme as LE for Little Endian.  

     Systems that support the IEEE graceful underflow are
     designated with G; those that support only flush-to-zero
     underflow are designated with F.  We assume that G systems
     also support the IEEE infinity and NaN arithmetic, but
     that F systems support neither.

 */

#ifdef BEG
#define DNAN dnan1
#define RMinx RMinDeNorm
#define HighDMinx DZero
#define LowDMinx DMinDeNorm
#define highpart(x) *((long *) &x)
#define lowpart(x) *((long *) &x + 1)
#define VER 1
#endif

#ifdef BEF
#define DNAN dnan1
#define RMinx RMinNorm
#define HighDMinx DMinNorm
#define LowDMinx DZero
#define highpart(x) *((long *) &x)
#define lowpart(x) *((long *) &x + 1)
#define VER 2
#endif

#ifdef LEG
#define DNAN dnan2
#define RMinx RMinDeNorm
#define HighDMinx DZero
#define LowDMinx DMinDeNorm
#define highpart(x) *((long *) &x + 1)
#define lowpart(x) *((long *) &x)
#define VER 1
#endif

#ifdef LEF
#define DNAN dnan2
#define RMinx RMinNorm
#define HighDMinx DMinNorm
#define LowDMinx DZero
#define highpart(x) *((long *) &x + 1)
#define lowpart(x) *((long *) &x)
#define VER 2
#endif

#define Zero  0.0e0
#define Half  0.5e0
#define One   1.0e0
#define RHalf ((float) 0.5)
#define ROne  ((float) 1.0)

#define RExpShift       23
#define RExpBias        ((long) 127)
#define RMinNormExp     ((long) -126)
#define RMinNormBiasExp ((long) 1)
#define RMaxNormBiasExp ((long) 254)
#define RExpMask        0x7f800000L
#define RMinNorm        0x00800000L
#define RMinDeNorm      0x00000001L
#define RMaxDenorm      0x007fffffL
#define MRMaxDenorm     0x807fffffL

#define DExpShift       20
#define DExpBias        ((long) 1023)
#define DMinNormExp     ((long) -1022)
#define DMinNormBiasExp ((long) 1)
#define DMaxNormBiasExp ((long) 2046)
#define DExpMask        0x7ff00000L
#define SignMask        0x80000000L
#define DMinNorm        0x00100000L
#define DZero           0x00000000L
#define DMinDeNorm      0x00000001L
#define DMaxDenorm      0x000fffffL
#define MDMaxDenorm     0x800fffffL

#define allof(x) *((long *) &x)
#define ABS(xxx) ((xxx>Zero)?(xxx):(-xxx))

static long int dnan1[2]  = {0x7ff00004L, 0 } ;
static long int dnan2[2]  = {0, 0x7ff00004L } ;

double  copysign_d ( double argVal, double argSign );
float   copysign_f ( float argVal, float argSign );
int     finite_d ( double arg);
int     finite_f ( float arg );
int     isnan_d ( double arg );
int     isnan_f ( float arg );
double  logb_d ( double arg );
float   logb_f ( float arg );
void    machar_d ( long int *ibeta, long int *it, long int *irnd,
          long int *ngrd, long int *machep, long int *negep, long int *iexp,
          long int *minexp, long int *maxexp, double *eps, double *epsneg,
          double *xmin, double *xmax );
void    machar_s ( long int *ibeta, long int *it, long int *irnd, 
          long int *ngrd, long int *machep, long int *negep, long int *iexp,
          long int *minexp, long int *maxexp, float *eps, float *epsneg, 
          float *xmin, float *xmax );
double  nextafter_d ( double argx, double argy );
float   nextafter_f ( float argx, float argy ); 
double  scalb_d ( double arg, long n );
float   scalb_f ( float arg, long n );
double  test_d ( double argx );
float   test_f ( float argx );
