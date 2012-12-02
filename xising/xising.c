/*
 * xising
 * 
 * Xwindow program to simulate the 2-d Ising model. The underlying
 * approach uses many "demons" circulating around the lattice trying
 * to flip spins.
 * 
 * The latest version of this program can be found at
 * http://penguin.phy.bnl.gov/www/xtoys/xtoys.html
 * 
 * Michael Creutz   creutz@wind.phy.bnl.gov
 * 
   November 1999, version 2.11: corrects problems on alphas pointed out
   by Kari Rummukainen

 * 
 * compiling:  cc -O -o xising -L/usr/X11R6/lib xising.c -lm -lX11
 * 
 * On suns:
 * 
 * cc -O -I/usr/openwin/include -L/usr/openwin/lib newxising.c -lm -lX11
 * 
 * If you find some machine supporting X that this does not work on,
 * please let me know.
 * 
 */

/* the usual includes */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/*
 * initial lattice dimensions: ncols will be truncated to a multiple
 * of WORDSIZE nrows will be changed if a multiple of 29 since I let
 * the demons hop 29 sites to get a bit of randomness
 */
int nrows = 256, ncols = 256;
long **field = NULL;		/* to store the system; it will be
				 * dynamically allocated to a [ncols]
				 * by [nrows/WORDSIZE] array */

/* some convenient definitions */
#define WORDSIZE  (CHAR_BIT*sizeof(long))
#define LEFTBIT  (~LONG_MAX)
#define nextrow (row+1-nrows*(row>=nrows-1))
/* nshift is here to allow for different endianness */
#define nextcol (col+nshift*(1-hwords*(col==topcol)))
#define previousrow (row-1+nrows*(row<1))
#define previouscol (col-nshift*(1-hwords*(col==bottomcol)))

char stringbuffer[100];
static char *progname;
char bits[256];			/* for counting bits */
/* hwords will store the number of long words in a single row */
int hwords, volume;

/* bitcount counts the set bits in a word; call setup() before using */
#define bitcount8(i)  (bits[(i)&255])
#define bitcount16(i) (bitcount8(i)+bitcount8((i)>>8))
/* check if long int is 32 bits */
#if ((0x7FFFFFFF)==(LONG_MAX))	/* 32 bits */
#define bitcount(i)   (bitcount16(i)+bitcount16((i)>>16))

/* note, current linux glibc has a bug in mrand48, so it only gives 16+1 
   random bits.  When that is fixed, this can go back to mrand48() */
#define RND (mrand48()^lrand48())
#else				/* assume 64 bits */
#define bitcount32(i)   (bitcount16(i)+bitcount16((i)>>16))
#define bitcount(i)   (bitcount32(i)+bitcount32((i)>>32))
#define RND (mrand48()^(mrand48()<<32))
#endif

/* things to monitor temperature, etc. */
double beta;			/* the inverse temperature */
long energycount = 0, volumecount = 0;
int heater = 0, cooler = 0;
int canonical = 1;
int boundary = 0, algorithm = 0, paused = 0, iteration = 0;
/* variables to allow for peculiar bit orderings on some machines */
int nshift, topcol, bottomcol;

/* for local algorithm; use WORDSIZE two-bit demons */
long demon0, demon1, checkermask;
long *work0 = NULL, *work1 = NULL;

/*
 * bitprob represents the probability for the corresponding demon bit
 * to be set multiplied by 1<<16 and rounded to an integer; used to
 * speed up demon refreshing
 */
long bitprob0, bitprob1;

/* stuff for the cluster algorithm; use three demon bits */
#define DEMONBITS 3
int ndemons;
long **cluster = NULL;		/* to store cluster shape */
long **sadx = NULL, **sady = NULL;	/* labels sad demons */
int **unchecked = NULL;		/* to reduce redundant calculations */
int boxtop, boxbottom;		/* will enclose cluster */
long **demon = NULL;		/* for demons in cluster algorithm */
long **newdemon = NULL;
int *activerow = NULL;
long demonindex = 0;		/* for shuffling demons */
long clustergrowth = 0;
int direction = 1;		/* for vertical sweeps through
				 * cluster */

/* things for regulating the updating speed */
int speed = 0;			/* updating delay proportional to
				 * speed   */
int delay = 0;			/* counter for implementing speed
				 * control */
struct timeval timeout;		/* timer for speed control */

/* various window stuff */
Display *display;
int screen;
Window window, quitbutton, pausebutton, algobutton, heatwindow,
  boundwindow, speedbutton, canonbutton, makebutton();
GC gc, gcpen, gcclear;
XImage *spinimage = NULL, *clusterimage = NULL;
XFontStruct *font = NULL;
int font_height, font_width;
unsigned int windowwidth, windowheight;
XSizeHints size_hints;
int depth, darkcolor, lightcolor, black, white, spinup, spindown;
long event_mask;
/* all my functions */
void refreshdemons(int step);	/* refreshes the cluster demons, step
				 * steps over some of them for only a
				 * partial refresh */
void drawbutton(), openwindow(), setup(), makebuttons(), growcluster(),
  startcluster(), localupdate(), repaint(), setalg(),
  setheater(), setboundary(), fixthermometer(), check(),
  leftshift(), rightshift(), allocarray(), heatbath();
/* text for the buttons */
char *boundtext[4] = {"periodic", "antiperiodic",
"up boundary", "down"};
char *heatertext[3] = {"run free", "heat", "cool"};
char *algotext[2] = {"local", "cluster"};

/* general layout dimensions */
#define LEFTOFFSET 64
#define TOPOFFSET 124
#define HEATERHEIGHT 54
#define HEATERWIDTH 86
#define BOUNDHEIGHT 72
#define BOUNDWIDTH  116
#define BUTTONWIDTH 48
#define BUTTONHEIGHT 18
#define TTOP 128
#define THGT (nrows*2)
#define TBOTTOM (TTOP+THGT)
#define TLEFT 24
#define TWDTH 5

int main(argc, argv)
  int argc;
  char **argv;
{
  unsigned int width, height, i;
  XEvent report;
  progname = argv[0];
  openwindow(argc, argv);
  setup();
  makebuttons();
  /* loop forever, looking for events */
  while (1) {
    if (clustergrowth)		/* finish cluster regardless of
				 * algorithm */
      while (clustergrowth) {
	growcluster();
      }
    else if ((0 == XPending(display)) && (0 == paused)) {
      /* no events waiting */
      if (delay) {		/* don't update yet */
	delay--;
	/*
	 * this use of select() seems a kludge to me; why can't
	 * usleep() be more standard?
	 */
	timeout.tv_sec = 0;
	timeout.tv_usec = 50000;/* .05 sec per delay unit */
	select(0, NULL, NULL, NULL, &timeout);
      } else {
	delay = speed;		/* reset delay counter */
	if (1 == algorithm) {
	  startcluster();
	} else
	  localupdate();
      }
    }
    if (paused | XPending(display)) {
      XNextEvent(display, &report);	/* find out what X wants done */
      switch (report.type) {
      case Expose:		/* the window has been exposed and
				 * needs redrawing */
	if ((report.xexpose.window) != window)
	  break;		/* cuts down flashing, but you might
				 * remove this line if things aren't
				 * being redrawn */
	if (report.xexpose.count != 0)
	  break;		/* more in queue, wait for them */
	repaint();
	break;
      case ConfigureNotify:	/* the window has been resized */
	width = report.xconfigure.width;
	height = report.xconfigure.height;
	if ((width < size_hints.min_width) ||
	    (height < size_hints.min_height)) {
	  fprintf(stderr, "%s: window too small to proceed.\n", progname);
	  XUnloadFont(display, font->fid);	/* I'm still too neat */
	  XFreeGC(display, gc);	/* from my Amiga days */
	  XFreeGC(display, gcpen);
	  XFreeGC(display, gcclear);
	  XCloseDisplay(display);
	  exit(1);
	}
	if ((width != windowwidth) || (height != windowheight)) {
	  windowwidth = width;
	  windowheight = height;
	  /* new lattice dimensions */
	  ncols = width - 75;
	  nrows = (height - 164 - 24) / 2;
	  makebuttons();
	}
	break;
      case ButtonPress:
	if (report.xbutton.window == quitbutton) {
	  XUnloadFont(display, font->fid);
	  XFreeGC(display, gc);
	  XFreeGC(display, gcpen);
	  XFreeGC(display, gcclear);
	  XCloseDisplay(display);
	  exit(1);
	} else if (report.xbutton.window == pausebutton) {
	  paused = 1 - paused;
	  drawbutton(pausebutton, 0, 0, BUTTONWIDTH, BUTTONHEIGHT,
		     "pause", 1 - 2 * paused);
	} else if (report.xbutton.window == algobutton) {
	  /* toggle algorithm */
	  algorithm = 1 - algorithm;
	  /* 1 for cluster, 0 for local updating */
	  setalg();
	  for (i = 0; i < 2; i++)
	    drawbutton(algobutton, 0, i * 18, 68, 18,
		       algotext[i], 1 - 2 * (algorithm == i));
	} else if (report.xbutton.window == heatwindow) {
	  /* adjust heater */
	  setheater((report.xbutton.y * 3) / HEATERHEIGHT);
	} else if (report.xbutton.window == boundwindow)
	  /* set boundary conditions */
	{
	  setboundary((report.xbutton.y * 4) / BOUNDHEIGHT);
	} else if (report.xbutton.window == canonbutton) {
	  canonical = 1 - canonical;
	  if (canonical)
	    drawbutton(canonbutton, 0, 0, font_width * 15, BUTTONHEIGHT,
		       "canonical", -1);
	  else
	    drawbutton(canonbutton, 0, 0, font_width * 15, BUTTONHEIGHT,
		       "microcanonical", 1);
	} else if (report.xbutton.window == speedbutton) {
	  /*
	   * reset speed; speed=0 is the fastest, 10 the slowest
	   */
	  speed = 10 - (11 * (report.xbutton.x - 1)) / (BOUNDWIDTH - 2);
	  if (speed < 0)
	    speed = 0;
	  if (speed > 10)
	    speed = 10;
	  delay = speed;
	  drawbutton(speedbutton, 0, 0, BOUNDWIDTH, 18, "speed", -1);
	  drawbutton(speedbutton, 1 + ((10 - speed) * (BOUNDWIDTH - 2)) / 11,
		     1, (BOUNDWIDTH - 2) / 11, 16, "", 2);
	} else {		/* do an update on a random mouse
				 * press */
	  if (canonical &&
	      (report.xbutton.y > TOPOFFSET) &&
	      (report.xbutton.x < LEFTOFFSET)) {
	    beta = 0.5 / (0.6 + 1.0 * (TBOTTOM - report.xbutton.y) / THGT);
	    bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
	    bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));

	    fixthermometer();
	  } else {
	    if (1 == algorithm)
	      startcluster();
	    else
	      localupdate();
	  }
	}
	break;
      default:
	break;
      }
    }
  }
  return 1;
}

void startcluster()
{				/* flip the old cluster and start a
				 * new one */
  long xedge, yedge;
  int row, col, bit;
  long currentword;
  long index0, index1;
  /* flip previous cluster and adjust demons on boundary */
  for (row = 0; row < nrows; row++) {
    for (col = 0; col < hwords; col++) {
      index0 = (demonindex + 2 * (col + row * hwords));
      index1 = (index0 + 1);
      while (index0 >= ndemons)
	index0 -= ndemons;	/* faster than using % */
      while (index1 >= ndemons)
	index1 -= ndemons;
      if (0 == unchecked[row][col]) {	/* don't bother if unchecked */
	unchecked[row][col] = 1;/* reset unchecked */
	currentword = cluster[row][col];
	/* find cluster edges */
	xedge = currentword ^ ((currentword << 1)
		 | (1 & (cluster[row][nextcol] >> (WORDSIZE - 1))));
	yedge = currentword ^ cluster[nextrow][col];
	/* update demons and field */
	for (bit = 0; bit < DEMONBITS; bit++) {
	  demon[bit][index0]
	    = (newdemon[bit][index0] & xedge) |
	    ((~xedge) & demon[bit][index0]);
	  demon[bit][index1]
	    = (newdemon[bit][index1] & yedge) |
	    ((~yedge) & demon[bit][index1]);
	}
	field[row][col] ^= currentword;
      }				/* end of if checked */
      /* do any cooling or heating required */
      /*
       * This doesn't work very well in the cluster algorithm; the
       * demons hold the heat for a while before dumping it, thus the
       * thermomenter doesn't read reliably.  For major heating and
       * cooling it seems better to go to the local algorithm.
       */
      if ((!canonical) && (drand48() < 0.05)) {
	if (heater) {
	  demon[0][index0] ^= (RND & RND);
	  demon[0][index1] ^= (RND & RND);
	} else if (cooler) {
	  demon[0][index0] &= RND;
	  demon[0][index1] &= RND;
	} else {		/* rotate the demon's first bit to
				 * help scramble */
	  demon[0][index0] = (demon[0][index0] << 1)
	    | (1 & (demon[0][index0] >> (WORDSIZE - 1)));
	  demon[0][index1] = (demon[0][index1] << 1)
	    | (1 & (demon[0][index1] >> (WORDSIZE - 1)));
	}
      }
      /* accumulate counts */
      currentword = demon[0][index0];
      energycount += bitcount(currentword);
      currentword = demon[0][index1];
      energycount += bitcount(currentword);
      volumecount += WORDSIZE;
    }				/* end of col loop */
  }				/* end of row loop */

  if (canonical) {
    if (heater || cooler) {
      beta += (cooler - heater) * .0005;
      bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
      bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));
      fixthermometer();
    }
  }
  /* clear previous cluster */
  for (row = 0; row < nrows; row++) {
    activerow[row] = 0;
    for (col = 0; col < hwords; col++)
      cluster[row][col] = 0;
  }
  /* randomize starting demon location */
  demonindex = lrand48() % ndemons;
  /* pick random site to seed cluster */
  row = lrand48() % nrows;
  col = lrand48() % hwords;
  clustergrowth = 1;
  cluster[row][col] = 1 << (lrand48() % WORDSIZE);
  activerow[row] = activerow[previousrow] = activerow[nextrow] = 1;
  /* set up bounding box for cluster */
  boxbottom = boxtop = row;
  /* monitor temperature */
  if (canonical) {
    refreshdemons(5);		/* refresh 1/10th of the demons */
  } else if (volumecount >= volume) {	/* update beta and
					 * thermometer */
    beta = .5 * log((double) (-1. + volumecount / (0.5 * energycount)));
    bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
    bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));
    fixthermometer();
    energycount = 0;
    volumecount = 0;
  }				/* end of temperature monitoring */
  /* copy lattice to window */
  XPutImage(display, window, gc, spinimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET, ncols, nrows);
  /* display the new cluster */
  XPutImage(display, window, gc, clusterimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET + 28 + nrows, ncols, nrows);
  return;
}				/* end of startcluster */

void growcluster()
/* this new version (July `97) is supposed to be more readable */
{
  int rowloop, row, col, nrow, ncol, newstuff;
  long edge, currentword, newsites;
  clustergrowth = 0;
  direction = -direction;	/* change sweeping direction */
  /* start where last activity was */
  if (direction > 0)
    row = boxtop;
  else
    row = boxbottom;
  for (rowloop = 0; rowloop < 4 * nrows; rowloop++) {
    /* multiple passes seem to help */
    if (activerow[row]) {
      activerow[row] = 0;
      nrow = nextrow;
      for (col = 0; col < hwords; col++) {
	ncol = nextcol;
	currentword = cluster[row][col];
	/* look for new stuff if part of cluster is around */
	newstuff = (0 !=
             (currentword | cluster[nrow][col] | cluster[row][ncol]));
	while (newstuff) {
	  newstuff = 0;
	  /* make sure necessary demons calculated */
	  if (unchecked[row][col])
	    check(row, col);
	  /* now put things from/to next column */
	  if (1 & (sadx[row][col])) {
	    if (1 & (currentword ^
	      ((cluster[row][ncol]) >> (WORDSIZE - 1)))) {
	      newstuff = 1;
	      cluster[row][ncol] |= LEFTBIT;
	      currentword |= 1;
	    }
	  }
	  /* put things from/to next row */
	  edge = currentword ^ (cluster[nrow][col]);
	  if ((newsites = edge & sady[row][col])) {
	    /* single "=" intended */
	    currentword |= newsites;
	    cluster[nrow][col] |= newsites;
	    newstuff = 1;
	    activerow[nrow] = 1;
	  }
	  /* find cluster edge within current word */
	  edge = (currentword ^ (currentword << 1)) & (~1);
	  while ((newsites = edge & (sadx[row][col]))) {
	    /* need to grow */
	    currentword |= newsites | (newsites >> 1);
	    edge = (currentword ^ (currentword << 1)) & (~1);
	    newstuff = 1;
	  }
	  if (newstuff) {
	    cluster[row][col] = currentword;
	    activerow[row] = activerow[previousrow] = 1;
	    clustergrowth = 1;
	  }
	}			/* end of while newstuff */
      }				/* end of column loop */
    }				/* end of if activerow */
    if (direction > 0) {	/* alternate directions */
      if (activerow[row])	/* back up */
	boxtop = row = nextrow;
      else
	row = previousrow;
    } else {
      if (activerow[row])	/* back up */
	boxbottom = row = previousrow;
      else
	row = nextrow;
    }
    if (speed)			/* speed is not an issue; draw
				 * progress */
      XPutImage(display, window, gc, clusterimage,
	0, row, LEFTOFFSET, TOPOFFSET + 28 + nrows + row, ncols, 1);
  }				/* end of row loop */
  /* display new cluster */
  XPutImage(display, window, gc, clusterimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET + 28 + nrows, ncols, nrows);
  return;
}				/* end of growcluster */

void check(row, col)
  int row, col;
{				/* calculate sadx and sady when
				 * unchecked */
  /* demons are sad if they cannot flip a bond */
  long right, bottom, carry0, carry1;
  int bit;
  long index0, index1;
  right = field[row][col] ^ ((field[row][col] << 1)
		   | (1 & (field[row][nextcol] >> (WORDSIZE - 1))));
  bottom = field[row][col] ^ field[nextrow][col];
  if (boundary) {		/* correct for antiperiodic
				 * boundaries */
    if (row == nrows - 1)
      bottom ^= (-1);
    if (col == topcol)
      right ^= 1;
  }
  index0 = (demonindex + 2 * (col + row * hwords));	/* start of xdemons */
  index1 = (index0 + 1);	/* start of ydemons */
  while (index0 >= ndemons)
    index0 -= ndemons;
  while (index1 >= ndemons)
    index1 -= ndemons;
  /* copy demons to newdemon */
  for (bit = 0; bit < DEMONBITS; bit++) {
    newdemon[bit][index0] = demon[bit][index0];
    newdemon[bit][index1] = demon[bit][index1];
  }
  /* subtract 1 from newdemon */
  carry0 = carry1 = (-1);
  for (bit = 0; bit < DEMONBITS; bit++) {
    newdemon[bit][index0] ^= carry0;
    newdemon[bit][index1] ^= carry1;
    carry0 &= newdemon[bit][index0];
    carry1 &= newdemon[bit][index1];
  }
  /* add two if neighbor antiparallel */
  for (bit = 1; bit < DEMONBITS; bit++) {
    newdemon[bit][index0] ^= right;
    newdemon[bit][index1] ^= bottom;
    right &= (~newdemon[bit][index0]);
    bottom &= (~newdemon[bit][index1]);
  }
  carry0 ^= right;
  carry1 ^= bottom;
  /* make demons sad if overflow */
  sadx[row][col] = carry0;
  sady[row][col] = carry1;
  unchecked[row][col] = 0;
  return;
}

void setup()
{
  int i, j, count;
  /* initialize "bits" used for bitcounts */
  for (i = 0; i < 256; i++) {
    j = i;
    count = j & 1;
    while (j)
      count += ((j >>= 1) & 1);
    bits[i] = count;		/* the number of set bits in i */
  }
  /* make checkermask have alternate bits set */
  checkermask = 2;
  /*
   * Do this 32 times in case we are on a 64 bit machine.  The lost
   * time otherwise is highly insignificant.
   */
  for (i = 0; i < 32; i++)
    checkermask |= (checkermask << 2);
  return;
}

void repaint()
/* this fixes the window up whenever it is uncovered */
{
  int i;
  XSetPlaneMask(display, gc, AllPlanes);
  drawbutton(quitbutton, 0, 0, BUTTONWIDTH, BUTTONHEIGHT, "quit", 1);
  drawbutton(pausebutton, 0, 0, BUTTONWIDTH, BUTTONHEIGHT,
	     "pause", 1 - 2 * paused);
  drawbutton(speedbutton, 0, 0, BOUNDWIDTH, 18, "speed", -1);
  drawbutton(speedbutton, 1 + ((10 - speed) * (BOUNDWIDTH - 2)) / 11, 1,
	     (BOUNDWIDTH - 2) / 11, 16, "", 2);
  if (canonical)
    drawbutton(canonbutton, 0, 0, font_width * 15, BUTTONHEIGHT,
	       "canonical", -1);
  else
    drawbutton(canonbutton, 0, 0, font_width * 15, BUTTONHEIGHT,
	       "microcanonical", 1);
  /* draw the algorithm buttons */
  for (i = 0; i < 2; i++)
    drawbutton(algobutton, 0, i * 18, 68, 18,
	       algotext[i], 1 - 2 * (algorithm == i));
  /* draw the heater buttons */
  for (i = 0; i < 3; i++)
    drawbutton(heatwindow, 0, i * HEATERHEIGHT / 3,
	       HEATERWIDTH, HEATERHEIGHT / 3,
	       heatertext[i], 1 - 2 * (i == (heater + 2 * cooler)));
  /* draw the boundary condition buttons */
  for (i = 0; i < 4; i++)
    drawbutton(boundwindow, 0, i * BOUNDHEIGHT / 4,
	       BOUNDWIDTH, BOUNDHEIGHT / 4,
	       boundtext[i], 1 - 2 * (i == boundary));

  /* draw thermometer */
  drawbutton(window, TLEFT - 14, TTOP - 14,
	     TWDTH + 28, TBOTTOM - TTOP + 34, "", 2);
  drawbutton(window, TLEFT - 3, TTOP - 2,
	     TWDTH + 6, TBOTTOM - TTOP + 4, "", 3);
  XSetForeground(display, gcpen, black);
  XDrawLine(display, window, gcpen, TLEFT - 3, TTOP - 2,
	    TLEFT + TWDTH + 2, TTOP - 2);
  XDrawLine(display, window, gcpen, TLEFT - 3, TTOP - 2,
	    TLEFT - 3, TBOTTOM + 2);
  XDrawLine(display, window, gcpen, TLEFT + TWDTH + 2, TBOTTOM + 2,
	    TLEFT + TWDTH + 2, TTOP - 2);
  XFillArc(display, window, gcpen,
	   TLEFT - (3 * TWDTH) / 2 - 3, TBOTTOM - TWDTH - 3,
	   4 * TWDTH + 4, 4 * TWDTH + 4, 0, 64 * 360);
  XSetForeground(display, gcpen, white);
  XFillArc(display, window, gcpen,
	   TLEFT - (3 * TWDTH) / 2 - 2, TBOTTOM - TWDTH - 2,
	   4 * TWDTH, 4 * TWDTH, 0, 64 * 360);
  XSetForeground(display, gcpen, darkcolor);
  XFillArc(display, window, gcpen,
	   TLEFT - (3 * TWDTH) / 2, TBOTTOM - TWDTH,
	   4 * TWDTH - 2, 4 * TWDTH - 2, 0, 64 * 360);
  XSetForeground(display, gcpen, spinup);
  XFillArc(display, window, gcpen,
    TLEFT - TWDTH / 2, TBOTTOM, TWDTH + 2, TWDTH + 2, 10, 64 * 250);
  XSetForeground(display, gcpen, black);

  for (i = TBOTTOM - 10; i > TTOP; i -= (THGT / 20)) {
    XDrawLine(display, window, gcpen, TLEFT - 7, i, TLEFT - 4, i);
    XDrawLine(display, window, gcpen, TLEFT + TWDTH + 3, i,
	      TLEFT + TWDTH + 6, i);
  }
  i = TBOTTOM - THGT * (-0.6 + 1.0 / log(1. + sqrt(2.)));
  XDrawLine(display, window, gcpen, TLEFT - 15, i, TLEFT + TWDTH + 15, i);
  XDrawString(display, window, gcpen, TLEFT + TWDTH + 15, i + 4, "T", 1);
  XDrawString(display, window, gcpen, TLEFT + TWDTH + 23, i + 8, "c", 1);
  fixthermometer();

  /* write various strings */
  sprintf(stringbuffer, "%d by %d lattice", ncols, nrows);
  XDrawString(display, window, gcpen, 150, 92,
	      stringbuffer, strlen(stringbuffer));
  XDrawString(display, window, gcpen, 20, 92, "beta=", 5);
  XDrawString(display, window, gcpen, LEFTOFFSET, 114, "spins:", 6);
  XDrawString(display, window, gcpen, LEFTOFFSET, 142 + nrows, "changes:", 8);
  XDrawString(display, window, gcpen, TLEFT - 8, TBOTTOM + 40, "MJC", 3);

  /* draw border and redraw images */
  drawbutton(window, LEFTOFFSET - 4, TOPOFFSET - 4,
	     ncols + 8, nrows + 8, NULL, 4);
  drawbutton(window, LEFTOFFSET - 4, TOPOFFSET + 28 + nrows - 4,
	     ncols + 8, nrows + 8, NULL, 4);
  XPutImage(display, window, gc, spinimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET, ncols, nrows);
  XPutImage(display, window, gc, clusterimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET + 28 + nrows, ncols, nrows);

  /* this might speed things up a little bit */
  XSetPlaneMask(display, gc, spinup ^ spindown);
  return;
}

void refreshdemons(int step)
/*
 * refresh cluster demons at beta note: Calling this after each
 * cluster update gives the canonical algorithm.
 */

/*
 * To generate a bit set with probability p: write p as a binary
 * fraction.  Compare a random bit string with the digits of this
 * fraction and find the first place at which they are the same.
 * Output the bit  that occupies that place.  This is done here in
 * parallel on WORDSIZE demons at a time.
 */
{
  long morebits, randombits, bitprob, currentbit;
  double bfactor, prob;
  long i;
  int bit, j;
  bfactor = exp(-beta);
  for (bit = 0; bit < DEMONBITS; bit++) {
    bfactor *= bfactor;		/* gives exp(-2 beta) for the first
				 * bit, exp(-4 beta) for the next,
				 * ... */
    prob = bfactor / (1. + bfactor);	/* prob that current demon
					 * bit is 1 */
    bitprob = ((1 << 16) * prob);	/* convert to 16 bit integer */

    for (i = 0; i < ndemons; i += step) {
      demon[bit][i] = 0;
      morebits = (-1);
      for (j = 0; j < 16; j++) {
	randombits = RND;	/* new set of random bits */
	currentbit = -((bitprob >> (15 - j)) & 1);
	/*
	 * set demon to current bit if it equals randombits and has
	 * not been fixed yet
	 */
	demon[bit][i] |= (currentbit & randombits & morebits);
	/* shut off morebits where they are equal */
	morebits &= (currentbit ^ randombits);
	if (0 == morebits)
	  break;		/* all bits decided */
      }
    }
  }
  return;
}

void heatbath()
/*
 * refresh local demons at beta
 */
/*
 * To generate a bit set with probability p: write p as a binary
 * fraction.  Compare a random bit string with the digits of this
 * fraction and find the first place at which they are the same.
 * Output the bit  that occupies that place.  This is done here in
 * parallel on WORDSIZE demons at a time.
 */
{
  long morebits, randombits, currentbit;
  int j;
  demon0 = demon1 = 0;
  morebits = (-1);
  for (j = 15; (j >= 0) && (morebits != 0); j--) {
    randombits = RND;		/* new set of random bits */
    currentbit = -((bitprob0 >> j) & 1);
    /* current bit is a word filled with the j'th bit of bitprob0 */
    /*
     * set demon to current bit if it equals randombits and has not
     * been fixed yet
     */
    demon0 |= (currentbit & randombits & morebits);
    /* shut off morebits where they are equal */
    morebits &= (currentbit ^ randombits);
    if (0 == morebits)
      break;			/* all bits decided */
  }
  /* now do demon1 */
  morebits = (-1);
  for (j = 15; (j >= 0) && (morebits != 0); j--) {
    randombits = RND;
    currentbit = -((bitprob1 >> j) & 1);
    demon1 |= (currentbit & randombits & morebits);
    morebits &= (currentbit ^ randombits);
    if (0 == morebits)
      break;
  }
  return;
}

void localupdate()
/* the local updating routine */
{
  long *up, *down;		/* pointers to above and below rows */
  long temp0, temp1, reject, spins, left, right, top, bottom;
  int row = 0, col, rowloop;
  checkermask = (~checkermask);
  for (rowloop = 0; rowloop < nrows; rowloop++) {
    row = (row + 29);		/* jump ahead 29 rows; allows for
				 * energy to move long  distances */
    while (row >= nrows)
      row -= nrows;
    /* get neighbors */
    up = field[previousrow];
    down = field[nextrow];
    leftshift(field[row], work0);
    rightshift(field[row], work1);
    /* to keep speed up, only refresh occasionally */
    if ((canonical) && ((row & 3) == 3))
      heatbath();
    for (col = 0; col < hwords; col++) {
      spins = field[row][col];
      /* determine antiparallel neighbors */
      right = work0[col] ^ spins;
      left = work1[col] ^ spins;
      top = up[col] ^ spins;
      bottom = down[col] ^ spins;
      if (0 == row)
	switch (boundary) {
	case 1:
	  top ^= (~0);
	  break;
	case 2:
	  top = (~0) ^ spins;
	  break;
	case 3:
	  top = spins;
	  break;
	}
      else if (nrows - 1 == row)
	switch (boundary) {
	case 1:
	  bottom ^= (~0);
	  break;
	case 2:
	  bottom = (~0) ^ spins;
	  break;
	case 3:
	  bottom = spins;
	  break;
	}
      /* add up top, bottom, and left antiparallel spins */
      temp0 = top ^ bottom ^ left;
      temp1 = (top & (bottom | left)) | (bottom & left);
      /* add in right */
      reject = temp0 & temp1 & right;
      temp1 ^= (temp0 & right);
      temp0 ^= right;
      /* add  demon0 */
      reject ^= (temp0 & temp1 & demon0);
      temp1 ^= (temp0 & demon0);
      temp0 ^= demon0;
      /* add in demon1 */
      reject ^= (temp1 & demon1);
      temp1 ^= demon1;
      /* subtract two */
      temp1 = ~temp1;
      reject ^= temp1;
      /* only change one color sites */
      reject |= checkermask;
      /* make changes */
      demon0 = (demon0 & reject) | (temp0 & (~reject));
      demon1 = (demon1 & reject) | (temp1 & (~reject));
      field[row][col] = (spins ^ (~reject));
      /* put changes into cluster image */
      cluster[row][col] &= checkermask;
      cluster[row][col] |= (~reject);
      energycount += bitcount(demon0 & (~checkermask));
    }
    /* heat or cool */
    if (!canonical)
      if (heater)
	demon0 |= (RND & RND & RND & RND & RND & RND);
      else if (cooler)
	demon1 &= (RND | RND | RND | RND);
  }
  /* heat or cool for canonical case */
  if (canonical && (heater || cooler)) {
    beta += (cooler - heater) * .000005 * nrows;
    bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
    bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));
    fixthermometer();
  }
  /* monitor temperature every 50 "half" iterations */
  iteration++;
  if (iteration == 50) {	/* update beta and thermometer */
    if (!canonical) {
      beta = .25 * log((double) (-1. + 50 * volume / (2. * energycount)));
      bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
      bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));
      fixthermometer();
    }
    iteration = 0;
    energycount = 0;
  }
  /* copy lattice to window */
  XPutImage(display, window, gc, spinimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET, ncols, nrows);
  XPutImage(display, window, gc, clusterimage, 0, 0,
	    LEFTOFFSET, TOPOFFSET + 28 + nrows, ncols, nrows);
  return;
}

void fixthermometer()
/* fill the thermometer with fluid */
{
  int temp;
  sprintf(stringbuffer, "%.3f", beta);
  XDrawImageString(display, window, gcpen, LEFTOFFSET, 92,
		   stringbuffer, strlen(stringbuffer));
  temp = THGT * (-0.6 + 0.5 / beta);
  if (temp < 0)
    temp = 0;
  if (temp > THGT)
    temp = THGT;
  XSetForeground(display, gcpen, darkcolor);
  XFillRectangle(display, window, gcpen, TLEFT, TBOTTOM - temp, TWDTH, temp);
  XFillRectangle(display, window, gcclear, TLEFT, TTOP, TWDTH, THGT - temp);
  XSetForeground(display, gcpen, black);
  return;
}

void rightshift(source, destination)
/* shifts a row of cells to the right by one bit */
  long *source, *destination;
{
  int col;
  long carry = 0;
  for (col = bottomcol; col != (topcol + nshift); col += nshift) {
    destination[col] = (~LEFTBIT) & (source[col] >> 1);
    if (carry)
      destination[col] |= LEFTBIT;
    carry = source[col] & 1;
  }
  if (carry)
    destination[bottomcol] |= LEFTBIT;
  switch (boundary) {
  case 1:
    destination[bottomcol] ^= LEFTBIT;
    break;
  case 2:
    destination[bottomcol] |= LEFTBIT;
    break;
  case 3:
    destination[bottomcol] &= (~LEFTBIT);
    break;
  }
  return;
}

void leftshift(source, destination)
/* shifts a row of cells one bit to the left */
  long *source, *destination;
{
  int col;
  long carry = 0;
  for (col = topcol; col != (bottomcol - nshift); col -= nshift) {
    destination[col] = (source[col] << 1);
    if (carry)
      destination[col] |= 1;
    carry = source[col] >> (WORDSIZE - 1);
  }
  if (carry)
    destination[topcol] |= 1;
  switch (boundary) {
  case 1:
    destination[topcol] ^= 1;
    break;
  case 2:
    destination[topcol] |= 1;
    break;
  case 3:
    destination[topcol] &= (~1);
    break;
  }
  return;
}

void setalg()
/* take preparatory action if algorithm button hit */
{
  int row, col;
  energycount = volumecount = iteration = 0;
  if (1 == algorithm) {		/* reset cluster stuff */
    if (boundary > 1) {
      setboundary(0);
    }
    refreshdemons(1);
    clustergrowth = 0;
    for (col = 0; col < hwords; col++)
      for (row = 0; row < nrows; row++) {
	cluster[row][col] = 0;
	unchecked[row][col] = 1;
      }
  }
  return;
}

void setheater(value)
/* acts on heater buttons */
  int value;
{
  switch (value) {
  case 0:
    heater = cooler = 0;
    break;
  case 1:
    heater = 1 - heater;
    break;
  case 2:
    cooler = 1 - cooler;
  }
  drawbutton(heatwindow, 0, 0,
	     HEATERWIDTH, HEATERHEIGHT / 3, heatertext[0],
	     2 * (heater | cooler) - 1);
  drawbutton(heatwindow, 0, (HEATERHEIGHT) / 3,
      HEATERWIDTH, HEATERHEIGHT / 3, heatertext[1], 1 - 2 * heater);
  drawbutton(heatwindow, 0, (HEATERHEIGHT * 2) / 3,
      HEATERWIDTH, HEATERHEIGHT / 3, heatertext[2], 1 - 2 * cooler);
  return;
}

void setboundary(value)
/* take action if boundary buttons hit */
  int value;
{
  drawbutton(boundwindow, 0, boundary * BOUNDHEIGHT / 4,
	     BOUNDWIDTH, BOUNDHEIGHT / 4, boundtext[boundary], 1);
  /* don't allow fixed boundaries with cluster algorithm */
  if ((algorithm == 1) && (value > 1))
    fprintf(stderr,
	"fixed boundaries not implemented for cluster algorithm\n");
  else
    boundary = value;
  drawbutton(boundwindow, 0, boundary * BOUNDHEIGHT / 4,
	     BOUNDWIDTH, BOUNDHEIGHT / 4, boundtext[boundary], -1);
  return;
}

void openwindow(argc, argv)
/*
 * a lot of this is taken from the basicwin program in the Xlib
 * Programming Manual
 */
  int argc;
  char **argv;
{
  char *window_name = "2-d Ising Model";
  char *icon_name = "Ising";
  Pixmap icon_pixmap;
  char *display_name = NULL;
  XColor xcolor, colorcell;
  Colormap cmap;
  /* icon */
#define icon_bitmap_width 16
#define icon_bitmap_height 16
  static char icon_bitmap_bits[] = {
    0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0xf8,
    0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

  /* open up the display */
  if ((display = XOpenDisplay(display_name)) == NULL) {
    fprintf(stderr, "%s: cannot connect to X server %s\n",
	    progname, XDisplayName(display_name));
    exit(-1);
  }
  screen = DefaultScreen(display);
  depth = DefaultDepth(display, screen);
  cmap = DefaultColormap(display, screen);
  /* depth=1; *//* uncomment to test monochrome display */
  spindown = darkcolor = black = BlackPixel(display, screen);
  spinup = lightcolor = white = WhitePixel(display, screen);
  if (depth > 1) {		/* color? This is not the right way
				 * to do it, but ... */
    if (XAllocNamedColor(display, cmap, "coral", &colorcell, &xcolor))
      darkcolor = colorcell.pixel;
    if (XAllocNamedColor(display, cmap, "turquoise", &colorcell, &xcolor))
      lightcolor = colorcell.pixel;
    if (XAllocNamedColor(display, cmap, "grey", &colorcell, &xcolor))
      spinup = colorcell.pixel;
    if (XAllocNamedColor(display, cmap, "black", &colorcell, &xcolor))
      spindown = colorcell.pixel;
  }
  /* make the main window */
  windowwidth = ncols + 75;
  windowheight = 2 * nrows + 164 + 24;
  window = XCreateSimpleWindow(display, RootWindow(display, screen),
    0, 0, windowwidth, windowheight, 4, BlackPixel(display, screen),
			       lightcolor);

  /* make the icon */
  icon_pixmap = XCreateBitmapFromData(display, window,
				icon_bitmap_bits, icon_bitmap_width,
				      icon_bitmap_height);

  size_hints.flags = PPosition | PSize | PMinSize;
  size_hints.min_width = windowwidth;
  size_hints.min_height = 210;
#ifdef X11R3
  size_hints.x = 0;
  size_hints.y = 0;
  size_hints.width = windowwidth;
  size_hints.height = windowheight;
  XSetStandardProperties(display, window, window_name, icon_name,
			 icon_pixmap, argv, argc, &size_hints);
#else
  {
    XWMHints wm_hints;
    XClassHint class_hints;
    XTextProperty windowName, iconName;
    if (XStringListToTextProperty(&window_name, 1, &windowName) == 0) {
      fprintf(stderr, "%s: structure allocation for windowName failed.\n"
	      ,progname);
      exit(-1);
    }
    if (XStringListToTextProperty(&icon_name, 1, &iconName) == 0) {
      fprintf(stderr, "%s: structure allocation for iconName failed.\n"
	      ,progname);
      exit(-1);
    }
    wm_hints.initial_state = NormalState;
    wm_hints.input = True;
    wm_hints.icon_pixmap = icon_pixmap;
    wm_hints.flags = StateHint | IconPixmapHint | InputHint;
    class_hints.res_name = progname;
    class_hints.res_class = "Basicwin";
    XSetWMProperties(display, window, &windowName, &iconName,
		  argv, argc, &size_hints, &wm_hints, &class_hints);
  }
#endif
  /* make a default cursor */
  XDefineCursor(display,window, XCreateFontCursor(display,XC_hand2));

  /* pick the events to look for */
  event_mask = ExposureMask | ButtonPressMask | StructureNotifyMask;
  XSelectInput(display, window, event_mask);
  /* pick font: 9x15 is supposed to almost always be there */
  if ((font = XLoadQueryFont(display, "9x15")) == NULL) {
    fprintf(stderr, "%s: Cannot open 9x15 font\n", progname);
    exit(-1);
  }
  font_height = font->ascent + font->descent;
  font_width = font->max_bounds.width;
  /*
   * make graphics contexts: gc for black on white, gcpen for
   * various, defaults to black on lightcolor, gclear does the
   * reverse, clears thermometer top
   */

  gc = XCreateGC(display, window, 0, NULL);
  XSetFont(display, gc, font->fid);
  XSetForeground(display, gc, spinup);
  XSetBackground(display, gc, spindown);

  gcpen = XCreateGC(display, window, 0, NULL);
  XSetFont(display, gcpen, font->fid);
  XSetForeground(display, gcpen, black);
  XSetBackground(display, gcpen, lightcolor);

  gcclear = XCreateGC(display, window, 0, NULL);
  XSetFont(display, gcclear, font->fid);
  XSetForeground(display, gcclear, spinup);
  XSetBackground(display, gcclear, black);

  /* show the window */
  XMapWindow(display, window);
  return;
}

void allocarray(array, xsize, ysize, datasize)
  char ***array;
  int xsize, ysize, datasize;
/*
 * dynamically allocates a two dimensional (*array)[xize][ysize] with
 * datasize bytes per element
 */
{
  int i;
  if ((*array) != ((char **) NULL))	/* free previous allocation */
    free(*array);
  /* allocate memory for array of xsize pointers plus data */
  *array = (char **) malloc(sizeof(char *) * xsize + datasize * xsize * ysize);
  if (*array == (char **) NULL) {
    fprintf(stderr, "%s: memory allocation problems\n", progname);
    exit(-1);
  }
  /* initialize the pointers */
  for (i = 0; i < xsize; i++)
    *((*array) + i) = ((char *) (*array + xsize) + (i * ysize * datasize));
  return;
}

void makebuttons()
/*
 * creates the buttons and allocates arrays; call after any resizing
 */
{
  int row, col;
  XDestroySubwindows(display, window);
  XSetForeground(display, gcpen, lightcolor);
  XFillRectangle(display, window, gcpen, 0, 0, windowwidth, windowheight);
  quitbutton = makebutton(4, 4, BUTTONWIDTH, BUTTONHEIGHT);
  pausebutton = makebutton(58, 4, BUTTONWIDTH, BUTTONHEIGHT);
  boundwindow = makebutton(windowwidth - BOUNDWIDTH - 4, 4,
			   BOUNDWIDTH, BOUNDHEIGHT);
  heatwindow = makebutton(windowwidth - HEATERWIDTH - BOUNDWIDTH - 12, 4,
			  HEATERWIDTH, HEATERHEIGHT);
  algobutton = makebutton(18, 36, 68, 36);
  speedbutton = makebutton(windowwidth - BOUNDWIDTH - 4,
			   BOUNDHEIGHT + 26, BOUNDWIDTH, 18);
  canonbutton = makebutton(LEFTOFFSET + 4 + (ncols - font_width * 15) / 2,
		  windowheight - 26, font_width * 15, BUTTONHEIGHT);

  /* reset various lattice dimensions */
  hwords = ncols / WORDSIZE;
  ncols = hwords * WORDSIZE;	/* truncate ncols to a multiple of
				 * wordsize */
  if (0 == (nrows % 29))
    nrows--;			/* so demon jump of 29 rows works ok */
  volume = nrows * ncols;
  ndemons = 2 * nrows * hwords;

  /* allocate various arrays */
  if (NULL != spinimage) {
    (*spinimage).data = NULL;
    XDestroyImage(spinimage);
  }
  if (NULL != clusterimage) {
    (*clusterimage).data = NULL;
    XDestroyImage(clusterimage);
  }
  allocarray(&field, nrows, hwords, sizeof(long));
  allocarray(&cluster, nrows, hwords, sizeof(long));
  allocarray(&sadx, nrows, hwords, sizeof(long));
  allocarray(&sady, nrows, hwords, sizeof(long));
  allocarray(&unchecked, nrows, hwords, sizeof(long));
  allocarray(&demon, DEMONBITS, ndemons, sizeof(long));
  allocarray(&newdemon, DEMONBITS, ndemons, sizeof(long));
  if (NULL != work0)
    free(work0);
  if (NULL != work1)
    free(work1);
  if (NULL != activerow)
    free(activerow);
  work0 = (long *) malloc(sizeof(long) * hwords);
  work1 = (long *) malloc(sizeof(long) * hwords);
  activerow = (int *) malloc(sizeof(int) * nrows);
  if ((NULL == work0) || (NULL == work1)) {
    fprintf(stderr, "%s: memory allocation problems\n", progname);
    exit(-1);
  }
  energycount = volumecount = iteration = 0;

  /* make image structures */
  spinimage = XCreateImage(display, (Visual *) & window, 1, XYBitmap, 0,
			   (char *) *field, ncols, nrows, 32, 0);
  clusterimage = XCreateImage(display, (Visual *) & window, 1, XYBitmap, 0,
			    (char *) *cluster, ncols, nrows, 32, 0);
  if ((NULL == spinimage) || (NULL == clusterimage)) {
    fprintf(stderr, "%s: memory allocation problems\n",
	    progname);
    exit(-1);
  }
  (*spinimage).bitmap_unit = WORDSIZE;
  (*clusterimage).bitmap_unit = WORDSIZE;

  /* test for byte order on client */
  **field = 1;
  if (*(char *) (*field)) {
    (*spinimage).byte_order = LSBFirst;
    (*clusterimage).byte_order = LSBFirst;
    /* printf("byte_order=LSBFirst\n"); */
  } else {
    (*spinimage).byte_order = MSBFirst;
    (*clusterimage).byte_order = MSBFirst;
    /* printf("byte_order=MSBFirst\n"); */
  }
  /* find out in what order bits are put on screen */
  **field = 0xff;
  if (XGetPixel(spinimage, 0, 0)) {	/* bits are backward */
    nshift = (-1);
    topcol = 0;
    bottomcol = hwords - 1;
    /* printf("bitmap_bit_order=LSBFirst\n"); */
  } else {
    nshift = 1;
    topcol = hwords - 1;
    bottomcol = 0;
  }

  /* set initial state partially random and clear cluster array */
  for (col = 0; col < hwords; col++)
    for (row = 0; row < nrows; row++) {
      field[row][col] = ~(RND * (drand48() > .77));
      cluster[row][col] = 0;
      unchecked[row][col] = 1;
    }
  /* equilibrate initial demons at critical temperature */
  beta = 0.5 * log(1.0 + sqrt(2.0));	/* critical coupling */
  bitprob0 = (1 << 16) * exp(-4 * beta) / (1. + exp(-4 * beta));
  bitprob1 = (1 << 16) * exp(-8 * beta) / (1. + exp(-8 * beta));
  refreshdemons(1);
  heatbath();
  return;
}

Window makebutton(xoffset, yoffset, xsize, ysize)
  int xoffset, yoffset, xsize, ysize;
/*
 * Puts button of specified dimensions on window.  Nothing is drawn.
 */
{
  Window buttonwindow;
  long event_mask;
  Cursor cursor;
  buttonwindow = XCreateSimpleWindow(display, window, xoffset, yoffset,
				xsize, ysize, 0, black, lightcolor);
  event_mask = ButtonPressMask | ExposureMask;	/* look for mouse-button
						 * presses */
  XSelectInput(display, buttonwindow, event_mask);
  /* use a hand cursor for buttons */
  cursor=XCreateFontCursor(display,XC_hand2);
  XDefineCursor(display,buttonwindow,cursor);
  XMapWindow(display, buttonwindow);
  return buttonwindow;
}

void drawbutton(buttonwindow, xoffset, yoffset, xsize, ysize, text, state)
  Window buttonwindow;
  int xoffset, yoffset, xsize, ysize, state;
  char *text;
/*
 * Draw a button in buttonwindow of specified dimensions with text
 * centered.  Color of button determined by sign of "state." size of
 * border determined by magnitude.
 */
{
  int textlength, i, j;
  int cdark, clight, cup, cdown;
  int cleft, cright, cbutton, ctext;
  cup = lightcolor;
  cdown = darkcolor;
  cdark = black;
  clight = white;
  if (state < 0) {
    cbutton = cdown;
    ctext = clight;
    cleft = cdark;
    cright = clight;
  } else {
    cbutton = cup;
    ctext = cdark;
    cleft = clight;
    cright = cdark;
  }
  if (1 == depth)
    clight = cleft = cright = black;
  j = abs(state);
  XSetForeground(display, gcpen, cbutton);
  XFillRectangle(display, buttonwindow, gcpen, xoffset + j, yoffset + j,
		 xsize - 2 * j, ysize - 2 * j);
  XSetForeground(display, gcpen, cleft);
  XFillRectangle(display, buttonwindow, gcpen, xoffset, yoffset, xsize, j);
  XFillRectangle(display, buttonwindow, gcpen, xoffset, yoffset, j, ysize);
  XSetForeground(display, gcpen, cright);
  for (i = 0; i < j; i++) {
    XDrawLine(display, buttonwindow, gcpen,
	      xoffset + i, yoffset + ysize - i - 1,
	      xoffset + xsize - i - 1, yoffset + ysize - i - 1);
    XDrawLine(display, buttonwindow, gcpen,
	      xoffset + xsize - i - 1, yoffset + i,
	      xoffset + xsize - i - 1, yoffset + ysize - i - 1);
  }
  if (NULL != text) {
    textlength = strlen(text);
    XSetForeground(display, gcpen, ctext);
    XDrawString(display, buttonwindow, gcpen,
		xoffset + (xsize - textlength * font_width) / 2,
		yoffset + (ysize + font_height) / 2 - 2,
		text, textlength);
  }
  XSetForeground(display, gcpen, darkcolor);
  return;
}
