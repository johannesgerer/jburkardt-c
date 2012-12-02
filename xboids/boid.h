#ifndef __BOID_H__
#define __BOID_H__

#include <X11/Xlib.h>

#include "vec.h"

typedef struct {
   Vec pos, vel;
   int X,Y;
   int tail_lX, tail_lY, tail_rX, tail_rY, tail_X, tail_Y;
   XPoint shadow[4];
   int wing_level;
   int onscreen; /* boolean */
   int upstroke; /* boolean */
   int perching; /* boolean */
   int perch_timer;
} _Boid, *Boid;

Boid new_boid(int W, int H);
void boid_perspective(Boid boid, int W, int H);
int boid_isonscreen(Boid boid, int W, int H);
void boid_move(Boid boid, Boid allboids[], int numboids,
	       Vec real_center, Vec real_avgvel,
	       int W, int H);
Vec boid_perceive_center(Boid boid, Vec real_cent, int numboids);
Vec boid_av_vel(Boid boid, Vec real_avgvel, int numboids);
Vec boid_chill_out(Boid boid, Boid boids[], int numboids);
   
#endif  /*  __BOID_H__  */
