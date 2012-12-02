#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <X11/Xlib.h>

#include "boid.h"
#include "vec.h"

Display *disp;
Window Root;
Window win;
Visual *vis;
int scr;
unsigned int depth;
GC fg,bg;
int blk,wht;
int width,height;

Colormap cmap;
XColor darkgray;
XColor yellow;
XColor blue;
XColor outlines[32];

#define GRATUITOUS

main(int argc, char **argv)
{
   Pixmap background, freshmap;
   int numboids=30;
   Boid *boids;
   int i;
   Vec center, avg_velocity;
   int root=0;
   char *filename;
   int rx,ry,rw,rh;
   int minx,miny,maxx,maxy;

   srandom(time(NULL));

   SetupDisplay();

   for(i=1; i<argc; i++) {
      if(!strcmp(argv[i], (char *)"-root")) {
	 root=1;
      }else if(!strcmp(argv[i], (char *)"-win")) {
	 root=0;
      }else if(!strcmp(argv[i], (char *)"-numboids")) {
	 numboids=atoi(argv[++i]);
      }
   }

   if(!root) {
      width=640;
      height=480;
      OpenWindow();
   }else {
      width=1152;
      height=900;
      win=Root;
   }

   SetupColormap();

   background=XCreatePixmap(disp, win, width, height, depth);
   freshmap=XCreatePixmap(disp, win, width, height, depth);

   boids = (Boid *)calloc(numboids, sizeof(_Boid));
   for(i=0; i<numboids; i++) {
      boids[i] = new_boid(width, height);
   }

   center = zero_vec();
   avg_velocity = zero_vec();

   XSetWindowBackgroundPixmap(disp,win,freshmap);
   for(;;) {

#if 1
     XCopyArea(disp, background, freshmap, fg, 0, 0, width, height, 0, 0);
#else
     XClearArea(disp, freshmap, 0, 0, width, height, False);
#endif

      minx=INT_MAX;miny=INT_MAX;maxx=0;maxy=0;
      vec_clear(center);
      vec_clear(avg_velocity);

      for(i=0; i<numboids; i++) {
	 vec_add(center, boids[i]->pos);
	 vec_add(avg_velocity, boids[i]->vel);
      }

      for(i=0; i<numboids; i++) {
#ifndef GRATUITOUS
	 if(boids[i]->onscreen) {
	    if (boids[i]->X<minx) minx=boids[i]->X;
	    if (boids[i]->Y<miny) miny=boids[i]->Y;
	    if (boids[i]->X>maxx) maxx=boids[i]->X;
	    if (boids[i]->Y>maxy) maxy=boids[i]->Y;
	 }
#endif

	 boid_move(boids[i], boids, numboids,
		   center, avg_velocity,
		   width, height);

	 if(boids[i]->onscreen) {
	    draw_boid(boids[i], freshmap);
#ifndef GRATUITOUS
	    if (boids[i]->X<minx) minx=boids[i]->X;
	    if (boids[i]->Y<miny) miny=boids[i]->Y;
	    if (boids[i]->X>maxx) maxx=boids[i]->X;
	    if (boids[i]->Y>maxy) maxy=boids[i]->Y;
#endif
	 }
      }

#ifdef GRATUITOUS
      XClearWindow(disp,win);
#else
      rx=minx;ry=miny;rw=maxx-minx;rh=maxy-miny;
      if (rx<0) rx=0;
      if (ry<0) ry=0;
      if (rw>width) rw=width;
      if (rh>height) rh=height;
      XClearArea(disp,win,rx,ry,rw,rh,False);
      XSync(disp,False);
      XCopyArea(disp, background, freshmap, fg, rx, ry, rw, rh, rx, ry);
#endif

      XFlush(disp);

      usleep(1000000/30);
   }

}
