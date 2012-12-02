#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <X11/Xlib.h>

#include "boid.h"
#include "vec.h"

extern Display *disp;
extern Window Root; 
extern Window win;
extern Visual *vis; 
extern int scr; 
extern unsigned int depth;
extern GC fg,bg;
extern int blk,wht;
extern int width,height;

extern Colormap cmap;
extern XColor darkgray;
extern XColor yellow;
extern XColor blue;
extern XColor outlines[32];

void draw_boid(Boid boid, Pixmap freshmap)
{
   int i;
   XPoint lwing[3], rwing[3];
   int outline_index;
   
   XSetForeground(disp, fg, darkgray.pixel);
   XFillPolygon(disp, freshmap, fg, boid->shadow, 4,
		Nonconvex, CoordModeOrigin);
   
   
   lwing[0].x = boid->X;
   lwing[1].x = boid->tail_lX;
   lwing[2].x = boid->tail_X;
   
   lwing[0].y = boid->Y;
   lwing[1].y = boid->tail_lY;
   lwing[2].y = boid->tail_Y;
   
   rwing[0].x = boid->X;
   rwing[1].x = boid->tail_rX;
   rwing[2].x = boid->tail_X;
   
   rwing[0].y = boid->Y;
   rwing[1].y = boid->tail_rY;
   rwing[2].y = boid->tail_Y;
   
   outline_index = ((int)boid->pos->z) >> 10;
   
   /* if moving right => lwing behind rwing */
   if(boid->vel->x > 0) {
      if(boid->tail_lY < boid->Y) {
	 XSetForeground(disp, fg, yellow.pixel);
      }
      else {
	 XSetForeground(disp, fg, blue.pixel);
      }
      XFillPolygon(disp, freshmap, fg, lwing, 3,
		   Convex, CoordModeOrigin);
      XSetForeground(disp, fg, outlines[outline_index].pixel);
      XDrawLines(disp, freshmap, fg, lwing, 3,
		 CoordModeOrigin);
      
      if(boid->tail_rY < boid->Y) {
	 XSetForeground(disp, fg, blue.pixel);
      }
      else {
	 XSetForeground(disp, fg, yellow.pixel);
      }
      XFillPolygon(disp, freshmap, fg, rwing, 3,
		   Convex, CoordModeOrigin);       
      XSetForeground(disp, fg, outlines[outline_index].pixel);
      XDrawLines(disp, freshmap, fg, rwing, 3,
		 CoordModeOrigin);
   }
   else {
      if(boid->tail_rY < boid->Y) {
	 XSetForeground(disp, fg, yellow.pixel);
      }
      else {
	 XSetForeground(disp, fg, blue.pixel);
      }
      XFillPolygon(disp, freshmap, fg, rwing, 3,
		   Convex, CoordModeOrigin);
      XSetForeground(disp, fg, outlines[outline_index].pixel);
      XDrawLines(disp, freshmap, fg, rwing, 3,
		 CoordModeOrigin);
      
      if(boid->tail_lY < boid->Y) {
	 XSetForeground(disp, fg, blue.pixel);
      }
      else {
	 XSetForeground(disp, fg, yellow.pixel);
      }
      XFillPolygon(disp, freshmap, fg, lwing, 3,
		   Convex, CoordModeOrigin);
      XSetForeground(disp, fg, outlines[outline_index].pixel);
      XDrawLines(disp, freshmap, fg, lwing, 3,
		 CoordModeOrigin);
      
   }
   
}   

