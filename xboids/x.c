#include <stdio.h>
#include <stdlib.h>
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

void SetupDisplay(void)
{
   XGCValues gcvals;

   if(!(disp=XOpenDisplay(NULL))) {
      perror("Cannot open display\n");
      exit(1);
   }

   Root=DefaultRootWindow(disp);
   scr=DefaultScreen(disp);
   depth=DefaultDepth(disp,scr);
   blk=BlackPixel(disp,scr);
   wht=WhitePixel(disp,scr);
   vis=DefaultVisual(disp,scr);
   fg=XCreateGC(disp,Root,(unsigned long)0,&gcvals);
   bg=XCreateGC(disp,Root,(unsigned long)0,&gcvals);
   XSetForeground(disp,fg,wht);
   XSetBackground(disp,bg,blk);
}

void OpenWindow(void)
{
   win=XCreateSimpleWindow(disp,Root,0,0,width+2,height+2,0,blk,blk);
   XMapWindow(disp,win);
   XSync(disp,False);
}

void SetupColormap(void)
{
   int i;

   cmap = DefaultColormap(disp, scr);
   darkgray.red = 0x4000;
   darkgray.green=0x4000;
   darkgray.blue=0x4000;
   XAllocColor(disp,cmap,&darkgray);
   yellow.red = 0xf000;
   yellow.green=0xf000;
   yellow.blue=0x0;
   XAllocColor(disp,cmap,&yellow);
   blue.red=0x0;
   blue.green=0x0;
   blue.blue=0xf000;
   XAllocColor(disp,cmap,&blue);
   for(i=0;i<32;i++) {
      outlines[i].red=0xfa00 - (i<<10);
      outlines[i].green=0xa000 - (i<<10);
      outlines[i].blue=0x7f00 - (i<<10);
      XAllocColor(disp,cmap,&outlines[i]);
   }
}
