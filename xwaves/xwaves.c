/* Michael Creutz          */
/* creutz@wind.phy.bnl.gov */
/* to compile:             */
/* cc -O fourier.c -lX11 -lm     */
/* version of August, 1995 */
/* January 2000: corrected array bounds on xpoints */
 
# include <X11/Xlib.h>
# include <X11/Xutil.h>
# include <X11/Xos.h>
# include <X11/Xatom.h> 
# include <X11/cursorfont.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <unistd.h>
# include <math.h>
      
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
struct timeval timeout;

/* some colors to use */
#define WHITE 0
#define BLACK 1
#define RED 2
#define BLUE 3
#define GREEN 4
#define YELLOW 5
#define MAGENTA 6
#define CYAN 7
#define ORANGE 8
#define PURPLE 9
#define PINK 10
#define BROWN 11
#define GREY 12
#define TURQUOISE 13
#define GOLD 14
#define NAVY 15
#define TAN 16
#define VIOLET 17

/* NSITES  = number of lattice sites */
# define NSITES 256

double phi[NSITES],cpi[NSITES];    /* position space field and momentum     */
double fphi[NSITES], fcpi[NSITES]; /* fourier space momentum                */
double cn[NSITES],sn[NSITES];      /* for cosine and sin of angles          */
double omega[NSITES];              /* frequency for given component         */
double comega[NSITES],somega[NSITES];/* cos and sin of dt*omega             */
double pi;                         /* 3.14....                              */
int xcenter,ycenter,vx,vy;         /* center of drifting blob               */
XPoint xpoints[NSITES+7];          /* for drawing polygons                  */
# define dt (.02*(1-.8*slow))

void drawbutton(),openwindow(),makebuttons(),update(),repaint(),
     cleanup(),drawit(),drawwave(),drawblob(),fastfourier(),
     xtof(),ftox(),makedispersion();

/* dimensions for control boxes */
# define BUTTONWIDTH 72
# define BUTTONHEIGHT 18

/* initial button settings */
int paused=0;             /* is the system paused?    */
int damped=0;             /* is the system damped?    */
int reversed=0;           /* are we going backwards?  */
int mass=1;               /* does the wave have mass? (0 or 1 only) */
int blob=0;               /* display as blob?         */
int drift=0;              /* let blob drift           */
int slow=0;               /* run at slow speed        */

int equation=0;           /* 0 for light, 1 for mesons, 2 for water */
# define WATER 2
char * equtext[3]={"light","mesons","water"};

char stringbuffer[256];   /* generally useful         */
 /* to simplify typing */
# define line(x1,y1,x2,y2) XDrawLine(display,window,gc1,x1,y1,x2,y2)
# define line2(x1,y1,x2,y2) XDrawLine(display,window,gc2,x1,y1,x2,y2)
# define color(i) XSetForeground(display,gc1,translate[i])
# define color2(i) XSetForeground(display,gc2,translate[i])

/* various window stuff */
Display *display;
int screen;
static char *progname;
Window window,quitbutton,pausebutton,vacuumbutton,rightbutton,dampbutton,
   reversebutton,blobbutton,driftbutton,slowbutton,packetbutton,
   equbutton,makebutton();
XColor xcolor,colorcell;
Colormap cmap;
GC gc1,gc2;
int windowwidth,windowheight,playwidth,playheight;
XFontStruct *font=NULL;
int font_height,font_width;
XSizeHints size_hints;
int darkcolor,lightcolor,black,white;
long translate[256];    /* for converting colors   */
double norm=1/(1.*NSITES);

int main ( int argc, char **argv )
{
  int i,j,k;
  double phipeg,w,theta,theta0;
  unsigned int width, height;
  XEvent report;
  progname=argv[0];
  pi=4*atan(1.0);
  for (i=0;i<NSITES;i++)
  {
    cn[i]=cos(2*pi*i/(1.*NSITES));
    sn[i]=sin(2*pi*i/(1.*NSITES));
  }

/* initial configuration */
  for (i=0;i<NSITES;i++)
  {
    fphi[i]=fcpi[i]=0.;
  }
  fphi[15]=NSITES/8.;
  fphi[17]=NSITES/8.;
  ftox();
  makedispersion();

  openwindow(argc,argv);
  makebuttons();
  xcenter=playwidth/2;
  ycenter=playheight/2;
  vx=0;
  vy=0;

/* loop forever, looking for events */

  while(1)
  {if ((0==XPending(display))&&(0==paused))
     update();
   else
    {XNextEvent(display,&report); 
     switch (report.type)
      {case Expose:
        if ((report.xexpose.window)!=window) break; /* cuts down flashing,
           but you might remove this line if things aren't being redrawn */
        if (report.xexpose.count!=0) break; /* more in queue, wait for them */
        repaint();  
        break;
       case ConfigureNotify:
        width=report.xconfigure.width;
        height=report.xconfigure.height;
        if ((width<size_hints.min_width)||(height<size_hints.min_height))
            {fprintf(stderr,"%s: window too small to proceed.\n",progname);
             cleanup();
            } 
        else if ((width!=windowwidth)||(height!=windowheight))
            {windowwidth=width;
             windowheight=height;
             playwidth=windowwidth;
             playheight=windowheight-3*BUTTONHEIGHT-1;
             makebuttons();
             xcenter=playwidth/2;
             ycenter=playheight/2;
            } 
        break; 
       case ButtonPress:
        if (report.xbutton.window==quitbutton)
             cleanup();
        else if (report.xbutton.window==pausebutton)
            {paused=1-paused;
             drawbutton(pausebutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "pause",1-2*paused);
            }
        else if (report.xbutton.window==dampbutton)
            {damped=1-damped;
             drawbutton(dampbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "damp",1-2*damped);
            } 
        else if (report.xbutton.window==blobbutton)
            {blob=1-blob;
             drawbutton(blobbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "blob",1-2*blob);
             drawit();
            } 
        else if (report.xbutton.window==rightbutton)
            {
              for (i=NSITES/2;i<NSITES;i++)
	       {fphi[i]=0;
                fcpi[i]=0;
               }  
             fphi[0]=fcpi[0]=0;
             ftox();
             drawit();
            } 
        else if (report.xbutton.window==vacuumbutton)
            {
             for (i=0;i<NSITES;i++)
              phi[i]=cpi[i]=0.0;
             xtof();
             drawit();
            } 
        else if (report.xbutton.window==reversebutton)
            {reversed=1-reversed;
             vx=-vx;
             vy=-vy;
             for (i=0;i<NSITES;i++)
              {cpi[i]*=-1;
	      }
             xtof();
             drawbutton(reversebutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "reverse",1-2*reversed);
             drawit();
            } 
        else if (report.xbutton.window==driftbutton)
            {drift=1-drift;
             if (drift)
	       {vx=1;
                vy=2;
	       }
             else
	       {vx=vy=0;
                xcenter=playwidth/2;
                ycenter=playheight/2;
	       }
             drawbutton(driftbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "drift",1-2*drift);
             drawit();
            } 
        else if (report.xbutton.window==slowbutton)
            {slow=1-slow;
             drawbutton(slowbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,
               "slow",1-2*slow);
             makedispersion();
            } 
        else if (report.xbutton.window==packetbutton)
            {
             for (i=0;i<NSITES;i++)
               fphi[i]=fcpi[i]=0.;
             fphi[15]=NSITES/8.;
             fphi[17]=NSITES/8.;
             ftox();
             drawit();
            } 
        else if (report.xbutton.window==equbutton)
	  {equation=report.xbutton.x/BUTTONWIDTH;
           for (i=0;i<3;i++)
             drawbutton(equbutton,i*BUTTONWIDTH,0,BUTTONWIDTH,BUTTONHEIGHT,
               equtext[i],1-2*(equation==i));
           makedispersion();
	  }
        else /* a mouse click in main window */
         {i=report.xbutton.x;
          j=report.xbutton.y;
          if (j>playheight) update(); /* outside of playfield */
          else                        /* perturb field */
	   {if (blob) /* do some trig */
             {phipeg=(i-xcenter)*(i-xcenter)
                    +(j-ycenter)*(j-ycenter);
              phipeg=4*sqrt(phipeg)/(1.*playheight)-.60;
              if (j!=ycenter)
                theta0=atan((xcenter-i)/(1.0*(ycenter-j)));
              else
                theta0=pi/2.0;
              if (j>ycenter) /* put theta0 between 0 and 2 pi */
                theta0+=pi;
              if (theta0<0.) theta0+=2*pi;
             }
            else  /* not displayed as a blob */
             {theta0=2*pi*i/(1.*playwidth);
              phipeg=(j-playheight/2)/(0.5*playheight);
             }
             /* now construct the perturbed field */
            for (k=0;k<NSITES;k++)
             {theta=2*pi*k/(NSITES-1);
              w=sin(.5*(theta0-theta));
              w*=(w*100);
              w=1/(1+w);
              phi[k]=phipeg*w+phi[k]*(1-w);
             }    
            xtof(); 
            drawit();
	   }
         }
        break;
       default:
        break;
      } /* end of switch */
    } /* end of if XPending */
  } /* end of while(1) */
} /* end of main */

void update()
{int i;
 double temp;
  /* a kludgey replacement for usleep(): */
/* if (slow) */
  {timeout.tv_sec=0;
   timeout.tv_usec=50000;
   select(0,NULL,NULL,NULL,&timeout); 
  }
 if (damped)
  for (i=0;i<NSITES;i++)
   fcpi[i]*=0.98;
 for (i=0;i<NSITES;i++)
   {temp   =fcpi[i]*comega[i]+fphi[i]*somega[i];
    fphi[i]=fphi[i]*comega[i]-fcpi[i]*somega[i];
    fcpi[i]=temp;
   } 
 ftox();
 drawit();
 return;
}

void drawit()
{
  if ( blob ) 
  {
    drawblob();
  }
  else 
  {
    drawwave();
  }
  return;
}

void drawwave()
{
  int i,iy;
  color(CYAN);
  color2(MAGENTA);
  for (i=0;i<NSITES;i++)
  {iy=(playheight/2)*(1+phi[i]);
   if (iy>=playheight) iy=playheight-1;
   else if (iy<0) iy=0;
   xpoints[i].x=(i*playwidth)/(NSITES-1);
   xpoints[i].y=iy;
  }
  xpoints[NSITES].x=playwidth;
  xpoints[NSITES].y=playheight-1;
  xpoints[NSITES+1].x=0;
  xpoints[NSITES+1].y=playheight-1;
  XFillPolygon(display,window,gc1,xpoints,NSITES+2,Nonconvex,CoordModeOrigin); 

  xpoints[NSITES].y=0;
  xpoints[NSITES+1].y=0;
  XFillPolygon(display,window,gc2,xpoints,NSITES+2,Nonconvex,CoordModeOrigin);
  return;
}

void drawblob()
{int i,ix,iy;
 double r;
 color(CYAN);
 color2(MAGENTA);
 xcenter+=vx;
 ycenter+=vy;

 for (i=0;i<NSITES;i++)
  {r=.15+.25*phi[i];
   iy=ycenter-playheight*r*cn[i];
   ix= xcenter-playheight*r*sn[i];
   if (iy>(playheight-1)) iy=playheight-1, vy=-abs(vy);
   if (iy< 0            ) iy=0           , vy= abs(vy);
   if (ix>playwidth     ) ix=playwidth   , vx=-abs(vx);
   if (ix<0             ) ix=0           , vx= abs(vx);
   xpoints[i].x=ix;
   xpoints[i].y=iy;
  }
  /* draw inside, allowing for boundary if imposed */
 XFillPolygon(display,window,gc2,xpoints,NSITES,Nonconvex,
   CoordModeOrigin);
  /* draw outside of region */ 
 xpoints[NSITES].x=xpoints[0].x;
 xpoints[NSITES].y=xpoints[0].y;
 xpoints[NSITES+1].x=xcenter;
 xpoints[NSITES+1].y=0;
 xpoints[NSITES+2].x=playwidth;
 xpoints[NSITES+2].y=0;
 xpoints[NSITES+3].x=playwidth;
 xpoints[NSITES+3].y=playheight-1;
 xpoints[NSITES+4].x=0;
 xpoints[NSITES+4].y=playheight-1;
 xpoints[NSITES+5].x=0;
 xpoints[NSITES+5].y=0;
 xpoints[NSITES+6].x=xcenter;
 xpoints[NSITES+6].y=0;
 XFillPolygon(display,window,gc1,xpoints,NSITES+7,Nonconvex,
   CoordModeOrigin);
 return;
}

void repaint()
   /* this fixes the window up whenever it is uncovered */
{int i;
 color(BROWN);
 XFillRectangle(display,window,gc1,0,playheight,
              playwidth,windowheight-playheight);
 drawbutton(pausebutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"pause",1-2*paused);
 drawbutton(dampbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"damp",1-2*damped);        
 drawbutton(quitbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"quit",1);
 drawbutton(vacuumbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"vacuum",1);
 drawbutton(reversebutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"reverse",1-2*reversed);
 drawbutton(blobbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"blob",1-2*blob);
 drawbutton(rightbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"right",1);
 drawbutton(driftbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"drift",1-2*drift);
 drawbutton(slowbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"slow",1-2*slow);
 drawbutton(packetbutton,0,0,BUTTONWIDTH,BUTTONHEIGHT,"packet",1);
 for (i=0;i<3;i++)
  drawbutton(equbutton,i*BUTTONWIDTH,0,BUTTONWIDTH,BUTTONHEIGHT,
       equtext[i],1-2*(equation==i));
 color(BLACK);
 line(0,playheight,windowwidth,playheight);
 drawit();
 return;
}

void openwindow(argc,argv)
/* a lot of this is taken from basicwin in the Xlib Programming Manual */
int argc;
char **argv;
{char *window_name="Wave tank";
 char *icon_name="Wave";
 long event_mask;
 Pixmap icon_pixmap;
 char *display_name=NULL;
 int i;
# define icon_bitmap_width 16
# define icon_bitmap_height 16
 static char icon_bitmap_bits[] = {
   0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0xf8,
   0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

/* open up the display */
 if ((display=XOpenDisplay(display_name))==NULL)
  {fprintf(stderr,"%s: cannot connect to X server %s\n",
    progname,XDisplayName(display_name));
   exit(-1);
  }
 screen=DefaultScreen(display);
 cmap=DefaultColormap(display,screen);

 darkcolor=black=BlackPixel(display,screen);
 lightcolor=white=WhitePixel(display,screen);
 translate[WHITE]=white;
 translate[BLACK]=black;
 translate[2]=lightcolor;
 translate[3]=darkcolor;
 for(i=4;i<256;i++)
   translate[i]=translate[i%4]; 

if (DefaultDepth(display,screen)>1)
{
 if (XAllocNamedColor(display,cmap,"salmon",&colorcell,&xcolor))
              darkcolor=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
              lightcolor=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"red",&colorcell,&xcolor))
              translate[RED]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"blue",&colorcell,&xcolor))
              translate[BLUE]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"green",&colorcell,&xcolor))
              translate[GREEN]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"cyan",&colorcell,&xcolor))
              translate[CYAN]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"orange",&colorcell,&xcolor))
              translate[ORANGE]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"purple",&colorcell,&xcolor))
              translate[PURPLE]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"yellow",&colorcell,&xcolor))
              translate[YELLOW]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"pink",&colorcell,&xcolor))
              translate[PINK]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"brown",&colorcell,&xcolor))
              translate[BROWN]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"grey",&colorcell,&xcolor))
              translate[GREY]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"turquoise",&colorcell,&xcolor))
              translate[TURQUOISE]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"gold",&colorcell,&xcolor))
              translate[GOLD]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"magenta",&colorcell,&xcolor))
              translate[MAGENTA]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"navy",&colorcell,&xcolor))
              translate[NAVY]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"tan",&colorcell,&xcolor))
              translate[TAN]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"violet",&colorcell,&xcolor))
              translate[VIOLET]=colorcell.pixel;
}
  /* feel free to type in more colors, I got bored */
 for(i=18;i<256;i++)
   translate[i]=translate[i%18];

    /* make the main window */
 windowwidth=6*BUTTONWIDTH;
 windowheight=600;
 playwidth=windowwidth;
 playheight=windowheight-3*BUTTONHEIGHT-1;
 window=XCreateSimpleWindow(display,RootWindow(display,screen),
   0,0,windowwidth,windowheight,4,translate[14],lightcolor);
    /* make the icon */
 icon_pixmap=XCreateBitmapFromData(display,window,
   icon_bitmap_bits,icon_bitmap_width,icon_bitmap_height);

 size_hints.flags=PPosition | PSize | PMinSize;
 size_hints.min_width=windowwidth;
 size_hints.min_height=100;
#ifdef X11R3
 size_hints.x=x;
 size_hints.y=y;
 size_hints.width=windowwidth;
 size_hints.height=windowheight;
 XSetStandardProperties(display,window,window_name,icon_name,
    icon_pixmap,argv,argc,&size_hints);
#else
 {XWMHints wm_hints;
  XClassHint class_hints;
  XTextProperty windowName, iconName;
  if (XStringListToTextProperty(&window_name,1,&windowName)==0)
   {fprintf(stderr,"%s: structure allocation for windowName failed.\n"
      ,progname);
    exit(-1);
   }
  if (XStringListToTextProperty(&icon_name,1,&iconName)==0)
   {fprintf(stderr,"%s: structure allocation for iconName failed.\n"
       ,progname);
    exit(-1);
   }
  wm_hints.initial_state=NormalState;
  wm_hints.input=True;
  wm_hints.icon_pixmap=icon_pixmap;
  wm_hints.flags=StateHint|IconPixmapHint|InputHint;
  class_hints.res_name=progname;
  class_hints.res_class="Basicwin";
  XSetWMProperties(display,window,&windowName,&iconName,
       argv,argc,&size_hints,&wm_hints,&class_hints);
 }
#endif

   /* pick the events to look for */
 event_mask=ExposureMask|ButtonPressMask|StructureNotifyMask;
 XSelectInput(display,window,event_mask);
   /* pick font: 9x15 is supposed to almost always be there */
 if ((font=XLoadQueryFont(display,"9x15"))==NULL) 
   if ((font=XLoadQueryFont(display,"fixed"))==NULL)
     {fprintf(stderr,"%s: Sorry, having font problems.\n",progname);
      exit(-1);
     }
 font_height=font->ascent+font->descent;
 font_width=font->max_bounds.width;

  /* make graphics context: */
 gc1=XCreateGC(display,window,0,NULL);
 XSetFont(display,gc1,font->fid);
 XSetForeground(display,gc1,black);
 XSetBackground(display,gc1,lightcolor); 

 gc2=XCreateGC(display,window,0,NULL);
 XSetFont(display,gc2,font->fid);
 XSetForeground(display,gc2,black);
 XSetBackground(display,gc2,lightcolor); 

   /* show the window */
 XMapWindow(display,window);
 return;
}

void makebuttons()
{
     /* first destroy any old buttons */
 XDestroySubwindows(display,window);
     /* now make the new buttons */
 quitbutton=makebutton(0,windowheight-BUTTONHEIGHT,BUTTONWIDTH,BUTTONHEIGHT);
 vacuumbutton=makebutton(3*BUTTONWIDTH,windowheight-BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 pausebutton=makebutton(BUTTONWIDTH,windowheight-BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 dampbutton=makebutton(2*BUTTONWIDTH,windowheight-BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 reversebutton=makebutton(windowwidth-2*BUTTONWIDTH,windowheight-3*BUTTONHEIGHT
                       ,BUTTONWIDTH,BUTTONHEIGHT);
 blobbutton=makebutton(windowwidth-BUTTONWIDTH,windowheight-3*BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 rightbutton=makebutton(windowwidth-2*BUTTONWIDTH,windowheight-2*BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 driftbutton=makebutton(windowwidth-BUTTONWIDTH,windowheight-2*BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 slowbutton=makebutton(windowwidth-BUTTONWIDTH,windowheight-BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 equbutton=makebutton(0,windowheight-3*BUTTONHEIGHT,
                        3*BUTTONWIDTH,BUTTONHEIGHT);
 packetbutton=makebutton(windowwidth-2*BUTTONWIDTH,windowheight-BUTTONHEIGHT,
                        BUTTONWIDTH,BUTTONHEIGHT);
 return;
}

void cleanup()
{XUnloadFont(display,font->fid);
 XFreeGC(display,gc1); 
 XCloseDisplay(display);
 exit(1);
} 

Window makebutton(xoffset,yoffset,xsize,ysize)
int xoffset,yoffset,xsize,ysize;
/* Puts button of specified dimensions on window.  Nothing is drawn. */
{Window buttonwindow;
 long event_mask;
 buttonwindow=XCreateSimpleWindow(display,window,xoffset,yoffset,
        xsize,ysize,0,black,lightcolor);
 event_mask=ButtonPressMask|ExposureMask; /* look for mouse-button presses */
 XSelectInput(display,buttonwindow,event_mask);
 XMapWindow(display,buttonwindow);
 return buttonwindow;
}

void drawbutton(buttonwindow,xoffset,yoffset,xsize,ysize,text,state)
Window buttonwindow;
int xoffset,yoffset,xsize,ysize,state;
char * text;
/* Draw a button in buttonwindow of specified dimensions with text
   centered.  Color of button determined by sign of "state." 
   size of border determined by magnitude. */
{int textlength,i,j;
 int cdark,clight,cup,cdown;
 int cleft,cright,cbutton,ctext;
 cup=lightcolor;
 cdown=darkcolor;
 cdark=black;
 clight=white;
 if (state<0)
  {cbutton=cdown;
   ctext=clight;
   cleft=cdark;
   cright=clight;
  } 
 else
  {cbutton=cup;
   ctext=cdark;
   cleft=clight;
   cright=cdark;
  }
 j=abs(state);
 XSetForeground(display,gc1,cbutton);  
 XFillRectangle(display,buttonwindow,gc1,xoffset+j,yoffset+j,
        xsize-2*j,ysize-2*j);
 XSetForeground(display,gc1,cleft);
 XFillRectangle(display,buttonwindow,gc1,xoffset,yoffset,xsize,j);
 XFillRectangle(display,buttonwindow,gc1,xoffset,yoffset,j,ysize);
 XSetForeground(display,gc1,cright);
 for (i=0;i<j;i++) 
  {XDrawLine(display,buttonwindow,gc1,
     xoffset+i,yoffset+ysize-i-1,xoffset+xsize-i-1,yoffset+ysize-i-1);
   XDrawLine(display,buttonwindow,gc1,
     xoffset+xsize-i-1,yoffset+i,xoffset+xsize-i-1,yoffset+ysize-i-1);
  }
 if (NULL!=text)
   {textlength=strlen(text);
    XSetForeground(display,gc1,ctext);
    XDrawString(display,buttonwindow,gc1,
          xoffset+(xsize-textlength*font_width)/2,
          yoffset+(ysize+font_height)/2-2,
          text,textlength);
   }
 XSetForeground(display,gc1,black);
 return;
}

void xtof()
{fastfourier(phi,cpi,fphi,fcpi,NSITES,1,1,1);
 return;
}

void ftox()
{int i;
 fastfourier(fphi,fcpi,phi,cpi,NSITES,1,1,-1);
 for (i=0;i<NSITES;i++)
   {phi[i]*=norm;
    cpi[i]*=norm;
   }
 return;
}

void makedispersion()
/* set up dispersion relations for the models */
{int j;
 double k;
 for (j=0;j<NSITES;j++)
    {if (2*j<NSITES) k=j; 
     else k=NSITES-j;
     if (equation==WATER)
      omega[j]=4*sqrt(k);
     else
      omega[j]=sqrt(k*k+NSITES*NSITES*equation/200);
     comega[j]=cos(dt*omega[j]);
     somega[j]=sin(dt*omega[j]);
    }
 return;
}

void fastfourier(double *rex,double *imx, double *ref,double *imf,
                 int n,int d1,int d2,int direction)
/* recursively calls itself on two halves of the vector */
/* (rex[],imx[]) is the input vector
   (ref[],imf[]) is the output
   n=vector length
   d1=spacing between elements of input vector
   d2= same for output
   direction= +/-1 for forward and reverse fourier transforms
   vector is not normalized */  
{int j,Ndn,nd2,i1,i2;
 double c,s,temp1,temp2,temp3;
 if (n==2) /* do two elements case by hand */
   {*ref=*rex+*(rex+d1);
    *imf=*imx+*(imx+d1);
    *(ref+d2)=*rex-*(rex+d1);
    *(imf+d2)=*imx-*(imx+d1);
    return;
   }
 nd2=n>>1; /* half the vector length */
  /* recursively call fastfourier */
 fastfourier(rex   ,imx   ,ref       ,imf       ,nd2,2*d1,d2,direction);
 fastfourier(rex+d1,imx+d1,ref+nd2*d2,imf+nd2*d2,nd2,2*d1,d2,direction);
  /* combine the results */
 Ndn=NSITES/n;
 for (j=0;j<nd2;j++)
  {c=cn[j*Ndn];
   s=direction*sn[j*Ndn];
   i1=d2*j;
   i2=d2*(j+nd2);
   temp1  =ref[i1]+c*ref[i2]-s*imf[i2];
   temp2  =imf[i1]+c*imf[i2]+s*ref[i2];
   temp3  =ref[i1]-c*ref[i2]+s*imf[i2];
   imf[i2]=imf[i1]-c*imf[i2]-s*ref[i2];
   ref[i1]=temp1;
   imf[i1]=temp2;
   ref[i2]=temp3;
  }
 return;
}








