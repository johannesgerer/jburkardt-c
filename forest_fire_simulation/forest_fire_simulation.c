 /*        Xfires, a forest fire simulator for X windows.
           
                       Michael Creutz   
                   creutz@wind.phy.bnl.gov 

 compiling:  cc -O -o xfires -L/usr/X11R6/lib xfires.c -lX11 

 version of 24 Nov 95
 LZW disabled, 10 October 1998

 minor modifications to eliminate gcc -Wall warnings, March 1999

The latest version is kept at "http://thy.phy.bnl.gov/www/xtoys/xtoys.html"
A text description of this program is there as well.
 */

/* because of the unisys patent on LZW, it may be illegal to
   distribute this program with the following "if" modified to "while"*/
#define LZW if

# include <X11/Xlib.h>
# include <X11/Xutil.h>
# include <X11/Xos.h>
# include <X11/Xatom.h> 
# include <X11/cursorfont.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

/* size and position parameters for buttons, etc. */
# define PLAYTOP 116
# define PLAYLEFT 14
# define BUTTONWIDTH 68
# define BBWIDTH 96
# define SRHEIGHT (2*18)
# define BPW (sizeof(long))

     /* lattice dimensions; ncols will be truncated to 
        a multiple of bytes per long word (BPW),
        and then two units are lost for boundaries    */
int nrows,ncols,volume,mask; 
int block,           /* size of cells displayed on the screen          */
    blockshift;      /* log_2(block)                                   */
unsigned char *field[2]={NULL,NULL}; /* pointers to old and new system */
/* 
  which field is current?            
*/
int old = 0;
int new = 1;
                 
# define barren 0
# define tree 1
# define fire 2

int paused=0,        /* is the system paused?                  */
    iteration=0;     /* counting sweeps                        */
long mrand48(),lrand48(); /* generates a random word, lrand is non-negative */
double drand48(); /* generates a random double */

/* 
  for fast but crude random number generation 
*/
# define RSIZE 127
long births[RSIZE]; 
int randomizer=57,birthindex=0;

char stringbuffer[256];              /* generally useful       */
  /* things for regulating the updating speed */
int speed=0; /* updating delay proportional to speed   */
int delay=0; /* counter for implementing speed control */
struct timeval timeout; /* timer for speed control */

/* various window stuff */
Display *display;
int screen;
static char *progname;
Window window,quitbutton,pausebutton,playground,blockbutton,srbutton,
    speedbutton,makebutton();
XColor xcolor,colorcell;
Colormap cmap;
GC gc,gcpen;
int windowwidth,windowheight;
XFontStruct *font=NULL;
int font_height,font_width;
XSizeHints size_hints;
int darkcolor,lightcolor,black,white;
XImage *spinimage=NULL;
long translate[256];                /* for converting colors   */
char * srtext[2]={"save","restore"};
void drawbutton();
void openwindow();
void makebuttons();
void update();
void repaint(),
     cleanup(),showpic(),fixboundary(),loadpic(),savepic(),compress(),
     decompress(),inittable();

int main ( int argc, char *argv[] )
{
  unsigned int height;
  int i;
  XEvent report;
  unsigned int width;
  unsigned int x;
  unsigned int y;

  progname = argv[0];
/* 
  Set the initial lattice size.
*/
  block = 1;
  nrows = 200 / block;
  ncols = (~(BPW-1)) & (200 / block );
  volume = nrows * ncols;
/* 
  Set up array for fast random number generation.
*/
  for (i=0;i<RSIZE;i++)
  {
    births[i]=(lrand48()>>2);
  }
/* 
  Open the window and make the buttons.
*/
  openwindow ( argc, argv );

  makebuttons();

/* 
  Loop forever, looking for events.
*/
  while(1)
    {if (0==paused)
      {if (delay) /* don't update yet */
       {delay--;
        /* this use of select() seems a kludge to me; 
           why can't usleep() be more standard? */
        timeout.tv_sec=0;
        timeout.tv_usec=100000; /* .1 sec per delay unit */
        select(0,NULL,NULL,NULL,&timeout); 
       }
     else
      {delay=speed;
       update();
      }
    }
   if (paused|XPending(display))
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
             XSetForeground(display,gcpen,lightcolor);  
             XFillRectangle(display,window,gcpen,0,0,windowwidth,windowheight);
             nrows=(height-PLAYTOP-10)/block;
             ncols=(~(BPW-1))&((width-PLAYLEFT-16)/block);
             volume=nrows*ncols;
             makebuttons();
             showpic();
            } 
        break; 
       case ButtonPress:
        if (report.xbutton.window==quitbutton)
             cleanup();
        else if (report.xbutton.window==pausebutton)
            {paused=1-paused;
             drawbutton(pausebutton,0,0,BUTTONWIDTH,18,"pause",1-2*paused);
            } 
        else if (report.xbutton.window==blockbutton)
          {i=(1<<(4*report.xbutton.x/BBWIDTH));
           if (i!=block) /* reset for new block size */
             {block=i;
              nrows=(windowheight-PLAYTOP-10)/block;
              ncols=(~(BPW-1))&((windowwidth-PLAYLEFT-16)/block);
              volume=nrows*ncols;
              XSetForeground(display,gcpen,lightcolor);  
              XFillRectangle(display,window,gcpen,0,0,
                             windowwidth,windowheight);
              makebuttons();
              showpic();
             }
          }
        else if (report.xbutton.window==playground) /* start a fire */
            {x=report.xbutton.x/block;
             y=report.xbutton.y/block;
             field[old][x+ncols*y]=field[new][x+ncols*y]=fire;
             fixboundary();
             showpic();
            }
        else if (report.xbutton.window==srbutton) /* save or restore */
            {y=(report.xbutton.y*2)/SRHEIGHT;
             drawbutton(srbutton,0,y*(SRHEIGHT/2),BUTTONWIDTH,SRHEIGHT/2,
                  srtext[y],-1);
             if (0==y)
                savepic(field[old],ncols,nrows);
             else if (1==y) /* load in gif image */
              {loadpic(field[old],ncols,nrows);
               for (i=0;i<volume;i++)
                  field[new][i]=1&(field[old][i]&=3);
               showpic();
              }
             drawbutton(srbutton,0,y*(SRHEIGHT/2),BUTTONWIDTH,SRHEIGHT/2,
                  srtext[y],1);
            }
        else if (report.xbutton.window==speedbutton)
          {  /* reset speed */
           speed=10-(11*report.xbutton.x)/BBWIDTH;
           if (speed<0) speed=0;
           if (speed>10) speed=10;
           delay=speed;
           drawbutton(speedbutton,0,0,BBWIDTH,18,"speed",-1);
           drawbutton(speedbutton,1+((10-speed)*(BBWIDTH-2))/11,1,
                      (BBWIDTH-2)/11,16,"",2);
          }
        else
/* 
  do a sweep when mouse clicked 
*/
         update(); 
        break;
       default:
        break;
      } /* end of switch */
    } /* end of if XPending */
  } /* end of while(1) */
  return 0;
}

void update ( void )
{
  int i;
  int nfires;
  int ntrees;
  int newtrees = volume / 32;
  unsigned char *pold;
  unsigned char *pnew;
  unsigned char *pup;
  unsigned char *pdown;
  unsigned char *pleft;
  unsigned char *pright;
  unsigned char *ptop;
/* 
  On each cycle, we add volume/32 trees.
*/
  while ( newtrees )
  {

    newtrees--;
    i = volume;

    while ( volume <= i )
    {
      i = mask & births[birthindex];
      births[birthindex]^=births[randomizer];
      if ( RSIZE <= (++birthindex) )
      {
        birthindex = 0;
      }
      if ( RSIZE <= (++randomizer) ) 
      {
        randomizer = 0;
      }
    }
 
    if ( field[old][i] != fire ) 
    {
      field[new][i]=field[old][i]=tree;
    }

  } 
/* 
  Spread fires. 
*/
  pnew = field[new]+ncols;
  ptop = field[old]+volume-ncols;
  for (pold=field[old]+ncols; pold<ptop; pold++,pnew++)
  {
    if ( fire == *pold )
    {
      pright=pnew+1;
      pleft=pnew-1;
      pup=pnew-ncols;
      pdown=pnew+ncols;
      *pright=fire & (*pright | (*pright<<1)); 
      *pleft =fire & (*pleft  | (*pleft <<1));
      *pup   =fire & (*pup    | (*pup   <<1));
      *pdown =fire & (*pdown  | (*pdown <<1));
      *pold=*pnew=barren; /* old fire dies */
    }
  }
/* 
  spread fires from top and bottom boundaries 
*/
 for (i=0;i<ncols;i++)
  if (fire==field[old][i]) 
   {if (tree==field[new][i+ncols]) field[new][i+ncols]=fire;
    field[old][i]=field[new][i]=barren;
   }

  for (i=volume-ncols;i<volume;i++)
  {
    if (fire==field[old][i]) 
    {
      if ( tree == field[new][i-ncols] ) 
      {
        field[new][i-ncols]=fire;
      }
      field[old][i]=field[new][i]=barren;
    }
  }
/*
  Swap screens.
*/
  old = new;
  new = 1 - new;

  fixboundary();

  showpic(); 

  iteration++;

  if ( 25 <= iteration )
  {
    ntrees=0;
    nfires=0;
    for (i=ncols;i<volume-ncols;i++)
    {
      ntrees+=(field[old][i]==tree);
      nfires+=(field[old][i]==fire);
    }
    iteration=0;
    sprintf(stringbuffer,"%d fires, %d trees  ",nfires,ntrees);
    XDrawImageString(display,window,gc,8,82,stringbuffer,strlen(stringbuffer));
    if (0==nfires) field[old][(int) (volume*drand48())]=fire;
  }

  return;
}

void showpic ( void ) 
/* 
  display the field 
*/
{
  int row,col,i1,i2,color,j,j1,j2,blocktop=block;
  char *picture=(*spinimage).data;
  if (block>4) blocktop=block-1;
  if (8==(*spinimage).depth) 
   {if (block>1) /* I wish I knew how to do this faster */
     for (row=0;row<nrows;row++)
      for (col=0;col<ncols;col++)
       {color=translate[field[old][row*ncols+col]];
        j=block*(col+block*ncols*row);
        if (color!=picture[j])
         for (i1=0;i1<blocktop;i1++)
          {j1=i1*block*ncols+j;
           for (i2=0;i2<blocktop;i2++)
            picture[j1+i2]=color;
          }
       }
   else
    {for (j=0;j<volume;j++)
      picture[j]=translate[field[old][j]];
    }
  }
 else /* depth is not 8, use xputpixel (this is really ugly) */
  {if (block>1) /* I wish I knew how to do this faster */
    for (row=0;row<nrows;row++)
     for (col=0;col<ncols;col++)
      {color=translate[field[old][row*ncols+col]];
       if (color!=XGetPixel(spinimage,j1=block*col,j2=block*row))
         for (i2=0;i2<blocktop;i2++)
           for (i1=0;i1<blocktop;i1++)
            XPutPixel(spinimage,j1+i1,j2+i2,color);
      }
   else
    for (row=0;row<nrows;row++)
     {j1=row*ncols; 
      for (col=0;col<ncols;col++)
        XPutPixel(spinimage,col,row,translate[field[old][j1+col]]);
     }
  }
 XPutImage(display,playground,gc,spinimage,0,0,0,0,block*ncols,block*nrows); 
 return;
}

void fixboundary ( void )
/* copies edges for periodicity */
{
  int i;
  for (i=0;i<ncols;i++)
  {
    field[old][i]=field[old][volume-2*ncols+i];
    field[old][volume-ncols+i]=field[old][ncols+i];
  }
  for (i=0;i<volume;i+=ncols)
  {
    field[old][i]=field[old][i+ncols-2];
    field[old][i+ncols-1]=field[old][i+1];
  }
  return;
} 

void repaint ( void )
/* this fixes the window up whenever it is uncovered */
{
  int i;

 drawbutton(pausebutton,0,0,BUTTONWIDTH,18,"pause",1-2*paused);
 drawbutton(quitbutton,0,0,BUTTONWIDTH,18,"quit",1);
 drawbutton(window,PLAYLEFT-4,PLAYTOP-4,block*ncols+8,block*nrows+8,NULL,4);
 drawbutton(speedbutton,0,0,BBWIDTH,18,"speed",-1);
 drawbutton(speedbutton,1+((10-speed)*(BBWIDTH-2))/11,1,
               (BBWIDTH-2)/11,16,"",2);
 for (i=0;i<2;i++)
  drawbutton(srbutton,0,i*(SRHEIGHT/2),BUTTONWIDTH,SRHEIGHT/2,srtext[i],1);
  /* fix the block buttons */
 for (i=1;i<=4;i++)
  {sprintf(stringbuffer,"%d",(1<<i)/2);
   drawbutton(blockbutton,(i-1)*BBWIDTH/4,0,BBWIDTH/4,18,
         stringbuffer,1-2*(i==(blockshift+1)));
  } 
  /* write various strings */
 XDrawString(display,window,gc,24,62,"cell size",9);
 sprintf(stringbuffer,"%d by %d lattice  ",ncols-2,nrows-2);
 XDrawImageString(display,window,gc,14,100,
       stringbuffer,strlen(stringbuffer));
 XDrawString(display,window,gc,windowwidth-40,108
       ,"MJC",3); 
 showpic();  
 return;
}
 
void openwindow ( int argc, char *argv[] )
/* a lot of this is taken from basicwin in the Xlib Programming Manual */
{
  char *window_name="Forest fires";
  char *icon_name="fires";
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
/*
  open up the display 
*/
 if ((display=XOpenDisplay(display_name))==NULL)
  {fprintf(stderr,"%s: cannot connect to X server %s\n",
    progname,XDisplayName(display_name));
   exit(-1);
  }
 screen=DefaultScreen(display);
 cmap=DefaultColormap(display,screen);
 darkcolor=black=BlackPixel(display,screen);
 lightcolor=white=WhitePixel(display,screen);
 translate[0]=white;
 translate[1]=black;
 translate[2]=lightcolor;
 translate[3]=darkcolor;
 if (XAllocNamedColor(display,cmap,"firebrick",&colorcell,&xcolor))
              darkcolor=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
              lightcolor=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"black",&colorcell,&xcolor))
              translate[barren]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"forest green",&colorcell,&xcolor))
              translate[tree]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"yellow",&colorcell,&xcolor))
              translate[fire]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"grey",&colorcell,&xcolor))
              translate[3]=colorcell.pixel;

   /* fill out the color table for future uses */
 for(i=4;i<256;i++)
   translate[i]=translate[i%4]; 
   /* make the main window */
 windowwidth=(block*ncols+PLAYLEFT+16);
 windowheight=(PLAYTOP+block*nrows+10);
 window=XCreateSimpleWindow(display,RootWindow(display,screen),
   0,0,windowwidth,windowheight,4,black,lightcolor);
/* 
  Make the icon. 
*/
 icon_pixmap = XCreateBitmapFromData ( display, window,
   icon_bitmap_bits, (unsigned int ) icon_bitmap_width, ( unsigned int ) icon_bitmap_height );

   /* hints for window manager */
 size_hints.flags=PPosition | PSize | PMinSize;
 size_hints.min_width=windowwidth;
 size_hints.min_height=windowheight-160;
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

/* make graphics contexts: 
      gc for black on white
      gcpen for varying purposes */
   
 gc=XCreateGC(display,window,0,NULL);
 XSetFont(display,gc,font->fid);
 XSetForeground(display,gc,black);
 XSetBackground(display,gc,lightcolor); 

 gcpen=XCreateGC(display,window,0,NULL); 
 XSetFont(display,gcpen,font->fid);
 XSetForeground(display,gcpen,darkcolor);
 XSetBackground(display,gcpen,lightcolor); 

/* 
  show the window 
*/
  XMapWindow(display,window);
  return;
}

void makebuttons ( void )
{
  int i;
  long event_mask;
  XEvent report;
  Cursor cursor;
/*
  first destroy any old buttons 
*/
  XDestroySubwindows(display,window);
/* 
  now make the new buttons 
*/
  quitbutton  =makebutton(4,4,BUTTONWIDTH,18);
  pausebutton =makebutton(BUTTONWIDTH+8,4,BUTTONWIDTH,18);
  blockbutton =makebutton(16,32,BBWIDTH,18);
  srbutton    =makebutton(windowwidth-BUTTONWIDTH-4,4,BUTTONWIDTH,18*2);
  speedbutton=makebutton(windowwidth-BBWIDTH-8,
           46,
           BBWIDTH,18);
  playground=XCreateSimpleWindow(display,window,
    PLAYLEFT,PLAYTOP,block*ncols,block*nrows,0,black,translate[3]);
  event_mask=ExposureMask|ButtonReleaseMask|ButtonPressMask|
            PointerMotionHintMask|ButtonMotionMask; 
  XSelectInput(display,playground,event_mask);
  XMapWindow(display,playground);
/* 
  wait for playground to be displayed befor proceeding 
*/
  i=1; /* a flag */
  while (i)
    {XNextEvent(display,&report); 
     switch (report.type)
      {case Expose:
          if (report.xexpose.window!=playground) i=0;
       default:
        break;
      }
    }
/* 
  make image structure 
*/
 if (NULL!=spinimage) 
  {XDestroyImage(spinimage); 
   spinimage=NULL;
  }
 spinimage=XGetImage((Display *) display, (Drawable) playground,
           0,0,block*ncols,block*nrows,AllPlanes,ZPixmap);
 if (NULL==spinimage)
       {fprintf(stderr,"trouble creating image structure\n");
        exit(-1);
       } 
/* 
  make special cursors to be cute 
*/
 cursor=XCreateFontCursor(display,XC_sb_up_arrow); 
 XDefineCursor(display,playground,cursor);
 cursor=XCreateFontCursor(display,XC_hand2);
 XDefineCursor(display,blockbutton,cursor);
 XDefineCursor(display,quitbutton,cursor);
 XDefineCursor(display,pausebutton,cursor);
 XDefineCursor(display,srbutton,cursor);
 XDefineCursor(display,speedbutton,cursor);

/* 
  reallocate various arrays 
*/
 for (i=0;i<2;i++)
  {if (NULL!=field[i])  free((char *) field[i]);
   if (NULL==( field[i]= (unsigned char *) malloc(volume)))
     {fprintf(stderr,"allocation problems\n");
      cleanup();
     } 
  } 
/* 
  set blockshift to log_2 of block 
*/
 blockshift=0;
 i=block;
 while(i>>=1)
  blockshift++; 
/* 
  mask used for random site selection 
*/
 mask=1;
 while (mask<volume) mask=1|(mask<<1);
/* 
  set initial state random 
*/
 for (i=0;i<volume;i++)
  field[old][i]=field[new][i]=lrand48()&tree;

/* 
  draw everything 
*/
  repaint();
  return;
}

Window makebutton ( int xoffset, int yoffset, int xsize, int ysize )

/* Puts button of specified dimensions on window.  Nothing is drawn. */
{
  Window buttonwindow;
  long event_mask;
  buttonwindow=XCreateSimpleWindow(display,window,xoffset,yoffset,
        xsize,ysize,0,black,lightcolor);
 event_mask=ButtonPressMask|ExposureMask; /* look for mouse-button presses */
 XSelectInput(display,buttonwindow,event_mask);
 XMapWindow(display,buttonwindow);
 return buttonwindow;
}

void drawbutton ( Window buttonwindow, int xoffset, int yoffset, int xsize, 
  int ysize, char *text, int state )

/* 
  Draw a button in buttonwindow of specified dimensions with text
  centered.  Color of button determined by sign of "state." 
  size of border determined by magnitude. 
*/
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
 XSetForeground(display,gcpen,cbutton);  
 XFillRectangle(display,buttonwindow,gcpen,xoffset+j,yoffset+j,
        xsize-2*j,ysize-2*j);
 XSetForeground(display,gcpen,cleft);
 XFillRectangle(display,buttonwindow,gcpen,xoffset,yoffset,xsize,j);
 XFillRectangle(display,buttonwindow,gcpen,xoffset,yoffset,j,ysize);
 XSetForeground(display,gcpen,cright);
 for (i=0;i<j;i++) 
  {XDrawLine(display,buttonwindow,gcpen,
     xoffset+i,yoffset+ysize-i-1,xoffset+xsize-i-1,yoffset+ysize-i-1);
   XDrawLine(display,buttonwindow,gcpen,
     xoffset+xsize-i-1,yoffset+i,xoffset+xsize-i-1,yoffset+ysize-i-1);
  }
 if (NULL!=text)
   {textlength=strlen(text);
    XSetForeground(display,gcpen,ctext);
    XDrawString(display,buttonwindow,gcpen,
          xoffset+(xsize-textlength*font_width)/2,
          yoffset+(ysize+font_height)/2-2,
          text,textlength);
   }
 return;
}

void cleanup ( void )
{
  XUnloadFont(display,font->fid);
  XFreeGC(display,gc); 
  XFreeGC(display,gcpen); 
  XCloseDisplay(display);
  XDestroyImage(spinimage);

  if ( field[0] != NULL ) 
  {
    free((char *) field[0]);
  }

  if ( field[1] != NULL ) 
  {
    free((char *) field[1]);
  }

  exit ( 1 );
} 

/* from here on is the stuff for saving and restoring from a gif file */
/* for info on how gif works, see ftp://network.ucsd.edu/graphics/GIF.shar.Z */

char *picturename="xfires.gif";

void loadpic ( unsigned char *data, int xsize, int ysize )  
/* 
  load GIF image to field  
*/
{
  int i,j,filesize,gwidth,gheight,gvolume;
 unsigned char *ptr, *ptr1, *rawgif;
 int colorbits,codesize;
 FILE *infile;
 if (NULL==(infile=fopen(picturename,"r")))
   {fprintf(stderr,"couldn't open input file\n");
    return;
   }
    /* find the file size */
 fseek(infile, 0L, 2);
 filesize = ftell(infile);
 fseek(infile, 0L, 0);
  /* make a place in memory for the file */
 if (NULL==(rawgif= (unsigned char *) malloc(filesize) ))
    {fprintf(stderr, "not enough memory to read gif file\n");
     return;
    }
 ptr=rawgif;
  /* read in the file */ 
 if (fread(ptr, filesize, 1, infile) != 1)
    {fprintf(stderr, "read failed\n");
     free((char *) rawgif);
     return;
    }
 fclose(infile);
  /* check for GIF signature */
 if (strncmp((char *) ptr,"GIF87a", 6))
     {fprintf(stderr, "not a GIF87a file\n");
      free((char *) rawgif);
      return;
     }
 ptr+=6;
 ptr+=4;                    /* skip over screen size  */
 colorbits=1+((*ptr)&0xf);    /* how many bits of color */
 ptr+=2;                    /* skip over background   */ 
 if (*ptr++)                /* should be zero         */
   {fprintf(stderr, "corrupt GIF file\n");
    free((char *) rawgif);
    return;
   }
 ptr+=(3*(1<<colorbits));   /* skip over colormap     */
 if (','!=(*ptr++))         /* image separator        */
   {fprintf(stderr, "corrupt GIF file\n");
    free((char *) rawgif);
    return;
   }
 ptr+=4;                    /* skip over image offset */
 gwidth=(*ptr)+0x100*(*(ptr+1));
 ptr+=2;
 gheight=(*ptr)+0x100*(*(ptr+1));
 ptr+=2;
 if (0x80&(*ptr++))           /* skip over local color map */
   ptr+=(3*(1<<colorbits));
      /* should now be at start of data */
 codesize=(*ptr++);
  /* make a place for the decompressed file */
 gvolume=gwidth*gheight;
 ptr1=(unsigned char *) malloc(gvolume);
 decompress(codesize,ptr,ptr1,gvolume);
 free((char *) rawgif);
  /* map picture into data, allowing for different dimensions */ 
 for (j=0;j<ysize;j++)
  {if (j>=gheight) break;
   for (i=0;i<xsize;i++)
    {if (i>=gwidth) break;
     data[i+j*xsize]=ptr1[i+j*gwidth];
    }
  }
 free((char *) ptr1);
 fixboundary();
 return;
}

void savepic ( unsigned char *data, int xsize, int ysize )  
/* 
  save the field as a GIF image 
*/
{
  int i;
  int colorbits=5,codesize=5; /* assume a 32 color image */
  FILE *outfile;
  if (NULL==(outfile=fopen(picturename,"w")))
   {fprintf(stderr,"couldn't open output file\n");
    return;
   }
/* 
  GIF signature 
*/
  fwrite("GIF87a",6,1,outfile);
/* 
  screen descriptor 
*/
  stringbuffer[0]=xsize&0xff;   /* screen width  */
  stringbuffer[1]=(xsize>>8)&0xff;
  stringbuffer[2]=ysize&0xff;   /* screen height */
  stringbuffer[3]=(ysize>>8)&0xff;
  stringbuffer[4]=(0x80)             /* M=1; global color map follows */
                |((colorbits-1)<<4) /* -1+ bits of color reslution   */
                |(colorbits-1);     /* -1+bits per pixel in image    */
  stringbuffer[5]=0; /* background color */
  stringbuffer[6]=0; /* should be zero   */

  fwrite(stringbuffer,7,1,outfile);
/* 
  global color map 
*/
  for (i=0;i<(1<<colorbits);i++)
  {colorcell.pixel=translate[i];
   XQueryColor(display,cmap,&colorcell);
   fputc(colorcell.red>>8,outfile);
   fputc(colorcell.green>>8,outfile);
   fputc(colorcell.blue>>8,outfile);
  }
/* 
  image descriptor 
*/
  stringbuffer[0]=','; /* image descriptor separator */
  stringbuffer[1]=0;   /* image offset               */
  stringbuffer[2]=0;
  stringbuffer[3]=0;
  stringbuffer[4]=0;
  stringbuffer[5]=xsize&0xff;         /* image width  */
  stringbuffer[6]=(xsize>>8)&0xff;
  stringbuffer[7]=ysize&0xff;         /* image height */
  stringbuffer[8]=(ysize>>8)&0xff;
  stringbuffer[9]=0; /* use global color map, no interlace  */ 
  fwrite(stringbuffer,10,1,outfile);
/* 
  start of image data 
*/
  fputc(codesize,outfile); 
  compress(codesize,data,outfile,volume); 
/* 
  gif terminator 
*/
  fputc(';',outfile);
  fclose(outfile);
  return;
}

/* LZW compression */
/* hash function assumes TABLELENGTH is a power of 2 */
# define TABLELENGTH (1<<13)
char **addresses=NULL;          /* where to find the string         */
int  *codes=NULL,               /* the code value                   */
     *linktonext=NULL,          /* the next index in the hash chain */
     *lengths=NULL,             /* the length of the coded string   */
     *codeindex=NULL;           /* the index for a given code       */ 
int nextcode;                   /* the next unused code             */

/* hashit is supposed to give a unique fairly random number in the table for
each length a and string b */

# define hashit(a,b) (51*a+53*(57*b[0]+59*(61*b[a-1]+b[a>>1])))&(TABLELENGTH-1)

void compress ( initcodesize, ptr, outfile, size )
int initcodesize; /* the initial compression bits */
char * ptr;       /* where the data comes from    */
FILE * outfile;   /* where the output goes        */
int size;         /* how much data                */
{int currentcode,prefixcode=0,codesize,maxbits=12,maxcode;
 int clearcode,eoicode,currentplace=0,length,blocksize=0,bitoffset;
 int findcode();
 unsigned long outputword;
 unsigned char blockbuffer[256]; /* to hold data blocks before writing */
  /* allocate space for hash tables */
 if (NULL==(codes=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentplace=size;
  }
 if (NULL==(linktonext=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentplace=size;
  }
 if (NULL==(lengths=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentplace=size;
  }
   /* need one extra place in codeindex for overflow before resetting: */
 if (NULL==(codeindex=(int *) malloc(sizeof(int)*4097)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentplace=size;
  }
 if (NULL==(addresses=(char **) malloc(sizeof(char *)*TABLELENGTH)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentplace=size;
  }
   /* set up initial code table */
 inittable(initcodesize);
 clearcode=(1<<initcodesize);
 eoicode=clearcode+1;
 codesize=initcodesize+1;
 maxcode=1<<codesize;
 nextcode=eoicode+1;
  /* start with a clear code */
 outputword=clearcode;
 bitoffset=codesize; 
   /* now do the compressing */
 while (currentplace<size)
  {   /* check if codesize needs increasing */
   if (nextcode>maxcode)
    {codesize++;
     maxcode=1<<codesize;
     /* if too big, then reset compressor */
     if (codesize>maxbits)
      {if (bitoffset) outputword|=(clearcode<<bitoffset);
       else outputword=clearcode;
       bitoffset+=maxbits;
       inittable(initcodesize);
       codesize=initcodesize+1;
       maxcode=1<<codesize;
       nextcode=eoicode+1;
      }
    }
    /* look for an unstored string */
   length=1;
   LZW (nextcode>
            (currentcode=findcode(length,(char *)(ptr+currentplace))))
    {prefixcode=currentcode;
     length++;
     if ((currentplace+length)>=size) break;
    }
   nextcode++;
   currentplace+=(length-1);
    /* output the prefix code */
   if (bitoffset)  outputword|=(prefixcode<<bitoffset);
   else outputword=prefixcode;
   bitoffset+=codesize;
    /* output finished bytes to blocks */
   while (bitoffset>=8)
     {blockbuffer[blocksize]=outputword&0xff;
      outputword>>=8;
      bitoffset-=8;
      blocksize++;
        /* output filled block */
      if (blocksize>=254)
        {fputc((char) blocksize, outfile);
         fwrite(blockbuffer,blocksize,1,outfile);
         blocksize=0;
        }
     }
  }
   /* output the end of information code */
 if (bitoffset) outputword|=(eoicode<<bitoffset);
 else outputword=eoicode;
 bitoffset+=codesize;
   /* finish outputting the data */
 while (bitoffset>=0) 
     {blockbuffer[blocksize]=(char) (outputword&0xff);
      outputword>>=8;
      bitoffset-=8;
      blocksize++;
      if (blocksize>=254)
        {fputc((char) blocksize, outfile);
         fwrite(blockbuffer,blocksize,1,outfile);
         blocksize=0;
        }
     }
   /* output the last block */
 if (blocksize) 
   {fputc((char) blocksize, outfile);
    fwrite(blockbuffer,blocksize,1,outfile);
   }
  /* a final zero block count */
 fputc(0, outfile);
  /* deallocate tables */
 if (NULL!=codes) free((char *) codes);
 if (NULL!=linktonext) free((char *) linktonext);
 if (NULL!=lengths) free((char *) lengths);
 if (NULL!=codeindex) free((char *) codeindex);
 if (NULL!=addresses) free((char *) addresses);
 codes=linktonext=lengths=codeindex=NULL;
 addresses=(char **) NULL;
 return;
}

void decompress ( int initcodesize, unsigned char *ptr, unsigned char *ptr1, int size )
/*
  Parameters:

    int initcodesize;

    unsigned char *ptr, *ptr1, compressed data from ptr go to ptr1

    int size, an upper limit purely as a check.
*/
{
  int i,currentcode,codesize=0,maxbits=12,blocksize;
 int clearcode,eoicode,codemask=0;
 int bitoffset=0,indx,oldindx=0;
 int currentplace=0,oldplace=0;
 int findcode();
 unsigned long inputword=0;
 unsigned char *p1, *p2;
   /* first deblock the data */
 p1=p2=ptr;
 blocksize=(*p1++);
 while (blocksize)
  {while (blocksize--)
    (*p2++)=(*p1++); /* a wonderful example of how abstruse C can be */
   blocksize=(*p1++);
  }
   /* set up initial code table */
 currentcode=clearcode=(1<<initcodesize);
 eoicode=clearcode+1;
  /* allocate space for hash table */
 if (NULL==(codes=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"decompress: trouble allocating tables\n");
   currentcode=eoicode;
  }
 if (NULL==(linktonext=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"decompress: trouble allocating tables\n");
   currentcode=eoicode;
  }
 if (NULL==(lengths=(int *) malloc(sizeof(int)*TABLELENGTH)))
  {fprintf(stderr,"decompress: trouble allocating tables\n");
   currentcode=eoicode;
  }
   /* need one extra place in codeindex for overflow before resetting: */
 if (NULL==(codeindex=(int *) malloc(sizeof(int)*4097)))
  {fprintf(stderr,"compress: trouble allocating tables\n");
   currentcode=eoicode;
  }
 if (NULL==(addresses=(char **) malloc(sizeof(char*)*TABLELENGTH)))
  {fprintf(stderr,"decompress: trouble allocating tables\n");
   currentplace=eoicode;
  }
 while (currentcode!=eoicode)
  {if (currentcode==clearcode) /* reset the decompressor */
      {inittable(initcodesize);
       codesize=initcodesize+1;
       nextcode=eoicode+1;
       codemask=(1<<codesize)-1;
       oldindx=(-1);
      }
   else /* code represents data */
    {indx=codeindex[currentcode]; /* where in table is currentcode */
     if (indx>=0) /* it is there */
      { /* put it into the output */
       for (i=0;i<lengths[indx];i++)
         ptr1[currentplace+i]=addresses[indx][i];
       if (oldindx>=0) /* first character treated differently */
        {findcode(lengths[oldindx]+1,(char *) (ptr1+oldplace));
         nextcode++;  /* add new code to table */
        }
       oldplace=currentplace;
       currentplace+=lengths[indx];
       oldindx=indx;
      }
     else /* not in table yet; must be old code plus last=first character */
      {for (i=0;i<lengths[oldindx];i++)
         ptr1[currentplace+i]=addresses[oldindx][i];
       ptr1[currentplace+lengths[oldindx]]=addresses[oldindx][0];
              /* store new code */
       findcode(lengths[oldindx]+1,(char *)(ptr1+currentplace)); 
       oldplace=currentplace;
       currentplace+=lengths[oldindx]+1;
       oldindx=codeindex[nextcode++];
      }
      /* crude error checking */
     if ((oldindx<0)||(currentplace>size))
      {fprintf(stderr,"gif file appears to be corrupt\n");
       break;
      }
    } 
     /* check if codesize needs increasing */     
   if (nextcode>codemask)
      if (codesize<maxbits)
        {codesize++;
         codemask=(1<<codesize)-1;
        }
   while (bitoffset<codesize) /* read some more data */
    {if (bitoffset) inputword|=(((int)(*ptr++))<<bitoffset);
     else inputword=(*ptr++);
     bitoffset+=8;
    } 
    /* strip off current code */
   currentcode=inputword&codemask;
   inputword>>=codesize;
   bitoffset-=codesize;
   if (currentcode>nextcode) 
      {fprintf(stderr,"gif file appears to be corrupt\n");
       break;
      }
  }
  /* deallocate tables */
 if (NULL!=codes) free((char *) codes);
 if (NULL!=linktonext) free((char *) linktonext);
 if (NULL!=lengths) free((char *) lengths);
 if (NULL!=codeindex) free((char *) codeindex);
 if (NULL!=addresses) free((char *) addresses);
 codes=linktonext=lengths=codeindex=NULL;
 addresses=(char **) NULL;
 return;
}

void inittable ( int size )

{
  int i,findcode();
  for (i=0;i<TABLELENGTH;i++)
  {linktonext[i]=(-1);
   codes[i]=(-1);
   lengths[i]=(-1);
  }
  for (i=0;i<4096;i++)
   codeindex[i]=(-1);
/* store initial codes for raw characters */
  nextcode=0;
  for (i=0;i<(1<<size);i++)
  {stringbuffer[i]=i;
   findcode(1,(char *) (stringbuffer+i));
   nextcode++;
  }
  return;
}

int findcode ( int length, char *string )
/* 
  return code for string of given length;  
  if not found, store it and return nextcode 
*/
{
  int i;
  int j;
  int indx;
  int previousindex;

  indx=hashit(length,string); 
/* 
  look for string in table 
*/
  previousindex = indx;

  while ( indx > 0 )
  {
/*
  Is the length right?
*/
    if (lengths[indx]==length)
    {
      for (j=0;j<length;j++)
      {
        if ((string[j]) != (addresses[indx][j])) break;
      }
 
      if ( j == length ) return codes[indx]; 
    }
    previousindex=indx;
    indx=linktonext[indx];
  } 
/* 
  not found, so store it 
*/
  indx=previousindex;
  i=indx;
/*
  Find an unused slot in the table.
*/
  while ( codes[i]>=0 )
  {
    i++;
    if (i>=TABLELENGTH) i-=TABLELENGTH;
  }
/*
  Link to it, and move to it.
*/
  if (i!=indx) 
  {
    linktonext[indx]=i;
    indx=i;
  }

  codes[indx]=nextcode; /* save the new code */
  lengths[indx]=length;
  addresses[indx]=string;
  codeindex[nextcode]=indx; 

  return nextcode;
}  
