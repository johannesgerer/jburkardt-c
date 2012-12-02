#
#gcc TEST.c -o TEST -lX11 -lXext -lIm -L/usr/X11/lib -I/usr/X11/include
#
#
# K 1997

XBOIDS_SRCS = main.c vec.c boid.c xboids.c x.c

XBOIDS_OBJS = $(XBOIDS_SRCS:.c=.o)

PROGS = xboids

#CC = gcc
CFLAGS = -g -O2
LFLAGS = -L/usr/X11/lib -lm -lX11 -lXext

all : xboids

xboids : $(XBOIDS_OBJS)
	$(CC) -o $@ $(XBOIDS_OBJS) $(LFLAGS)

clean :
	$(RM) $(PROGS) core *.o *~

