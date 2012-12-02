#
#  On Chili and Kachina, this is where the C compiler is:
#
CC = /usr/bin/gcc

c_simple.out: c_simple
	./c_simple > c_simple.out
c_simple: c_simple.o midpoint.o f.o
	$(CC) -o c_simple c_simple.o midpoint.o f.o -lm
c_simple.o: c_simple.c c_simple.h
	$(CC) -c c_simple.c
midpoint.o: midpoint.c c_simple.h
	$(CC) -c midpoint.c
f.o: f.c c_simple.h
	$(CC) -c f.c
clean:
	/bin/rm c_simple.out c_simple c_simple.o midpoint.o f.o



