## Makefile for the "c_simple" example.
#
#  Discussion:
#
#    This makefile is used to manipulate the "c_simple" code.
#
#    make run
#      asks the makefile to do what is necessary to run the code.
#    make c_simple
#      creates the executable program, but does not run it.
#    make c_simple.o
#      compiles the c_simple.c source code.
#    make midpoint.o
#      compiles the midpoint.c source code.
#    make f.o
#      compiles the f.c source code.
#    make clean
#      removes the output file, executable program and compiled files.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    10 May 2013
#
#  Author:
#
#    John Burkardt
#
CC = /usr/bin/gcc

run: c_simple
	./c_simple > c_simple_output.txt
c_simple: c_simple.o midpoint.o f.o
	$(CC) -o c_simple c_simple.o midpoint.o f.o -lm
c_simple.o: c_simple.c c_simple.h
	$(CC) -c c_simple.c
midpoint.o: midpoint.c c_simple.h
	$(CC) -c midpoint.c
f.o: f.c c_simple.h
	$(CC) -c f.c
clean:
	rm -f c_simple_output.txt
	rm -f c_simple
	rm -f c_simple.o
	rm -f midpoint.o
	rm -f f.o

