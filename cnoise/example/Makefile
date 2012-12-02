include ../Config.mk

all: test

test: test.o
	gcc -L../ $(GSL_LIB) test.o -o test -lcnoise -lgsl

test.o: test.c
	gcc -O3 -c test.c $(GSL_INC) -o test.o

clean:
	rm -f test
	rm -f test.o
