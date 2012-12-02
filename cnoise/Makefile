include Config.mk

all: libcnoise.a

libcnoise.a: cnoise.o
	ar rcs libcnoise.a cnoise.o

cnoise.o: cnoise.c cnoise.h
	$(CC) $(CFLAGS) $(GSL_INC) -c cnoise.c -o cnoise.o

clean:
	rm -f libcnoise.a
	rm -f cnoise.o
