CC = gcc
CFLAGS = -O3 -ffast-math -funroll-loops -ftree-vectorize

## Linux (depends on how you installed gsl, if you used your distro repository this will probably work)
GSL_INC =
GSL_LIB =

## MacOSX (assuming you used "sudo port install gsl")
# GSL_INC = -I/opt/local/include
# GSL_LIB = -L/opt/local/lib
