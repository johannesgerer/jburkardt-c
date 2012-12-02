#!/bin/bash
#
gcc -c -g -O2 main.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling main.c"
  exit
fi
rm compiler.txt
#
gcc -c -g -O2 boid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling boid.c"
  exit
fi
rm compiler.txt
#
gcc -c -g -O2 vec.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vec.c"
  exit
fi
rm compiler.txt
#
gcc -c -g -O2 x.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling x.c"
  exit
fi
rm compiler.txt
#
gcc -c -g -O2 xboids.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xboids.c"
  exit
fi
rm compiler.txt
#
gcc main.o boid.o vec.o x.o xboids.o -L/usr/X11/lib -lm -lX11 -lXext
if [ $? -ne 0 ]; then
  echo "Errors loading mri_to_pgm.c"
  exit
fi
#
rm *.o
mv a.out ~/binc/$ARCH/xboids
#
echo "Executable installed as ~/binc/$ARCH/xboids"
