#!/bin/bash
#
gcc -c -g rectangle.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rectangle.c"
  exit
fi
rm compiler.txt
#
gcc rectangle.o $HOME/libc/$ARCH/toms886.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rectangle.o"
  exit
fi
rm rectangle.o
#
mv a.out rectangle
./rectangle > rectangle_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rectangle"
  exit
fi
rm rectangle
#
echo "Test results written to rectangle_output.txt."
