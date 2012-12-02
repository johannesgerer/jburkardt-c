#!/bin/bash
#
gcc -c -g mandelbrot.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot.c"
  exit
fi
rm compiler.txt
#
gcc mandelbrot.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot.o"
  exit
fi
rm mandelbrot.o
#
mv a.out ~/binc/$ARCH/mandelbrot
#
echo "Executable installed as ~/binc/$ARCH/mandelbrot"
