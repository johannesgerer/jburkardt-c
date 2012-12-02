#!/bin/bash
#
gcc -c -g mandelbrot_ascii.c
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot_ascii.c"
  exit
fi
#
gcc mandelbrot_ascii.o
if [ $? -ne 0 ]; then
  echo "Errors loading mandelbrot_ascii.c"
  exit
fi
#
rm mandelbrot_ascii.o
mv a.out ~/binc/$ARCH/mandelbrot_ascii
#
echo "Executable installed as ~/binc/$ARCH/mandelbrot_ascii."
