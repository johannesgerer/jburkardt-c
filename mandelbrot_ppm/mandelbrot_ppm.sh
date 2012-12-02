#!/bin/bash
#
gcc -c -g mandelbrot_ppm.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mandelbrot_ppm.c"
  exit
fi
rm compiler.txt
#
gcc mandelbrot_ppm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mandelbrot_ppm.o"
  exit
fi
rm mandelbrot_ppm.o
#
mv a.out ~/binc/$ARCH/mandelbrot_ppm
#
echo "Executable installed as ~/binc/$ARCH/mandelbrot_ppm"
