#!/bin/bash
#
cp sphere_triangle_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include sphere_triangle_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_triangle_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv sphere_triangle_monte_carlo.o ~/libc/$ARCH/sphere_triangle_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/sphere_triangle_monte_carlo.o"
