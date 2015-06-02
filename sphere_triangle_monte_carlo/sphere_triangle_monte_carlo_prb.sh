#!/bin/bash
#
gcc -c -I/$HOME/include sphere_triangle_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_triangle_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc sphere_triangle_monte_carlo_prb.o /$HOME/libc/$ARCH/sphere_triangle_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_triangle_monte_carlo_prb.o."
  exit
fi
#
rm sphere_triangle_monte_carlo_prb.o
#
mv a.out sphere_triangle_monte_carlo_prb
./sphere_triangle_monte_carlo_prb > sphere_triangle_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_triangle_monte_carlo_prb."
  exit
fi
rm sphere_triangle_monte_carlo_prb
#
echo "Program output written to sphere_triangle_monte_carlo_prb_output.txt"
