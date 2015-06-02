#!/bin/bash
#
gcc -c -I/$HOME/include ellipse_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc ellipse_monte_carlo_prb.o /$HOME/libc/$ARCH/ellipse_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipse_monte_carlo_prb.o."
  exit
fi
#
rm ellipse_monte_carlo_prb.o
#
mv a.out ellipse_monte_carlo_prb
./ellipse_monte_carlo_prb > ellipse_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipse_monte_carlo_prb."
  exit
fi
rm ellipse_monte_carlo_prb
#
echo "Program output written to ellipse_monte_carlo_prb_output.txt"
