#!/bin/bash
#
gcc -c -g -I/$HOME/include line_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc line_monte_carlo_prb.o /$HOME/libc/$ARCH/line_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_monte_carlo_prb.o."
  exit
fi
#
rm line_monte_carlo_prb.o
#
mv a.out line_monte_carlo_prb
./line_monte_carlo_prb > line_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_monte_carlo_prb."
  exit
fi
rm line_monte_carlo_prb
#
echo "Program output written to line_monte_carlo_prb_output.txt"
