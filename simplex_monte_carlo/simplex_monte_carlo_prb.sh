#!/bin/bash
#
gcc -c -g -I/$HOME/include simplex_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc simplex_monte_carlo_prb.o /$HOME/libc/$ARCH/simplex_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simplex_monte_carlo_prb.o."
  exit
fi
#
rm simplex_monte_carlo_prb.o
#
mv a.out simplex_monte_carlo_prb
./simplex_monte_carlo_prb > simplex_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simplex_monte_carlo_prb."
  exit
fi
rm simplex_monte_carlo_prb
#
echo "Program output written to simplex_monte_carlo_prb_output.txt"
