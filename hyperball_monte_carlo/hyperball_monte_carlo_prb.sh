#!/bin/bash
#
gcc -c -g -I/$HOME/include hyperball_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc hyperball_monte_carlo_prb.o /$HOME/libc/$ARCH/hyperball_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hyperball_monte_carlo_prb.o."
  exit
fi
#
rm hyperball_monte_carlo_prb.o
#
mv a.out hyperball_monte_carlo_prb
./hyperball_monte_carlo_prb > hyperball_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hyperball_monte_carlo_prb."
  exit
fi
rm hyperball_monte_carlo_prb
#
echo "Program output written to hyperball_monte_carlo_prb_output.txt"
