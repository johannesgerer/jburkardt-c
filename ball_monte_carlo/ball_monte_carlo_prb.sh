#!/bin/bash
#
gcc -c -g -I/$HOME/include ball_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc ball_monte_carlo_prb.o /$HOME/libc/$ARCH/ball_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_monte_carlo_prb.o."
  exit
fi
#
rm ball_monte_carlo_prb.o
#
mv a.out ball_monte_carlo_prb
./ball_monte_carlo_prb > ball_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ball_monte_carlo_prb."
  exit
fi
rm ball_monte_carlo_prb
#
echo "Program output written to ball_monte_carlo_prb_output.txt"
