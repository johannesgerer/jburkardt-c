#!/bin/bash
#
gcc -c -g -I/$HOME/include cube_monte_carlo_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_monte_carlo_prb.c"
  exit
fi
rm compiler.txt
#
gcc cube_monte_carlo_prb.o /$HOME/libc/$ARCH/cube_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_monte_carlo_prb.o."
  exit
fi
#
rm cube_monte_carlo_prb.o
#
mv a.out cube_monte_carlo_prb
./cube_monte_carlo_prb > cube_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_monte_carlo_prb."
  exit
fi
rm cube_monte_carlo_prb
#
echo "Program output written to cube_monte_carlo_prb_output.txt"
