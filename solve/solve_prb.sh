#!/bin/bash
#
gcc -c -I/$HOME/include solve_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling solve_prb.c"
  exit
fi
rm compiler.txt
#
gcc solve_prb.o /$HOME/libc/$ARCH/solve.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading solve_prb.o."
  exit
fi
#
rm solve_prb.o
#
mv a.out solve_prb
./solve_prb > solve_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running solve_prb."
  exit
fi
rm solve_prb
#
echo "Program output written to solve_prb_output.txt"
