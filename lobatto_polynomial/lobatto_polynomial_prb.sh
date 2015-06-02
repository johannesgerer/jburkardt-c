#!/bin/bash
#
gcc -c -I/$HOME/include lobatto_polynomial_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling lobatto_polynomial_prb.c"
  exit
fi
#
gcc -o lobatto_polynomial_prb lobatto_polynomial_prb.o /$HOME/libc/$ARCH/lobatto_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lobatto_polynomial_prb.o."
  exit
fi
#
rm lobatto_polynomial_prb.o
#
./lobatto_polynomial_prb > lobatto_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lobatto_polynomial_prb."
  exit
fi
rm lobatto_polynomial_prb
#
echo "Program output written to lobatto_polynomial_prb_output.txt"
