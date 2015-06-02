#!/bin/bash
#
gcc -c -I/$HOME/include polynomial_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling polynomial_prb.c"
  exit
fi
#
gcc polynomial_prb.o /$HOME/libc/$ARCH/polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polynomial_prb.o."
  exit
fi
#
rm polynomial_prb.o
#
mv a.out polynomial_prb
./polynomial_prb > polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polynomial_prb."
  exit
fi
rm polynomial_prb
#
echo "Program output written to polynomial_prb_output.txt"
