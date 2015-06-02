#!/bin/bash
#
gcc -c -g -I/$HOME/include jacobi_polynomial_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_polynomial_prb.c"
  exit
fi
rm compiler.txt
#
gcc jacobi_polynomial_prb.o /$HOME/libc/$ARCH/jacobi_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_polynomial_prb.o."
  exit
fi
#
rm jacobi_polynomial_prb.o
#
mv a.out jacobi_polynomial_prb
./jacobi_polynomial_prb > jacobi_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running jacobi_polynomial_prb."
  exit
fi
rm jacobi_polynomial_prb
#
echo "Program output written to jacobi_polynomial_prb_output.txt"
