#!/bin/bash
#
gcc -c -I/$HOME/include legendre_polynomial_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_polynomial_prb.c"
  exit
fi
#
gcc legendre_polynomial_prb.o /$HOME/libc/$ARCH/legendre_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_polynomial_prb.o"
  exit
fi
#
rm legendre_polynomial_prb.o
#
mv a.out legendre_polynomial_prb
./legendre_polynomial_prb > legendre_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running legendre_polynomial_prb."
  exit
fi
rm legendre_polynomial_prb
#
echo "Program output written to legendre_polynomial_prb_output.txt"
