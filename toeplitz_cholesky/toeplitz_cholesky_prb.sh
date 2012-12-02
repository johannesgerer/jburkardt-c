#!/bin/bash
#
gcc -c -g -I/$HOME/include toeplitz_cholesky_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_cholesky_prb.c"
  exit
fi
rm compiler.txt
#
gcc toeplitz_cholesky_prb.o /$HOME/libc/$ARCH/toeplitz_cholesky.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toeplitz_cholesky_prb.o."
  exit
fi
#
rm toeplitz_cholesky_prb.o
#
mv a.out toeplitz_cholesky_prb
./toeplitz_cholesky_prb > toeplitz_cholesky_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toeplitz_cholesky_prb."
  exit
fi
rm toeplitz_cholesky_prb
#
echo "Program output written to toeplitz_cholesky_prb_output.txt"
