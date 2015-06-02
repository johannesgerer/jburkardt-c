#!/bin/bash
#
gcc -c -O2 -I$HOME/include blas3_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3_prb.c"
  exit
fi
rm compiler.txt
#
gcc blas3_prb.o $HOME/libc/$ARCH/blas3.o $HOME/libc/$ARCH/blas0.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas3_prb.o."
  exit
fi
#
rm blas3_prb.o
#
mv a.out blas3_prb
./blas3_prb > blas3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas3_prb."
  exit
fi
rm blas3_prb
#
echo "Program output written to blas3_prb_output.txt"
