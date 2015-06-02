#!/bin/bash
#
gcc -c -O2 -I$HOME/include blas2_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas2_prb.c"
  exit
fi
rm compiler.txt
#
gcc blas2_prb.o $HOME/libc/$ARCH/blas2.o $HOME/libc/$ARCH/blas0.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas2_prb.o."
  exit
fi
#
rm blas2_prb.o
#
mv a.out blas2_prb
./blas2_prb > blas2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas2_prb."
  exit
fi
rm blas2_prb
#
echo "Program output written to blas2_prb_output.txt"
