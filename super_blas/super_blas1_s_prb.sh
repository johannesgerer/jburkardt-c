#!/bin/bash
#
gcc -c -g -I/$HOME/include super_blas1_s_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling super_blas1_s_prb.c"
  exit
fi
rm compiler.txt
#
gcc super_blas1_s_prb.o /$HOME/libc/$ARCH/super_blas.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading super_blas1_s_prb.o"
  exit
fi
#
rm super_blas1_s_prb.o
#
mv a.out super_blas1_s_prb
./super_blas1_s_prb > super_blas1_s_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running super_blas1_s_prb"
  exit
fi
rm super_blas1_s_prb
#
echo "Program output written to super_blas1_s_prb_output.txt"
