#!/bin/bash
#
gcc -c -g fd1d_bvp_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_bvp_prb.c."
  exit
fi
rm compiler.txt
#
gcc fd1d_bvp_prb.o /$HOME/libc/$ARCH/fd1d_bvp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_bvp_prb.o."
  exit
fi
#
rm fd1d_bvp_prb.o
#
mv a.out fd1d_bvp_prb
./fd1d_bvp_prb > fd1d_bvp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fd1d_bvp_prb."
  exit
fi
rm fd1d_bvp_prb
#
echo "Program output written to fd1d_bvp_prb_output.txt"
