#!/bin/bash
#
gcc -c -g -I/$HOME/include sparse_interp_nd_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_interp_nd_prb.c"
  exit
fi
rm compiler.txt
#
gcc sparse_interp_nd_prb.o /$HOME/libc/$ARCH/sparse_interp_nd.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_interp_nd_prb.o."
  exit
fi
#
rm sparse_interp_nd_prb.o
#
mv a.out sparse_interp_nd_prb
./sparse_interp_nd_prb > sparse_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_interp_nd_prb."
  exit
fi
rm sparse_interp_nd_prb
#
echo "Program output written to sparse_interp_nd_prb_output.txt"
