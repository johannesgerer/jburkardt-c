#!/bin/bash
#
gcc -c bvec_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling bvec_prb.c."
  exit
fi
#
gcc bvec_prb.o /$HOME/libc/$ARCH/bvec.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bvec_prb.o."
  exit
fi
#
rm bvec_prb.o
#
mv a.out bvec_prb
./bvec_prb > bvec_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bvec_prb."
  exit
fi
rm bvec_prb
#
echo "Program output written to bvec_prb_output.txt"
