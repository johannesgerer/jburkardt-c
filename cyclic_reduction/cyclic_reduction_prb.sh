#!/bin/bash
#
gcc -c -g cyclic_reduction_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cyclic_reduction_prb.c."
  exit
fi
rm compiler.txt
#
gcc cyclic_reduction_prb.o /$HOME/libc/$ARCH/cyclic_reduction.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cyclic_reduction_prb.o."
  exit
fi
#
rm cyclic_reduction_prb.o
#
mv a.out cyclic_reduction_prb
./cyclic_reduction_prb > cyclic_reduction_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cyclic_reduction_prb."
  exit
fi
rm cyclic_reduction_prb
#
echo "Program output written to cyclic_reduction_prb_output.txt"
