#!/bin/bash
#
gcc -c -g -I/$HOME/include blend_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blend_prb.c"
  exit
fi
rm compiler.txt
#
gcc blend_prb.o /$HOME/libc/$ARCH/blend.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blend_prb.o."
  exit
fi
#
rm blend_prb.o
#
mv a.out blend_prb
./blend_prb > blend_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blend_prb."
  exit
fi
rm blend_prb
#
echo "Program output written to blend_prb_output.txt"
