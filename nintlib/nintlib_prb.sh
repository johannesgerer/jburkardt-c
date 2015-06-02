#!/bin/bash
#
gcc -c -g -I/$HOME/include nintlib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib_prb.c"
  exit
fi
rm compiler.txt
#
gcc nintlib_prb.o /$HOME/libc/$ARCH/nintlib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nintlib_prb.o."
  exit
fi
#
rm nintlib_prb.o
#
mv a.out nintlib_prb
./nintlib_prb > nintlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nintlib_prb."
  exit
fi
rm nintlib_prb
#
echo "Program output written to nintlib_prb_output.txt"
