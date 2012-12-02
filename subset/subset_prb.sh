#!/bin/bash
#
gcc -c -g subset_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_prb.c."
  exit
fi
rm compiler.txt
#
gcc subset_prb.o /$HOME/libc/$ARCH/subset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subset_prb.o."
  exit
fi
#
rm subset_prb.o
#
mv a.out subset_prb
./subset_prb > subset_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subset_prb."
  exit
fi
rm subset_prb
#
echo "Program output written to subset_prb_output.txt"
