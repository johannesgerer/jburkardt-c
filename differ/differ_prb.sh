#!/bin/bash
#
gcc -c -g -I/$HOME/include differ_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling differ_prb.c"
  exit
fi
rm compiler.txt
#
gcc differ_prb.o /$HOME/libc/$ARCH/differ.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading differ_prb.o."
  exit
fi
#
rm differ_prb.o
#
mv a.out differ_prb
./differ_prb > differ_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running differ_prb."
  exit
fi
rm differ_prb
#
echo "Program output written to differ_prb_output.txt"
