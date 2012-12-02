#!/bin/bash
#
gcc -c -g -I/$HOME/include latinize_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latinize_prb.c"
  exit
fi
rm compiler.txt
#
gcc latinize_prb.o /$HOME/libc/$ARCH/latinize.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latinize_prb.o."
  exit
fi
#
rm latinize_prb.o
#
mv a.out latinize_prb
./latinize_prb > latinize_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latinize_prb."
  exit
fi
rm latinize_prb
#
echo "Program output written to latinize_prb_output.txt"
