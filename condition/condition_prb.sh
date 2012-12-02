#!/bin/bash
#
gcc -c -g -I/$HOME/include condition_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling condition_prb.c"
  exit
fi
rm compiler.txt
#
gcc condition_prb.o /$HOME/libc/$ARCH/condition.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading condition_prb.o."
  exit
fi
#
rm condition_prb.o
#
mv a.out condition_prb
./condition_prb > condition_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running condition_prb."
  exit
fi
rm condition_prb
#
echo "Program output written to condition_prb_output.txt"
