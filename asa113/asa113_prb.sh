#!/bin/bash
#
gcc -c -g asa113_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa113_prb.o /$HOME/libc/$ARCH/asa113.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa113_prb.o."
  exit
fi
#
rm asa113_prb.o
#
mv a.out asa113_prb
./asa113_prb > asa113_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa113_prb."
  exit
fi
rm asa113_prb
#
echo "Program output written to asa113_prb_output.txt"
