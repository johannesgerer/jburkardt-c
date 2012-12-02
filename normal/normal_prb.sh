#!/bin/bash
#
gcc -c -g normal_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_prb.c."
  exit
fi
rm compiler.txt
#
gcc normal_prb.o /$HOME/libc/$ARCH/normal.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_prb.o."
  exit
fi
#
rm normal_prb.o
#
mv a.out normal_prb
./normal_prb > normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running normal_prb."
  exit
fi
rm normal_prb
#
echo "Program output written to normal_prb_output.txt"
