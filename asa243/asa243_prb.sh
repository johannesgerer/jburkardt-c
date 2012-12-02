#!/bin/bash
#
gcc -c -g asa243_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa243_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa243_prb.o /$HOME/libc/$ARCH/asa243.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa243_prb.o."
  exit
fi
#
rm asa243_prb.o
#
mv a.out asa243_prb
./asa243_prb > asa243_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa243_prb."
  exit
fi
rm asa243_prb
#
echo "Program output written to asa243_prb_output.txt"
