#!/bin/bash
#
gcc -c -g asa152_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa152_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa152_prb.o /$HOME/libc/$ARCH/asa152.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa152_prb.o."
  exit
fi
#
rm asa152_prb.o
#
mv a.out asa152_prb
./asa152_prb > asa152_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa152_prb."
  exit
fi
rm asa152_prb
#
echo "Program output written to asa152_prb_output.txt"
