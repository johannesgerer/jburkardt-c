#!/bin/bash
#
gcc -c -g asa091_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa091_prb.o /$HOME/libc/$ARCH/asa091.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa091_prb.o."
  exit
fi
#
rm asa091_prb.o
#
mv a.out asa091_prb
./asa091_prb > asa091_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa091_prb."
  exit
fi
rm asa091_prb
#
echo "Program output written to asa091_prb_output.txt"
