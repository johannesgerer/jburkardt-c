#!/bin/bash
#
gcc -c -g asa245_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa245_prb.o /$HOME/libc/$ARCH/asa245.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa245_prb.o."
  exit
fi
#
rm asa245_prb.o
#
mv a.out asa245_prb
./asa245_prb > asa245_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa245_prb."
  exit
fi
rm asa245_prb
#
echo "Program output written to asa245_prb_output.txt"
