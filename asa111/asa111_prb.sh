#!/bin/bash
#
gcc -c -g asa111_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa111_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa111_prb.o /$HOME/libc/$ARCH/asa111.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa111_prb.o."
  exit
fi
#
rm asa111_prb.o
#
mv a.out asa111_prb
./asa111_prb > asa111_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa111_prb."
  exit
fi
rm asa111_prb
#
echo "Program output written to asa111_prb_output.txt"
